#pragma once
// newton/simulate.hpp — the fixed-step time loop that turns a Particle plus
// a force into a Trajectory.
//
// Physics: this file does not add any new physics beyond the four
// integrators in newton/integrators/ — it is the loop that calls one of
// them n times, records samples, and (in simulate_until) stops precisely
// when a user-supplied scalar "event" function changes sign. Everything
// interesting here is about NUMERICS, not mechanics: getting the time
// bookkeeping and the event-crossing refinement right is what separates a
// simulation you can trust to 1e-9 from one that is merely plausible.
//
// C++ note: fixed-step only, deliberately (docs/research/numerical-
// integration.md §5.1) — the two symplectic integrators in this library
// (semi_implicit_euler, velocity_verlet) only get their bounded-energy
// guarantee under a CONSTANT dt; an adaptive step size would silently throw
// that guarantee away. Reproducibility (same input -> bit-identical output)
// and trivial convergence testing are the other two reasons. Adaptive
// stepping, if ever added, belongs in a separate RK-family-only code path.

#include <cmath>       // std::floor, std::isfinite
#include <cstddef>     // std::size_t
#include <limits>      // std::numeric_limits
#include <stdexcept>   // std::invalid_argument, std::logic_error

#include <newton/core/particle.hpp>
#include <newton/core/trajectory.hpp>
#include <newton/integrators/euler.hpp>
#include <newton/integrators/rk4.hpp>
#include <newton/integrators/semi_implicit_euler.hpp>
#include <newton/integrators/velocity_verlet.hpp>

namespace newton {

/// Which one-step integrator simulate()/simulate_until() should use. See
/// newton/integrators/*.hpp for what each one costs and how it behaves —
/// RK4 is the default because it is the safest general-purpose choice for a
/// student exploring many different force laws, not because it is always
/// the best (see velocity_verlet.hpp for when Verlet is the better pick).
enum class Method { Euler, SemiImplicitEuler, VelocityVerlet, RK4 };

/// Options for a single simulate()/simulate_until() run.
struct SimOptions {
    double t0{0.0};                  ///< start time [s]
    double t_end{1.0};                ///< end time [s]
    double dt{1e-3};                  ///< step size [s] — no sensible default exists across
                                       ///< problems; pick it from the fastest timescale in your
                                       ///< force (docs/research/numerical-integration.md §3.1),
                                       ///< e.g. period/100 for an oscillator.
    Method method{Method::RK4};
    int    record_every{1};           ///< record every Nth step; the first and last states are
                                       ///< always recorded regardless of this stride, with no
                                       ///< duplicate entry when the stride would have recorded
                                       ///< them anyway.
};

namespace detail {

/// A run consists of zero or more full nominal-dt steps and, when the
/// interval is not an exact multiple of dt, one shorter terminal step.  The
/// latter is deliberate: silently rounding the step count can end before or
/// after the caller's requested t_end.
struct SimulationStepPlan {
    long long full_steps{0};
    double terminal_dt{0.0};

    [[nodiscard]] bool has_terminal_step() const { return terminal_dt > 0.0; }
    [[nodiscard]] long long total_steps() const {
        return full_steps + (has_terminal_step() ? 1 : 0);
    }
};

/// Validate the part of SimOptions that controls the time interval before
/// deriving a step count.  In particular, a negative dt used to turn into a
/// negative count, then an enormous vector reserve after conversion to
/// size_t; reject it at the API boundary instead.
inline void validate_simulation_time_options(const SimOptions& opt) {
    if (!std::isfinite(opt.t0) || !std::isfinite(opt.t_end)) {
        throw std::invalid_argument("newton::SimOptions requires finite t0 and t_end");
    }
    if (opt.t_end < opt.t0) {
        throw std::invalid_argument("newton::SimOptions requires t_end >= t0");
    }
    if (!std::isfinite(opt.dt) || !(opt.dt > 0.0)) {
        throw std::invalid_argument("newton::SimOptions requires a finite dt > 0");
    }
}

/// Split [t0, t_end] into full dt steps plus at most one shorter final step.
/// Full-step times are still formed directly from t0 + i*dt by the callers;
/// no clock is accumulated with repeated t += dt.
[[nodiscard]] inline SimulationStepPlan make_step_plan(const SimOptions& opt) {
    validate_simulation_time_options(opt);

    const double span = opt.t_end - opt.t0;
    if (!std::isfinite(span)) {
        throw std::invalid_argument("newton::SimOptions interval is too large to simulate");
    }

    const double quotient = span / opt.dt;
    if (!std::isfinite(quotient)) {
        throw std::invalid_argument("newton::SimOptions would require too many steps");
    }

    const double full_steps_as_double = std::floor(quotient);
    if (full_steps_as_double >= static_cast<double>(std::numeric_limits<long long>::max())) {
        throw std::invalid_argument("newton::SimOptions would require too many steps");
    }

    long long full_steps = static_cast<long long>(full_steps_as_double);
    double full_end = opt.t0 + static_cast<double>(full_steps) * opt.dt;
    if (!std::isfinite(full_end)) {
        throw std::invalid_argument("newton::SimOptions interval cannot be represented at this dt");
    }

    // A quotient that rounded infinitesimally upward can otherwise put the
    // nominal last step past t_end.  Back up one full step and make that
    // portion part of the final shortened step instead.
    if (full_end > opt.t_end) {
        --full_steps;
        full_end = opt.t0 + static_cast<double>(full_steps) * opt.dt;
        if (!std::isfinite(full_end) || full_end > opt.t_end) {
            throw std::invalid_argument("newton::SimOptions interval cannot be represented at this dt");
        }
    }

    const double terminal_dt = opt.t_end - full_end;
    if (terminal_dt < 0.0) {
        throw std::invalid_argument("newton::SimOptions interval cannot be represented at this dt");
    }
    return {full_steps, terminal_dt};
}

[[nodiscard]] inline bool opposite_sign(double a, double b) {
    return (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0);
}

constexpr int kEventBisectionIterations = 40;

/// Treat two event times separated by only a few representable time units as
/// the same sampled endpoint.  This is deliberately a *time* tolerance used
/// only after a real sign-change bracket exists, not a generic tolerance on
/// arbitrary user event values.  It prevents roundoff from turning a root at
/// a nominal step boundary into a second, nearly coincident final sample.
[[nodiscard]] inline bool same_time_to_roundoff(double a, double b, double step_dt) {
    constexpr double kEndpointSnapUlps = 64.0;
    const double scale = std::fmax(std::fabs(step_dt), std::fmax(std::fabs(a), std::fabs(b)));
    const double roundoff = kEndpointSnapUlps * std::numeric_limits<double>::epsilon() * scale;
    const double bisection_resolution = std::ldexp(std::fabs(step_dt), -kEventBisectionIterations);
    return std::fabs(a - b) <= roundoff + bisection_resolution;
}

/// Internal: dispatches to the requested integrator's step function. Not
/// part of the public API — both simulate() and simulate_until() need this
/// exact switch (the latter also needs it to re-step at a bisected dt), so
/// it lives here once instead of being copy-pasted twice.
///
/// C++ note: the trailing throw is unreachable in practice — Method has
/// exactly four enumerators and every one of them has a case above — but
/// the compiler cannot prove a switch over an `enum class` is exhaustive
/// (nothing stops a bad static_cast from producing a fifth value), so
/// -Wreturn-type demands a value on every path. Throwing here fails loudly
/// if that ever happens instead of silently fabricating a State.
template <class Force>
[[nodiscard]] inline State dispatch_step(Method method, const State& s, double mass,
                                          double dt, Force&& force) {
    switch (method) {
        case Method::Euler:             return integrators::euler_step(s, mass, dt, force);
        case Method::SemiImplicitEuler: return integrators::semi_implicit_euler_step(s, mass, dt, force);
        case Method::VelocityVerlet:    return integrators::velocity_verlet_step(s, mass, dt, force);
        case Method::RK4:               return integrators::rk4_step(s, mass, dt, force);
    }
    throw std::logic_error("newton::detail::dispatch_step: unknown Method");
}

}  // namespace detail

/// Fixed-step time loop: integrates a Particle under `force` from opt.t0 to
/// opt.t_end and returns the recorded Trajectory.
///
/// C++ note (CRITICAL numerics detail, research §7.2): time is computed
/// fresh each iteration as `t = t0 + i*dt`, NEVER accumulated with
/// `t += dt`. `dt` is generally not exactly representable in binary
/// floating point (0.001 isn't), so every `+=` commits a fresh rounding
/// error, and those errors compound step after step — over n steps the
/// drift can grow like a random walk (~eps*|t|*sqrt(n)) or, worst case,
/// linearly (~eps*|t|*n). Each xxx_step() function only ever adds dt ONCE
/// to whatever t it's handed, so if this loop just let a step's output t
/// feed into the next step's input t, the SAME accumulation would happen —
/// just spread invisibly across many function calls instead of one visible
/// `+=` line. The fix is to recompute the exact t from t0 and the integer
/// step index after every step, which caps the rounding error at ~eps*|t|
/// for every single sample instead of letting it compound. This is also why
/// the full-step count is an integer, `floor((t_end - t0)/dt)`, rather than
/// a `while (t < t_end)` loop — a float-compared loop condition can take a
/// different number of iterations on different compilers/platforms depending
/// on the last-bit rounding of `t`, which would make trajectories (and tests)
/// non-reproducible. A single shortened final step lands exactly on t_end
/// when the interval is not a multiple of dt.
template <class Force>
[[nodiscard]] Trajectory simulate(const Particle& p, Force&& force, const SimOptions& opt) {
    const detail::SimulationStepPlan plan = detail::make_step_plan(opt);
    const long long n = plan.total_steps();
    // Guard record_every against 0/negative input: it feeds a `% ` below,
    // and `x % 0` is undefined behaviour. A stride of "1" (record everything)
    // is the sane fallback for a nonsensical stride.
    const int stride = opt.record_every > 0 ? opt.record_every : 1;

    Trajectory traj;
    traj.mass = p.mass;
    traj.states.reserve(static_cast<std::size_t>(n / stride) + 2);  // +2: slack for the initial/final samples

    State s = p.initial_state(opt.t0);
    traj.states.push_back(s);  // index 0 (the initial state) is always recorded

    for (long long i = 0; i < n; ++i) {
        const bool is_terminal_step = plan.has_terminal_step() && i == plan.full_steps;
        const double step_dt = is_terminal_step ? plan.terminal_dt : opt.dt;
        s = detail::dispatch_step(opt.method, s, p.mass, step_dt, force);
        s.t = is_terminal_step ? opt.t_end
                               : opt.t0 + static_cast<double>(i + 1) * opt.dt;  // see the CRITICAL note above

        const bool is_last_step = (i + 1 == n);
        // i+1 ranges over 1..n here (index 0 was pushed once, above), so
        // this can never duplicate the initial sample; is_last_step covers
        // the final sample even when stride doesn't naturally land on it.
        if (is_last_step || (i + 1) % stride == 0) {
            traj.states.push_back(s);
        }
    }
    return traj;
}

/// Event-terminated run: integrates like simulate(), but stops the instant
/// `event(state)` changes sign between two consecutive steps.
///
/// Physics/numerics: naively reporting "the state at the first step where
/// the sign changed" would put the crossing somewhere inside the last
/// interval — an O(dt) error in, say, a projectile's range, no matter how
/// accurate the integrator is per step. Instead, once a sign change is
/// detected between the previous state (known to be on the event's starting
/// side) and the freshly-stepped next state, this function BISECTS ON THE
/// STEP SIZE dt' in (0, dt]: it re-integrates from the previous state with
/// a trial sub-step, checks the event's sign there, and narrows the
/// bracket, ~40 times (each iteration halves the bracket, so 40 iterations
/// resolve dt to a fraction of it far below any floating-point-meaningful
/// step size). Because every trial re-uses the REAL integrator at a shorter
/// step rather than linearly interpolating between the two known states,
/// the final crossing lands on event ≈ 0 to the integrator's own
/// truncation-error accuracy, not to O(dt) — e.g. an RK4 run lands on the
/// event to ~O(dt'^4), not ~O(dt).
///
/// If the event function never changes sign, this behaves exactly like
/// simulate() and runs to opt.t_end.
/// @param event  callable: double operator()(const State&) const — its sign
///               change (not its magnitude) is what stops the run.
template <class Force, class EventFn>
[[nodiscard]] Trajectory simulate_until(const Particle& p, Force&& force, const SimOptions& opt,
                                         EventFn&& event) {
    const detail::SimulationStepPlan plan = detail::make_step_plan(opt);
    const long long n = plan.total_steps();
    const int stride = opt.record_every > 0 ? opt.record_every : 1;

    Trajectory traj;
    traj.mass = p.mass;
    traj.states.reserve(static_cast<std::size_t>(n / stride) + 2);

    State s = p.initial_state(opt.t0);
    traj.states.push_back(s);
    double e_prev = event(s);
    // Starting on an event surface is common (a projectile launched from the
    // ground).  Arm endpoint detection only after the trajectory has first
    // left that surface, so the initial state is not treated as an impact.
    bool event_has_left_initial_surface = e_prev != 0.0;

    for (long long i = 0; i < n; ++i) {
        const State s_prev = s;  // the bisection below always re-steps from here
        const bool is_terminal_step = plan.has_terminal_step() && i == plan.full_steps;
        const double step_dt = is_terminal_step ? plan.terminal_dt : opt.dt;
        State s_next = detail::dispatch_step(opt.method, s_prev, p.mass, step_dt, force);
        s_next.t = is_terminal_step ? opt.t_end
                                    : opt.t0 + static_cast<double>(i + 1) * opt.dt;
        const double e_next = event(s_next);

        // An exactly-zero later endpoint is an event too.  It is handled
        // directly instead of bisected, which also avoids manufacturing a
        // nearly-zero extra interval in the recorded trajectory.
        if (event_has_left_initial_surface && e_next == 0.0) {
            traj.states.push_back(s_next);
            return traj;
        }

        // A strict sign disagreement places a root inside this actual step.
        // `opposite_sign` avoids multiplying possibly large event values.
        if (detail::opposite_sign(e_prev, e_next)) {
            double lo = 0.0;       // event(step(s_prev, lo)) has e_prev's sign, by construction
            double hi = step_dt;   // event(step(s_prev, hi)) == e_next, the opposite sign
            State s_mid = s_next;
            for (int iter = 0; iter < detail::kEventBisectionIterations; ++iter) {
                const double mid = 0.5 * (lo + hi);
                s_mid = detail::dispatch_step(opt.method, s_prev, p.mass, mid, force);
                s_mid.t = s_prev.t + mid;
                const double e_mid = event(s_mid);
                if (e_mid == 0.0 || detail::opposite_sign(e_prev, e_mid)) {
                    hi = mid;  // s_mid is already on (or past) the far side: root is in (lo, mid]
                } else {
                    lo = mid;  // s_mid is still on the near side: root is in (mid, hi]
                }
            }

            // A mathematically exact root at a sampled endpoint can evaluate
            // a few ULPs to either side of zero.  The following sign change
            // then makes bisection return a state only roundoff away from the
            // preceding endpoint.  Reuse that existing endpoint instead of
            // adding a numerically meaningless tiny final interval.
            if (detail::same_time_to_roundoff(s_mid.t, s_prev.t, step_dt)) {
                if (traj.states.back().t != s_prev.t) traj.states.push_back(s_prev);
                return traj;
            }
            if (detail::same_time_to_roundoff(s_mid.t, s_next.t, step_dt)) {
                traj.states.push_back(s_next);
                return traj;
            }
            traj.states.push_back(s_mid);
            return traj;
        }

        s = s_next;
        e_prev = e_next;
        if (e_prev != 0.0) event_has_left_initial_surface = true;

        const bool is_last_step = (i + 1 == n);
        if (is_last_step || (i + 1) % stride == 0) {
            traj.states.push_back(s);
        }
    }
    return traj;  // event never fired: identical behaviour to simulate()
}

namespace events {

/// Event function for simulate_until(): height above y = y_level. Its sign
/// change (positive to negative) marks the instant a falling particle
/// crosses that height — the canonical "stop at ground impact" event.
/// @param y_level ground height [m]
[[nodiscard]] inline auto ground_impact(double y_level = 0.0) {
    return [y_level](const State& s) { return s.r.y - y_level; };
}

/// Event function for simulate_until(): vertical velocity. Its sign change
/// (positive to negative) marks the apex of a vertical throw, where v.y
/// passes through zero.
[[nodiscard]] inline auto apex() {
    return [](const State& s) { return s.v.y; };
}

}  // namespace events

}  // namespace newton
