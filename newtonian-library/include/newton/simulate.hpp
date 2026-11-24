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

#include <cmath>       // std::llround
#include <cstddef>     // std::size_t
#include <stdexcept>   // std::logic_error

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
/// the step count is an integer, `n = llround((t_end - t0)/dt)`, rather
/// than a `while (t < t_end)` loop — a float-compared loop condition can
/// take a different number of iterations on different compilers/platforms
/// depending on the last-bit rounding of `t`, which would make trajectories
/// (and tests) non-reproducible.
template <class Force>
[[nodiscard]] Trajectory simulate(const Particle& p, Force&& force, const SimOptions& opt) {
    const long long n = std::llround((opt.t_end - opt.t0) / opt.dt);
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
        s = detail::dispatch_step(opt.method, s, p.mass, opt.dt, force);
        s.t = opt.t0 + static_cast<double>(i + 1) * opt.dt;  // see the CRITICAL note above

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
    const long long n = std::llround((opt.t_end - opt.t0) / opt.dt);
    const int stride = opt.record_every > 0 ? opt.record_every : 1;

    Trajectory traj;
    traj.mass = p.mass;
    traj.states.reserve(static_cast<std::size_t>(n / stride) + 2);

    State s = p.initial_state(opt.t0);
    traj.states.push_back(s);
    double e_prev = event(s);

    for (long long i = 0; i < n; ++i) {
        const State s_prev = s;  // the bisection below always re-steps from here
        State s_next = detail::dispatch_step(opt.method, s_prev, p.mass, opt.dt, force);
        s_next.t = opt.t0 + static_cast<double>(i + 1) * opt.dt;
        const double e_next = event(s_next);

        // Strict sign change: e_prev*e_next < 0 means the two endpoints
        // disagree on sign, so a root lies strictly between them. Using a
        // strict "<" (not "<=") matters at i == 0: a particle launched from
        // exactly the event surface (e.g. ground_impact() for a projectile
        // that starts at y = 0) has e_prev == 0 at the very first state, and
        // e_prev*e_next == 0 for that first step regardless of e_next's
        // sign — that must NOT be treated as a crossing, or the run would
        // terminate immediately at t0.
        if (e_prev * e_next < 0.0) {
            double lo = 0.0;       // event(step(s_prev, lo)) has e_prev's sign, by construction
            double hi = opt.dt;    // event(step(s_prev, hi)) == e_next, the opposite sign
            constexpr int kBisectionIterations = 40;  // halves the bracket 40x: dt/2^40 is far below any meaningful step size
            State s_mid = s_next;
            for (int iter = 0; iter < kBisectionIterations; ++iter) {
                const double mid = 0.5 * (lo + hi);
                s_mid = detail::dispatch_step(opt.method, s_prev, p.mass, mid, force);
                s_mid.t = s_prev.t + mid;
                const double e_mid = event(s_mid);
                if (e_prev * e_mid <= 0.0) {
                    hi = mid;  // s_mid is already on (or past) the far side: root is in (lo, mid]
                } else {
                    lo = mid;  // s_mid is still on the near side: root is in (mid, hi]
                }
            }
            traj.states.push_back(s_mid);
            return traj;
        }

        s = s_next;
        e_prev = e_next;

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
