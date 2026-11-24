#pragma once
// newton/systems/rocket.hpp — vertical variable-mass rocket flight (Taylor §3.2).
//
// Physics: WHY THIS CANNOT REUSE newton::Particle / newton::simulate(). Every
// integrator in this library advances a *constant-mass* system,
// dr/dt = v, dv/dt = F(t,r,v)/m, with m fixed for the whole run — that is
// Newton's second law applied to one unchanging lump of matter. A burning
// rocket is not that: propellant is continuously ejected, so the mass in
// "F = ma" is itself a function of time. The tempting shortcut — integrate
// d(mv)/dt = F_ext, on the theory that F = dp/dt — is simply wrong here,
// because
//
//     d(mv)/dt = m dv/dt + v dm/dt ,
//
// and dm/dt != 0 for a burning rocket, so that extra v*dm/dt term smuggles
// in a fictitious force with no physical source (research
// numerical-integration.md §7.8: "F = ma is not d(mv)/dt"). The correct
// equation of motion treats the ejected propellant properly: momentum
// conservation of (rocket + expelled gas) gives, in the rocket's own frame,
//
//     m(t) dv/dt = F_ext + v_ex * |dm/dt|      (Taylor eq. 3.14/3.16)
//
// i.e. "F = ma" *does* hold, but only once the thrust term v_ex * burn_rate
// (the reaction to ejecting mass at rate burn_rate [kg/s] with exhaust speed
// v_ex relative to the rocket) is included as a genuine force on the
// right-hand side, divided by the *current* mass m(t) — not the constant
// mass Particle assumes. For vertical flight in uniform gravity (F_ext =
// -m(t) g) this becomes the scalar ODE this file integrates:
//
//     dv/dt = v_ex * burn_rate / m(t)  -  g       (while m(t) > m_dry)
//     dv/dt =                         -  g        (after burnout: thrust = 0)
//     dy/dt = v
//
// with m(t) = m0 - burn_rate * t. Since Particle/State/Trajectory have no
// slot for a time-varying mass — and don't need one for anything else in
// this v0.1 library — the rocket gets its own small scalar (y, v) state and
// its own RK4 loop instead of going through simulate().
//
// C++ note: kept fully self-contained (own RocketSample, own RK4 stages)
// rather than forcing the Vec3-based State/Trajectory types to carry a
// 1-DOF problem with y and z permanently pinned at zero — a plain
// {double t, m, v, y;} is simpler to read, plot, and test.

#include <algorithm>
#include <cmath>
#include <vector>

#include <newton/core/constants.hpp>

namespace newton::systems {

/// Inputs to a vertical single-stage rocket burn (Taylor §3.2).
struct RocketParams {
    double m0;                      ///< initial (wet) mass [kg]
    double m_dry;                   ///< mass at burnout [kg]; must be < m0
    double v_ex;                    ///< exhaust speed relative to the rocket [m/s]
    double burn_rate;               ///< |dm/dt| while burning [kg/s]; must be > 0
    double g{constants::g_earth};   ///< uniform gravitational acceleration [m/s^2]
};

/// One sampled instant of the simulated flight.
struct RocketSample {
    double t;  ///< time since launch [s]
    double m;  ///< instantaneous mass [kg]
    double v;  ///< vertical velocity, positive = up [m/s]
    double y;  ///< altitude above the launch point [m]
};

namespace detail {

/// Mass at time t, computed fresh from the closed form m0 - burn_rate*t on
/// every call — NEVER by accumulating m -= burn_rate*dt step by step, which
/// would drift from the true m(t) by an amount that compounds with the
/// number of steps taken (exactly the trap DESIGN.md and research §7.8 call
/// out: "m(t) exactly, never accumulated"). `burning` selects the formula;
/// see rocket_accel()'s C++ note for why the CALLER decides this rather than
/// this function re-deriving it from t vs. t_burnout.
[[nodiscard]] inline double rocket_mass(const RocketParams& p, bool burning, double t) {
    return burning ? p.m0 - p.burn_rate * t : p.m_dry;
}

/// dv/dt at time t: thrust/mass minus gravity while burning, else free fall
/// under gravity alone once the propellant is spent.
///
/// C++ note: `burning` is passed in by the caller — which phase of
/// simulate_rocket()'s burn/coast loop is currently stepping — rather than
/// re-derived here by comparing t to a t_burnout threshold. A threshold
/// comparison (`t >= t_burnout`) is ambiguous at exactly the one instant
/// that matters most: t == t_burnout is simultaneously the LAST t+h stage
/// (k4) of the burn phase's final step and the FIRST t stage (k1) of the
/// coast phase's first step. Comparing against t_burnout there silently
/// picks the coast formula (thrust = 0) for the burn phase's own k4 stage —
/// but the physically correct value approaching burnout from below is the
/// one-sided limit WITH thrust still firing (mass has only just reached
/// m_dry, not before). That swap corrupts the burn phase's last RK4 step by
/// an O(dt) amount — the thrust acceleration right before cutoff is
/// v_ex*burn_rate/m_dry, tens to hundreds of m/s^2, not a small
/// correction — right at the burnout instant the file header promises is an
/// exact sample. Threading `burning` through from simulate_rocket()'s own
/// phase bookkeeping sidesteps the ambiguity instead of trying to break a
/// tie at an instant that is genuinely shared by both phases.
[[nodiscard]] inline double rocket_accel(const RocketParams& p, bool burning, double t) {
    if (!burning) return -p.g;
    return p.v_ex * p.burn_rate / rocket_mass(p, true, t) - p.g;
}

}  // namespace detail

/// Integrates vertical flight from t=0 to t=t_end with classic RK4 on the
/// coupled scalar system dy/dt = v, dv/dt = a(t) (see the file header for
/// why a(t) switches off thrust at burnout rather than folding variable mass
/// into a constant-mass Particle). a(t) here does not depend on y or v — the
/// forces in this model, thrust and gravity, are both pure functions of time
/// — so all four RK4 stages are still written out explicitly (matching the
/// standard form in newton/integrators/rk4.hpp) even though k2v and k3v end
/// up evaluating the same midpoint.
///
/// The run is split into a burn phase, [0, min(t_burnout, t_end)], and a
/// coast phase, [that, t_end], each re-indexed on its own t = phase_start +
/// i*h grid (never t += h) so floating-point error never accumulates across
/// steps — the same rule newton::simulate() follows (research §7.2) — and,
/// by construction, the *last* step of the burn phase lands exactly on
/// t_burnout. That guarantees the burnout instant is itself a sample: it is
/// where thrust cuts off and v(t)'s slope has a corner, so a fixed dt-grid
/// that merely straddled it would blur the kink and bias any burnout-time
/// reading taken from the output.
///
/// @param params  rocket parameters (masses, exhaust speed, burn rate, g)
/// @param dt      nominal step size [s]; each phase subdivides its own span
///                into round(span/dt) equal steps, so both the burnout
///                instant and t_end land exactly on a sample
/// @param t_end   run duration [s], measured from t=0
/// @returns one RocketSample per step, including t=0, the burnout instant
///          (if it falls within [0, t_end]), and t_end
[[nodiscard]] inline std::vector<RocketSample>
simulate_rocket(const RocketParams& params, double dt, double t_end) {
    const double t_burnout = (params.m0 - params.m_dry) / params.burn_rate;

    std::vector<RocketSample> samples;
    samples.push_back({0.0, params.m0, 0.0, 0.0});

    double y = 0.0;
    double v = 0.0;

    double phase_start = 0.0;
    // phase 0 is the burn, [0, min(t_burnout, t_end)]; phase 1 is the coast,
    // [that, t_end]. `burning` is fixed per phase (see rocket_accel()'s C++
    // note above for why this — not a t-vs-t_burnout comparison — is what
    // decides the formula at the shared boundary instant).
    const double phase_ends[2] = {std::min(t_burnout, t_end), t_end};
    const bool phase_burning[2] = {true, false};
    for (int phase = 0; phase < 2; ++phase) {
        const double phase_end = phase_ends[phase];
        const bool burning = phase_burning[phase];
        const double span = phase_end - phase_start;
        if (span <= 0.0) continue;  // burnout at/before t=0, or already past t_end

        const long long n = std::max<long long>(1, std::llround(span / dt));
        const double h = span / static_cast<double>(n);

        for (long long i = 1; i <= n; ++i) {
            const double t = phase_start + static_cast<double>(i - 1) * h;

            // Classic RK4 stages for the pair (y, v).
            const double k1v = detail::rocket_accel(params, burning, t);
            const double k1y = v;

            const double k2v = detail::rocket_accel(params, burning, t + 0.5 * h);
            const double k2y = v + 0.5 * h * k1v;

            const double k3v = detail::rocket_accel(params, burning, t + 0.5 * h);
            const double k3y = v + 0.5 * h * k2v;

            const double k4v = detail::rocket_accel(params, burning, t + h);
            const double k4y = v + h * k3v;

            v = v + (h / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
            y = y + (h / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);

            const double t_next = phase_start + static_cast<double>(i) * h;
            samples.push_back({t_next, detail::rocket_mass(params, burning, t_next), v, y});
        }

        phase_start = phase_end;
    }

    return samples;
}

}  // namespace newton::systems
