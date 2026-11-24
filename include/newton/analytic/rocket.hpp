#pragma once
// newton/analytic/rocket.hpp — the rocket equation (Taylor §3.2).
//
// Physics: a rocket's equation of motion m dv/dt = -ṁ v_ex + F_ext follows
// from momentum conservation applied to the rocket-plus-exhaust system, NOT
// from F = ma with a fixed mass (see systems/rocket.hpp and
// docs/research/numerical-integration.md §7.8 for why "F = d(mv)/dt" is the
// wrong equation here). In free space (F_ext = 0) this integrates exactly
// to the Tsiolkovsky rocket equation; adding a constant external gravity
// g adds one more exact term.
//
// C++ note: both functions call std::log, which is not constexpr in C++17,
// so both are plain inline.

#include <cmath>

namespace newton::analytic {

/// Tsiolkovsky rocket equation: Δv = v_ex ln(m0/m) (Taylor §3.2) — the
/// speed gained in free space (no external force) burning from initial
/// mass m0 down to current mass m, at exhaust speed v_ex.
///
/// Physics: only mass ratio matters, not absolute mass or burn rate/time —
/// this is why staging helps (each stage restarts the ratio from a smaller
/// m0) and why chemical rockets (v_ex a few km/s) need large mass ratios to
/// reach orbital speeds (~8 km/s).
/// @param v_ex exhaust speed relative to the rocket [m/s]
/// @param m0   initial (wet) mass [kg]
/// @param m    current mass [kg], 0 < m <= m0
[[nodiscard]] inline double tsiolkovsky_delta_v(double v_ex, double m0, double m) {
    return v_ex * std::log(m0 / m);
}

/// Rocket velocity under constant thrust and constant external gravity g:
/// v(t) = v_ex ln(m0/(m0 - burn_rate·t)) - g t (Taylor §3.2), integrating
/// dv/dt = v_ex·burn_rate/m(t) - g with m(t) = m0 - burn_rate·t (constant
/// burn rate). Reduces to tsiolkovsky_delta_v() when g = 0.
///
/// Physics: only valid while burning, i.e. burn_rate·t < m0 (mass must stay
/// positive) — beyond burnout the rocket coasts under gravity alone and
/// this closed form no longer applies (see systems::simulate_rocket for the
/// full burn-then-coast numeric profile). At t == m0/burn_rate the
/// logarithm's argument is 0 and this returns -infinity; past that it's
/// negative and this returns NaN.
/// @param v_ex      exhaust speed relative to the rocket [m/s]
/// @param m0        initial (wet) mass [kg]
/// @param burn_rate mass flow rate |dm/dt| while burning [kg/s], > 0
/// @param g         gravitational acceleration [m/s^2]
/// @param t         elapsed burn time [s], must satisfy burn_rate*t < m0
[[nodiscard]] inline double rocket_velocity_gravity(double v_ex, double m0, double burn_rate, double g,
                                                     double t) {
    return v_ex * std::log(m0 / (m0 - burn_rate * t)) - g * t;
}

}  // namespace newton::analytic
