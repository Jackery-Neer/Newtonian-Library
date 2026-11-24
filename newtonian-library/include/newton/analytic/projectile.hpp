#pragma once
// newton/analytic/projectile.hpp — vacuum projectile motion (Taylor §1.6).
//
// Physics: with gravity the only force, x and y decouple into two
// independent 1D constant-acceleration problems (a_x = 0, a_y = -g), giving
// the textbook parabola. This is the library's simplest validation oracle:
// a degree-2 polynomial in t that RK4 and velocity Verlet must reproduce to
// machine precision (see docs/research/numerical-integration.md §4.1(a) and
// the analytic table in docs/research/taylor-ch1-7-map.md §8, row 2).
//
// C++ note: std::sin/std::cos are not constexpr in C++17, so every function
// here is plain `inline` rather than `constexpr` — the same trade-off
// core/vec3.hpp makes for norm() (see its C++ note).

#include <cmath>

#include <newton/core/constants.hpp>
#include <newton/core/vec3.hpp>

namespace newton::analytic {

/// Horizontal range on level ground: R = v0² sin(2θ) / g (Taylor §1.6).
/// Maximized at θ = 45° (pi/4 rad), where R = v0²/g — a handy sanity check.
/// @param v0    launch speed [m/s]
/// @param theta launch angle above horizontal [rad]
/// @param g     gravitational acceleration [m/s^2], defaults to Earth's
[[nodiscard]] inline double range(double v0, double theta, double g = constants::g_earth) {
    return v0 * v0 * std::sin(2.0 * theta) / g;
}

/// Apex (maximum) height above the launch point: h = v0² sin²θ / (2g)
/// (Taylor §1.6).
/// @param v0    launch speed [m/s]
/// @param theta launch angle above horizontal [rad]
/// @param g     gravitational acceleration [m/s^2], defaults to Earth's
[[nodiscard]] inline double apex_height(double v0, double theta, double g = constants::g_earth) {
    const double s = std::sin(theta);
    return v0 * v0 * s * s / (2.0 * g);
}

/// Total time of flight back to launch height: T = 2 v0 sinθ / g
/// (Taylor §1.6).
/// @param v0    launch speed [m/s]
/// @param theta launch angle above horizontal [rad]
/// @param g     gravitational acceleration [m/s^2], defaults to Earth's
[[nodiscard]] inline double time_of_flight(double v0, double theta, double g = constants::g_earth) {
    return 2.0 * v0 * std::sin(theta) / g;
}

/// Position at time t, launched from the origin: x = v0 cosθ · t,
/// y = v0 sinθ · t - ½ g t² (Taylor §1.6). z is always 0 — this is a planar
/// trajectory embedded in R^3, matching the rest of the library's "2D
/// problems just use z = 0" convention (see core/vec3.hpp).
/// @param t     elapsed time since launch [s]
/// @param v0    launch speed [m/s]
/// @param theta launch angle above horizontal [rad]
/// @param g     gravitational acceleration [m/s^2], defaults to Earth's
[[nodiscard]] inline Vec3 position(double t, double v0, double theta, double g = constants::g_earth) {
    return {v0 * std::cos(theta) * t, v0 * std::sin(theta) * t - 0.5 * g * t * t, 0.0};
}

}  // namespace newton::analytic
