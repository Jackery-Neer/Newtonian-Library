#pragma once
// newton/analytic/drag.hpp — closed-form drag solutions (Taylor §2.2-2.4).
//
// Physics: linear drag (F = -b v) is a *linear* ODE, so the full vector
// motion — including the coupling to gravity — has an exact solution; it is
// the single richest 2D validation oracle in the book (research map §8, row
// 6). Quadratic drag (F = -c|v|v) is only integrable in the two decoupled
// 1D cases below (pure horizontal coast, pure vertical fall from rest) —
// the general 2D quadratic-drag trajectory has NO closed form (that's the
// whole reason this library's numeric integrators exist; see
// forces/drag.hpp and docs/research/taylor-ch1-7-map.md §2).
//
// C++ note: every function here calls std::exp/log/tanh/cosh, none of which
// are constexpr in C++17, so nothing in this file is constexpr (matches
// core/vec3.hpp's norm()/normalized() precedent). LinearDragParams::tau()
// and v_ter() are pure arithmetic and ARE constexpr.

#include <cmath>

#include <newton/core/constants.hpp>
#include <newton/core/vec3.hpp>

namespace newton::analytic {

/// Parameters for the linear-drag closed forms below: m ẍ = -m g ŷ - b v,
/// i.e. gravity along -y at magnitude g (matching forces::UniformGravity's
/// default 2D convention, "y is up") plus linear drag -b v (Taylor §2.2).
struct LinearDragParams {
    double m{1.0};                 ///< mass [kg]
    double b{0.0};                 ///< linear drag coefficient [kg/s]
    double g{constants::g_earth};  ///< gravitational acceleration [m/s^2]

    /// Time constant τ = m/b [s] — the e-folding time for velocity
    /// relaxation toward the drift solution (Taylor §2.2).
    [[nodiscard]] constexpr double tau() const { return m / b; }

    /// Terminal (vertical, free-fall) speed v_ter = m g / b [m/s]
    /// (Taylor §2.2) — the speed at which drag exactly balances gravity.
    [[nodiscard]] constexpr double v_ter() const { return m * g / b; }
};

/// Velocity at time t under gravity + linear drag, from initial velocity v0
/// (Taylor §2.2-2.3). Vector form of both 1D results at once: writing
/// a_grav = (0, -g, 0) (matching LinearDragParams's -y gravity convention),
/// the ODE dv/dt = a_grav - v/τ integrates to
///   v(t) = v0 e^(-t/τ) + a_grav τ (1 - e^(-t/τ)).
/// The x-component (a_grav_x = 0) reduces to Taylor's vx(t) = vx0 e^(-t/τ);
/// the y-component reduces to vy(t) = -v_ter + (vy0+v_ter) e^(-t/τ), using
/// v_ter = g τ.
/// @param p  drag/gravity parameters
/// @param v0 initial velocity [m/s]
/// @param t  elapsed time [s]
[[nodiscard]] inline Vec3 linear_drag_velocity(const LinearDragParams& p, Vec3 v0, double t) {
    const double tau = p.tau();
    const double decay = std::exp(-t / tau);
    const Vec3 a_grav{0.0, -p.g, 0.0};
    return v0 * decay + a_grav * tau * (1.0 - decay);
}

/// Position at time t under gravity + linear drag, from initial position r0
/// and velocity v0 (Taylor §2.2-2.3). Obtained by integrating
/// linear_drag_velocity() in closed form:
///   r(t) = r0 + (v0 - a_grav τ) τ (1 - e^(-t/τ)) + a_grav τ t.
/// Components reduce to Taylor's x(t) = vx0 τ (1-e^(-t/τ)) and
/// y(t) = -v_ter t + (vy0+v_ter) τ (1-e^(-t/τ)), exactly as with
/// linear_drag_velocity().
/// @param p  drag/gravity parameters
/// @param r0 initial position [m]
/// @param v0 initial velocity [m/s]
/// @param t  elapsed time [s]
[[nodiscard]] inline Vec3 linear_drag_position(const LinearDragParams& p, Vec3 r0, Vec3 v0, double t) {
    const double tau = p.tau();
    const double decay = std::exp(-t / tau);
    const Vec3 a_grav{0.0, -p.g, 0.0};
    return r0 + (v0 - a_grav * tau) * tau * (1.0 - decay) + a_grav * tau * t;
}

/// The exact 2D trajectory y(x) for gravity + linear drag (Taylor §2.3):
///   y(x) = ((vy0+v_ter)/vx0) x + v_ter τ ln(1 - x/(vx0 τ))
/// eliminating time entirely — a rare case where a 2D drag problem has a
/// closed-form shape, and the best 2D integrator oracle in the book
/// (research map §8, row 6).
///
/// Physics: only valid for x < vx0·τ. As t -> infinity, x(t) -> vx0·τ (see
/// linear_drag_position) — the horizontal position has a genuine finite
/// asymptote, so y(x) has a real vertical asymptote there, not a numerical
/// artifact. Calling this with x >= vx0·τ evaluates ln of a non-positive
/// number: x == vx0·τ gives -infinity, x > vx0·τ gives NaN. Callers must
/// keep x strictly inside (-infinity, vx0·τ).
/// @param p   drag/gravity parameters
/// @param vx0 initial horizontal speed [m/s], assumed > 0
/// @param vy0 initial vertical speed [m/s]
/// @param x   horizontal position [m], must be < vx0 * p.tau()
[[nodiscard]] inline double linear_drag_trajectory_y(const LinearDragParams& p, double vx0, double vy0,
                                                      double x) {
    const double tau = p.tau();
    const double v_ter = p.v_ter();
    return (vy0 + v_ter) / vx0 * x + v_ter * tau * std::log(1.0 - x / (vx0 * tau));
}

/// Quadratic-drag terminal speed for vertical fall: v_ter = sqrt(m g / c)
/// (Taylor §2.4) — the speed at which quadratic drag c v² balances gravity.
/// @param m mass [kg]
/// @param c quadratic drag coefficient [kg/m]
/// @param g gravitational acceleration [m/s^2]
[[nodiscard]] inline double quadratic_drag_terminal_speed(double m, double c, double g) {
    return std::sqrt(m * g / c);
}

/// Fall speed under quadratic drag from rest: v(t) = v_ter tanh(g t / v_ter)
/// (Taylor §2.4) — the exact solution to m v̇ = m g - c v²; monotonically
/// approaches v_ter as t -> infinity (tanh saturates to 1).
/// @param v_ter terminal speed [m/s], from quadratic_drag_terminal_speed()
/// @param g     gravitational acceleration [m/s^2]
/// @param t     elapsed time since release [s]
[[nodiscard]] inline double quadratic_drag_fall_speed(double v_ter, double g, double t) {
    return v_ter * std::tanh(g * t / v_ter);
}

/// Fall distance under quadratic drag from rest:
/// y(t) = (v_ter²/g) ln cosh(g t / v_ter) (Taylor §2.4) — the time integral
/// of quadratic_drag_fall_speed().
/// @param v_ter terminal speed [m/s], from quadratic_drag_terminal_speed()
/// @param g     gravitational acceleration [m/s^2]
/// @param t     elapsed time since release [s]
[[nodiscard]] inline double quadratic_drag_fall_distance(double v_ter, double g, double t) {
    return (v_ter * v_ter / g) * std::log(std::cosh(g * t / v_ter));
}

/// Horizontal speed under quadratic drag alone (no gravity component along
/// the motion): v(t) = v0 / (1 + c v0 t / m) (Taylor §2.4) — power-law
/// decay, in contrast to linear drag's exponential decay. Never reaches
/// zero in finite time, but the associated distance x(t) = (m/c) ln(1 + c v0
/// t/m) grows without bound (logarithmically) rather than approaching a
/// finite asymptote the way linear drag's x_∞ = vx0 τ does.
/// @param m  mass [kg]
/// @param c  quadratic drag coefficient [kg/m]
/// @param v0 initial speed [m/s]
/// @param t  elapsed time [s]
[[nodiscard]] inline double quadratic_drag_horizontal_speed(double m, double c, double v0, double t) {
    return v0 / (1.0 + c * v0 * t / m);
}

}  // namespace newton::analytic
