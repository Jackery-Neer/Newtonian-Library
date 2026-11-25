#pragma once
// newton/forces/drag.hpp — fluid resistance (Taylor §2.1-2.4).
//
// Physics: Taylor models air resistance as a magnitude f(v) = b|v| + c|v|^2
// opposing the velocity, i.e. f = -b v - c |v| v. The linear term dominates
// for small, slow objects (Stokes drag, e.g. a dust mote); the quadratic term
// dominates for everyday projectiles (a baseball, a raindrop) because it's
// caused by turbulent momentum transfer to the air rather than viscosity.
// Chapter 2's whole point (see the research map) is that quadratic drag in
// 2D has no closed-form solution — this is the library's founding
// justification for numerical integration.
//
// C++ note: both forces depend on v but not on r or t, so those two
// parameters are unused and commented out per the Force concept convention
// (see gravity.hpp).

#include <newton/core/constants.hpp>
#include <newton/core/vec3.hpp>

namespace newton::forces {

/// Linear (viscous/Stokes) drag: F = -b v (Taylor §2.2-2.3).
///
/// Physics: proportional to speed, always opposing v. Dominant for small
/// objects at low Reynolds number (a smoke particle, an oil drop in
/// Millikan's experiment). Admits exact solutions — see
/// newton/analytic/drag.hpp.
///
/// Physics: NOT conservative. It depends on v, so it fails Taylor §4.2
/// condition (i) (F must depend on position only), and no potential U(r)
/// exists. It always does negative work, at the rate P = F.v = -b |v|^2, so
/// mechanical energy only decreases (use analysis::work_done / energy_balance,
/// passing it as F_nc, to keep the ledger).
struct LinearDrag {
    double b{0.0};  // linear drag coefficient [kg/s]

    [[nodiscard]] Vec3 operator()(double /*t*/, const Vec3& /*r*/, const Vec3& v) const {
        return -b * v;
    }
};

/// Quadratic (turbulent/form) drag: F = -c |v| v (Taylor §2.4).
///
/// Physics: NOT conservative: it depends on v, so it fails Taylor §4.2
/// condition (i), and no potential U(r) exists. It always does negative work,
/// at the rate P = F.v = -c |v|^3, so mechanical energy only decreases.
///
/// Physics: proportional to speed squared, opposing v. Dominant for
/// everyday-sized objects moving through air (baseballs, cars, skydivers).
///
/// C++ note: written as -c|v|v rather than the equivalent -c|v|^2 v_hat
/// (v_hat = v/|v|) precisely because it is singularity-free: at v = 0 this
/// evaluates to exactly the zero vector with no division, whereas the
/// v_hat form needs normalized(v), which has to special-case v = 0 to avoid
/// 0/0. Same physics, one fewer edge case (this is the algebraic identity
/// |v|^2 v_hat = |v|^2 * v/|v| = |v| v).
struct QuadraticDrag {
    double c{0.0};  // quadratic drag coefficient [kg/m]

    [[nodiscard]] Vec3 operator()(double /*t*/, const Vec3& /*r*/, const Vec3& v) const {
        return -c * norm(v) * v;
    }
};

/// Linear drag coefficient for a smooth sphere in air at STP (Taylor §2.1),
/// b = sphere_air_beta * D, from the diameter D [m]. Handy for
/// LinearDrag{linear_drag_sphere_air(D)}.
[[nodiscard]] inline double linear_drag_sphere_air(double D) {
    return constants::sphere_air_beta * D;
}

/// Quadratic drag coefficient for a smooth sphere in air at STP
/// (Taylor §2.1), c = sphere_air_gamma * D^2, from the diameter D [m].
/// Handy for QuadraticDrag{quadratic_drag_sphere_air(D)}.
[[nodiscard]] inline double quadratic_drag_sphere_air(double D) {
    return constants::sphere_air_gamma * D * D;
}

}  // namespace newton::forces
