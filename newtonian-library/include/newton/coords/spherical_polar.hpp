#pragma once
// newton/coords/spherical_polar.hpp — spherical polar coordinates
// (r, theta, phi) and the conversion into them from Cartesian.
//
// Physics: this is the PHYSICS convention, the one Taylor uses and the one
// his spherical gradient, divergence and curl formulas assume (§4.8):
//   r      distance from the origin                     [0, inf)
//   theta  polar angle, measured DOWN from the +z axis  [0, pi]
//   phi    azimuthal angle around z, from +x toward +y  [-pi, pi]
//
//   x = r sin(theta) cos(phi),   y = r sin(theta) sin(phi),   z = r cos(theta)
//
// Many math courses swap the names theta and phi. A formula copied from a
// math text can silently give wrong answers here, so check which angle is
// measured from the z axis before using one. phi here is the same angle as
// the phi in cylindrical_polar.hpp.
//
// C++ note: SphericalPolar is a plain aggregate, deliberately NOT a Vec3;
// see cylindrical_polar.hpp for why. This header includes nothing from
// coords/. The conversion back to Cartesian lives in cartesian.hpp.

#include <cmath>

#include <newton/core/vec3.hpp>

namespace newton::coords {

/// A point in spherical polar coordinates (physics convention).
struct SphericalPolar {
    double r{0.0};      ///< distance from the origin [m], >= 0
    double theta{0.0};  ///< polar angle down from +z [rad], in [0, pi]
    double phi{0.0};    ///< azimuthal angle around z, from +x toward +y [rad], in [-pi, pi]
};

/// Cartesian -> spherical:
///   r = sqrt(rho^2 + z^2),  theta = atan2(rho, z),  phi = atan2(y, x),
/// where rho = sqrt(x^2 + y^2) is the distance from the z axis.
/// Physics: textbooks write theta = acos(z / r). That is the same angle,
/// since tan(theta) = rho / z, but acos has two problems atan2 avoids.
/// (1) It divides by r, so the origin gives 0/0 = NaN. (2) It loses
/// precision near the poles: for a point 1e-10 rad off the +z axis, z/r
/// rounds to exactly 1.0 and acos returns 0, snapping the point onto the
/// axis. atan2(rho, z) has no division and keeps full precision at every
/// angle. At the origin and on the z axis the angles are undefined; atan2
/// still returns finite values there, and the point converts back correctly.
/// @param v position [m]
[[nodiscard]] inline SphericalPolar to_spherical(const Vec3& v) {
    const double rho = std::hypot(v.x, v.y);  // distance from the z axis
    return {std::hypot(rho, v.z), std::atan2(rho, v.z), std::atan2(v.y, v.x)};
}

}  // namespace newton::coords
