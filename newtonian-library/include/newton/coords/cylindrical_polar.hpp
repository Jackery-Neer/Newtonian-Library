#pragma once
// newton/coords/cylindrical_polar.hpp — cylindrical polar coordinates
// (rho, phi, z) and the conversion into them from Cartesian.
//
// Physics: cylindrical polars are plane polars (Taylor §1.7) with the z axis
// left alone (Taylor problems 1.47-1.49). rho is the distance from the z
// AXIS, not from the origin (that is spherical r); phi is the angle around
// the z axis, measured from +x toward +y; z is unchanged. This phi is the
// same angle as the phi in spherical_polar.hpp, so converting between the
// two systems never touches it.
//
// C++ note: CylindricalPolar is a plain aggregate, deliberately NOT a Vec3.
// Vec3's API (+, dot, cross, norm) is Cartesian: adding two (rho, phi, z)
// triples component by component gives a meaningless answer, so this type
// offers no arithmetic at all. It is a label for a location. Do the math in
// Vec3 and convert at the edges. The conversion back to Cartesian lives in
// cartesian.hpp, the only coords header that sees every coordinate type.
// This header includes nothing from coords/, so the coords headers can never
// include each other in a circle.

#include <cmath>

#include <newton/core/vec3.hpp>

namespace newton::coords {

/// A point in cylindrical polar coordinates.
struct CylindricalPolar {
    double rho{0.0};  ///< distance from the z axis [m], >= 0
    double phi{0.0};  ///< angle around z, from +x toward +y [rad], in [-pi, pi]
    double z{0.0};    ///< height along the z axis [m]
};

/// Cartesian -> cylindrical: rho = sqrt(x^2 + y^2), phi = atan2(y, x), z = z.
/// Physics: phi comes from std::atan2(y, x), not std::atan(y / x). The ratio
/// y/x is the same for (1, 1) and (-1, -1), so atan alone cannot tell which
/// quadrant the point is in (Taylor problem 1.42). atan2 sees the sign of
/// each component and returns the correct angle in [-pi, pi]. On the z axis
/// (x = y = 0) phi is undefined: atan2 still returns a finite angle (0 or
/// ±pi, depending on the signs of the zeros), and every one of those converts
/// back to the same point because rho = 0 multiplies it away.
/// C++ note: std::hypot(x, y) is sqrt(x^2 + y^2) computed without
/// overflowing when x^2 alone would exceed the largest double.
/// @param v position [m]
[[nodiscard]] inline CylindricalPolar to_cylindrical(const Vec3& v) {
    return {std::hypot(v.x, v.y), std::atan2(v.y, v.x), v.z};
}

}  // namespace newton::coords
