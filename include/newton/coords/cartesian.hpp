#pragma once
// newton/coords/cartesian.hpp — conversions from cylindrical and spherical
// polar coordinates back to Cartesian (Vec3), plus the direct conversions
// between cylindrical and spherical.
//
// Physics: the library always simulates in Cartesian (Taylor §1.6). Forces,
// integrators and every Vec3 operation assume x, y, z components, because
// the x, y, z unit vectors are the same everywhere in space. Cylindrical and
// spherical coordinates are for setting up a problem and reading off
// results: convert at the edges, compute in Vec3.
//
// C++ note: this is the only coords header that includes the others. The
// two coordinate headers include only vec3.hpp, so the include graph has no
// cycle. A cycle (A includes B, B includes A) leaves one struct undefined
// while the other is being compiled, even with #pragma once.
//
// Every function is named after the type it RETURNS, and overloads are told
// apart by their argument type: to_cartesian(CylindricalPolar) and
// to_cartesian(SphericalPolar) can share a name because their parameters
// differ. Naming by the input instead (from_spherical) breaks as soon as two
// functions take the same input and return different types, because C++
// cannot overload on return type alone.
//
// One consequence: write the type out when you call these with literal
// numbers, as in to_cartesian(SphericalPolar{2.0, pi / 2, 0.0}). A bare
// to_cartesian({2.0, pi / 2, 0.0}) is ambiguous and will not compile. That
// is a feature: the same three numbers name different points in the two
// systems, so the compiler makes you say which one you mean.

#include <cmath>

#include <newton/coords/cylindrical_polar.hpp>
#include <newton/coords/spherical_polar.hpp>
#include <newton/core/vec3.hpp>

namespace newton::coords {

/// Cylindrical -> Cartesian: (rho cos(phi), rho sin(phi), z)
/// (Taylor problems 1.47-1.49).
/// @param c point in cylindrical polars: rho [m], phi [rad], z [m]
[[nodiscard]] inline Vec3 to_cartesian(const CylindricalPolar& c) {
    return {c.rho * std::cos(c.phi), c.rho * std::sin(c.phi), c.z};
}

/// Spherical -> Cartesian:
/// (r sin(theta) cos(phi), r sin(theta) sin(phi), r cos(theta)).
/// @param s point in spherical polars: r [m], theta [rad], phi [rad]
[[nodiscard]] inline Vec3 to_cartesian(const SphericalPolar& s) {
    const double rho = s.r * std::sin(s.theta);  // distance from the z axis
    return {rho * std::cos(s.phi), rho * std::sin(s.phi), s.r * std::cos(s.theta)};
}

/// Spherical -> cylindrical, without going through Cartesian:
/// rho = r sin(theta), phi unchanged, z = r cos(theta).
/// Physics: the two systems share phi, so only the half-plane at fixed phi
/// changes coordinates, and inside it (r, theta) are just plane polars for
/// the point (rho, z), with theta measured from the z axis.
/// @param s point in spherical polars: r [m], theta [rad], phi [rad]
[[nodiscard]] inline CylindricalPolar to_cylindrical(const SphericalPolar& s) {
    return {s.r * std::sin(s.theta), s.phi, s.r * std::cos(s.theta)};
}

/// Cylindrical -> spherical, without going through Cartesian:
/// r = sqrt(rho^2 + z^2), theta = atan2(rho, z), phi unchanged.
/// Physics: atan2 rather than acos(z / r), for the reasons given at
/// to_spherical(const Vec3&) in spherical_polar.hpp.
/// @param c point in cylindrical polars: rho [m], phi [rad], z [m]
[[nodiscard]] inline SphericalPolar to_spherical(const CylindricalPolar& c) {
    return {std::hypot(c.rho, c.z), std::atan2(c.rho, c.z), c.phi};
}

}  // namespace newton::coords
