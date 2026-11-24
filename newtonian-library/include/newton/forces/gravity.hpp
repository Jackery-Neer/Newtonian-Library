#pragma once
// newton/forces/gravity.hpp — uniform and point-mass gravity.
//
// Physics: two very different regimes both called "gravity" in intro
// mechanics. Near a planet's surface the field is uniform (Taylor §2.1 sets
// this up as the baseline projectile problem); at planetary/orbital scale it
// is the inverse-square law of Taylor §4.8 (central forces). Both are
// position-only forces — velocity v is unused by either.
//
// C++ note: both structs are aggregates (see vec3.hpp) with all-public
// members and default initializers, so `UniformGravity{}` and
// `UniformGravity{2.0}` both work. `t` and `v` are unused parameters in every
// operator() here; they're commented out (`/*t*/`) rather than removed so the
// Force concept's full signature stays visible at the call site while
// `-Wunused-parameter` (from -Wextra) stays quiet.

#include <cmath>

#include <newton/core/constants.hpp>
#include <newton/core/vec3.hpp>

namespace newton::forces {

/// Uniform gravitational field: F = m g, constant everywhere.
///
/// Physics: valid near a planet's surface where the field doesn't vary
/// appreciably over the trajectory's extent (Taylor's projectile setup,
/// §2.1). Default g points in -y with magnitude g_earth — this library's 2D
/// convention is "y is up" (see constants.hpp).
struct UniformGravity {
    double mass{1.0};                                  // [kg]
    Vec3   g{0.0, -constants::g_earth, 0.0};           // [m/s^2]

    [[nodiscard]] Vec3 operator()(double /*t*/, const Vec3& /*r*/, const Vec3& /*v*/) const {
        return mass * g;
    }
};

/// Newtonian point-mass gravity: F = -G M m (r - center) / |r - center|^3
/// (Taylor §4.8, central forces). Attractive, along the line joining the
/// orbiting particle to the central body.
///
/// Physics: GM is stored as one product (not G and M separately) because
/// that's the combination every orbital-mechanics formula actually uses —
/// it's also how you'd plug in a tabulated "standard gravitational
/// parameter" (e.g. Earth's GM) without needing the body's mass in kg.
///
/// Optional Plummer softening `eps` replaces |d|^2 with (|d|^2 + eps^2) in
/// the denominator (raised to the 3/2 power), which caps the force as
/// d -> 0 instead of letting it diverge. This is a numerical convenience
/// (keeps close-encounter integrations finite) that changes the physics at
/// short range — with eps > 0 this is no longer exact inverse-square gravity
/// inside r ~ eps. Leave eps = 0 for the textbook law.
struct PointGravity {
    double GM{1.0};        // G * (central mass) [m^3/s^2]
    double mass{1.0};      // orbiting particle's mass [kg]
    Vec3   center{};       // position of the attracting body [m]
    double softening{0.0}; // Plummer eps [m]

    [[nodiscard]] Vec3 operator()(double /*t*/, const Vec3& r, const Vec3& /*v*/) const {
        const Vec3 d = r - center;                    // vector from center to particle
        const double denom = std::pow(norm2(d) + softening * softening, 1.5);
        if (denom <= 0.0) return {};                  // guard: particle exactly at an unsoftened center
        return (-GM * mass / denom) * d;
    }
};

}  // namespace newton::forces
