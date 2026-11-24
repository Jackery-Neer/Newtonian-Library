#pragma once
// newton/forces/spring.hpp — Hooke's law forces (Taylor §5.1, §5.3).
//
// Physics: Hooke's law F = -k x is the leading-order term of the force near
// any smooth potential minimum (Taylor's motivation for treating it as
// nearly universal, §5.1) — it's why "everything oscillates like a spring"
// for small displacements. This file has the isotropic form (restoring force
// always points straight at an anchor point) and the anisotropic form (an
// independent spring constant per axis, used to build Lissajous figures,
// §5.3).

#include <newton/core/vec3.hpp>

namespace newton::forces {

/// Isotropic Hooke's-law spring toward an anchor point (Taylor §5.1):
/// F = -k (|d| - natural_length) d_hat, with d = r - anchor.
///
/// Physics: with natural_length = 0 this is F = -k d, the simple harmonic
/// oscillator about `anchor` in any direction — displace it and it's pulled
/// straight back. A nonzero natural_length models a real spring that pushes
/// back once compressed past its rest length and pulls back once stretched
/// past it (e.g. a mass hanging from a spring of length natural_length at
/// equilibrium).
///
/// C++ note: d_hat = normalized(d) is zero-guarded (see vec3.hpp), so a
/// particle sitting exactly on the anchor with natural_length > 0 (an
/// unstable equilibrium) produces F = 0 rather than NaN.
struct Spring {
    double k{1.0};              // spring constant [N/m]
    Vec3   anchor{};            // anchor point [m]
    double natural_length{0.0}; // rest length [m]

    [[nodiscard]] Vec3 operator()(double /*t*/, const Vec3& r, const Vec3& /*v*/) const {
        const Vec3 d = r - anchor;
        return -k * (norm(d) - natural_length) * normalized(d);
    }
};

/// Anisotropic spring about the origin, independent constant per axis
/// (Taylor §5.3): F = (-kx x, -ky y, -kz z). Setting kx = ky (= kz) recovers
/// the isotropic oscillator; unequal constants decouple the axes into
/// independent 1D oscillators whose phase-space combination traces a
/// Lissajous figure — closed only when the frequency ratios are rational.
struct AnisotropicSpring {
    double kx{1.0}, ky{1.0}, kz{1.0};  // spring constants [N/m], one per axis

    [[nodiscard]] Vec3 operator()(double /*t*/, const Vec3& r, const Vec3& /*v*/) const {
        return {-kx * r.x, -ky * r.y, -kz * r.z};
    }
};

}  // namespace newton::forces
