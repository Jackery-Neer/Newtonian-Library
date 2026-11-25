#pragma once
// newton/forces/electromagnetic.hpp — magnetic and Lorentz forces
// (Taylor §2.5, §2.7).
//
// Physics: a charge q in electric field E and magnetic field B feels the
// Lorentz force F = q(E + v x B). The magnetic part alone (F = q v x B) is
// always perpendicular to v, so it does no work — |v| is exactly conserved
// and the motion is a helix along B (circular in the plane perpendicular to
// B, uniform drift along B).
//
// The magnetic force depends on v, so it is NOT conservative by Taylor §4.2
// condition (i) (no potential U(r) exists), yet it does no work: the WORK it
// books, W_mag = integral of F.v dt, is zero to roundoff for any integrator
// (~1e-17 J measured). That is not the same as T staying constant: a
// low-order integrator still changes |v| through its own truncation error
// (explicit and semi-implicit Euler grow T by ~170% over 100 s at dt = 0.01,
// Verlet by ~2e-5, RK4 by ~1e-10), so the kinetic energy is only as good as
// the integrator, while W_mag = 0 is exact.

#include <newton/core/vec3.hpp>

namespace newton::forces {

/// Magnetic force F = q v x B (Taylor §2.5). Does no work: because F is
/// always perpendicular to v (a x b is perpendicular to both a and b), the
/// power delivered F . v is identically zero, so |v| is exactly constant
/// along the true trajectory — a useful sanity check on any integrator (see
/// analysis/energy.hpp and analytic/magnetic.hpp for the cyclotron
/// formulas).
struct MagneticForce {
    double q{1.0};                    // charge [C]
    Vec3   B{0.0, 0.0, 1.0};          // magnetic field [T]

    [[nodiscard]] Vec3 operator()(double /*t*/, const Vec3& /*r*/, const Vec3& v) const {
        return q * cross(v, B);
    }
};

/// Full Lorentz force F = q (E + v x B) (Taylor §2.7).
///
/// Physics: the electric part q E does work (it's parallel to F, generally
/// not to v) and can accelerate or decelerate the charge; the magnetic part
/// only steers it. Crossed uniform E and B fields (E perpendicular to B)
/// produce the classic E x B drift: a cycloidal path whose guiding center
/// drifts at velocity (E x B)/|B|^2, independent of q and m.
struct LorentzForce {
    double q{1.0};  // charge [C]
    Vec3   E{};     // electric field [V/m] = [N/C]
    Vec3   B{};     // magnetic field [T]

    [[nodiscard]] Vec3 operator()(double /*t*/, const Vec3& /*r*/, const Vec3& v) const {
        return q * (E + cross(v, B));
    }
};

}  // namespace newton::forces
