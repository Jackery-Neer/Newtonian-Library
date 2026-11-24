#pragma once
// newton/analysis/momentum.hpp — momentum and angular-momentum diagnostics.
//
// Physics: for an isolated system, internal forces cancel in equal-and-
// opposite pairs (Newton's third law), so total linear momentum
// P = sum m_i v_i is conserved (Taylor §3.1); angular momentum about a fixed
// point, L = sum r_i x p_i, is conserved for a single free particle about any
// point, and for a system whose net external torque vanishes — in
// particular any central force about its own center (Taylor §3.4). These are
// standard "is my simulation obeying physics" checks, independent of and
// complementary to the energy checks in energy.hpp.
//
// C++ note: center_of_mass()/total_momentum() take std::vector<Particle>
// rather than std::vector<State>, because computing R = (1/M) sum m_i r_i
// needs a mass alongside every position, and Particle is the type that
// carries mass with (r, v). This v0.1 library has no plural counterpart to
// Trajectory (one shared mass, many states) for multiple bodies — that's the
// v0.2 System type (see notes/04-extending.md); here the caller assembles
// the vector<Particle> snapshot itself, one entry per body, at whatever
// instant they want to check.

#include <vector>

#include <newton/core/particle.hpp>
#include <newton/core/vec3.hpp>

namespace newton::analysis {

/// Angular momentum about a reference point: ell = (r - about) x (m v)
/// (Taylor §3.4). Conserved for a central force (F = f(r) r-hat) taken about
/// that force's center, or trivially for any free (F = 0) particle about any
/// fixed point.
/// @param mass  particle mass [kg]
/// @param s     state (r, v, t) to evaluate
/// @param about reference point [m]; defaults to the origin
[[nodiscard]] inline Vec3 angular_momentum(double mass, const State& s, const Vec3& about = {}) {
    return cross(s.r - about, mass * s.v);
}

/// Center of mass R = (1/M) sum m_i r_i, M = sum m_i (Taylor §3.3). Returns
/// the zero vector for an empty list or non-positive total mass — the
/// defining formula divides by M, so M <= 0 is guarded the same way
/// Vec3::normalized() guards |v| == 0.
[[nodiscard]] inline Vec3 center_of_mass(const std::vector<Particle>& particles) {
    Vec3 weighted_sum{};
    double total_mass = 0.0;
    for (const Particle& p : particles) {
        weighted_sum += p.mass * p.r;
        total_mass += p.mass;
    }
    if (total_mass <= 0.0) return {};
    return weighted_sum / total_mass;
}

/// Total (linear) momentum P = sum m_i v_i (Taylor §3.1). For an isolated
/// system this is exactly conserved: internal forces exert equal and
/// opposite impulses in pairs (Newton III), so only a net external force can
/// change P — the standard multi-body structural check (research
/// numerical-integration.md §4.3).
[[nodiscard]] inline Vec3 total_momentum(const std::vector<Particle>& particles) {
    Vec3 sum{};
    for (const Particle& p : particles) sum += p.mass * p.v;
    return sum;
}

}  // namespace newton::analysis
