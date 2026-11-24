#pragma once
// newton/analysis/energy.hpp — energy diagnostics for a simulated trajectory.
//
// Physics: mechanical energy E = T + U is conserved for any force derivable
// from a time-independent potential, F = -grad(U) (Taylor §4.1-4.4). No
// integrator step function in this library ever evaluates U — a numerical
// integrator only ever sees forces — so watching E(t) computed *after the
// fact* from a recorded Trajectory is the standard, integrator-agnostic way
// to grade how well a run conserved energy without touching the stepper's
// internals. This is the single most useful integrator diagnostic (research
// numerical-integration.md §4.3): a symplectic method (semi-implicit Euler,
// velocity Verlet) shows bounded oscillation around 0 drift with no trend;
// RK4 shows a small monotone drift; explicit Euler on an oscillator shows
// unbounded exponential growth. See energy_drift() below.
//
// C++ note: PotentialFn is a template callable (`double U(const Vec3& r)`),
// the same "template, not std::function" contract the Force concept uses
// (DESIGN.md). These functions typically run once per sample in a
// post-processing pass rather than inside a hot integration loop, but
// matching the rest of the library's calling convention costs nothing and
// keeps the two contracts easy to remember together.

#include <cmath>
#include <vector>

#include <newton/core/particle.hpp>
#include <newton/core/trajectory.hpp>
#include <newton/core/vec3.hpp>

namespace newton::analysis {

/// Kinetic energy T = 1/2 m |v|^2 (Taylor §4.1).
/// @param mass particle mass [kg]
/// @param v    velocity [m/s]
[[nodiscard]] inline double kinetic_energy(double mass, const Vec3& v) {
    return 0.5 * mass * norm2(v);
}

/// Total mechanical energy E = T + U(r) at one instant.
/// @param mass particle mass [kg]
/// @param s    state (r, v, t) to evaluate
/// @param U    potential energy function: double U(const Vec3& r) [J]
template <class PotentialFn>
[[nodiscard]] double total_energy(double mass, const State& s, PotentialFn&& U) {
    return kinetic_energy(mass, s.v) + U(s.r);
}

/// Relative energy drift over a trajectory: (E(t) - E(0)) / |E(0)|, one entry
/// per sample in traj.states, same length and order — so drift[i] lines up
/// with traj.states[i] for plotting drift(t).
/// @param traj trajectory to diagnose; traj.mass is used as the particle mass
/// @param U    potential energy function: double U(const Vec3& r) [J]
template <class PotentialFn>
[[nodiscard]] std::vector<double> energy_drift(const Trajectory& traj, PotentialFn&& U) {
    std::vector<double> drift;
    drift.reserve(traj.size());
    if (traj.empty()) return drift;

    const double E0 = total_energy(traj.mass, traj.states.front(), U);
    // Guard the (legitimate) E0 == 0 case — e.g. a particle released from
    // rest where U was zeroed at the release point — so a genuinely-zero
    // reference energy can't turn every later sample's drift into +-inf/NaN
    // by dividing by zero. Mirrors the zero-guard pattern in
    // Vec3::normalized() (core/vec3.hpp).
    const double scale = (std::fabs(E0) > 0.0) ? std::fabs(E0) : 1.0;

    for (const State& s : traj.states) {
        drift.push_back((total_energy(traj.mass, s, U) - E0) / scale);
    }
    return drift;
}

/// U(r) = m g r.y: potential energy for newton::forces::UniformGravity's
/// default convention (g points in -y — this library's 2D "y is up" rule,
/// see core/constants.hpp). Pass g as a positive magnitude (e.g.
/// constants::g_earth), not a signed vector component: F = -dU/dr then comes
/// out to (0, -m g, 0), matching that force.
/// @param mass particle mass [kg]
/// @param g    gravitational acceleration magnitude [m/s^2]
[[nodiscard]] inline auto uniform_gravity_U(double mass, double g) {
    return [mass, g](const Vec3& r) { return mass * g * r.y; };
}

/// U(r) = 1/2 k |r|^2: potential energy for an isotropic spring anchored at
/// the origin (newton::forces::Spring with anchor = {} and
/// natural_length = 0 — the reduced form F = -k r, Taylor §5.1).
/// @param k spring constant [N/m]
[[nodiscard]] inline auto spring_U(double k) {
    return [k](const Vec3& r) { return 0.5 * k * norm2(r); };
}

}  // namespace newton::analysis
