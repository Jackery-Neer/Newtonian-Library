#pragma once
// newton/core/particle.hpp — the state of a point mass.
//
// Physics: Newton's second law m r̈ = F(t, r, ṙ) is a second-order ODE. Every
// integrator in this library first recasts it as a first-order system in the
// pair (r, v):  dr/dt = v,  dv/dt = F(t, r, v)/m  (Taylor §1.6). `State` is
// exactly that pair, stamped with the time it belongs to.
//
// C++ note: both types are aggregates (see vec3.hpp for why). `mass` defaults
// to 1.0, not 0.0 — a zero default would make a = F/m infinite the moment you
// forget to set it, and 1.0 makes "natural units" problems read cleanly.

#include <newton/core/vec3.hpp>

namespace newton {

/// A snapshot of a particle's motion: position, velocity, and the time at
/// which they hold. Integrator steps map State -> State.
struct State {
    Vec3   r{};      ///< position [m]
    Vec3   v{};      ///< velocity [m/s]
    double t{0.0};   ///< time [s]
};

/// A point mass with initial conditions. The simulate() driver takes one of
/// these and produces a Trajectory.
struct Particle {
    double mass{1.0};  ///< mass [kg], must be > 0
    Vec3   r{};        ///< initial position [m]
    Vec3   v{};        ///< initial velocity [m/s]

    /// The particle's initial conditions as a State at t = t0.
    [[nodiscard]] constexpr State initial_state(double t0 = 0.0) const {
        return {r, v, t0};
    }
};

/// Linear momentum p = m v (Taylor §1.5).
[[nodiscard]] constexpr Vec3 momentum(const Particle& p) { return p.mass * p.v; }

}  // namespace newton
