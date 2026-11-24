#pragma once
// newton/integrators/euler.hpp — explicit (forward) Euler.
//
// Physics: every integrator in this library solves the first-order system
// r' = v, v' = a = F(t, r, v)/m that m r'' = F becomes once you carry v as a
// second unknown (Taylor §1.6, and see particle.hpp). Forward Euler is the
// most literal possible discretisation: evaluate the derivative at the OLD
// state and take one linear step along it.
//
//   a_n     = F(t_n, r_n, v_n) / m
//   r_{n+1} = r_n + dt * v_n        // <-- uses v_n, NOT v_{n+1}
//   v_{n+1} = v_n + dt * a_n
//
// That one substitution — v_n here, v_{n+1} in semi_implicit_euler.hpp — is
// the ONLY difference between this file and the next one, and it is the
// difference between an integrator that is useless for oscillatory motion
// and one that is a perfectly reasonable default.
//
// Physics: on the undamped oscillator x'' = -omega^2 x, this step multiplies
// the energy by exactly (1 + (omega*dt)^2) every single step (derivation in
// docs/research/numerical-integration.md §2.1) — so E_n = E_0 (1+k^2)^n with
// k = omega*dt. There is NO stable dt: shrinking dt shrinks the growth RATE
// but the growth is still unconditional and exponential in simulated time.
// Geometrically, in the phase plane (omega*x, v) the true motion is a circle
// and forward Euler always steps along the TANGENT to that circle — a
// tangent line never re-enters the circle, so every step lands strictly
// farther from the origin than the one before. This is why the method ships
// here at all: it is the sharpest, cheapest demonstration of why the other
// three integrators in this directory look the way they do. Use it as a
// pedagogical foil (see tests/test_integrators.cpp) or for purely
// dissipative relaxation where you only need a rough answer — never as a
// simulation default.
//
// C++ note: `Force&&` is a forwarding reference on a template parameter
// (not an rvalue-reference overload) — it binds to lambdas, functors, and
// std::function alike, and because the concrete callable type is known at
// compile time the compiler can fully inline `force(...)` into the step
// instead of going through indirect/virtual dispatch. See
// docs/research/cpp-library-design.md §5 "std::function vs template
// callables" for the measured cost of getting this wrong.

#include <newton/core/particle.hpp>

namespace newton::integrators {

/// One step of explicit (forward) Euler.
///
/// Physics: first-order accurate (global error O(dt)); energy on any
/// oscillator grows every step (see file header) — do not use for anything
/// that should conserve or bound energy. It is exact for zero-force,
/// constant-velocity motion only.
/// @param s      current state (r, v, t)
/// @param mass   particle mass [kg]
/// @param dt     step size [s]
/// @param force  callable: Vec3 operator()(double t, const Vec3& r, const Vec3& v) const, returns force [N]
/// @returns the state at s.t + dt
template <class Force>
[[nodiscard]] State euler_step(const State& s, double mass, double dt, Force&& force) {
    const Vec3 a = force(s.t, s.r, s.v) / mass;  // evaluated once, at the OLD (r, v)
    return State{s.r + dt * s.v, s.v + dt * a, s.t + dt};
}

}  // namespace newton::integrators
