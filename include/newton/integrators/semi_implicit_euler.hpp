#pragma once
// newton/integrators/semi_implicit_euler.hpp — symplectic Euler (Euler–Cromer).
//
// Physics: same cost as forward Euler — one force evaluation, fully
// explicit, no linear solve anywhere — but velocity is advanced FIRST, and
// position then uses the NEW velocity instead of the old one:
//
//   a_n     = F(t_n, r_n, v_n) / m
//   v_{n+1} = v_n + dt * a_n
//   r_{n+1} = r_n + dt * v_{n+1}     // <-- v_{n+1}, contrast euler.hpp's v_n
//
// "Semi-implicit" is a slight misnomer: nothing is solved implicitly here,
// it is exactly as explicit as forward Euler. The name refers to r_{n+1}
// looking like it depends on itself (it's defined using v_{n+1}) even though
// v_{n+1} is already known by the time you compute it.
//
// Physics: that one reordering is enough to make the per-step map on phase
// space (omega*x, v) AREA-PRESERVING — det = 1 exactly, versus det = 1+k^2
// for forward Euler (k = omega*dt; full derivation in
// docs/research/numerical-integration.md §2.2). A map that exactly preserves
// phase-space volume is called SYMPLECTIC, and symplectic integrators have a
// remarkable property (backward error analysis): they don't approximately
// solve your Hamiltonian, they solve a nearby "shadow" Hamiltonian exactly.
// Energy of the shadow system never drifts, so the true energy is trapped
// within a band of relative half-width ~ omega*dt/2 around a constant,
// forever — it oscillates but never trends. This argument needs dt to stay
// FIXED for the whole run; change dt mid-simulation and you jump shadow
// Hamiltonians, and drift comes back.
//
// Physics: this is still only first-order accurate PER STEP (same local
// truncation error as forward Euler) — symplecticity buys bounded long-run
// energy, not pointwise accuracy. For a central force (F parallel to r) the
// angular momentum r x v is conserved EXACTLY, to machine roundoff, by this
// scheme (see the research doc for the one-line algebraic proof), which
// makes it a surprisingly good cheap orbit integrator despite being "only"
// order 1. This is the recommended default when you need something cheap,
// safe, and qualitatively right over a long run — reach for velocity_verlet
// or rk4 when you need pointwise accuracy instead.
//
// C++ note: see euler.hpp for why `Force&&` beats std::function here.

#include <newton/core/particle.hpp>

namespace newton::integrators {

/// One step of semi-implicit (symplectic) Euler, a.k.a. Euler–Cromer.
///
/// Physics: first-order accurate, same as forward Euler — but energy stays
/// bounded (oscillates, never drifts) instead of growing without limit. The
/// default "cheap and safe" choice for position-dependent forces F(r).
/// @param s      current state (r, v, t)
/// @param mass   particle mass [kg]
/// @param dt     step size [s]
/// @param force  callable: Vec3 operator()(double t, const Vec3& r, const Vec3& v) const, returns force [N]
/// @returns the state at s.t + dt
template <class Force>
[[nodiscard]] State semi_implicit_euler_step(const State& s, double mass, double dt, Force&& force) {
    const Vec3 a = force(s.t, s.r, s.v) / mass;  // evaluated at the OLD (r, v), same as forward Euler
    const Vec3 v_new = s.v + dt * a;              // velocity first...
    const Vec3 r_new = s.r + dt * v_new;          // ...then position, using the NEW velocity
    return State{r_new, v_new, s.t + dt};
}

}  // namespace newton::integrators
