#pragma once
// newton/integrators/rk4.hpp — classic 4th-order Runge-Kutta.
//
// Physics: RK4 is the generic 4-stage, 4th-order Runge-Kutta method applied
// to the coupled first-order system y' = f(t, y) with y = (r, v),
// f = (v, a) and a(t, r, v) = F(t, r, v)/m. Each stage evaluates the
// derivative at a trial point built from the previous stages, then a
// weighted average of all four gives the step:
//
//   k1r = v_n
//   k1v = a( t_n,          r_n,               v_n              )
//
//   k2r = v_n + 0.5*dt*k1v
//   k2v = a( t_n + 0.5*dt, r_n + 0.5*dt*k1r,   v_n + 0.5*dt*k1v )
//
//   k3r = v_n + 0.5*dt*k2v
//   k3v = a( t_n + 0.5*dt, r_n + 0.5*dt*k2r,   v_n + 0.5*dt*k2v )
//
//   k4r = v_n + dt*k3v
//   k4v = a( t_n + dt,     r_n + dt*k3r,       v_n + dt*k3v     )
//
//   r_{n+1} = r_n + (dt/6)*(k1r + 2*k2r + 2*k3r + k4r)
//   v_{n+1} = v_n + (dt/6)*(k1v + 2*k2v + 2*k3v + k4v)
//
// STAGE ORDER MATTERS: at each stage you must compute k_ir before k_iv,
// because k_iv's position argument is r_n + (fraction)*dt*k_ir — the very
// k_ir you just computed. Concretely: k2v needs k1r, k3v needs k2r, k4v
// needs k3r. If you get this wrong — for instance, reuse k1r where k2r
// belongs, because the mistake still type-checks and still runs — the
// method silently degrades from 4th order to 2nd order. It does not crash,
// it does not warn, it just quietly gets less accurate; nothing catches this
// except a convergence-order test (tests/test_integrators.cpp fits the
// observed order by halving dt and asserts it lands near 4 — that test
// exists specifically to catch this class of bug).
//
// Physics: RK4 is 4th-order accurate (global error O(dt^4)) but pays for it
// with FOUR force evaluations per step (no "first same as last" reuse in
// the classic tableau, unlike Verlet's one). It is neither symplectic nor
// time-reversible: on an undamped oscillator its amplification factor is
// slightly less than 1, so it numerically DAMPS oscillators — a slow,
// monotone, secular energy decay (ΔE/E ≈ -(omega*dt)^6/72 per step) rather
// than Verlet's bounded oscillation. At sane step sizes this decay is
// minuscule (~1e-9 per step at dt = T/100), but over a very long run (say,
// 10^7 periods of a planetary orbit) it accumulates into a real, visible
// spiral-in that a symplectic method would never produce — this is the
// "phase-perfect but energy-wrong over aeons" mirror image of Verlet's
// "energy-perfect but phase-wrong over aeons." See
// docs/research/numerical-integration.md §2.4 for the numbers.
//
// When to reach for this file instead of velocity_verlet.hpp: any
// velocity-dependent force (drag, magnetic, damping) where you want the
// full 4th-order accuracy rather than Verlet's predictor–corrector
// second-order patch; any problem where you care about hitting a specific
// TIME accurately (a targeted rendezvous, a projectile's range) rather than
// preserving long-run structure; moderate horizons in general. Prefer
// Verlet for millions of periods of genuinely conservative dynamics, where
// its energy-boundedness matters more than pointwise accuracy.
//
// C++ note: written out stage-by-stage rather than as a generic Butcher-
// tableau loop — slower to extend to other RK methods, but every line is a
// direct transcription of the formulas above, which is the point for a
// library meant to be read alongside the maths.

#include <newton/core/particle.hpp>

namespace newton::integrators {

/// One step of classic (4-stage) Runge-Kutta, RK4.
///
/// Physics: fourth-order accurate (global error O(dt^4)); costs 4 force
/// evaluations per step; not symplectic — introduces a slow, monotone
/// energy decay on oscillators rather than Verlet's bounded wobble (see
/// file header). The accuracy workhorse for velocity-dependent forces and
/// for problems where pointwise accuracy in time matters more than
/// long-run structure.
/// @param s      current state (r, v, t)
/// @param mass   particle mass [kg]
/// @param dt     step size [s]
/// @param force  callable: Vec3 operator()(double t, const Vec3& r, const Vec3& v) const, returns force [N]
/// @returns the state at s.t + dt
template <class Force>
[[nodiscard]] State rk4_step(const State& s, double mass, double dt, Force&& force) {
    // Stage 1: at the state itself. k1r before k1v — k1v does not actually
    // need k1r, but computing "position slope, then velocity slope" at every
    // stage is the pattern that keeps the ordering rule uniform and easy to
    // audit stage-to-stage.
    const Vec3 k1r = s.v;
    const Vec3 k1v = force(s.t, s.r, s.v) / mass;

    // Stage 2: at the midpoint, predicted using stage 1.
    const Vec3 k2r = s.v + 0.5 * dt * k1v;
    const Vec3 k2v = force(s.t + 0.5 * dt, s.r + 0.5 * dt * k1r, s.v + 0.5 * dt * k1v) / mass;

    // Stage 3: at the midpoint again, but now predicted using stage 2 —
    // a second, refined estimate of the midpoint state.
    const Vec3 k3r = s.v + 0.5 * dt * k2v;
    const Vec3 k3v = force(s.t + 0.5 * dt, s.r + 0.5 * dt * k2r, s.v + 0.5 * dt * k2v) / mass;

    // Stage 4: at the endpoint, predicted using stage 3.
    // (k4r = s.v + dt*k3v is never named as its own variable: the collapsed
    // position update below absorbs it algebraically — see the comment
    // there. k4v itself still needs r_n + dt*k3r as its position argument.)
    const Vec3 k4v = force(s.t + dt, s.r + dt * k3r, s.v + dt * k3v) / mass;

    // Collapsed position update. Because k1r, k2r, k3r, k4r are each just
    // v_n plus a multiple of dt times a k_*v, the weighted sum
    // (dt/6)*(k1r + 2*k2r + 2*k3r + k4r) telescopes:
    //   k1r + 2*k2r + 2*k3r + k4r
    //     = v_n + 2(v_n + 0.5dt*k1v) + 2(v_n + 0.5dt*k2v) + (v_n + dt*k3v)
    //     = 6*v_n + dt*(k1v + k2v + k3v)
    // so (dt/6)*(...) = dt*v_n + (dt^2/6)*(k1v + k2v + k3v) exactly — no k4v
    // term survives, and no k_ir needs to be kept around after it has been
    // used inside a force call. This both saves arithmetic and makes the
    // "position advances by v*dt plus an acceleration correction" structure
    // visible in the code, not just on paper.
    const Vec3 r_new = s.r + dt * s.v + (dt * dt / 6.0) * (k1v + k2v + k3v);

    // Velocity update keeps the full 4-term Simpson-like weighting — nothing
    // collapses here because the k_*v are the primary unknowns, not derived
    // quantities.
    const Vec3 v_new = s.v + (dt / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);

    return State{r_new, v_new, s.t + dt};
}

}  // namespace newton::integrators
