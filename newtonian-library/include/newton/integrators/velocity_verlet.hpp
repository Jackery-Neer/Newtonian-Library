#pragma once
// newton/integrators/velocity_verlet.hpp — kick-drift-kick Verlet.
//
// Physics: for a force that depends only on (t, r) — no velocity — velocity
// Verlet is the textbook second-order symplectic workhorse. Written in
// "kick-drift-kick" (KDK) form:
//
//   a_n     = F(t_n, r_n, v_n) / m
//   v_half  = v_n + 0.5*dt*a_n            // half kick
//   r_{n+1} = r_n + dt*v_half             // full drift
//   a_{n+1} = F(t_{n+1}, r_{n+1}, ...) / m
//   v_{n+1} = v_half + 0.5*dt*a_{n+1}     // half kick
//
// It is exact (to roundoff) for constant acceleration at ANY dt — that is
// the very first thing to unit-test, and it instantly catches sign errors
// and dropped-1/2 typos. Its conserved-quantity story: bounded energy
// oscillation of relative amplitude ~ (omega*dt)^2/4 with NO secular drift
// (it is symplectic and time-reversible, like semi_implicit_euler, but one
// order more accurate), while its numerical frequency is slightly off
// (omega_num = omega*(1 + (omega*dt)^2/24 + ...)), so PHASE drifts linearly
// in time even though energy does not — "gets structure and statistics
// right forever, gets long-run phase wrong." See
// docs/research/numerical-integration.md §2.3 for the full derivation; this
// is the trade you are making every time you pick Verlet over RK4.
//
// ---------------------------------------------------------------------------
// THE SUBTLETY: what if the force also depends on velocity? (drag, damping)
// ---------------------------------------------------------------------------
//
// Physics: the derivation above silently assumes F = F(t, r). The moment F
// also depends on v, the final half-kick becomes IMPLICIT:
//
//   v_{n+1} = v_half + 0.5*dt * F(t_{n+1}, r_{n+1}, v_{n+1}) / m
//                                                     ^^^^^^^ appears on both sides
//
// You cannot literally write that line in code — v_{n+1} is what you're
// trying to compute. Four ways to resolve this are catalogued in the
// research doc (§2.3); this file implements resolution (2), the
// PREDICTOR–CORRECTOR, because it is the right default for a student
// library: it costs nothing extra, needs no special-casing per force, and
// stays second order. The other three (closed-form solve for pure linear
// drag, fixed-point iteration, operator splitting / exponential drag) are
// each better in some narrower situation but need either a force-specific
// closed form or extra iteration — see the research doc if you outgrow this.
//
// The fix: PREDICT v_{n+1} with an explicit, Euler-quality estimate, use the
// prediction (only) to evaluate the force, then correct:
//
//   v_pred  = v_half + 0.5*dt*a_n                  // predictor (== v_n + dt*a_n, algebraically)
//   a_{n+1} = F(t_{n+1}, r_{n+1}, v_pred) / m       // the ONLY place v_pred is used
//   v_{n+1} = v_half + 0.5*dt*a_{n+1}               // corrector: real half-kick, using a_{n+1}
//
// WHY THIS STAYS SECOND ORDER (the part worth understanding, not just
// copying): the predictor v_pred is only a first-order (Euler-quality)
// estimate of v_{n+1}, so it carries a local error of O(dt^2). But v_pred is
// used only to evaluate a force, and that force's contribution to v_{n+1} is
// itself scaled by 0.5*dt in the corrector step. An O(dt^2) error in the
// FORCE ARGUMENT produces, after being multiplied by dt in the update, an
// O(dt^3) local error in v_{n+1} — which is exactly the same order as
// velocity Verlet's own local truncation error. In other words: an
// approximate velocity is good enough here precisely because it only ever
// feeds into a term that is already being damped by one power of dt. Get
// this wrong — e.g. skip the predictor and just reuse v_half or v_n
// directly in place of v_pred — and you silently lose an order of accuracy
// (this is exactly what tests/test_integrators.cpp's convergence-order test
// is designed to catch: watch the fitted order for a drag problem drop from
// ~2 to ~1). Symplecticity itself is NOT preserved once F depends on v — but
// a velocity-dependent force is generally non-conservative anyway (drag
// dissipates energy), so there was no exact conservation law to lose in the
// first place; what you keep is the accuracy order, which is the thing that
// actually mattered.
//
// C++ note: this is a STANDALONE, stateless step function — it evaluates F
// TWICE per call (once for a_n, once for a_{n+1} using v_pred), because it
// has nowhere to cache a_n from the previous call. A production integrator
// would keep an integrator OBJECT that carries a_n forward and gets exactly
// one force evaluation per step (see research doc §5.4); that caching
// version is a deliberate v0.2 optimization, left out here so the KDK
// structure and the predictor–corrector fix stay the whole story, with
// nothing tucked away in mutable state.

#include <newton/core/particle.hpp>

namespace newton::integrators {

/// One step of velocity Verlet (kick-drift-kick), with predictor-corrector
/// velocity handling so velocity-dependent forces (drag, magnetic, damping)
/// stay second order — see the file header for why that predictor step is
/// necessary and why it doesn't cost you an order of accuracy.
///
/// Physics: second-order accurate (global error O(dt^2)) in both r and v;
/// exact for constant acceleration at any dt; energy-bounded and
/// time-reversible for F(t, r) forces, merely second-order (not symplectic)
/// once F depends on v.
/// @param s      current state (r, v, t)
/// @param mass   particle mass [kg]
/// @param dt     step size [s]
/// @param force  callable: Vec3 operator()(double t, const Vec3& r, const Vec3& v) const, returns force [N]
/// @returns the state at s.t + dt
template <class Force>
[[nodiscard]] State velocity_verlet_step(const State& s, double mass, double dt, Force&& force) {
    const Vec3 a_n = force(s.t, s.r, s.v) / mass;      // force eval #1, at the OLD state

    const Vec3 v_half = s.v + 0.5 * dt * a_n;          // half kick
    const Vec3 r_new = s.r + dt * v_half;               // full drift

    const Vec3 v_pred = v_half + 0.5 * dt * a_n;        // predictor: explicit estimate of v_{n+1}
    const Vec3 a_new = force(s.t + dt, r_new, v_pred) / mass;  // force eval #2, at the NEW (r, v_pred)

    const Vec3 v_new = v_half + 0.5 * dt * a_new;       // corrector: real half kick with a_{n+1}

    return State{r_new, v_new, s.t + dt};
}

}  // namespace newton::integrators
