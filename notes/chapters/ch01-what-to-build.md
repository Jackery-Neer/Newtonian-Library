# Chapter 1 → what you can build now

A build plan for `newton` v0.1 at `/Users/jackery/Github/newtonian-library`, scoped to exactly what Taylor ch. 1 gives you. Nothing here needs physics you haven't read.

> **How this document was made (2026-09-10).** Two Opus agents researched — one the chapter itself (from public problem indexes and course pages; the book's text is not available online and nothing from it is reproduced here), one the actual repo. A third Opus agent turned both into proposals. Fable agents then verified each piece against the code and the physics, sent every inconsistency back to Opus for revision, and re-verified — three rounds — before the parent Fable model did a final pass and applied the last corrections by hand. The full research inputs are in `docs/research/chapter-01-research.md` and `chapter-01-build-proposals.md` if you want the long form. The third input, a September snapshot of the repo (`chapter-01-repo-audit.md`), was deleted on 2026-11-26 as superseded; the current coverage picture is `notes/progress/2026-11-24-coverage-audit.md`. Every "measured" number below came from prototype programs compiled against the real headers.

---

## 1. TL;DR

- **Chapter 1 demands exactly one numerical thing from a simulator.** Problems 1.50 and 1.51 are the chapter's only two `[Computer]` problems: integrate `phi_ddot = -(g/R) sin(phi)` from rest at 20 deg and at 90 deg, and plot it against `phi_0 cos(sqrt(g/R) t)`. That's it. Everything else in the chapter is algebra, closed forms, and frame bookkeeping.
- **Chapter 1 also hands you ~20 free unit tests.** Problems 1.1–1.9 and 1.13–1.23 are closed-form vector identities you can assert to machine epsilon. That is the cheapest confidence you will ever buy.
- **The repo's first half of ch. 1 is done and well tested.** `Vec3`, the `(r, v)` state recast, the four integrators, the Cartesian force path, and the projectile oracle all exist and are exercised.
- **The repo's second half of ch. 1 is missing entirely.** No polar coordinates (§1.7), no frames (§1.4), no pairwise forces or system momentum (§1.5). Polar is the biggest single gap — and `docs/research/taylor-ch1-7-map.md` §1.4 already specified the module, so `DESIGN.md` is the document that's out of date, not your plan.
- **Top 3 to build next, in order: (1) `newton::coords` polar, (2) `analytic::position_const_accel`, (3) the half-pipe `systems::CircularTrack`.** Polar unlocks §1.7; the constant-acceleration closed form takes ten minutes and makes the library's sharpest integrator test callable; the half-pipe is Taylor's actual homework.
- **Use momentum conservation and turning-point drift as your integrator metrics, not energy.** Chapter 1 has no concept of energy, so "symplectic integrators conserve energy" is a sentence you are not yet entitled to. Turning points of a periodic motion must not wander — that says the same thing with ch. 1 vocabulary.
- **The turning-point metric is sharp, but only if you chain the legs correctly.** Measured over a genuine 6 periods at `dt = 1e-3`: RK4 drifts **9.4e-15 rad**, Euler drifts **9.3e-3 rad** — twelve orders of magnitude. Proposal 6 below spells out the one-step "nudge" that makes the chained measurement actually cover 6 periods instead of silently stopping at 4.
- **One decision to make now while it's free:** add a `depends_on_velocity` trait to your force functors while every answer is still `false`. Velocity Verlet assumes `F = F(r)`; ch. 2's drag breaks that.

---

## 2. Chapter 1 in one page

Seven sections, a summary page, then problems 1.1–1.51.

| § | Title | What it gives your library |
|---|---|---|
| 1.1 | Classical Mechanics | Domain of validity. Nothing to implement. |
| 1.2 | Space and Time | The whole vector toolkit. Your `Vec3` spec and ~20 tests. |
| 1.3 | Mass and Force | Forces **sum**. Net force is a vector sum, never a replacement. |
| 1.4 | Newton's First and Second Laws; Inertial Frames | `F = dp/dt`, Galilean transformations, what makes a frame inertial. |
| 1.5 | The Third Law and Conservation of Momentum | `F_12 = -F_21`, strong (central) form, `dP/dt = F_ext`. |
| 1.6 | Newton's Second Law in Cartesian Coordinates | Why Cartesian is the right internal representation. Tilted axes. |
| 1.7 | Two-Dimensional Polar Coordinates | The rotating basis, the half-pipe, the chapter's only numerics. |

### The equations you will actually implement

```
§1.2  vectors
  r.s     = r s cos(theta)            (1.6)   and   sum_i r_i s_i   (1.7)
  r x s   = (r_y s_z - r_z s_y, r_z s_x - r_x s_z, r_x s_y - r_y s_x)   (1.9)
  |r x s| = r s sin(theta), right-hand rule
  d/dt (r.s)   = r.s_dot + r_dot.s
  d/dt (r x s) = r x s_dot + r_dot x s        ORDER MATTERS
  a x (b x c)  = (a.c) b - (a.b) c             "BAC minus CAB"

§1.3-1.4  dynamics
  F_net = sum of all forces acting     w = m g
  p = m v      F = dp/dt      F = m a (constant mass)
  F = 0  =>  v = const                         Newton I
  boost:       r' = r - V t,   v' = v - V,    a' = a      (S' inertial)
  accelerate:  r'' = r - (1/2) A t^2,         a'' = a - A (S'' NOT inertial)

§1.5  systems
  F_12 = -F_21                         weak third law
  F_12 parallel to (r_1 - r_2)         strong (central) form
  P = sum p_alpha     dP/dt = F_ext    F_ext = 0 => P const
  |F_mag_12| <= (v_1 v_2 / c^2) |F_el_12|      where the third law breaks

§1.6  Cartesian
  m x_ddot = F_x,  m y_ddot = F_y,  m z_ddot = F_z        (1.35)
  tilted axes on a slope at phi:  W = (-m g sin phi, -m g cos phi, 0)

§1.7  polar
  x = r cos phi,  y = r sin phi,  r = sqrt(x^2+y^2),  phi = atan2(y,x)   (1.37)
  r_hat   = ( cos phi, sin phi, 0)                                        (1.59)
  phi_hat = (-sin phi, cos phi, 0)
  d r_hat/dt = phi_dot phi_hat     (1.42)    d phi_hat/dt = -phi_dot r_hat   (1.46)
  v = r_dot r_hat + r phi_dot phi_hat
  a = (r_ddot - r phi_dot^2) r_hat + (r phi_ddot + 2 r_dot phi_dot) phi_hat
  F_r = m (r_ddot - r phi_dot^2),  F_phi = m (r phi_ddot + 2 r_dot phi_dot)   (1.48)
  phi_ddot = -(g/R) sin phi                                               (1.51)
  small angle: phi_ddot = -w^2 phi, w = sqrt(g/R)                          (1.55)
               phi(t) = A sin(wt) + B cos(wt)                              (1.56)
               phi(t) = phi_0 cos(sqrt(g/R) t),  T = 2 pi sqrt(R/g)        (1.57)
```

### Worked examples, paraphrased

**Example 1.2 — the oscillating skateboard (§1.7).** A point mass is released from rest inside a frictionless circular trough of radius `R` and slides back and forth. Work in polar coordinates centred on the trough's centre of curvature: `r = R` is pinned, so `r_dot = r_ddot = 0`, and eq. (1.48) collapses. The tangential component gives the whole equation of motion, `phi_ddot = -(g/R) sin(phi)` — eq. (1.51). The radial component is not thrown away: it determines the normal force. Taylor states outright that (1.51) has no elementary closed-form solution, which is the book's first explicit argument for numerical methods. Then he linearises it and gets SHM with period `2 pi sqrt(R/g)`, amplitude-independent. Note what this problem *is*: a constraint problem, not a string-tension problem — identical mathematics, different physical story.

**Example 1.1 — a block on an incline (§1.6).** *(Placement is from memory, not web-confirmed — see §6.)* A block slides on an inclined plane. The lesson is the axis choice: one axis along the slope, one normal to it, so gravity resolves as `mg sin(theta)` down-slope and `mg cos(theta)` into the surface. The normal force is found not from a force law but from the kinematic condition that the perpendicular acceleration is zero. That one move — constraint force from a kinematic condition — is the technique you'll reuse in the half-pipe and in Atwood.

**Two set-piece derivations in the running text, not numbered examples.** The N-particle momentum-conservation derivation, eqs. (1.25)–(1.29), where internal forces cancel pairwise. And a counterexample (Figure 1.8) where two charges moving with perpendicular velocities exert magnetic forces that violate `F_12 = -F_21`; the violation is of relative order `v1 v2 / c^2`, and momentum is rescued only by giving the field its own momentum.

### Problem archetypes

| Problems | Roughly | Archetype |
|---|---|---|
| 1.1–1.9, 1.13–1.23 | §1.2 | Vector identities, dot/cross proofs, triple products, law of cosines/sines, triangle and parallelepiped geometry. **All closed-form: a ready-made test suite.** |
| 1.10–1.12 | §1.2 | Given a trajectory, differentiate it: uniform circular motion, a prescribed ellipse, an elliptical helix. No dynamics at all. |
| 1.24, 1.25, 1.44 | — | ODE warm-ups. Taylor assumes **zero** prior ODE knowledge. |
| 1.26, 1.27, 1.46 | §1.4 | Same motion described from a boosted frame, an accelerating frame, and a rotating turntable. Which are inertial? |
| 1.28–1.34 | §1.5 | Momentum derivation written out for N = 3 and N = 4; inelastic collision; the converse (momentum conservation implies the third law); where the third law fails for moving charges; angular momentum given central forces. |
| 1.35–1.40 | §1.6 | Six constant-force problems that differ in axis choice: four in-flight projectiles (1.35 golf ball, 1.36 bundle dropped from a plane, 1.39 ball thrown up a slope, 1.40 the monotonic-range condition) and two pucks sliding *on* an incline with a normal force present (1.37 up the slope, 1.38 across a tilted board). |
| 1.41–1.49 | §1.7 | Twirled mass on a string in gravity-free space; the `atan2` quadrant fix; deriving the rotating basis; the turntable quantitatively; cylindrical coordinates; puck between coaxial cylinders. |
| **1.50, 1.51** | §1.7 | **The only two `[Computer]` problems.** The half-pipe at 20 deg and at 90 deg. |

---

## 3. Coverage matrix

| Ch. 1 concept | Status | Where in the repo |
|---|---|---|
| Vector algebra: dot / cross / norm (§1.2) | **Implemented** | `include/newton/core/vec3.hpp` — `dot`, `cross`, `norm2`, `norm`, `normalized` (zero-guarded), operators, doc comments citing §1.2. Tested in `tests/test_vec3.cpp`, including `dot(cross(a,b), a) == 0`. Missing: `angle`, triple products. → **Proposal 1** |
| Newton II as a first-order system (§1.6) | **Implemented** | `include/newton/core/particle.hpp` — `State {r, v, t}` is the `(r, v)` pair, the file header derives it. All four steppers are `State -> State`. Driver: `newton::simulate` / `simulate_until` in `include/newton/simulate.hpp`. |
| Cartesian component ODEs (§1.6) | **Implemented** | Carried componentwise by `Vec3` throughout; e.g. `include/newton/integrators/rk4.hpp` advances `r` and `v` as whole vectors. The decoupling case is stated in `include/newton/analytic/projectile.hpp`'s header. |
| Superposition of forces (§1.3) | **Implemented** | `include/newton/forces/compose.hpp` — `forces::sum(f1, f2, ...)`. |
| Uniform gravity / constant force (§1.6) | **Implemented** | `include/newton/forces/gravity.hpp::UniformGravity`, `include/newton/forces/drives.hpp::ConstantForce`. |
| General constant-acceleration closed form (§1.6) | **Partial** | Only the specialised case exists: `include/newton/analytic/projectile.hpp::position(t, v0, theta, g)`, hardcoded to launch from the origin with `a = (0,-g,0)`. `tests/test_integrators.cpp` writes `r0 + v0*T + 0.5*T*T*a` inline rather than calling the library. → **Proposal 2** |
| Tilted axes / normal force (§1.6) | **Missing** (no new code needed) | Nothing in `examples/` does an incline, though it is expressible today with a lambda satisfying the existing `Force` concept. → **Proposal 3** |
| Polar conversion (§1.7) | **Missing** | No `coords` module, no `to_polar`/`from_polar`. `docs/research/taylor-ch1-7-map.md` §1.4 specifies one; `DESIGN.md` and `include/newton/newton.hpp` silently omit it. → **Proposal 4** |
| Rotating unit vectors `r_hat`, `phi_hat` and their derivatives (§1.7) | **Missing** | Nothing computes them. This is the conceptual heart of §1.7 and has no representation in code or notes. → **Proposal 4** |
| Polar velocity / acceleration decomposition (§1.7) | **Missing** | No `polar_acceleration`. `notes/02-physics-to-code.md`'s ch. 1 section is titled "(Taylor §1.1–1.7)" but its key-equations block stops at `F = dp/dt = m a`; §1.7 contributes nothing. → **Proposal 4** |
| Uniform circular motion (§1.7) | **Partial — simulable, no ch. 1 API** | Producible today (`forces::Spring` with `natural_length = 0`, or `PointGravity`); `tests/test_analysis.cpp` hand-builds a circular fixture as scaffolding. No `a = v^2/r` helper in `newton::analytic`. Nearest equivalents are `include/newton/analytic/magnetic.hpp`'s `cyclotron_*`, which are §2.5 and unreadable for a ch. 1 student. → **Proposal 5** |
| Half-pipe / pendulum, problems 1.50–1.51 (§1.7) | **Missing** | No constraint machinery, no example. Doubly blocked: needs polar reasoning and a constraint. → **Proposal 6** |
| Inertial frames / Galilean transforms (§1.4) | **Missing** | No frame-transformation or inertial-frame concept appears anywhere in `include/`, `examples/`, or `tests/` (the only physics use of the word "frame" is `systems/rocket.hpp:19`'s rest frame in the momentum argument). Absent from `notes/02-physics-to-code.md` and from the v0.2 roadmap in `notes/04-extending.md` — the one ch. 1 gap not even acknowledged as future work. → **Proposal 7** |
| Third law & two-body momentum (§1.5) | **Partial — diagnostic only** | `include/newton/analysis/momentum.hpp::total_momentum(const std::vector<Particle>&)` computes `P = sum m_i v_i`, tested in `tests/test_analysis.cpp`. But Newton III appears nowhere in executable code: the Force concept is `Vec3 operator()(double t, const Vec3& r, const Vec3& v) const` and has no way to name a second particle. → **Proposal 9** |
| Centre of mass (§1.5 momentum; the `R_cm` *definition* is Taylor §3.3 — used here only as the one-line shorthand `P = M V_cm`) | **Partial** | `include/newton/analysis/momentum.hpp::center_of_mass` with an `M <= 0` guard, tested. Missing the dynamical half: nothing evolves more than one particle, so `R(t) = R_0 + V_0 t` cannot be demonstrated. No `total_mass` / `velocity_of_cm`. → **Proposal 9** |
| Multi-particle systems | **Missing** | `DESIGN.md` line 15 scopes the library as "single-particle Newtonian mechanics". `Trajectory` holds one `double mass`; `simulate()` takes `const Particle&`. `analysis/momentum.hpp`'s own C++ note documents this as a v0.1 limitation. → **Proposal 9** |
| Constraint problems (Atwood, rigid pendulum) | **Missing** | No constraint machinery. Atwood is buildable by hand as one DOF; a general solver is the `notes/04-extending.md` Lagrangian hook, correctly v0.2. → **Proposal 8** |

Of the 27 headers under `include/` (25 plus the two umbrellas, `newton.hpp` and `forces/forces.hpp`), 8 are fully ch.-1-readable, 7 partly, 10 later-chapter — §5 spells out which. Of the six examples, **one** (`examples/01_free_fall.cpp`) is ch. 1. That's the case for building.

---

## 4. Build-now projects, in recommended order

Ten proposals, ordered so each one's dependencies are already done. Every tolerance and "measured" figure below came out of prototype programs compiled against the real headers with `clang++ -std=c++17 -Wall -Wextra -Wpedantic -O2 -I include` and run.

---

### Proposal 1 — Finish §1.2 on `Vec3`: `angle`, triple products, twenty free tests

**Pitch.** Add the three vector functions ch. 1 uses and turn problems 1.4–1.23 into a test suite. Cheapest high-confidence win in the repo.

**Physics (§1.2).**
```
a.b     = |a||b| cos(theta)        (1.6)   and   sum_i a_i b_i   (1.7)
|a x b| = |a||b| sin(theta), right-hand rule       (1.9)
a.(b x c)                 scalar triple = signed parallelepiped volume
a x (b x c) = (a.c) b - (a.b) c                   "BAC minus CAB"
d/dt (r.s)   = r.s_dot + r_dot.s
d/dt (r x s) = r x s_dot + r_dot x s              ORDER MATTERS
```

**Files.**
- modify `/Users/jackery/Github/newtonian-library/include/newton/core/vec3.hpp`
- modify `/Users/jackery/Github/newtonian-library/tests/test_vec3.cpp`
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md`
- Umbrella line: **none** — `#include <newton/core/vec3.hpp>` is already line 19 of `include/newton/newton.hpp`.

Spec friction to resolve deliberately: `DESIGN.md:57` heads the core section "Already written (do not modify): core". That was an instruction to the v0.1 implementation agents, not a permanent freeze, but `notes/00-start-here.md` makes `DESIGN.md` authoritative — so amend `DESIGN.md:59-60` to list the three new functions in the same commit. (If you'd rather keep the freeze literal, put them in a new `include/newton/core/vec3_algebra.hpp`; I recommend against it — they belong next to `dot` and `cross`.)

**Signatures.**
```cpp
namespace newton {

/// Angle between two vectors, in [0, pi] (Taylor §1.2, problems 1.4–1.5).
/// @param a,b  any two nonzero vectors; returns 0 if either is zero
/// C++ note: atan2(|a x b|, a.b), NOT acos(a.b/(|a||b|)). acos loses
/// precision catastrophically near 0 and pi (its derivative is infinite
/// there) and can be handed an argument a hair outside [-1, 1] by rounding,
/// returning NaN. The atan2 form is accurate over the whole range and needs
/// no clamping. Not constexpr: std::atan2 isn't, in C++17 (see norm()).
[[nodiscard]] inline double angle(const Vec3& a, const Vec3& b);

/// Scalar triple product a . (b x c) — the signed volume of the
/// parallelepiped spanned by a, b, c (Taylor problem 1.21). Cyclic:
/// a.(bxc) == b.(cxa) == c.(axb), and zero iff the three are coplanar.
[[nodiscard]] constexpr double scalar_triple(const Vec3& a, const Vec3& b, const Vec3& c);

/// Vector triple product a x (b x c) = (a.c) b - (a.b) c — the "BAC minus
/// CAB" identity (Taylor problem 1.33 expects you to have it). Note a x (b x c)
/// != (a x b) x c: the cross product is NOT associative.
[[nodiscard]] constexpr Vec3 vector_triple(const Vec3& a, const Vec3& b, const Vec3& c);

}  // namespace newton
```
`scalar_triple` and `vector_triple` can be `constexpr` (pure arithmetic); `angle` cannot. Worth a sentence in the doc comments.

**Validation test.** All measured on `a = (1,2,3)`, `b = (-2,0.5,4)`, `c = (3,-1,0.25)`.

| Taylor | Assertion | Tol | Measured |
|---|---|---|---|
| 1.4 | `angle(a,b) == acos(dot(a,b)/(norm(a)*norm(b)))` | 1e-12 | — |
| 1.5 | `angle({1,1,1},{1,1,0}) == acos(sqrt(2/3))` = 35.2643896828 deg | 1e-12 | agrees to 12 digits |
| 1.7 | `dot(a,b) == norm(a)*norm(b)*cos(angle(a,b))` | 1e-12 | — |
| 1.9 | `norm2(a+b) - (norm2(a)+norm2(b)+2*dot(a,b)) == 0` (law of cosines) | 1e-12 | exactly 0 |
| 1.13 | `norm2(b) == dot(u,b)^2 + norm2(cross(u,b))`, `u = normalized(a)` | 1e-12 | exactly 0 |
| 1.14 | `norm(a+b) <= norm(a)+norm(b)` | exact `<=` | 7.5 <= 8.24166 |
| 1.15 | `norm(cross(a,b)) == norm(a)*norm(b)*sin(angle(a,b))` | 1e-12 | — |
| 1.18 | sides head-to-tail, so build `c = -(a + b)`: `0.5*norm(cross(a,b))`, `0.5*norm(cross(b,c))`, `0.5*norm(cross(c,a))` all equal (checked against the book, 2026-09-15 — an earlier draft mislabeled a different triangle problem here) | 1e-12 | — |
| (vectors-from-origin triangle; a different problem number — check your copy) | `0.5*norm(cross(b-a, c-a))` equals `0.5*norm(cross(b,c)+cross(c,a)+cross(a,b))` | 1e-12 | — |
| 1.21 | `scalar_triple` cyclic under a→b→c→a | 1e-12 | exactly 0, both rotations |
| 1.22 | unit vectors at th1=0.7, th2=1.9: `dot(e1,e2) == cos(th2-th1)`, `cross(e1,e2).z == sin(th2-th1)` | 1e-15 | 5.6e-17, 1.1e-16 |
| 1.33 | `vector_triple(a,b,c) == dot(a,c)*b - dot(a,b)*c` | 1e-12 | exactly 0 |
| 1.8(b), 1.17(b) | product rules by **central finite difference at h = 1e-5** on `r(t) = (cos t, 2t, e^{0.3t})`, `s(t) = (t^2, sin 2t, 1-t)`, swept over t in [0.1, 1.9] | **1e-8** | dot-rule 5.8e-10, cross-rule 6.2e-10 |

Two things about that last row. The FD tolerance is step-size-dependent and has a floor: at h = 1e-4 the error is 6.6e-8 (truncation-dominated, so 1e-8 would fail); at h = 1e-6 it is 2.7e-9 (roundoff-dominated). h = 1e-5 is the sweet spot — write that reasoning into the test as a comment so nobody "tightens" it later. And the test genuinely catches the operand-order bug: the wrong order `cross(s_dot,r) + cross(s,r_dot)` differs from the right answer by **16.36** at t ~ 0.95, nine orders above the tolerance.

**Example.** None — test-only.

**Difficulty.** S, 2–3 h.

**What you learn.** *Physics:* that the geometric and component definitions of the dot product being the same thing is a provable claim, not a convention; that `a x (b x c) != (a x b) x c`; that a finite-difference check of a product rule is a real test of your code. *C++:* `constexpr` vs `inline` and why `std::atan2`'s non-`constexpr`-ness in C++17 forces the distinction (the same trade-off `norm()` already makes); numerical conditioning of `acos` vs `atan2`; how the harness's `NT_CHECK(...)` variadic macro lets you write `NT_CHECK(v == Vec3{1,2,3})` with top-level commas inside braces — read `tests/test_util.hpp:87-94` before you need it.

---

### Proposal 2 — `analytic::position_const_accel`: ch. 1's most basic result, currently uncallable

**Pitch.** A general `r(t) = r0 + v0 t + (1/2) a t^2` in `newton::analytic`, so the sharpest integrator test in the library stops being an inline formula inside a test file.

**Physics (§1.6).**
```
F constant  =>  a = F/m constant  =>
    r(t) = r0 + v0 t + (1/2) a t^2      degree-2 polynomial in t
    v(t) = v0 + a t
```
That is the whole content of §1.6 when `F` does not depend on `r`, `v`, or `t`. `analytic/projectile.hpp` covers only `r0 = 0`, `a = (0,-g,0)`.

**Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/analytic/kinematics.hpp`
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — add **`#include <newton/analytic/kinematics.hpp>`** as the first line of the analytic block (immediately before line 36's `projectile.hpp`)
- modify `/Users/jackery/Github/newtonian-library/tests/test_integrators.cpp` (replace the inline `r0 + v0*T + 0.5*T*T*a` in block (a) with the library call, then strengthen the assertions)
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — new bullet in the Agent A3 / `newton::analytic` section

**Signatures.**
```cpp
namespace newton::analytic {

/// Position under constant acceleration: r(t) = r0 + v0 t + ½ a t²
/// (Taylor §1.6). The exact solution of m r̈ = F for constant F, and the
/// library's sharpest integrator test: a degree-2 polynomial in t, which
/// RK4 and velocity Verlet reproduce to roundoff at ANY step size while both
/// Euler variants do not.
/// @param t  elapsed time [s]    @param r0 initial position [m]
/// @param v0 initial velocity [m/s]   @param a constant acceleration [m/s^2]
[[nodiscard]] constexpr Vec3 position_const_accel(double t, const Vec3& r0,
                                                 const Vec3& v0, const Vec3& a);

/// Velocity under constant acceleration: v(t) = v0 + a t (Taylor §1.6).
[[nodiscard]] constexpr Vec3 velocity_const_accel(double t, const Vec3& v0, const Vec3& a);

/// Both at once, as a State ready to compare against a Trajectory sample.
/// @param s0 the state at t = s0.t; the returned State carries time s0.t + dt
[[nodiscard]] constexpr State state_const_accel(const State& s0, double dt, const Vec3& a);

}  // namespace newton::analytic
```
All three are `constexpr` — no `sqrt`, no trig. Worth a C++ note contrasting with `analytic/projectile.hpp`, which cannot be.

**Validation test.** Fixture: `m = 2.5 kg`, `F = (2, -7, 0.5) N`, `r0 = (1,2,-3) m`, `v0 = (4,-5,6) m/s`, `dt = 0.1 s`, `T = 5 s` (so `|a| = 2.9189 m/s^2`, `|r(T)| = 72.08 m`).

| Method | Assertion | Measured |
|---|---|---|
| RK4 | `close_vec3(final.r, position_const_accel(T,...), 1e-12)` | `|dr|` = 4.9e-14 m, `|dv|` = 2.2e-14 |
| VelocityVerlet | same, 1e-12 | `|dr|` = 7.8e-14 m, `|dv|` = 4.2e-14 |
| Euler | `NT_CLOSE_VEC3(final.r - exact, -0.5*a*T*dt, 1e-9)` and `NT_CHECK(norm(final.r - exact) > 0.1)` | error exactly `(-0.2, 0.7, -0.05)`; matches `-½aTdt` to 6.0e-14; `|err|` = 0.729726 m |
| SemiImplicitEuler | `NT_CLOSE_VEC3(final.r - exact, +0.5*a*T*dt, 1e-9)` | error exactly `(+0.2, -0.7, +0.05)`; matches to 6.6e-14 |

The Euler rows are better than "Euler is not exact": the two first-order methods **bracket** the exact answer with errors that are equal and opposite and exactly `∓(1/2) a T dt`, because Euler drifts the position with the *old* velocity and Euler–Cromer with the *new* one. Both get velocity exactly right (2.2e-14) because `a` is constant. Assert all of that. Also keep `DESIGN.md`'s existing check (d): final `state.t == t_end` to 1e-12.

**Example.** None — fold a one-line use into `examples/07_incline.cpp`.

**Difficulty.** S, 1–2 h.

**What you learn.** *Physics:* why "exact for constant acceleration" is the *first* thing to test in any integrator, and why it discriminates first-order from second-order methods with no convergence study at all. *C++:* `constexpr` functions on an aggregate type; why a test that re-derives a formula inline is a test of two things at once and therefore a worse test.

---

### Proposal 3 — `examples/07_incline.cpp`: tilted axes and the normal force (1.37–1.39)

**Pitch.** The §1.6 lesson that you may *choose* your axes, worked three ways against three closed forms — and zero new library code.

**Physics (§1.6).**
```
W = (-m g sin(phi), -m g cos(phi), 0)      tilted components of m g

puck ON the slope (problem 1.37): the surface supplies N = m g cos(phi),
   which exactly cancels the normal component (a_normal = 0 is a KINEMATIC
   condition, not a force law — that is how constraint forces are found),
   leaving  m x_ddot = -m g sin(phi),   t_return = 2 v0 / (g sin phi)

ball IN FLIGHT above the slope (problem 1.39): there is NO normal force —
   the ball is in the air. The tilted axes are pure bookkeeping.
   R     = 2 v0^2 sin(theta) cos(theta + phi) / (g cos^2(phi))
   R_max = v0^2 / [g (1 + sin(phi))]        at  2 theta + phi = pi/2
```
The "no normal force in flight" trap deserves a loud comment; it is the entire subject of the public confusion around problem 1.39.

**Files.** Create `/Users/jackery/Github/newtonian-library/examples/07_incline.cpp`. **No header changes, no umbrella line, no Makefile edit** — `EXAMPLE_SRCS := $(wildcard examples/*.cpp)` (Makefile line 33) picks it up, and `HEADERS := $(shell find include -name '*.hpp')` is a blanket prerequisite.

**Signatures.** None public. Inside the example, the force is a lambda satisfying the existing `Force` concept:
```cpp
// Tilted axes: x up the slope, y normal to the slope, z across it.
const auto tilted_gravity = [m, g, slope](double, const Vec3&, const Vec3&) {
    return Vec3{-m * g * std::sin(slope), -m * g * std::cos(slope), 0.0};
};
// Puck constrained to the surface: the normal force cancels the y component.
const auto on_slope = [m, g, slope](double, const Vec3&, const Vec3&) {
    return Vec3{-m * g * std::sin(slope), 0.0, 0.0};
};
const double normal_force = m * g * std::cos(slope);   // [N], from a_y = 0
```
`events::ground_impact(0.0)` already means "back to the slope", because y *is* the slope-normal axis. That reuse is the point.

**Validation test** (`v0 = 20 m/s`, `g = 9.80665`, RK4 `dt = 1e-4`, `simulate_until`):

| Case | Analytic | Measured | Assert |
|---|---|---|---|
| 1.39, theta = 25 deg above a phi = 20 deg slope | R = 27.607745329 m | 27.607745329 | rel 5.6e-13 → tol 1e-9 |
| 1.39, theta = theta_opt = 35 deg (2theta+phi = 90 deg) | R_max = 30.3934696673 m | 30.3934696674 | rel 8.4e-13 → tol 1e-9 |
| degenerate phi = 0 collapse | `v0^2 sin2theta/g` = 38.3287920252 m | 38.3287920252 | rel 4.9e-13 → tol 1e-9 |
| 1.37, puck kicked up the slope | t_return = 11.9258030017 s | 11.9258030017 | rel 9.2e-13 → tol 1e-9 |
| 1.37, constraint force | N/m = g cos phi = 9.21523664 m/s^2 | algebraic | 1e-12 |

The phi = 0 collapse is the best test in the set: it checks that your tilted-frame bookkeeping reduces to the already-trusted `analytic::range`. These belong in `tests/test_numeric_vs_analytic.cpp` too, which already has the `simulate_until`-range-to-1e-9 idiom.

**Example sketch.** Three labelled blocks in the style of `examples/01_free_fall.cpp`'s `row()` helper: (1) puck up the slope — numeric vs analytic return time and the normal force, with a note that N came from `a_y = 0`, not from a force law; (2) ball thrown at 25 deg above a 20 deg slope — numeric vs analytic range; (3) a sweep of theta from 5 deg to 85 deg in 1 deg steps, printing the maximising theta against `theta_opt = (90 deg - phi)/2` and `R_max`, written to `build/incline_range_sweep.csv` via the generic `io::write_csv(path, {"theta_deg","range_m","analytic_m"}, rows)` overload. Plot with `python3 scripts/plot_csv.py build/incline_range_sweep.csv`.

**Difficulty.** S, 2–3 h.

**What you learn.** *Physics:* that a free-body diagram's axis choice is a *decision* with consequences; that a constraint force is determined by a kinematic condition rather than a force law; that a normal force exists only while there is contact. *C++:* lambdas as `Force`-concept models and why template callables (not `std::function`) are the library's contract — `notes/01-cpp-for-this-library.md:219` is the note to read; `const` lambda captures by value; unnamed unused parameters to keep `-Wunused-parameter` quiet.

---

### Proposal 4 — `newton::coords`: the polar module (§1.7). The largest single ch. 1 gap

**Pitch.** Cartesian ↔ polar, the rotating basis `r_hat`/`phi_hat` and its time derivatives, the velocity/acceleration decomposition — plus the test that a Cartesian `simulate()` run actually reproduces `a_r = r_ddot - r phi_dot^2`.

**Physics (§1.7).**
```
x = r cos(phi),  y = r sin(phi)                             (1.37)
r = sqrt(x^2 + y^2),  phi = atan2(y, x)      <-- Taylor's arctan(y/x) is
                                                 INCOMPLETE; problem 1.42
r_hat   = ( cos phi, sin phi, 0)                            (1.59)
phi_hat = (-sin phi, cos phi, 0)
d r_hat/dt   =  phi_dot phi_hat                             (1.42)
d phi_hat/dt = -phi_dot r_hat                               (1.46)
v = r_dot r_hat + r phi_dot phi_hat
a = (r_ddot - r phi_dot^2) r_hat + (r phi_ddot + 2 r_dot phi_dot) phi_hat
F_r   = m (r_ddot - r phi_dot^2)                            (1.48)
F_phi = m (r phi_ddot + 2 r_dot phi_dot)
cylindrical (problems 1.47-1.49):  rho = sqrt(x^2+y^2), phi, z
```
The whole reason the acceleration has "extra" terms is that the basis rotates with the particle — `d r_hat/dt != 0`. Those terms are **kinematics, not forces**. Do not let them become forces; that is ch. 9.

**Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/coords/polar.hpp`
- create `/Users/jackery/Github/newtonian-library/tests/test_coords.cpp`
- create `/Users/jackery/Github/newtonian-library/examples/08_polar_coordinates.cpp`
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — new block right after the core block (after line 21), since `coords` depends only on `Vec3`:
  ```cpp
  // coords — coordinate systems other than Cartesian (Taylor §1.7)
  #include <newton/coords/polar.hpp>
  ```
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — **hard rule 4 (lines 29-31) must gain `newton::coords` to its namespace list**, plus a new file-by-file section. This corrects a spec omission rather than expanding scope: `docs/research/taylor-ch1-7-map.md` §1.4 already specifies the module.
- modify `/Users/jackery/Github/newtonian-library/notes/02-physics-to-code.md` — its ch. 1 "Key equations" block stops at `F = dp/dt = m a`; §1.7 contributes nothing today.

**Signatures.**
```cpp
namespace newton::coords {

/// Plane-polar coordinates (Taylor §1.7).
/// C++ note: deliberately NOT a Vec3. Vec3's entire API — dot, cross, norm,
/// normalized — is Cartesian, so a (r, phi) pair wearing a Vec3 would compile
/// into silently meaningless calls (what is norm({r, phi, 0})?). The library's
/// precedent for a named non-Cartesian tuple is analytic::SHM and
/// systems::RocketSample. Small aggregate, zero-initialized members.
struct Polar { double r{0.0};    ///< radius [m]
               double phi{0.0}; };  ///< polar angle [rad], in [-pi, pi] — std::atan2's CLOSED range; both ±pi are reachable

/// Cartesian -> plane polar, taking x and y and ignoring z (Taylor eq. 1.37).
/// Physics: phi comes from std::atan2(y, x), not std::atan(y/x). Taylor's own
/// problem 1.42 is exactly this point: arctan(y/x) cannot tell (-1,-1) from
/// (1,1), because the ratio is the same. atan2 takes the two components
/// separately and resolves the quadrant, returning a value in the closed
/// interval [-pi, pi] (see the branch-cut test for why it is closed at BOTH ends).
[[nodiscard]] inline Polar to_polar(const Vec3& r);

/// Plane polar -> Cartesian: (r cos phi, r sin phi, 0).
[[nodiscard]] inline Vec3 from_polar(const Polar& p);

/// The rotating basis (Taylor eq. 1.59). Both are unit vectors, mutually
/// perpendicular, with cross(r_hat, phi_hat) = +z_hat — a right-handed pair.
/// @param phi polar angle [rad]
[[nodiscard]] inline Vec3 r_hat(double phi);
[[nodiscard]] inline Vec3 phi_hat(double phi);

/// Time derivatives of the rotating basis (Taylor eqs. 1.42 and 1.46):
/// d r_hat/dt = phi_dot phi_hat,  d phi_hat/dt = -phi_dot r_hat. These two
/// lines are the entire reason polar acceleration has terms Cartesian
/// acceleration does not.
[[nodiscard]] inline Vec3 r_hat_dot(double phi, double phi_dot);
[[nodiscard]] inline Vec3 phi_hat_dot(double phi, double phi_dot);

/// Radial and angular rates read off a Cartesian state: r_dot = v . r_hat,
/// phi_dot = (v . phi_hat) / r.
/// Guarded: returns {0, 0} when r <= 0, where phi and phi_dot are undefined —
/// the same guard convention as Vec3::normalized().
struct PolarVel { double r_dot{0.0};     ///< [m/s]
                  double phi_dot{0.0}; };  ///< [rad/s]
[[nodiscard]] inline PolarVel polar_velocity(const Vec3& r, const Vec3& v);

/// Polar components of acceleration — the left-hand sides of Taylor eq. 1.48
/// divided by m:  a_r = r_ddot - r phi_dot^2,  a_phi = r phi_ddot + 2 r_dot phi_dot.
/// Physics: -r phi_dot^2 is the centripetal term and 2 r_dot phi_dot the
/// Coriolis-like term. Both are KINEMATIC consequences of the rotating basis,
/// not forces. Nothing in this library ever adds them to a force.
struct PolarAccel { double a_r{0.0};      ///< [m/s^2]
                    double a_phi{0.0}; }; ///< [m/s^2]
[[nodiscard]] constexpr PolarAccel polar_acceleration(double r, double r_dot, double r_ddot,
                                                      double phi_dot, double phi_ddot);

/// Rebuild a Cartesian vector from polar components at angle phi:
/// w = w_r r_hat(phi) + w_phi phi_hat(phi). The inverse of projecting with dot().
[[nodiscard]] inline Vec3 from_polar_components(double phi, double w_r, double w_phi);

/// Cylindrical polars (Taylor problems 1.47-1.49). Note the notational hazard
/// Taylor flags: many books write r for rho; this library never does.
struct Cylindrical { double rho{0.0};  ///< [m]
                     double phi{0.0};  ///< [rad]
                     double z{0.0}; };  ///< [m]
[[nodiscard]] inline Cylindrical to_cylindrical(const Vec3& r);
[[nodiscard]] inline Vec3 from_cylindrical(const Cylindrical& c);

}  // namespace newton::coords
```
`polar_acceleration` is the one `constexpr` function here; everything else calls `std::sin`/`std::cos`/`std::atan2`/`std::hypot`, none of which are `constexpr` in C++17.

**Validation test** (`tests/test_coords.cpp`):

1. **Round trip.** `to_polar(from_polar({r, phi}))` for r in {0.5, 1, 1000}, phi swept over the *open* interval (-pi, pi) → 1e-12 relative. Compare angles by wrapped difference, or keep ±pi out of the sweep: at the seam a raw comparison is meaningless. Then the branch cut, asserted where it actually lives — **on the sign of zero**: `NT_CHECK(std::atan2(-0.0, -1.0) == -constants::pi)` and `NT_CHECK(std::atan2(+0.0, -1.0) == +constants::pi)` (exact `==` is legitimate here: both sides are the correctly rounded double for ∓pi). And a trap worth its own comment: `to_polar(from_polar({1, -pi})).phi` comes back as **-pi**, not +pi — `std::sin(-pi_double)` is `-1.22e-16`, a *negative* y, because `pi_double` sits just below the true pi, so `atan2` lands on the -pi side. Both ±pi are fixed points to roundoff; the mathematical branch cut is only visible with a signed zero. Assert `NT_CLOSE(to_polar(from_polar({1, -pi})).phi, -constants::pi, 1e-12)` so nobody later "fixes" it to +pi. (An earlier draft of this plan asserted `== +pi`; the verifier compiled it on this toolchain and it fails — measured `-3.1415926535897931`.)
2. **Problem 1.42, the quadrant bug.** `to_polar({-1,-1,0}).phi == -3*pi/4` to 1e-15, *and* `NT_CHECK(std::atan(-1.0/-1.0) != to_polar({-1,-1,0}).phi)` — the naive formula gives +pi/4, off by pi. This is the memorable test.
3. **Basis geometry.** `dot(r_hat(phi), phi_hat(phi)) == 0`, both norms 1, `cross(r_hat, phi_hat) == {0,0,1}` — swept over phi, tol 1e-15.
4. **Basis derivatives vs finite difference** (eqs. 1.42, 1.46), h = 1e-5, phi_dot = 2.7, phi swept over [-3, 3]: **tol 1e-8**; measured 3.3e-10 for both (at h = 1e-6 it improves to 6.6e-11, so the tolerance has headroom either way).
5. **Uniform circular motion, analytic, no integrator.** R = 2 m, w = 3 rad/s: feed `polar_acceleration(R, 0, 0, w, 0)` and compare to `dot(-w^2 r, r_hat)` and `dot(-w^2 r, phi_hat)`. Tol 1e-12; measured `a_r` err **3.6e-15**, `a_phi` err **1.7e-15**, `phi_dot` and velocity reconstruction **2.2e-15**.
6. **The headline cross-check: a Cartesian simulation must reproduce the polar decomposition.** Fixture is ch.-1-only — a vacuum projectile under `forces::UniformGravity`, v0 = 50 m/s at 70 deg, launched from the origin, polar origin at the launch point, RK4, `dt = 1e-4`, `record_every = 1`, t_end = 9 s. For each sample: compute `r`, `phi`, and exactly `r_dot`, `phi_dot` via `polar_velocity`; get `r_ddot` and `phi_ddot` by **central finite difference of those two sequences**; assert `polar_acceleration(...)` matches `dot(g, r_hat)` and `dot(g, phi_hat)`. Skip samples with `r < 1 m` (near the polar origin the angular derivatives blow up, and phi is meaningless at r = 0). **Tolerance 1e-7 absolute in m/s^2**; measured over 89 799 usable samples: worst `|a_r|` error **3.1e-9**, worst `|a_phi|` error **8.1e-9** (3.1e-10 relative to g). The tolerance is tied to `dt` and you must say so in the test: the same run at `dt = 1e-3` gives **3.1e-7** — exactly 100x worse, the O(dt^2) signature of a central difference. Do not copy this test and change `dt` without changing the tolerance.
7. **Velocity reconstruction along that same trajectory is FD-free**, so it holds to roundoff: `from_polar_components(phi, r_dot, r*phi_dot) == s.v` to 1e-12 relative.
8. **Cylindrical** round trip + `to_cylindrical({x,y,z}).rho == to_polar({x,y,0}).r` to 1e-15.

**Example sketch** (`examples/08_polar_coordinates.cpp`):
- *Block 1:* the projectile of test 6. A 10-row table of `t, r, phi(deg), r_dot, phi_dot, a_r, a_phi, resid_r, resid_phi`, plus the max residuals over the whole run with the O(dt^2) note beside them. Writes `build/polar_projectile.csv` via the generic `io::write_csv` overload (the `Trajectory` overload's fixed 7-column header cannot carry polar columns).
- *Block 2:* **problem 1.40** as a live demo. Sweep the launch angle and find the largest theta for which `d|r|/dt` (use `dot(r, v)`, same sign, no division) never dips negative during the flight. Bisect on theta to 40 iterations. Measured with `UniformGravity`, v0 = 20 m/s, RK4 `dt = 1e-5`: **theta_max = 70.5287794 deg**, against `arcsin(sqrt(8/9))` = 70.5287794 deg — **7.0e-13 relative**. The most satisfying number in the chapter: a non-obvious closed form your simulator finds by search.
- *Block 3:* `r_hat`/`phi_hat` printed at four angles, with `d r_hat/dt` next to `phi_dot phi_hat`, so the reader sees eq. 1.42 as numbers.

**Difficulty.** M, 5–8 h.

**What you learn.** *Physics:* why a rotating basis produces acceleration terms that are *not* forces — the single most common conceptual error in the rest of the course; why Cartesian is the right internal representation and polar the right *reporting* representation; that `atan2` is a physics bug fix, not a C detail. *C++:* designing a small aggregate struct instead of overloading a vector type with a second meaning (and *why* — `norm()` of a polar pair is nonsense); `constexpr` where arithmetic allows and `inline` where `<cmath>` forbids; zero-guarding a function whose mathematical domain excludes a point, the way `normalized()` does; `std::hypot` over `sqrt(x*x+y*y)` for overflow safety; the harness gap — `NT_CLOSE_VEC3` is templated on `.x/.y/.z` (`tests/test_util.hpp:51`), so `PolarAccel` needs two `NT_CLOSE` calls, and bind an `expected` local first because the 3-argument macro splits on top-level commas.

---

### Proposal 5 — `analytic::UniformCircle`: the circular-motion oracle ch. 1 needs

**Pitch.** A ch.-1-framed closed form for uniform circular motion, so `a = v^2/R`, `T = 2 pi R/v` and problem 1.41's `T_string = m R w^2` stop being reachable only through the §2.5 cyclotron functions.

**Physics (§1.7, problems 1.10, 1.41, 1.45).**
```
r(t)  = R (cos(w t + phi0), sin(w t + phi0), 0)
v(t)  = R w phi_hat,          |v| = R |w|               constant speed
a(t)  = -R w^2 r_hat,         |a| = R w^2 = |v|^2 / R   toward the centre
T     = 2 pi / |w| = 2 pi R / |v|
Whatever force holds the circle has magnitude m R w^2 = m v^2 / R; for a mass
on a string in gravity-free space that force IS the tension (problem 1.41),
and the phi-component of Newton II (1.48) is identically satisfied.
Problem 1.45: |r| is constant  <=>  r . v = 0 at all times.
```

**Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/analytic/circular.hpp`
- create `/Users/jackery/Github/newtonian-library/examples/09_circular_motion.cpp`
- modify `/Users/jackery/Github/newtonian-library/tests/test_numeric_vs_analytic.cpp` (new ch.-1 block) and `/Users/jackery/Github/newtonian-library/tests/test_analytic.cpp` (internal consistency)
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — add **`#include <newton/analytic/circular.hpp>`** to the analytic block, after `kinematics.hpp`
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — bullet in the `newton::analytic` section

**Signatures.**
```cpp
namespace newton::analytic {

/// Uniform circular motion in the xy-plane (Taylor §1.7, problems 1.10, 1.41).
/// Physics: the speed is constant but the VELOCITY is not — the acceleration
/// is nonzero, constant in magnitude, and always directed at the centre.
/// Nothing here is oscillator theory: this is prescribed kinematics plus the
/// one force magnitude Newton II demands of it.
struct UniformCircle {
    double radius{1.0};  ///< R [m], must be > 0
    double omega{1.0};   ///< phi_dot [rad/s]; negative = clockwise
    double phi0{0.0};    ///< phase at t = 0 [rad]

    [[nodiscard]] double phi(double t) const;             ///< omega*t + phi0 [rad]
    [[nodiscard]] Vec3 position(double t) const;          ///< [m]
    [[nodiscard]] Vec3 velocity(double t) const;          ///< [m/s], = R*omega*phi_hat
    [[nodiscard]] Vec3 acceleration(double t) const;      ///< [m/s^2], = -R*omega^2*r_hat
    [[nodiscard]] double speed() const;                   ///< R*|omega| [m/s]
    [[nodiscard]] double centripetal_accel() const;       ///< R*omega^2 = speed^2/R [m/s^2]
    [[nodiscard]] double period() const;                  ///< 2*pi/|omega| [s]; guarded at omega == 0
    /// Magnitude of the central force required: m R omega^2 = m v^2 / R.
    /// For problem 1.41's mass on a string this is the tension [N].
    [[nodiscard]] double central_force(double mass) const;
};

/// a = v^2 / r, for when you have speed and radius rather than omega
/// (Taylor §1.7). Guarded: returns 0 for radius <= 0.
[[nodiscard]] constexpr double centripetal_accel_from_speed(double speed, double radius);

/// T = 2 pi r / v. Guarded: returns 0 for speed <= 0.
[[nodiscard]] inline double circular_period_from_speed(double speed, double radius);

}  // namespace newton::analytic
```
Optional add-on, 20 minutes: problems 1.11/1.12's prescribed ellipse and elliptical helix, `struct EllipticalPath { double b, c, omega, vz; ... }`, which exercise `coords::to_polar` on a non-circular path and satisfy `x^2/b^2 + y^2/c^2 = 1` to roundoff.

**Validation test.** Internal consistency (`test_analytic.cpp`, no integrator): `centripetal_accel() == speed()^2/radius` and `== centripetal_accel_from_speed(speed(), radius)`; `period() == circular_period_from_speed(speed(), radius)`; `dot(position(t), velocity(t)) == 0` and `norm(position(t)) == radius` over a swept t — all 1e-12.

Numeric-vs-analytic (`test_numeric_vs_analytic.cpp`), R = 3 m, w = 2 rad/s, m = 1.5 kg. The force is a **string**: constant magnitude `T = m R w^2 = 18 N` directed at the centre, `F = -T * normalized(r)` — problem 1.41's actual force, not a spring. One full period, RK4, `dt = 1e-4`:

| Assertion | Tol | Measured |
|---|---|---|
| `norm(s.r) == R` every sample (the constraint holds) | 1e-11 | 3.7e-14 |
| `norm(s.v) == analytic.speed()` every sample | 1e-11 | 5.0e-14 |
| `dot(s.r, s.v) == 0` every sample (**problem 1.45**) | 1e-11 | 2.0e-13 |
| `norm(a) == norm(v)^2/norm(r)` every sample | 1e-11 | 2.3e-13 |
| period: `2 * simulate_until(..., events::turning_point_x()).back().t` vs `2 pi/w` = 3.14159265358979 | 1e-12 rel | **3.0e-15** rel |
| tension magnitude `== analytic.central_force(m)` = 18 N | 1e-12 | algebraic |

Note the period trick, which needs Proposal 6's tiny event addition: starting at `(R,0,0)` with `v = (0,Rw,0)`, `v.x` starts at 0 and next vanishes a **half** period later, so the measured full period is `2 x t_event`. (This is a single-leg measurement, so the chaining discipline of Proposal 6 does not apply here.) With a linear central force `-m w^2 r` instead of the string, the same run gives 4.5e-14 / 9.2e-14 / 4.5e-13 / 5.4e-13 — indistinguishably good, because for exact circular initial conditions the two force laws agree. Use the string; it is the one ch. 1 poses.

**Example sketch** (`examples/09_circular_motion.cpp`): problem 1.41's astronaut. Prints R, w, speed, centripetal acceleration, period, and tension, then the numeric/analytic pair for each, plus the two running invariants (`max ||r|-R|` and `max |r.v|`) over one period with a line explaining that `r.v = 0` *is* problem 1.45's theorem. Writes `build/circle.csv` (the `Trajectory` overload works — ordinary Cartesian motion) and a second generic CSV of `t, r, phi_deg, r_dot, phi_dot` from `coords::polar_velocity`, demonstrating that `r_dot = 0` and `phi_dot = w` to roundoff — which is what makes the polar decomposition collapse to `a_r = -R w^2`.

**Difficulty.** S–M, 3–4 h (depends on Proposal 4 for the polar read-out and Proposal 6's event for the period).

**What you learn.** *Physics:* that constant *speed* is not constant *velocity*; that the force holding a circle is whatever the constraint supplies, with magnitude fixed by kinematics; problem 1.45's `|r|` const <=> `r` perpendicular to `v` as a runnable invariant rather than an exercise. *C++:* a struct of pure accessors as an analytic oracle (the `analytic::SHM` pattern without any of SHM's ch.-5 content); guarding a divide inside a `const` member function; `[[nodiscard]]` on accessors and why discarding `period()` is always a bug.

---

### Proposal 6 — The half-pipe: `systems::CircularTrack` + `events::turning_point_x`. Ch. 1's only `[Computer]` problems

**Pitch.** Integrate `phi_ddot = -(g/R) sin(phi)` from rest at 20 deg and at 90 deg, plot against the small-angle curve, and read the normal force off the radial equation — Taylor problems 1.50 and 1.51, exactly as posed. This is the headline deliverable of the whole plan.

**Physics (§1.7, Example 1.2).** A point mass on a rigid frictionless circular track of radius R, phi measured from the lowest point. Polar coordinates centred on the centre of curvature: `r = R` fixed, so `r_dot = r_ddot = 0`, and the two components of eq. (1.48) become
```
tangential:   m (R phi_ddot) = -m g sin(phi)
              multiply by R:   (m R^2) phi_ddot = -(m g R) sin(phi)
              =>  phi_ddot = -(g/R) sin(phi)                       (1.51)

radial:       m (0 - R phi_dot^2) = m g cos(phi) - N
              =>  N = m ( g cos(phi) + R phi_dot^2 )               the constraint force
                  at a turning point (phi_dot = 0):  N = m g cos(phi_0)

small-angle limit only:  phi_ddot = -w^2 phi,  w = sqrt(g/R)       (1.55)
                         phi(t)   = phi_0 cos(sqrt(g/R) t)         (1.57)
                         T_lin    = 2 pi sqrt(R/g)                 amplitude-independent
```
Taylor says outright that (1.51) has no elementary closed-form solution. That sentence is the book's first argument for numerical methods, and it is the reason this proposal exists.

**Scope discipline.** (1.55)–(1.57) are ch. 1 material (§1.7, problem 1.44), but they are the small-angle limit *of this one problem* — not oscillator theory, which is ch. 5. The doc comment must say so, and nothing here may use energy: the integrator-quality metric below is **turning-point drift**, which ch. 1 licenses, not energy drift, which it does not.

**Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/systems/circular_track.hpp`
- create `/Users/jackery/Github/newtonian-library/examples/10_half_pipe.cpp`
- create `/Users/jackery/Github/newtonian-library/tests/test_circular_track.cpp`
- modify `/Users/jackery/Github/newtonian-library/include/newton/simulate.hpp` — one new event in the existing `events` namespace (beside `ground_impact` and `apex`, lines 222–231)
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — add **`#include <newton/systems/circular_track.hpp>`** to the systems block (beside line 48's `rocket.hpp`)
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — a `newton::systems` file-by-file entry and the new `events` function in the `simulate.hpp` spec block (lines 222–226)

No new namespace: `newton::systems` is already on rule 4's list, and "self-contained special-purpose simulator" is exactly what this is.

**Signatures.**
```cpp
namespace newton {
namespace events {
/// Event function for simulate_until(): the x-component of velocity. Its sign
/// change marks a turning point of x — and, for a one-degree-of-freedom
/// problem packed into the x slot (see systems/circular_track.hpp), a turning
/// point of the generalized coordinate. apex() is the same idea on y.
[[nodiscard]] inline auto turning_point_x() {
    return [](const State& s) { return s.v.x; };
}
}  // namespace events
}  // namespace newton

namespace newton::systems {

/// Taylor's Example 1.2 — a point mass on a rigid, frictionless circular
/// track of radius R in a uniform field, reduced to ONE degree of freedom.
///
/// Physics: with r = R pinned, eq. (1.48)'s tangential component is the whole
/// equation of motion, phi_ddot = -(g/R) sin(phi) — nonlinear, with no
/// elementary closed-form solution (Taylor eq. 1.51). The radial component is
/// not thrown away: it DETERMINES the constraint force, N = m(g cos phi +
/// R phi_dot^2). phi is measured from the lowest point.
///
/// C++ note: this struct satisfies the Force concept with a DELIBERATE
/// reinterpretation of every slot, and that reinterpretation is the one thing
/// to understand before using it. The integrator does not care what a State's
/// numbers mean; it only ever computes a = F/mass and advances (r, v). So we
/// store phi [rad] in r.x, phi_dot [rad/s] in v.x, the "generalized force"
/// (-m g R sin phi [N m]; in ch. 3 you will learn to call this a torque) in
/// the returned x component, and the "generalized mass" (m R^2 [kg m^2]) in
/// the integrator's `mass` slot. Those two names are borrowed from ch. 7 and
/// you do not need them yet: "generalized mass" just means the coefficient
/// multiplying phi_ddot, exactly as m multiplies x_ddot, and "generalized
/// force" means the right-hand side of the phi equation. The idea is pure
/// algebra — multiply the tangential equation by R. Then
/// F/mass is exactly phi_ddot, and simulate(), simulate_until(), the event
/// bisection and io::write_csv all work unmodified. The hazard is real: the
/// y and z slots are meaningless, r.x is radians and not metres, and nothing
/// in the type system says so. Always go through pack()/phi()/phi_dot()
/// rather than touching .r.x by hand, and never hand one of these States to
/// analysis/ or coords/, which assume Cartesian metres.
struct CircularTrack {
    double R{1.0};                        ///< track radius [m]
    double g{constants::g_earth};         ///< field strength [m/s^2]
    double mass{1.0};                     ///< particle mass [kg]

    /// The coefficient of phi_ddot: m R^2 [kg m^2]. Pass as Particle::mass.
    [[nodiscard]] constexpr double generalized_mass() const { return mass * R * R; }

    /// Generalized force -m g R sin(phi) [N m] — the right-hand side of the
    /// phi equation — in the x slot.
    [[nodiscard]] Vec3 operator()(double /*t*/, const Vec3& q, const Vec3& /*qdot*/) const;

    /// Magnitude of the track's normal force, N = m(g cos phi + R phi_dot^2) [N],
    /// from the RADIAL component of eq. (1.48). Always non-negative for a
    /// track that can only push: if this ever came out negative the particle
    /// would have left the surface, which this model cannot represent.
    [[nodiscard]] double normal_force(double phi, double phi_dot) const;

    /// Slot packing/unpacking — the only sanctioned way in or out.
    [[nodiscard]] constexpr State pack(double phi, double phi_dot, double t = 0.0) const;
    [[nodiscard]] constexpr Particle particle(double phi, double phi_dot) const;  // uses generalized_mass()
    [[nodiscard]] static constexpr double phi(const State& s)     { return s.r.x; }
    [[nodiscard]] static constexpr double phi_dot(const State& s) { return s.v.x; }

    /// Where the mass actually is, relative to the centre of curvature:
    /// (R sin phi, -R cos phi, 0) [m]. Use this for plotting and for handing
    /// a genuine Cartesian position to coords/ or io/.
    [[nodiscard]] Vec3 cartesian(double phi) const;

    /// The SMALL-ANGLE reference curve, Taylor eq. (1.57):
    /// phi_0 cos(sqrt(g/R) t). A LIMIT of this problem, not a prediction of
    /// it — valid only as phi_0 -> 0, and the whole point of problems 1.50
    /// and 1.51 is to watch it fail. This is not general oscillator theory
    /// (Taylor ch. 5), and nothing in this library treats it as such.
    [[nodiscard]] double small_angle_phi(double phi0, double t) const;
    [[nodiscard]] double small_angle_period() const;   ///< 2 pi sqrt(R/g) [s]
};

}  // namespace newton::systems
```

#### The method for every turning-point quantity — and the trap in it

Taylor's numbers: `R = 5 m`, `g = 9.8 m/s^2`, so `w = sqrt(9.8/5) = 1.4 rad/s` exactly and `T_lin = 2 pi/1.4 = 4.48798951 s`; `phi_ddot = -1.96 sin phi`.

`simulate_until` refines a crossing by bisecting on the **step size**, re-stepping from the previous state 40 times
(`include/newton/simulate.hpp`), so the returned state sits on `phi_dot = 0` to the integrator's own accuracy rather
than to O(dt). It returns on the **first** crossing, so each run yields exactly **one** turning point. To get a series you
chain legs — and **chaining naively does not work.** Here is the correct loop:

```cpp
const systems::CircularTrack tk{5.0, 9.8, 1.0};
using Tk = systems::CircularTrack;

Particle q = tk.particle(phi0, 0.0);
double t0 = 0.0;
double phi_prev = phi0;
const double dt = 1e-3;

for (int leg = 0; leg < 12; ++leg) {   // 12 half-periods = 6 FULL periods
    SimOptions o;  o.t0 = t0;  o.t_end = t0 + 6.0;  o.dt = dt;  o.method = Method::RK4;
    const Trajectory tr = simulate_until(q, tk, o, events::turning_point_x());
    const State tp = tr.back();

    // Guard A: the event really fired. A leg whose t_end was too short would
    // silently return the t_end state instead.
    NT_CHECK(std::fabs(Tk::phi_dot(tp)) < 1e-10);
    // Guard B: this is a NEW turning point, not a duplicate of the last one.
    // Guard A CANNOT see a duplicate — a stalled leg also has phi_dot ~ 0.
    NT_CHECK(tp.t - t0 > 0.5 * dt);
    // Guard C: turning points must alternate in sign.
    NT_CHECK(leg == 0 || Tk::phi(tp) * phi_prev < 0.0);

    /* ... record phi and N at tp ... */
    phi_prev = Tk::phi(tp);

    // THE NUDGE. Advance one full step past the crossing before seeding the
    // next leg, so that leg's e_prev is unambiguously on the far side of zero.
    SimOptions one;  one.t0 = tp.t;  one.t_end = tp.t + dt;
    one.dt = dt;     one.method = o.method;
    const State s1 = simulate(tk.particle(Tk::phi(tp), Tk::phi_dot(tp)), tk, one)
                         .states.back();
    q  = tk.particle(Tk::phi(s1), Tk::phi_dot(s1));
    t0 = s1.t;
}
```

**Why the nudge is required — read this before you trust any turning-point number.** Only **leg 0** starts with
`phi_dot = 0` exactly, because you pass the literal `0.0`. For every later leg the seed `phi_dot` is the *bisection
residual*: the 40-iteration loop returns its last midpoint, which lands on either side of the root, so the seed is roughly
±1e-16 to ±1e-15 rather than 0. `simulate_until` deliberately ignores an **exact** initial zero, but a residual is not
exactly zero. Whenever it has the pre-crossing sign, the first step of the next leg still makes the sign test fire and the
endpoint-rounding rule returns a **duplicate of the previous turning point** (now exactly at `t0`). Guard A passes on
these stalled legs (`phi_dot` is ~0 there too), so without guard B the failure is completely silent.

The exact number of stalled legs is sensitive to floating-point details and the endpoint-rounding rule, so do not use an
old count as a correctness criterion. The invariant is the three guards: without the nudge, Guard B can detect a duplicate;
with it, the first state of each new leg is unambiguously on the far side of zero and the 12 guards must show twelve
strictly alternating turning points.

The nudge costs one extra step per leg and starts each leg `dt` past the turning point. The recorded turning points are the refined `tp` states, so nothing you measure is shifted.

**Validation test.** Measured with the nudged loop, `dt = 1e-3`, RK4, 12 legs = **6 full periods**:

| Quantity | phi0 = 20 deg | phi0 = 90 deg | Assert |
|---|---|---|---|
| legs covering a genuine 6 periods | 12 legs, 0 stalls, final t = 27.1344 s | 12 legs, 0 stalls, final t = 31.7841 s | guards A/B/C above |
| max \| \|phi_turning\| - phi0 \| over 6 periods | **9.4e-15 rad** | **2.2e-14 rad** | 1e-9 |
| max \| N/m - g cos phi0 \| at turning points, 6 periods | **3.1e-14** | **2.2e-13** | 1e-9 |
| reference `g cos phi0` | 9.20898768 m/s^2 | 6.0e-16, i.e. **0 to 1e-9** (not "exactly 0") | — |
| measured full period (2 x first leg) | **4.52240801 s** | **5.29735622 s** | see below |
| ratio to `T_lin` | 1.00766903 | 1.18034060 | see below |

**Pin `dt <= 1e-3` in the test and say why.** Same nudged method at `dt = 1e-2`: amplitude drift 4.8e-11 / 1.9e-10 and normal-force error 1.6e-10 / **1.9e-9** — the 90 deg case would fail a 1e-9 assertion. At `dt = 1e-4` nothing improves (7.0e-15 / 1.7e-13 drift; 2.3e-14 / 1.7e-12 normal force): `dt = 1e-3` is already at the roundoff floor.

**Do not scan recorded samples for a `phi_dot` sign change as a route to any of these assertions.** Measured over the identical 6-period run, taking the nearest recorded sample to each of the 12 crossings: amplitude error up to **6.7e-8** (20 deg) and **2.1e-7** (90 deg); normal-force error up to **6.8e-7** and **6.1e-6** — three to four orders outside 1e-9, because the nearest sample misses the turning point by O(phi_ddot dt^2/8) in phi and O(phi_ddot dt/2) in phi_dot, and `N = m(g cos phi + R phi_dot^2)` is first-order in that phi_dot error. If you ever need a turning-point series without chaining legs, *scan then refine*: re-step from the sample **before** the crossing with a bisected sub-step. If you keep a raw scan, loosen the assertion to ~1e-5 and write the O(dt^2) quantization floor into the comment.

**Period assertions, ch.-1-safe.** Do not assert an amplitude-correction formula — that is not in ch. 1. Assert instead: (i) the measured period **exceeds** `T_lin` in both cases, which follows from ch. 1 alone since `sin phi < phi` makes the true restoring torque weaker than the linear one; (ii) at 20 deg the measured period is within 1% of `T_lin` (measured +0.767%); (iii) at 90 deg it is more than 15% above `T_lin` (measured +18.03%). That trio captures "the small-angle limit is good at 20 deg and bad at 90 deg" without importing ch. 5.

**Small-angle curve comparison — assert the failure as well as the success.** Over 3 `T_lin` at `dt = 1e-3`, max `|phi_num(t) - phi0 cos(w t)|`:
- phi0 = 20 deg (= 0.34907 rad): **0.046138 rad (2.644 deg)**, i.e. 13% of the amplitude → assert `< 0.05`;
- phi0 = 90 deg (= 1.57080 rad): **3.1125 rad (178.3 deg)** — the two curves end up essentially in antiphase → assert `> 1.0`. Problem 1.51's expected answer *is* the divergence, so the test must demand it, not tolerate it.

**Integrator quality, energy-free.** Same refined-turning-point drift metric, nudged loop, `dt = 1e-3`, over a genuine **6 full periods**:

| Method | drift at 20 deg | drift at 90 deg |
|---|---|---|
| `Method::Euler` | **9.3e-3 rad** (2.7% of phi0) | **4.2e-2 rad** |
| `Method::SemiImplicitEuler` | 6.5e-7 rad | 2.6e-6 rad |
| `Method::VelocityVerlet` | 6.6e-14 rad | 2.4e-14 rad |
| `Method::RK4` | **9.4e-15 rad** | **2.2e-14 rad** |

Assert Euler drifts more than 1e-3 rad and RK4 less than 1e-9 rad. This is the ch.-1-licensed replacement for `examples/06_integrator_comparison.cpp`'s energy-drift chart: the turning points of a periodic motion must not wander, and you do not need the word "energy" to say so.

**Example sketch** (`examples/10_half_pipe.cpp`), reproducing both `[Computer]` problems:
- Header comment: the physical situation, §1.7 / Example 1.2, the analytic expectation (none — that is the point), and what to look at.
- Problem 1.50: phi0 = 20 deg, RK4, `dt = 1e-3`, `t_end = 14 s` (just over three periods, matching the window the problem asks for). Problem 1.51: phi0 = pi/2, `t_end = 21 s`.
- Per amplitude, print: measured period, `T_lin = 4.48799 s`, the ratio, the max deviation from `phi0 cos(w t)`, the turning-point drift for Euler / Euler–Cromer / RK4 **with the number of periods the drift was measured over printed next to it**, and `N/m` at a turning point beside `g cos phi0`.
- Print `N/m` at the lowest point too — measured **10.9820246 m/s^2 at phi0 = 20 deg** and **29.4000000 m/s^2 at 90 deg** — labelled honestly: *you cannot predict these two numbers yet. Ch. 4's energy methods will let you; write them down now and check them then.*
- `build/half_pipe.csv` via the generic overload, columns `t, phi20, small20, phi90, small90` on a common 0–21 s grid, so one `python3 scripts/plot_csv.py build/half_pipe.csv t phi20 small20` gives problem 1.50's figure and swapping the columns gives 1.51's.

**Difficulty.** M, 6–9 h (depends on Proposal 4 for the polar reasoning, though not on its code).

**What you learn.** *Physics:* that a constraint turns a 2-D vector problem into one scalar ODE; that the "unused" radial equation is where the constraint force comes from; that a nonlinear equation can be perfectly well-behaved and still have no closed form; exactly how far a linearisation can be trusted, measured rather than asserted. *C++:* reusing a general-purpose integrator for a generalized coordinate by reinterpreting a struct's slots — and the discipline that makes such a trick safe rather than a landmine (named `pack`/`phi`/`phi_dot` accessors, a doc comment that says what the slots mean, a rule that these States never reach Cartesian code); writing a lambda-returning factory in the established `events` idiom; **chaining stateful calls where each output seeds the next input, and why a seed that is "almost" a boundary value is a different input from one that is exactly the boundary value**; writing post-conditions that can actually detect the failure mode you have (guard A could not see the stall; guard B can).

---

### Proposal 7 — `newton::frames`: Galilean boosts and the turntable as *observation* (§1.4)

**Pitch.** Re-express a finished trajectory in another frame, so problems 1.26, 1.27 and 1.46 become three lines of post-processing — the one ch. 1 concept missing from both the code and the roadmap.

**Physics (§1.4, §1.5).**
```
Newton I defines the inertial frames: F = 0 => v = const.
A frame moving at constant V relative to an inertial frame is inertial:
   r' = r - V t,   v' = v - V,   a' = a      =>  F = m a is Galilean invariant
A uniformly accelerating frame is NOT inertial:
   r'' = r - (1/2) A t^2,   v'' = v - A t,   a'' = a - A
A rotating frame is NOT inertial. What the turntable observer SEES:
   r' = Rot(-w t) r,    v' = Rot(-w t) (v - w z_hat x r)
   |r'| = |r|, and a puck crossing at constant v0 appears to move at
   sqrt(v0^2 + rho^2 w^2) at radius rho.
```
Chapter 1 never adds a centrifugal or Coriolis **force**; it only asks what the moving observer measures. Fictitious forces are ch. 9. This module transforms *recorded output* and must never touch an equation of motion.

**Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/frames/galilean.hpp`
- create `/Users/jackery/Github/newtonian-library/tests/test_frames.cpp`
- create `/Users/jackery/Github/newtonian-library/examples/11_frames.cpp`
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — new block after the coords block:
  ```cpp
  // frames — re-expressing a recorded trajectory in another frame (Taylor §1.4)
  #include <newton/frames/galilean.hpp>
  ```
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — **hard rule 4's namespace list gains `newton::frames`** (bundle this with Proposal 4's `newton::coords` amendment if you do them together), plus a file-by-file section

**Signatures.**
```cpp
namespace newton::frames {

/// A second frame S' with axes parallel to the simulation frame's, whose
/// origin moves as  O(t) = origin0 + velocity*t + ½*accel*t²  (Taylor §1.4,
/// problem 1.26). With accel == 0 this is a Galilean boost and S' is inertial;
/// with accel != 0 it is not, and a free particle will appear to accelerate.
struct LinearFrame {
    Vec3 origin0{};   ///< S' origin at t = 0 [m]
    Vec3 velocity{};  ///< V, constant [m/s]
    Vec3 accel{};     ///< A [m/s^2]; nonzero makes S' NON-inertial
    /// True iff this frame is inertial given that the simulation frame was
    /// (Taylor §1.4: unaccelerated and non-rotating relative to an inertial
    /// frame => inertial).
    [[nodiscard]] constexpr bool is_inertial() const { return accel == Vec3{}; }
};

/// Re-express one recorded state as S' measures it:
/// r' = r - origin(t),  v' = v - V - A t,  t' = t.
///
/// Physics: this transforms an ALREADY-COMPUTED trajectory. It does not
/// change any equation of motion and adds no fictitious forces — Taylor's
/// ch. 1 asks only what a moving observer sees; centrifugal and Coriolis
/// FORCES are ch. 9.
///
/// Caveat worth knowing before you rely on a boost: "same physics in the
/// boosted frame" holds only when the force law is itself Galilean invariant.
/// forces::UniformGravity, forces::ConstantForce and any F(r) are. A
/// velocity-dependent force is not, because it is defined relative to
/// something physical (a medium, a field) rather than to your choice of
/// frame — boosting such a simulation changes the physics unless you boost
/// that something too.
[[nodiscard]] inline State observe_from(const State& s, const LinearFrame& f);
[[nodiscard]] inline Trajectory observe_from(const Trajectory& tr, const LinearFrame& f);

/// A frame rotating at constant rate about +z through the origin — Taylor's
/// turntable (problems 1.27, 1.46). NOT inertial at any nonzero omega.
struct RotatingFrame {
    double omega{0.0};  ///< [rad/s], counter-clockwise about +z
    [[nodiscard]] constexpr bool is_inertial() const { return omega == 0.0; }
};

/// What the turntable rider measures: r' = Rot(-omega t) r, and
/// v' = Rot(-omega t) (v - omega z_hat x r). In polar terms r' = r,
/// phi' = phi - omega t — the straight line across the turntable becomes a
/// spiral in and back out.
[[nodiscard]] inline State observe_from(const State& s, const RotatingFrame& f);
[[nodiscard]] inline Trajectory observe_from(const Trajectory& tr, const RotatingFrame& f);

}  // namespace newton::frames
```
Overloading `observe_from` on the frame type (rather than `boost`/`rotate_into`) is the shape that scales: a third frame type later adds two overloads and no new vocabulary.

**Validation test** (`tests/test_frames.cpp`):

1. **Galilean invariance, the headline.** Simulate `forces::UniformGravity{m}` from `(r0, v0) = ((0,10,0), (5,2,0))`, m = 2 kg, RK4 `dt = 1e-3`, 1.5 s. Simulate again from `(r0, v0 - V)` with `V = (3,4,-1)`. Then `observe_from(traj1, LinearFrame{{}, V, {}})` must equal `traj2` sample by sample. Tol **1e-11**; measured worst `|dr| + |dv|` = **4.0e-13**. This works precisely because `UniformGravity::operator()` ignores `v` — say so in the test comment.
2. **The free-fall frame (problem 1.26's third frame).** `observe_from(traj1, LinearFrame{{}, {}, g})` must turn the parabola into the straight line `r0 + v0 t`. Tol 1e-11; measured worst deviation **6.0e-13**. A parabola becoming a line is the most vivid possible statement of "this frame is not inertial".
3. **The turntable (problems 1.27, 1.46).** Puck crossing the centre at constant velocity, `v0 = 2 m/s` along y_hat, `w = 1.3 rad/s`, t in [-2, 2]: a rotation preserves length, so `norm(observed.r) == norm(s.r)` to **1e-14** (measured **4.4e-16**), and the observed speed equals `sqrt(v0^2 + rho^2 w^2)` with `rho = |r|` to **1e-12** (measured **1.8e-15**) — the analytic result the published solutions give for 1.27. Also assert the observed path is a spiral in and back out: `norm(observed.r)` decreases monotonically to 0 then increases.
4. **`is_inertial()`**: true for `{{}, V, {}}`, false for any nonzero `accel`, false for any nonzero `omega`.

**Example sketch** (`examples/11_frames.cpp`): one kicked puck, three descriptions (problem 1.26). A short table of `t, x, y` in the original frame, the boosted frame, and the uniformly-accelerating frame, with a one-line verdict under each ("inertial: straight line", "inertial: straight line, different slope", "non-inertial: parabola, yet no force acts"). Then the turntable (1.46): `rho, phi_ground_deg, phi_rotating_deg, observed_speed, sqrt(v0^2+rho^2 w^2)` at a few times, written to `build/frames.csv` with a `series` tag column (the idiom `examples/02_projectile_drag.cpp` already uses) so all four paths plot together.

**Difficulty.** M, 4–6 h.

**What you learn.** *Physics:* that Newton I is a *definition of a class of frames*, not a special case of Newton II; that acceleration is the Galilean invariant and therefore `F = ma` is frame-independent among inertial frames, while position and velocity are not; that a non-inertial frame shows you motion without force, which is the entire motivation for ch. 9. *C++:* overloading a free function on a tag-like parameter type as an extension point; `constexpr bool` predicates on an aggregate; transforming a container into a new one without mutating the input (`[[nodiscard]]`, value semantics, and why `Trajectory` being cheap to copy makes this honest); writing a doc comment whose job is to record a *limitation* so that a future you does not over-trust the function.

---

### Proposal 8 — `examples/12_atwood.cpp`: one constraint, two bodies, Newton III for the tension

**Pitch.** The Atwood machine from Newton II per mass plus the third law plus `a1 = -a2`, validated against two exact formulas — with no new library code.

**Physics (§1.4, §1.5).** Two masses on a light string over a frictionless pulley:
```
m1 a = m1 g - T
m2 a = T - m2 g           (taking a > 0 as m1 descending)
=>  a = (m1 - m2) g / (m1 + m2)
    T = 2 m1 m2 g / (m1 + m2)

Limits worth asserting:  m2 -> 0  =>  a -> g, T -> 0     (free fall)
                         m1 == m2 =>  a = 0,  T = m g     (balanced)
```
As one degree of freedom: coordinate `x` = how far `m1` has descended; the coefficient of `x_ddot` is `m1 + m2` and the right-hand side is `(m1 - m2) g` (Proposal 6's borrowed "generalized mass / force" vocabulary — same algebra, no ch. 7 physics) — constant, so the motion is constant-acceleration and `analytic::position_const_accel` from Proposal 2 is the exact oracle.

**Honesty note:** this is *not* a Taylor ch. 1 problem. A search of all 51 ch.-1 problem statements finds no Atwood machine and no pulley; Taylor's Atwood machine is a ch. 7 Lagrangian showcase. The physics is entirely ch. 1, so build it — but don't tell yourself you're following the book here, and say so in the file header.

**Files.** Create `/Users/jackery/Github/newtonian-library/examples/12_atwood.cpp`, plus a small block in `/Users/jackery/Github/newtonian-library/tests/test_numeric_vs_analytic.cpp`. No header changes, no umbrella line.

**Signatures.** None public. In the example:
```cpp
// One degree of freedom: x = descent of m1. Coefficient of x_ddot: m1 + m2;
// right-hand side: (m1 - m2) g — the two Newton II equations with T
// eliminated by the third law and the constraint a1 = -a2.
const auto atwood = [m1, m2, g](double, const Vec3&, const Vec3&) {
    return Vec3{(m1 - m2) * g, 0.0, 0.0};
};
const Particle dof{m1 + m2, Vec3{}, Vec3{}};
```
The same slot-reinterpretation discipline as Proposal 6 applies, though here it is milder — `x` really is metres.

**Validation test.** `m1 = 3 kg`, `m2 = 2 kg`, `g = 9.80665`:

| Assertion | Expected | Measured |
|---|---|---|
| `a == (m1-m2)g/(m1+m2)` | 1.96133 m/s^2 | — |
| `x(2 s)` numeric vs `analytic::position_const_accel(2, {}, {}, a)` | 3.9226599999999996 m | numeric 3.9226600000000134, abs err **1.4e-14** at RK4 `dt = 1e-2` → tol 1e-12 |
| Newton III consistency: `T` from body 1 (`m1(g-a)`) == `T` from body 2 (`m2(g+a)`) | 23.53596 N both | identical to printed precision → tol 1e-12 |
| `T == 2 m1 m2 g/(m1+m2)` | 23.53596 N | → tol 1e-12 |
| limit `m2 = 0`: `a == g`, `T == 0` | — | tol 1e-12 |
| limit `m1 == m2`: `a == 0`, `T == m g` | — | tol 1e-12 |

The exactness at `dt = 1e-2` is the point: the motion is constant-acceleration, so RK4 is exact and any error you see is a bug in your setup, not truncation.

**Example sketch.** The free-body reasoning as comments, then `a` and `T` three ways (from `m1`, from `m2`, from the closed form) so third-law consistency is visible as three equal numbers; then the numeric-vs-analytic descent at 2 s; then a small table sweeping `m2/m1` from 0 to 1 showing `a -> g` and `a -> 0`, written to `build/atwood.csv`.

**Difficulty.** S, 2 h.

**What you learn.** *Physics:* that a constraint is an *equation* you impose and a constraint force is what you solve for; that the third law is what lets one symbol `T` appear in two equations; that eliminating the constraint force gives you a single-DOF problem. *C++:* reducing two coupled bodies to one `Particle` by choosing a generalized coordinate, and when that is legitimate; reusing an analytic oracle written for a different purpose; sweeping a parameter into a `std::vector<std::vector<double>>` for the generic `io::write_csv`.

---

### Proposal 9 — `systems/nbody.hpp`: a System whose pairwise forces obey Newton III by construction

**Pitch.** A multi-body state, a pairwise-force concept that makes `F_ij = -F_ji` structurally true, and momentum conservation to roundoff as the sharpest invariant ch. 1 offers.

**Physics (§1.5).**
```
Newton III, weak form:   F_12 = -F_21
Newton III, strong form: F_12 parallel to (r_1 - r_2)      (central)
P = sum_alpha m_alpha v_alpha
dP/dt = F_external       internal forces cancel pairwise      (1.25)-(1.29)
F_external = 0  =>  P = const
R_cm = (1/M) sum m_alpha r_alpha,   M R_cm_ddot = F_external
    (the R_cm definition is Taylor §3.3; ch. 1 gives you P = M V_cm, which is all this uses)
Perfectly inelastic collision (problem 1.30):  v' = m1 v / (m1 + m2)
```
Momentum conservation is a *structural* consequence of how you apply a pair force, not a numerical accident — which is why it should hold to roundoff rather than to truncation error, and why it is a sharper test than anything that depends on step size.

**Boundary to document now:** ch. 1 also shows the third law *failing* for magnetic forces between moving charges (problem 1.32), with a relative violation of order `v1 v2 / c^2`. So if you ever add a velocity-dependent pair force, the momentum guarantee below is no longer automatic. Write that into the header while everything is still velocity-independent.

**Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/systems/nbody.hpp`
- create `/Users/jackery/Github/newtonian-library/tests/test_system.cpp`
- create `/Users/jackery/Github/newtonian-library/examples/13_two_body.cpp`
- modify `/Users/jackery/Github/newtonian-library/include/newton/analysis/momentum.hpp` — add `total_mass` and `velocity_of_cm` (its own C++ note already flags the missing multi-body `System` that these would serve)
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — add **`#include <newton/systems/nbody.hpp>`** to the systems block
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — **line 15's scope sentence ("single-particle Newtonian mechanics") must change**, both new concepts must be specified beside "The Force concept" (lines 41–51), and a new file-by-file section added. `notes/04-extending.md:167` already names this as the v0.2 headline item; this moves it from roadmap to spec.

**Signatures — including the two new concepts.**
```cpp
namespace newton::systems {

/// A snapshot of a whole system: one Particle per body, stamped with a time.
/// C++ note: std::vector<Particle> rather than a parallel mass array, because
/// analysis::center_of_mass() and analysis::total_momentum() already take
/// exactly that type — s.bodies passes to both verbatim, with no adapter.
struct SystemState {
    std::vector<Particle> bodies;  ///< one entry per body (mass, r, v)
    double t{0.0};                 ///< [s]
};

// ===================== The PairForce concept ==============================
//
//   Vec3 operator()(double t, const Vec3& ri, const Vec3& vi,
//                             const Vec3& rj, const Vec3& vj) const;
//
// Returns the force ON body i DUE TO body j [N]. The stepper evaluates it
// ONCE per unordered pair {i, j} and applies +f to i and -f to j, so
// F_ij = -F_ji holds bit for bit — Newton III BY CONSTRUCTION (Taylor §1.5),
// not because the force's author remembered to be symmetric. A pair force
// physically cannot break the weak third law in this library; the only way
// to break it is to not use this interface.
//
// Note for a future you: this does NOT satisfy the single-particle Force
// concept, and forces::sum() cannot mix the two — compose.hpp's fold calls
// every pack element as fs(t, r, v), which a 5-argument callable will not
// accept. A system takes two force arguments, not one list.
//
// ===================== The BodyForce concept ==============================
//
//   Vec3 operator()(double t, const Vec3& r, const Vec3& v, double m) const;
//
// Returns the EXTERNAL force on ONE body [N]. The trailing mass is not
// decoration, and this is deliberately NOT the single-particle Force concept
// (DESIGN.md "The Force concept"), which takes (t, r, v) and has no way to
// know which body is asking or how heavy it is. The consequence is concrete:
// forces::UniformGravity stores ONE fixed `mass` and returns `mass * g` to
// every caller, so used across bodies of different mass it applies the wrong
// force to all but one of them — with masses (1, 3) kg, UniformGravity{1.0}
// makes the total external force 2g instead of M g = 4g, and
// UniformGravity{2.0} makes dP/dt come out right while body 1 falls at 2g and
// body 2 at 2g/3. forces::UniformGravity therefore does NOT satisfy this
// concept (wrong arity, see the static_assert in tests/test_system.cpp), and
// must NEVER be adapted to it by dropping the mass argument. Use
// UniformField below, or write the field you mean:
//     [g](double, const Vec3&, const Vec3&, double m) { return m * g; }

/// External uniform field: F_i = m_i g. Satisfies BodyForce.
/// Physics: the mass cancels in a_i = F_i/m_i = g, so every body falls at the
/// same rate whatever its mass (Galileo), and dP/dt = (sum m_i) g = M g
/// exactly — the closed form the momentum test below checks against.
struct UniformField {
    Vec3 g{0.0, -constants::g_earth, 0.0};  ///< [m/s^2]
    [[nodiscard]] Vec3 operator()(double, const Vec3&, const Vec3&, double m) const { return m * g; }
};

/// No external force at all — the isolated system of Taylor §1.5, where P is
/// exactly conserved. Satisfies BodyForce.
struct NoExternalForce {
    [[nodiscard]] constexpr Vec3 operator()(double, const Vec3&, const Vec3&, double) const { return {}; }
};

/// A central pair interaction of constant magnitude f0, attractive for
/// f0 > 0:  f_ij = -f0 * normalized(ri - rj).  Satisfies PairForce.
/// Physics: central (parallel to ri - rj), so it obeys the STRONG form of
/// Newton III as well as the weak one (Taylor §1.5; problem 1.34 shows the
/// strong form is the extra assumption angular momentum needs). Deliberately
/// not an inverse-square or Hooke law: those belong to Taylor §4.8 and §5.1,
/// and nothing here needs a force law with a closed-form orbit.
/// Zero-guarded at coincident positions via normalized().
struct CentralPair { double f0{1.0}; /* [N] */ /* operator() ... */ };

/// Total force on every body: its BodyForce, plus the PairForce accumulated
/// once per unordered pair and applied +/-. Exposed separately from the
/// stepper because it is the piece worth reading and the piece worth testing.
template <class PairForce, class BodyForce>
[[nodiscard]] std::vector<Vec3> total_forces(const SystemState& s,
                                             PairForce&& pair, BodyForce&& body);

/// One velocity-Verlet step of the whole system. Mirrors
/// integrators::velocity_verlet_step's kick-drift-kick with the
/// predictor-corrector velocity handling, one body at a time, so there is one
/// update-equation idiom in the library rather than two.
template <class PairForce, class BodyForce>
[[nodiscard]] SystemState nbody_verlet_step(const SystemState& s, double dt,
                                            PairForce&& pair, BodyForce&& body);

struct SystemTrajectory { std::vector<SystemState> samples; /* size/empty/front/back */ };

/// Fixed-step loop over a system, with the same time bookkeeping as
/// newton::simulate (t = t0 + i*dt, never accumulated) and the same
/// record_every stride.
/// LIMITATION, stated rather than hidden: opt.method is IGNORED — only
/// velocity Verlet is implemented for systems. The honest fix is to template
/// the four single-body steppers on the state type so State and SystemState
/// both satisfy them and a single particle becomes the N = 1 case; that
/// touches four files DESIGN.md specifies signature by signature, so it is a
/// spec change, not a drop-in. See notes/04-extending.md's warning against
/// ever maintaining two integrator code paths.
template <class PairForce, class BodyForce>
[[nodiscard]] SystemTrajectory simulate_system(const SystemState& s0, const SimOptions& opt,
                                               PairForce&& pair, BodyForce&& body);

}  // namespace newton::systems

namespace newton::analysis {
/// M = sum m_i [kg].
[[nodiscard]] inline double total_mass(const std::vector<Particle>&);
/// V_cm = P / M [m/s]; zero-guarded for M <= 0 like center_of_mass().
[[nodiscard]] inline Vec3 velocity_of_cm(const std::vector<Particle>&);
}  // namespace newton::analysis
```

**Validation test** (`tests/test_system.cpp`). Fixture: two bodies, `m = (1, 3) kg`, `r = ((-1,0,0), (2,0.3,0))`, `v = ((0,0.5,0), (0.1,-0.2,0))`, so `M = 4 kg`, `|P0| = 0.3162 kg m/s`; `CentralPair{2.0}`; `dt = 1e-3`; 10^4 steps (10 s); minimum separation 0.58 m (no singularity).

1. **The concepts are enforced at compile time.**
   ```cpp
   template <class F, class = void> struct is_body_force : std::false_type {};
   template <class F>
   struct is_body_force<F, std::void_t<decltype(std::declval<const F&>()(
           0.0, std::declval<const Vec3&>(), std::declval<const Vec3&>(), 1.0))>>
       : std::true_type {};

   static_assert(!is_body_force<forces::UniformGravity>::value,
                 "forces::UniformGravity is single-mass and must not be usable as a BodyForce");
   static_assert(is_body_force<systems::UniformField>::value, "UniformField must satisfy BodyForce");
   ```
   Both `static_assert`s were compiled against the real headers: they pass.
2. **Newton III structurally.** With `NoExternalForce`, `total_forces` must sum to zero. For **N = 2** assert **exact** `== Vec3{}` — measured `(0, 0, 0)`, bit-exact, because `F[0] = f`, `F[1] = -f`, and IEEE addition of a value and its negation is exact. For **N = 3** assert 1e-15 instead, not exact: three pairs cancel pairwise but the *order* of the additions differs per body, so `(f01+f02) + (-f01+f12) + (-f02-f12)` need not round to zero — measured **2.2e-16**. That contrast is one of the best floating-point lessons in the repo; put it in a comment.
3. **Strong third law (problem 1.34).** `norm(cross(ri - rj, f_ij)) == 0` for `CentralPair` — tol 1e-15, measured **1.3e-16** with `|f| = 2 N`. Geometry only: angular momentum itself is ch. 3.
4. **Isolated system, momentum exactly conserved.** `max |P(t) - P0|` over 10^4 steps: tol **1e-12**, measured **2.4e-14**. Then the assertion that makes it meaningful: re-run at `dt = 5e-4` for 2x10^4 steps and assert the drift is *still* below 1e-12. Truncation error would fall by 4x; roundoff will not — and if your drift instead *grows* with step size, your force application is asymmetric, which is the bug this test exists to catch.
5. **Centre of mass moves uniformly.** `max |R_cm(t) - (R0 + V_cm t)|` over the same run: tol **1e-12**, measured **2.6e-14**. This is the dynamical half of §1.5 that the existing `center_of_mass` could never demonstrate.
6. **External field.** Same pair force, `body = systems::UniformField{{0, -9.80665, 0}}` (a `BodyForce`, *not* `forces::UniformGravity`). Then `P(t) = P0 + M g t` and `R_cm(t) = R0 + V_cm t + (1/2) g t^2`, tol **1e-9 absolute** on both. Measured over 10^4 steps: `max |P - (P0 + M g t)|` = **1.25e-12** and `max |R_cm - ...|` = **8.5e-13**, with `|P|` reaching ~392 kg m/s by t = 10 s (so 1.25e-12 is ~3e-15 relative — roundoff, as it should be, since the field is constant and Verlet is exact for constant acceleration).
7. **Inelastic collision (problem 1.30), algebra only.** Build the before state (`m1 = 2` at 7 m/s, `m2 = 3` at rest) and the after state (both at `v'`), and assert `total_momentum` is equal: `v' = m1 v/(m1+m2)` = **2.8 m/s**, tol 1e-12. No integration — collisions are not a force law, and ch. 1 treats this as momentum bookkeeping.

**Example sketch** (`examples/13_two_body.cpp`): the two-body fixture under `CentralPair` alone. Prints `|P(t) - P0|` and `|R_cm(t) - (R0 + V_cm t)|` at ten sample times, with a line explaining that these stay at 1e-14 because the cancellation is structural, not because the integrator is accurate. Then the same run with `UniformField` added, printing `P(t) - P0` beside `M g t`. Writes `build/two_body.csv` through the generic overload with a `body` tag column (0, 1, and 2 for the centre of mass) so all three paths plot together and the CoM's straight line / clean parabola is visible against the two wiggling bodies.

**Difficulty.** L, 12–20 h. Real design work, correctly scoped as v0.2 in `notes/04-extending.md`. Do the other nine first.

**What you learn.** *Physics:* that the third law is what makes momentum conservation *derivable* rather than empirical; the difference between the weak and strong forms and why the strong one is extra; that only external forces move the centre of mass; that a conservation law can be built into the structure of a program so that violating it requires bypassing the interface. *C++:* defining a *concept by documentation plus a detection idiom* (`std::void_t`, `std::declval`, partial specialization, `static_assert`) — the C++17 way to get a compile-time error instead of a physics bug; two template callable parameters with distinct arities in one function; why `forces::sum`'s fold cannot mix them; designing a new aggregate (`SystemState`) to be compatible with existing free functions instead of requiring an adapter; floating-point cancellation that *is* exact versus cancellation that merely rounds small, and how to write a test that distinguishes roundoff from truncation by changing `dt`.

---

### Proposal 10 — `notes/05-chapter-1-problems.md`: the study page that closes the loop

**Pitch.** A table mapping problems 1.1–1.51 to the library call, test, or example that answers each — so the library becomes a way of *doing the problem set*, not a parallel activity.

**Physics.** All of ch. 1; no new physics.

**Files.**
- create `/Users/jackery/Github/newtonian-library/notes/05-chapter-1-problems.md`
- modify `/Users/jackery/Github/newtonian-library/notes/02-physics-to-code.md` — add §1.7's polar set, §1.4's Galilean transformations, and §1.5's system momentum to its ch. 1 key-equations block
- modify `/Users/jackery/Github/newtonian-library/notes/00-start-here.md` — add the new note to the reading order at line 328
- modify `/Users/jackery/Github/newtonian-library/docs/research/taylor-ch1-7-map.md` — **fix three errors in its ch. 1 section** (lines 38–122), each verified against all 51 problem statements: (i) Atwood's machine is **not** a ch. 1 problem (zero hits for "Atwood" or "pulley" across all 51); (ii) a **conical** pendulum is not in ch. 1 — problem 1.41's twirling mass is explicitly in gravity-free space, so there is no cone; (iii) cylindrical coordinates are **ch. 1** (problems 1.47–1.49), not ch. 4 as the map defers them. Also soften the map's line-96 claim that ch. 1 is the chapter where the integrators "should all appear": ch. 1 justifies exactly one, and the `[Computer]` tag appears on exactly two problems.

**Signatures.** None. One table: *problem · what it asks (your own words) · library surface · file*. Group by the SECTION headers the book prints in its problem set (the boundaries here are from memory — in particular whether the ODE warm-ups 1.24–1.25 sit under §1.2 or §1.4 — so take them from your copy): 1.1–1.25 → §1.2 vectors and ODE warm-ups; 1.26–1.34 → §1.4–1.5 frames, third law, momentum; 1.35–1.40 → §1.6 Cartesian Newton II; 1.41–1.51 → §1.7 polar and cylindrical. Paraphrase every problem in your own words — do not transcribe Taylor's text.

**Validation.** A checklist at the top: every row naming a test must name a test that exists and passes. Run `make check` and paste the summary line.

**Difficulty.** S, 3–4 h (do it last; it is the index to everything above).

**What you learn.** *Physics:* the chapter as a *structure* rather than 51 unrelated exercises — you'll notice that 20 of them are the same vector-identity skill and that only two require a computer. *C++/tooling:* documentation as a first-class deliverable; keeping a spec (`DESIGN.md`), a manifest (`include/newton/newton.hpp`) and a research doc in sync, and what it feels like when they have silently drifted apart.

---

## 5. Already in the repo — read these next

Read in this order. Apart from the one comment flagged on `01_free_fall.cpp`, every line is ch. 1 physics or pure C++/numerics.

| File | Why now |
|---|---|
| `include/newton/core/vec3.hpp` | Is §1.2, end to end. Skim past the one-line §3.4 / §2.5 motivations in `cross`'s comment. |
| `include/newton/core/particle.hpp` | 41 lines that *are* §1.6: the header derives `dr/dt = v`, `dv/dt = F/m` from `m r_ddot = F`, and `momentum()` is `p = m v`. The single most valuable file in the repo for you — where the paper equation becomes a data structure. |
| `include/newton/core/trajectory.hpp` | Pure data: `mass` + `std::vector<State>`. |
| `include/newton/core/constants.hpp` | `pi`, `deg`, `g_earth`, `G` are yours. `sphere_air_beta`/`sphere_air_gamma` are ch. 2 — inert, ignore. |
| `include/newton/forces/compose.hpp` | `forces::sum(...)` is §1.3 superposition. Its C++17 fold expression is the only hard part, and it is C++, not physics. |
| `include/newton/analytic/projectile.hpp` | §1.6 verbatim: `range`, `apex_height`, `time_of_flight`, `position`. The library's first validation oracle. |
| `include/newton/simulate.hpp` | No new mechanics — time bookkeeping (`t = t0 + i*dt`, never `t += dt`) and event bisection. Read it as numerics. You'll use `simulate_until` constantly. |
| `examples/01_free_fall.cpp` | The only existing example whose *code* is all ch. 1. One comment (lines 109–111) justifies the impact-speed check in ch. 4's energy language — read the check, defer the justification. Study it beside `analytic/projectile.hpp`. |
| `tests/test_vec3.cpp` | Every check is a §1.2 identity, and the gentlest introduction to the harness you're about to extend. |
| `include/newton/io/csv.hpp` | No physics at all; you'll need both overloads. |

**Read now, but only in part.**
- `include/newton/forces/gravity.hpp` — `UniformGravity{mass, g}` is ch. 1. `PointGravity` is §4.8 with Plummer softening; defer.
- `include/newton/forces/drives.hpp` — `ConstantForce{F}` is ch. 1 and the natural partner for constant-acceleration exercises. `SinusoidalDrive` is §5.5; defer.
- `include/newton/analysis/momentum.hpp` — `center_of_mass` and `total_momentum` are §1.5 material. `angular_momentum` is §3.4; defer.
- `include/newton/integrators/euler.hpp` and `semi_implicit_euler.hpp` — the update equations are §1.6's first-order pair discretized, and the one-line difference between the two files (`v_n` vs `v_{n+1}`) is the clearest lesson in the repo. But both file headers justify themselves with oscillator energy growth, symplecticity, and shadow Hamiltonians — ch. 4–5 vocabulary. Read the code, defer the justification.
- `include/newton/integrators/velocity_verlet.hpp` and `rk4.hpp` — same split, more so. The one claim you can verify now: "exact for constant acceleration at any `dt`", since the true solution is a degree-2 polynomial (stated in `velocity_verlet.hpp` lines 14 and 92 and in `examples/01_free_fall.cpp`'s header; asserted for both methods in `tests/test_integrators.cpp` block (a)).
- `tests/test_forces.cpp` — the `UniformGravity`, `ConstantForce`, and `sum()` blocks are readable now; the drag / spring / drive / magnetic blocks are not.
- `tests/test_analytic.cpp` — its final block (vacuum projectile, range at 45 deg = `v0^2/g`) is the one ch.-1-readable section.

## Skip for now — needs ch. X

| File(s) | Chapter |
|---|---|
| `include/newton/forces/drag.hpp`, `include/newton/analytic/drag.hpp`, `examples/02_projectile_drag.cpp` | **Ch. 2** §2.1–2.4 |
| `include/newton/forces/electromagnetic.hpp`, `include/newton/analytic/magnetic.hpp`, `examples/03_charged_particle.cpp` | **Ch. 2** §2.5, §2.7 |
| `include/newton/systems/rocket.hpp`, `include/newton/analytic/rocket.hpp`, `examples/04_rocket.cpp` | **Ch. 3** §3.2 (variable mass) |
| `include/newton/analysis/momentum.hpp::angular_momentum` | **Ch. 3** §3.4 |
| `include/newton/analysis/energy.hpp`, `examples/06_integrator_comparison.cpp` | **Ch. 4** §4.1–4.4 |
| `include/newton/forces/spring.hpp`, `drives.hpp::SinusoidalDrive`, `include/newton/analytic/oscillator.hpp`, `include/newton/analysis/oscillation.hpp`, `examples/05_oscillations.cpp` | **Ch. 5** §5.1–5.6 |
| `tests/test_analytic.cpp`, `tests/test_numeric_vs_analytic.cpp`, `tests/test_analysis.cpp` | Mixed; each is dominated by ch. 2/4/5 fixtures |

**And the features not to build yet**, however tempting: drag and terminal velocity (ch. 2); the Lorentz force (ch. 2); variable mass (ch. 3); angular momentum as a diagnostic (ch. 3); **energy drift as an integrator metric** (ch. 4 — use Proposal 6's turning-point drift and Proposal 9's momentum conservation instead, both energy-free and both sharper); the closed form for `N` at the bottom of the half-pipe (ch. 4); damped/driven oscillators and any amplitude-dependent period formula, including the `1 + phi0^2/16` correction (ch. 5 and standard pendulum theory — **not** Taylor ch. 1); a general holonomic constraint solver (ch. 6–7); Kepler orbits (ch. 8); centrifugal and Coriolis **forces** (ch. 9 — ch. 1's `2 r_dot phi_dot` term is kinematic, and problems 1.27/1.46 ask only what the rotating observer *sees*).

**One design decision to make now, while it's free.** Velocity Verlet and the symplectic family assume `F = F(r)`. Every force ch. 1 actually *integrates* satisfies that (the magnetic forces of problems 1.32, 1.33 and 1.45 are velocity-dependent, but they appear only in the third-law discussion, never in an equation of motion); ch. 2's drag and ch. 5's damping do not. Add a trait — `static constexpr bool depends_on_velocity = false;` on each force functor, or a `depends_on_velocity()` predicate — now, while every answer is `false`, and have the integrator choice consult it. Discovering this in ch. 2 costs far more.

---

## 6. Source honesty

### Web-confirmed this session (URLs)

- **Problems 1.1–1.51 exist; all 51 statements, every equation-number binding, and every closed-form answer quoted above** — all 51 Chapter-1 solution PDFs from `https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%201/TaylorCMCh1p<N>.pdf`, downloaded and text-extracted. Index: https://stemjock.com/taylorcm.htm. Independently corroborated for the 51-problem count and the 1.1–1.23 topics by https://studysoup.com/tsg/physics/90/classical-mechanics/chapter/1765/1
- **The `[Computer]` tag appears on exactly two problems, 1.50 and 1.51** — verified exhaustively across all 51 extracted statements, not by sampling. (A naive grep also hits every file, but only inside embedded "Computer Modern" font-licence strings — a false positive, checked and excluded.)
- **Equation-number bindings:** (1.6) and (1.7) dot product, (1.9) cross product, (1.35) Newton II in Cartesian components, (1.37) the rectangular↔polar transformations and the incomplete `arctan`, (1.42) and (1.46) the basis derivatives, (1.48) Newton II in polar, (1.51) the half-pipe equation **and that it is on page 31**, (1.55)/(1.56)/(1.57) the small-angle solution, (1.59) `r_hat` — all cited verbatim inside the fetched statements of problems 1.7, 1.15, 1.35, 1.41, 1.42, 1.43, 1.44, 1.50, 1.51.
- **Example 1.2 is the skateboard in a curved trough, and is the source of (1.51)**; `R = 5 m`, `g = 9.8 m/s^2` — problem 1.50's fetched statement names it. Polar velocity/acceleration and `d(r_hat)/dt = phi_dot phi_hat` independently confirmed at https://dandavison.github.io/taylor-classical-mechanics.html
- **Inertial frames as problems 1.26/1.27 pose them** (constant-velocity frame inertial, constant-acceleration frame not, rotating turntable not): the fetched 1.26 solution works it in components.
- **The third law's strong/central form is an extra assumption, forward-referenced to ch. 3** — problem 1.34's statement. **The momentum derivation is eqs. (1.25)–(1.29)** — problem 1.28's statement. **The magnetic counterexample and Figure 1.8** — problem 1.32's statement.
- **Tilted axes with `W_x = -mg sin phi`, `W_y = -mg cos phi`** — the fetched 1.39 solution. **The "no normal force in flight" confusion** — https://www.physicsforums.com/threads/classical-mechanics-taylor-1-39.432209/
- **Closed-form answers** for 1.35 (`T = 2v0 sin θ/g`, `R = v0^2 sin 2θ/g`), 1.36 (lead `v0 sqrt(2h/g)` ≈ 224 m, ±0.2 s), 1.37 (`t_return = 2v0/(g sin θ)`), 1.39 (`R`, `R_max`), 1.40 (`θ_max = arcsin sqrt(8/9)` ≈ 70.5 deg), 1.41 (`T = mRω^2`), 1.49 (helix, `N = mRω^2`), 1.30 (`v' = m1 v/(m1+m2)`) — all read out of the fetched solution PDFs.
- **Plot windows 0–14 s and 0–21 s, and `ω = 1.4 rad/s` exactly for Taylor's numbers** — the fetched 1.50 and 1.51 solutions. **Problem 1.51's published verdict that the small-angle curve fails badly at 90 deg** — same source.
- **Which topics live in which later chapters** — https://courses.physics.ucsd.edu/2010/Fall/physics110a/HOMEWORK/hwsolutions1_F10.pdf (Taylor 1.46, 2.14, 3.11, 3.13, 4.2, 4.23, 4.36). Rotating-frame fictitious forces belong to the later chapter, corroborated by https://homepages.uc.edu/~brodjm/teaching/lecture_notes/mechanics.pdf
- **Ch. 1 taught as one week, sections 1–7** — https://www.usna.edu/Users/physics/cmorgan/SP333.php; ch. 1 as Newton II + polar with numerical solutions — https://rhettallain.com/classical-mechanics-course-material-and-videos/
- **Absences, verified by searching all 51 statements:** no Atwood machine, no pulley, no conical pendulum in ch. 1's problems.

No Taylor text, problem statement, or figure is reproduced anywhere in this report; every problem is paraphrased. The several unauthorised full-text copies that appeared in search results (archive.org full-text dump, pdfcoffee, silo.tips, slideshare, studylib, scribd, zlib.pub, sweetstudy, and an Internet Archive copy of the instructor's solution manual) were deliberately **not** opened.

### From memory, not web-confirmed

- The seven section titles and their page numbers as a table. Search summaries agreed, but `henry.pha.jhu.edu/chapters.html` returned 403 and the other TOC hits were pirated copies. Only "eq. (1.51) is on page 31" is web-confirmed.
- That **Example 1.1** is a block on an incline and sits in §1.6. A search summary agreed; no page was fetched. What *is* confirmed is that the tilted-axes technique is the backbone of problems 1.37–1.39.
- §1.3's treatment of mass, force, and superposition — no fetched page corroborates the specifics.
- Where the BAC−CAB identity is introduced (problem 1.33 assumes the reader has it; that much is confirmed).
- `h_apex = v0^2 sin^2 θ/(2g)`, `t_apex = v0 sin θ/g`, range maximal at 45 deg, and the Atwood formulas — elementary, no source needed.
- That drag is ch. 2, energy is ch. 4, oscillators are ch. 5, Lagrangians are ch. 6–7, orbits are ch. 8, and non-inertial frames are ch. 9 — corroborated by the UCSD problem set for 2/3/4 but memory for the rest.
- The `1 + phi0^2/16` pendulum-period correction. Standard physics, **not** Taylor ch. 1 — used above only as a tolerance guide, never cited to the book.

### Check against your copy

1. The seven section titles and page numbers (claimed: 1.2 → 4, 1.3 → 9, 1.5 → 17, 1.7 → 26, summary → 33, problems → 34).
2. That Example 1.1 is the block on an incline in §1.6 — and whether Taylor's version includes friction and solves for the normal force. **Proposal 3 is built on this assumption.**
3. How many numbered examples ch. 1 contains. Example 1.2 is confirmed to exist; whether the chapter stops there is unknown.
4. Whether an Atwood machine appears in ch. 1's **running text** (definitively absent from the 51 problems). If it's in the text, Proposal 8 earns a legitimate citation.
5. Whether a conical pendulum appears in the text (absent from the problems).
6. Taylor's difficulty-star convention — which ch. 1 problems carry which rating. No transcription preserves the stars, so none of the examples above can be graded by the book's own difficulty.
7. Whether `[Computer]` really appears on only 1.50 and 1.51. Strong evidence, but it depends on stemjock having transcribed the tag faithfully everywhere. A 30-second flip through pages 34–41 settles it.
8. The precise wording of eq. (1.35) — three scalar equations, one vector equation, or both.
9. Where §1.2 introduces BAC−CAB, if it does.
10. Whether Taylor measures `phi` from the bottom of the trough or from the vertical. **Proposal 6's sign convention for `N = m(g cos phi + R phi_dot^2)` assumes the former.**

### Photos that would help

Highest value first.

1. **Pages 34–41 — the full problem set, with `[Computer]` tags and difficulty stars.** Single highest value. Two things ride on it: whether `[Computer]` is really only on 1.50/1.51 (if a third problem carries it, Proposal 6 grows), and the star ratings, which no transcription preserves.
2. **Problems 1.50 and 1.51 specifically.** These are the only two problems Proposal 6 claims to reproduce. Confirm the plot windows (0–14 s and 0–21 s are second-hand), the exact `phi0` values, and whether Taylor asks for anything beyond `phi(t)` and the overlay.
3. **Pages 26–33 — all of §1.7, including Example 1.2 and eqs. (1.37)–(1.57).** The section Proposals 4 and 6 are built from. Specifically: the exact form and sign convention of eq. (1.48); whether `phi` is measured from the bottom of the trough or from the vertical; whether Example 1.2 computes the normal force explicitly.
4. **Pages 23–26 — §1.6, including Example 1.1.** Proposal 3 assumes this is the incline example. Does Taylor's version include friction? Does he solve for `N`? Is eq. (1.35) written as three scalars, one vector, or both?
5. **The chapter's table of contents and the "Principal Definitions and Equations of Chapter 1" page (p. 33).** One page confirms every section title, page number, and equation-number binding relied on second-hand, and settles how many numbered examples the chapter has.
6. **Pages 17–23 — §1.5, especially the N-particle derivation (1.25)–(1.29) and Figure 1.8's caption.** Proposal 9's doc comments paraphrase both; better to match Taylor's internal/external notation than invent one.
7. **Pages 9–16 — §1.3 and §1.4.** §1.3 is the only section with no fetched corroboration at all. Proposal 9's force split and Proposal 7's framing of Newton I would both be firmer with the text in view.
8. **Problems 1.47–1.49** (cylindrical coordinates). Proposal 4's cylindrical section rests largely on 1.49's closed-form helix. If 1.47's notation differs from `(rho, phi, z)`, the struct should match the book.

---

## 7. How to work

**One feature at a time, and always in this order: header → test → example → umbrella line → `make check`.**

1. **Header first, with the doc comment written before the body.** If you can't state the physics and the units in the comment, you don't yet know what the function does. Keep the style the repo already has: the equation in plain text, the Taylor section in parentheses, units in `///<` on every member, and a separate `C++ note:` paragraph whenever a language or numerics decision needed a reason (`atan2` over `acos`, `inline` over `constexpr`, a zero guard, a struct instead of a `Vec3`). Those notes are the most valuable thing in the repo for the person you'll be in six months.
2. **Test second, never last.** A header with no test is a hypothesis. Every proposal above names its assertions and their measured residuals so you can tell a real failure from a tolerance you set too tight. When a tolerance is tied to `dt` or to a finite-difference `h`, **write the dependence into the comment** — otherwise a future you will "tighten" it and break the build for no reason.
3. **Example third.** An example is where you find out the API is awkward. If calling your own function needs three lines of unpacking, fix the function, not the example.
4. **Umbrella line fourth.** `include/newton/newton.hpp` is the manifest. A header that isn't listed there effectively doesn't exist. And when `DESIGN.md` disagrees with your change — a new namespace, a widened scope sentence, a "do not modify" heading — amend `DESIGN.md` **in the same commit**. `notes/00-start-here.md` makes it the contract; a contract that silently drifts from the code is worse than no contract.
5. **`make check` fifth, every time.** `Makefile` discovers examples and tests by wildcard and treats every header as a prerequisite, so new files need no build edits. There is no excuse for a long gap between writing code and running it.
6. **Build them in the order above.** Proposals 1–3 are small, independent, and give you the habit. Proposal 4 unlocks §1.7 and is the real prerequisite for 5 and 6. Proposal 6 is the deliverable — it's Taylor's actual homework. Leave 9 for last; it is v0.2-sized and will change `DESIGN.md`'s scope sentence.
7. **When you assert a measurement, state what it was measured over.** The one bug in the previous version of this plan was a chained loop that reported six periods and silently covered four — every number in it was individually correct and the conclusion was still wrong. A guard that can't see the failure mode you actually have (`|phi_dot| < 1e-10` cannot detect a duplicate turning point) is worse than no guard, because it reports success. When you write a post-condition, ask what the wrong answer would look like and whether your check would notice.
