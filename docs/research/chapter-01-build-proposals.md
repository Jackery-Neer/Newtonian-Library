> **Research input — Build proposals, third revision.** Raw output of an Opus research agent (2026-09-10), kept for deep dives. It has been verified by Fable agents and condensed into `notes/chapters/ch01-what-to-build.md`, which is the authoritative, corrected version. KNOWN DEFECTS in this draft, fixed in the final report: (1) Proposal 4 asserts `to_polar(from_polar({1,-pi})).phi == +pi` — false in IEEE doubles, the round trip returns -pi; (2) Proposal 6 chains `simulate_until` legs without the one-step "nudge", so its "6 periods" figures were measured over ~4. Read the final report's Proposals 4 and 6 instead. No text from Taylor's book is reproduced; problems are paraphrased.

---

## Revision notes

Two verifier items, both resolved, plus one process change.

1. **[major] Body-force concept (Proposal 9).** The claim that external per-body forces can use the existing single-particle `Force` concept is withdrawn — it is false, for exactly the reason given: `forces::UniformGravity` carries one fixed `mass` and returns the same vector to every caller, so on masses (1, 3) kg it applies the wrong force to at least one body and the `P(t) = P(0) + M g t` check fails by O(1). Proposal 9 now defines a **separate `BodyForce` concept carrying the body's mass as a fourth argument** — `Vec3 operator()(double t, const Vec3& r, const Vec3& v, double m) const` — and ships `systems::UniformField{g}` (returns `m * g`) as the concrete object assertion 5 passes. The doc comment states that `forces::UniformGravity` is single-mass, does not satisfy `BodyForce` (wrong arity), and must never be adapted by dropping the mass argument; a `std::void_t` detection idiom turns that into a `static_assert` that I compiled and confirmed passes. Assertion 5 is restated with a measured number from a compiled prototype, not an estimate.
2. **[major] Half-pipe normal force (Proposal 6).** "Or scan for φ̇ sign changes" is removed as an acceptable route for any tight assertion. Proposal 6 now specifies **one** method — chained `simulate_until` legs on a new `events::turning_point_x()`, each leg re-launched from the previous refined crossing — and documents that `simulate_until` returns exactly one turning point per run (the strict `<` at `simulate.hpp:185` correctly skips the `t0` sample where φ̇ = 0, and the function returns on the first crossing, `simulate.hpp:201-202`). "Exactly 0" at 90° is replaced by "0 to 1e-9". The test now pins `dt ≤ 1e-3`, because at `dt = 1e-2` the refined 90° normal-force error is 1.4e-9 and a 1e-9 assertion would fail.
3. **Process.** Every tolerance and "measured" figure in this revision came out of three prototype programs compiled against the real headers with `clang++ -std=c++17 -Wall -Wextra -Wpedantic -O2 -I include` and run. Sources kept at `/private/tmp/claude-501/-Users-jackery-Github-newtonian-library/e1090fdc-4da5-488d-aae7-cac525791aa8/scratchpad/probe.cpp`, `probe2.cpp`, `probe3.cpp`. The contrast figures for rejected methods (scan-only turning points, Euler) were produced by the same program in the same run as the figures for the recommended method, so they are directly comparable.

---

# Chapter-1-only additions to `newton` — prioritized build sequence

You have finished Taylor ch. 1 and nothing further. Everything below is buildable and *validatable* with ch. 1 physics alone. Ten proposals, ordered so each one's dependencies are already done: quick wins first, the two genuinely large items last.

**The one-line summary of what ch. 1 actually demands of a simulator:** integrate `φ̈ = -(g/R) sin φ` from rest at 20° and at 90°, and plot it against `φ₀ cos(√(g/R) t)`. That is Taylor's only `[Computer]` work in the chapter (problems 1.50, 1.51). Proposal 6 is therefore the headline deliverable; Proposals 1–5 are what make it clean, and 7–10 are the rest of the chapter.

## Already there — study these first

Before writing anything, read these in this order. Every line is ch. 1 physics or pure C++/numerics.

| File | Why now |
|---|---|
| `/Users/jackery/Github/newtonian-library/include/newton/core/vec3.hpp` | Is §1.2, end to end. Skim past the one-line §3.4 / §2.5 motivations in `cross`'s comment. |
| `/Users/jackery/Github/newtonian-library/include/newton/core/particle.hpp` | 41 lines that *are* §1.6: the header derives `dr/dt = v`, `dv/dt = F/m` from `m r̈ = F`, and `momentum()` is `p = m v`. The single most valuable file in the repo for you. |
| `/Users/jackery/Github/newtonian-library/include/newton/core/trajectory.hpp` | Pure data. |
| `/Users/jackery/Github/newtonian-library/include/newton/core/constants.hpp` | `pi`, `deg`, `g_earth`, `G` are yours. `sphere_air_beta`/`sphere_air_gamma` are ch. 2 — inert, ignore. |
| `/Users/jackery/Github/newtonian-library/include/newton/forces/compose.hpp` | `forces::sum(...)` is §1.3 superposition of forces. Its C++17 fold expression is the only hard part, and it is C++, not physics. |
| `/Users/jackery/Github/newtonian-library/include/newton/analytic/projectile.hpp` | §1.6 verbatim; the library's first validation oracle. |
| `/Users/jackery/Github/newtonian-library/include/newton/simulate.hpp` | No new mechanics — time bookkeeping (`t = t0 + i*dt`, never `t += dt`) and event bisection. Read it as numerics. You will use `simulate_until` constantly below. |
| `/Users/jackery/Github/newtonian-library/examples/01_free_fall.cpp` | The only existing example with no post-ch.-1 physics. |
| `/Users/jackery/Github/newtonian-library/tests/test_vec3.cpp` | Every check is a §1.2 identity, and it is the gentlest introduction to the harness you are about to extend. |
| `/Users/jackery/Github/newtonian-library/include/newton/io/csv.hpp` | No physics at all; you need both overloads. |

Partly readable now: `forces/gravity.hpp` (`UniformGravity` yes, `PointGravity` is §4.8 — defer), `forces/drives.hpp` (`ConstantForce` yes, `SinusoidalDrive` is §5.5 — defer), `analysis/momentum.hpp` (`center_of_mass` / `total_momentum` yes, `angular_momentum` is §3.4 — defer), the four integrators (read the *update equations*; their file headers justify themselves with oscillator energy, which is ch. 4–5 vocabulary — skip the justification, keep the one claim you can verify: exact for constant acceleration at any `dt`).

---

## 1. Finish §1.2 on `Vec3`: `angle`, triple products, and nineteen free tests

**(a)** Add the three vector-algebra functions ch. 1 uses and turn Taylor's problems 1.4–1.23 into a test suite. Cheapest high-confidence win in the repo.

**(b) Physics (§1.2).**
```
a . b     = |a||b| cos(theta)                  (1.6)  and  sum_i a_i b_i   (1.7)
|a x b|   = |a||b| sin(theta)                          right-hand rule      (1.9)
a . (b x c)                                    scalar triple = signed parallelepiped volume
a x (b x c) = (a.c) b - (a.b) c                "BAC minus CAB"
d/dt (r.s) = r . s_dot + r_dot . s
d/dt (r x s) = r x s_dot + r_dot x s           ORDER MATTERS
```

**(c) Files.**
- modify `/Users/jackery/Github/newtonian-library/include/newton/core/vec3.hpp`
- modify `/Users/jackery/Github/newtonian-library/tests/test_vec3.cpp`
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md`

Umbrella line: **none** — `#include <newton/core/vec3.hpp>` is already line 19 of `include/newton/newton.hpp`.

Spec friction you must resolve deliberately: `DESIGN.md:57` heads the core section **"Already written (do not modify): core"**. That line was an instruction to the v0.1 implementation agents, not a permanent freeze, but `notes/00-start-here.md` makes `DESIGN.md` authoritative — so amend `DESIGN.md:59-60` to list the three new functions in the same commit. (If you would rather keep the freeze literal, put them in a new `include/newton/core/vec3_algebra.hpp`; I recommend against it — they belong next to `dot` and `cross`.)

**(d) Signatures.**
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
`scalar_triple` and `vector_triple` can be `constexpr` (pure arithmetic); `angle` cannot. That contrast is worth a sentence in the doc comments.

**(e) Closed forms and assertions.** All measured on `a = (1,2,3)`, `b = (-2,0.5,4)`, `c = (3,-1,0.25)`.

| Taylor | Assertion | Tolerance | Measured residual |
|---|---|---|---|
| 1.4 | `angle(a,b) == acos(dot(a,b)/(norm(a)*norm(b)))` | 1e-12 | — |
| 1.5 | `angle({1,1,1},{1,1,0}) == acos(sqrt(2/3))` = 35.2643896828° | 1e-12 | agrees to 12 digits |
| 1.7 | `dot(a,b) == norm(a)*norm(b)*cos(angle(a,b))` | 1e-12 | — |
| 1.9 | `norm2(a+b) - (norm2(a)+norm2(b)+2*dot(a,b)) == 0` (law of cosines) | 1e-12 | exactly 0 |
| 1.13 | `norm2(b) == dot(u,b)^2 + norm2(cross(u,b))`, `u = normalized(a)` | 1e-12 | exactly 0 |
| 1.14 | `norm(a+b) <= norm(a)+norm(b)` | exact `<=` | 7.5 ≤ 8.24166 |
| 1.15 | `norm(cross(a,b)) == norm(a)*norm(b)*sin(angle(a,b))` | 1e-12 | — |
| 1.18 | `0.5*norm(cross(b-a, c-a))` equals `0.5*norm(cross(b,c)+cross(c,a)+cross(a,b))` (problem 1.20's form) | 1e-12 | — |
| 1.21 | `scalar_triple` cyclic under a→b→c→a | 1e-12 | exactly 0, both rotations |
| 1.22 | unit vectors at θ₁=0.7, θ₂=1.9: `dot(e1,e2) == cos(θ₂-θ₁)`, `cross(e1,e2).z == sin(θ₂-θ₁)` | 1e-15 | 5.6e-17, 1.1e-16 |
| 1.33 | `vector_triple(a,b,c) == dot(a,c)*b - dot(a,b)*c` | 1e-12 | exactly 0 |
| 1.8(b), 1.17(b) | product rules by **central finite difference at h = 1e-5** on smooth `r(t) = (cos t, 2t, e^{0.3t})`, `s(t) = (t², sin 2t, 1-t)`, swept over t ∈ [0.1, 1.9] | **1e-8** | dot-rule 5.8e-10, cross-rule 6.2e-10 |

Two things about that last row. The FD tolerance is step-size-dependent and has a floor: at h = 1e-4 the error is 6.6e-8 (truncation-dominated, so 1e-8 would fail); at h = 1e-6 it is 2.7e-9 (roundoff-dominated). h = 1e-5 is the sweet spot — write that reasoning into the test as a comment so nobody "tightens" it later. And the test genuinely catches the operand-order bug: the wrong order `cross(s_dot,r) + cross(s,r_dot)` differs from the right answer by **16.36** at t ≈ 0.95, nine orders above the tolerance.

**(f) Example program.** None needed — this is test-only. (A three-line `angle` call belongs in Proposal 4's example instead.)

**(g)** S, 2–3 h.

**(h) Learning.** Physics: that the two definitions of the dot product (geometric and component) being the *same thing* is a provable claim, not a convention; that `a × (b × c) ≠ (a × b) × c`; that a finite-difference check of a product rule is a real test of your code. C++: `constexpr` vs `inline` and why `std::atan2`'s non-`constexpr`-ness in C++17 forces the distinction (the same trade-off `norm()` already makes); numerical conditioning of `acos` vs `atan2`; how the harness's `NT_CHECK(...)` variadic macro lets you write `NT_CHECK(v == Vec3{1,2,3})` with top-level commas inside braces — read the comment at `tests/test_util.hpp:87-94` before you need it.

---

## 2. `analytic::position_const_accel` — ch. 1's most basic result, currently uncallable

**(a)** A general `r(t) = r₀ + v₀t + ½at²` in `newton::analytic`, so the sharpest integrator test in the library stops being an inline formula inside a test file.

**(b) Physics (§1.6).**
```
F constant  =>  a = F/m constant  =>
    r(t) = r0 + v0 t + (1/2) a t^2          degree-2 polynomial in t
    v(t) = v0 + a t
```
This is the whole content of §1.6 when `F` does not depend on `r`, `v`, or `t`. `analytic/projectile.hpp` only covers the special case `r₀ = 0`, `a = (0,-g,0)`.

**(c) Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/analytic/kinematics.hpp`
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — add **`#include <newton/analytic/kinematics.hpp>`** as the first line of the analytic block (immediately before line 36's `projectile.hpp`)
- modify `/Users/jackery/Github/newtonian-library/tests/test_integrators.cpp` (replace the inline `r0 + v0*T + 0.5*T*T*a` in block (a) with the library call, then strengthen the assertions)
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — new bullet in the Agent A3 / `newton::analytic` section

**(d) Signatures.**
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

**(e) Closed form and assertions.** Fixture: `m = 2.5 kg`, `F = (2, -7, 0.5) N`, `r₀ = (1,2,-3) m`, `v₀ = (4,-5,6) m/s`, `dt = 0.1 s`, `T = 5 s` (so `|a| = 2.9189 m/s²`, `|r(T)| = 72.08 m`).

| Method | Assertion | Measured |
|---|---|---|
| RK4 | `close_vec3(final.r, position_const_accel(T,...), 1e-12)` | `|Δr|` = 4.9e-14 m, `|Δv|` = 2.2e-14 |
| VelocityVerlet | same, 1e-12 | `|Δr|` = 7.8e-14 m, `|Δv|` = 4.2e-14 |
| Euler | `NT_CLOSE_VEC3(final.r - exact, -0.5*a*T*dt, 1e-9)` and `NT_CHECK(norm(final.r - exact) > 0.1)` | error is exactly `(-0.2, 0.7, -0.05)`; matches `-½aTdt` to 6.0e-14; `|err|` = 0.729726 m |
| SemiImplicitEuler | `NT_CLOSE_VEC3(final.r - exact, +0.5*a*T*dt, 1e-9)` | error exactly `(+0.2, -0.7, +0.05)`; matches to 6.6e-14 |

The Euler rows are better than "Euler is not exact": the two first-order methods **bracket** the exact answer with errors that are equal and opposite and exactly `∓½ a T dt`, because Euler drifts the position with the *old* velocity and Euler–Cromer with the *new* one. Both get velocity exactly right (2.2e-14) because `a` is constant. Assert all of that. Also keep `DESIGN.md`'s existing check (d): final `state.t == t_end` to 1e-12.

**(f) Example.** None — fold a one-line use into `examples/07_incline.cpp` (Proposal 3).

**(g)** S, 1–2 h.

**(h) Learning.** Physics: why "exact for constant acceleration" is the *first* thing to test in any integrator, and why it discriminates first-order from second-order methods with no convergence study at all. C++: `constexpr` functions on an aggregate type; why a test that re-derives a formula inline is a test of two things at once and therefore a worse test.

---

## 3. `examples/07_incline.cpp` — tilted axes and the normal force (1.37–1.39)

**(a)** The §1.6 lesson that you may *choose* your axes, worked three ways against three closed forms — and zero new library code.

**(b) Physics (§1.6).** With one axis along a slope at angle φ above horizontal and one normal to it, gravity resolves as
```
W = (-m g sin(phi), -m g cos(phi), 0)          tilted components of m g

puck ON the slope (problem 1.37): the surface supplies N = m g cos(phi),
   which exactly cancels the normal component (a_normal = 0 is a KINEMATIC
   condition, not a force law — that is how constraint forces are found),
   leaving  m x_ddot = -m g sin(phi),   t_return = 2 v0 / (g sin phi)

ball IN FLIGHT above the slope (problem 1.39): there is NO normal force —
   the ball is in the air. The tilted axes are pure bookkeeping.
   R     = 2 v0^2 sin(theta) cos(theta + phi) / (g cos^2(phi))
   R_max = v0^2 / [g (1 + sin(phi))]            at  2 theta + phi = pi/2
```
The "no normal force in flight" trap is worth a loud comment; it is the entire subject of the public confusion around problem 1.39.

**(c) Files.** Create `/Users/jackery/Github/newtonian-library/examples/07_incline.cpp`. **No header changes, no umbrella line, no Makefile edit** — `EXAMPLE_SRCS := $(wildcard examples/*.cpp)` (Makefile line 35) picks it up automatically, and `HEADERS := $(shell find include -name '*.hpp')` is a blanket prerequisite.

**(d) Signatures.** None public. Inside the example, the force is a lambda satisfying the existing `Force` concept:
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

**(e) Closed forms and assertions** (`v0 = 20 m/s`, `g = 9.80665`, RK4 `dt = 1e-4`, `simulate_until`):

| Case | Analytic | Measured numeric | Assert |
|---|---|---|---|
| 1.39, θ = 25° above a φ = 20° slope | R = 27.607745329 m | 27.607745329 | rel 5.6e-13 → tol 1e-9 |
| 1.39, θ = θ_opt = 35° (so 2θ+φ = 90°) | R_max = 30.3934696673 m | 30.3934696674 | rel 8.4e-13 → tol 1e-9 |
| degenerate φ = 0 collapse | `v0² sin2θ/g` = 38.3287920252 m | 38.3287920252 | rel 4.9e-13 → tol 1e-9 |
| 1.37, puck kicked up the slope | t_return = 2v₀/(g sinφ) = 11.9258030017 s | 11.9258030017 | rel 9.2e-13 → tol 1e-9 |
| 1.37, constraint force | N/m = g cos φ = 9.21523664 m/s² | algebraic | 1e-12 |

The φ = 0 collapse is the best test in the set: it checks that your tilted-frame bookkeeping reduces to the already-trusted `analytic::range`. These belong in `tests/test_numeric_vs_analytic.cpp` too, which already has the `simulate_until`-range-to-1e-9 idiom.

**(f) Example sketch.** Prints three labelled blocks in the style of `examples/01_free_fall.cpp`'s `row()` helper: (1) puck up the slope — numeric vs analytic return time and the normal force, with a note that N came from `a_y = 0`, not from a force law; (2) ball thrown at 25° above a 20° slope — numeric vs analytic range; (3) a sweep of θ from 5° to 85° in 1° steps, printing the maximising θ and comparing with `θ_opt = (90° - φ)/2` and `R_max`, written to `build/incline_range_sweep.csv` via the generic `io::write_csv(path, {"theta_deg","range_m","analytic_m"}, rows)` overload. Plot with `python3 scripts/plot_csv.py build/incline_range_sweep.csv`.

**(g)** S, 2–3 h.

**(h) Learning.** Physics: that a free-body diagram's axis choice is a *decision* with consequences; that a constraint force is determined by a kinematic condition rather than a force law; that a normal force exists only while there is contact. C++: lambdas as `Force`-concept models and why template callables (not `std::function`) are the library's contract — `notes/01-cpp-for-this-library.md:219` is the note to read; `const` lambda captures by value; unnamed unused parameters in a lambda to keep `-Wunused-parameter` quiet.

---

## 4. `newton::coords` — the polar module (§1.7). The largest single ch. 1 gap

**(a)** Cartesian ↔ polar, the rotating basis `r̂`/`φ̂` and its time derivatives, and the velocity/acceleration decomposition — plus the test that a Cartesian `simulate()` run actually reproduces `a_r = r̈ - rφ̇²`.

**(b) Physics (§1.7).**
```
x = r cos(phi),  y = r sin(phi)                             (1.37)
r = sqrt(x^2 + y^2),  phi = atan2(y, x)       <-- Taylor's arctan(y/x) is
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
The whole reason the acceleration has "extra" terms is that the basis rotates with the particle — `d r̂/dt ≠ 0`. Those terms are **kinematics, not forces**. Do not let them become forces; that is ch. 9.

**(c) Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/coords/polar.hpp`
- create `/Users/jackery/Github/newtonian-library/tests/test_coords.cpp`
- create `/Users/jackery/Github/newtonian-library/examples/08_polar_coordinates.cpp`
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — new block right after the core block (after line 21), since `coords` depends only on `Vec3`:
  ```cpp
  // coords — coordinate systems other than Cartesian (Taylor §1.7)
  #include <newton/coords/polar.hpp>
  ```
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — **hard rule 4 (line 29-31) must gain `newton::coords` to its namespace list**, plus a new file-by-file section. This is correcting a spec omission, not expanding scope: `docs/research/taylor-ch1-7-map.md` §1.4 already specifies the module.
- modify `/Users/jackery/Github/newtonian-library/notes/02-physics-to-code.md` — its ch. 1 "Key equations" block stops at `F = dp/dt = m a`; §1.7 contributes nothing today.

**(d) Signatures.**
```cpp
namespace newton::coords {

/// Plane-polar coordinates (Taylor §1.7).
/// C++ note: deliberately NOT a Vec3. Vec3's entire API — dot, cross, norm,
/// normalized — is Cartesian, so a (r, phi) pair wearing a Vec3 would compile
/// into silently meaningless calls (what is norm({r, phi, 0})?). The library's
/// precedent for a named non-Cartesian tuple is analytic::SHM and
/// systems::RocketSample. Small aggregate, zero-initialized members.
struct Polar { double r{0.0};    ///< radius [m]
               double phi{0.0}; };  ///< polar angle [rad], in (-pi, pi]

/// Cartesian -> plane polar, taking x and y and ignoring z (Taylor eq. 1.37).
/// Physics: phi comes from std::atan2(y, x), not std::atan(y/x). Taylor's own
/// problem 1.42 is exactly this point: arctan(y/x) cannot tell (-1,-1) from
/// (1,1), because the ratio is the same. atan2 takes the two components
/// separately and resolves the quadrant, returning (-pi, pi].
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

**(e) Closed forms and assertions** in `tests/test_coords.cpp`:

1. **Round trip.** `to_polar(from_polar({r, phi}))` for r ∈ {0.5, 1, 1000}, φ swept over (-π, π] → 1e-12 relative. Then the branch cut, asserted on purpose: `to_polar(from_polar({1, -pi})).phi == +pi` — `atan2` returns (-π, π], so -π is not a fixed point. Document it rather than hiding it.
2. **Problem 1.42, the quadrant bug.** `to_polar({-1,-1,0}).phi == -3*pi/4` to 1e-15, *and* `NT_CHECK(std::atan(-1.0/-1.0) != to_polar({-1,-1,0}).phi)` — the naive formula gives +π/4, off by π. This is the memorable test.
3. **Basis geometry.** `dot(r_hat(phi), phi_hat(phi)) == 0`, both norms 1, `cross(r_hat, phi_hat) == {0,0,1}` — swept over φ, tol 1e-15.
4. **Basis derivatives vs finite difference** (eqs. 1.42, 1.46), h = 1e-5, φ̇ = 2.7, φ swept over [-3, 3]: **tol 1e-8**; measured 3.3e-10 for both (at h = 1e-6 it improves to 6.6e-11, so the tolerance has headroom either way).
5. **Uniform circular motion, analytic, no integrator.** R = 2 m, ω = 3 rad/s: feed `polar_acceleration(R, 0, 0, ω, 0)` and compare to `dot(-ω²r, r̂)` and `dot(-ω²r, φ̂)`. Tol 1e-12; measured `a_r` err **3.6e-15**, `a_φ` err **1.7e-15**, `φ̇` and velocity reconstruction **2.2e-15**.
6. **The headline cross-check: a Cartesian simulation must reproduce the polar decomposition.** Fixture is ch.-1-only — a vacuum projectile under `forces::UniformGravity`, v₀ = 50 m/s at 70°, launched from the origin, polar origin at the launch point, RK4, `dt = 1e-4`, `record_every = 1`, t_end = 9 s. For each sample: compute `r`, `φ`, and exactly `ṙ`, `φ̇` via `polar_velocity`; get `r̈` and `φ̈` by **central finite difference of those two sequences**; assert `polar_acceleration(...)` matches `dot(g, r̂)` and `dot(g, φ̂)`. Skip samples with `r < 1 m` (near the polar origin the angular derivatives blow up, and φ is meaningless at r = 0). **Tolerance 1e-7 absolute in m/s²**; measured over 89 799 usable samples: worst `|a_r|` error **3.1e-9**, worst `|a_φ|` error **8.1e-9** (3.1e-10 relative to g).
   The tolerance is tied to `dt` and you must say so in the test: the same run at `dt = 1e-3` gives **3.1e-7** — exactly 100× worse, the O(dt²) signature of a central difference. Do not copy this test and change `dt` without changing the tolerance.
7. **Velocity reconstruction along that same trajectory is FD-free**, so it holds to roundoff: `from_polar_components(φ, ṙ, r*φ̇) == s.v` to 1e-12 relative.
8. **Cylindrical** round trip + `to_cylindrical({x,y,z}).rho == to_polar({x,y,0}).r` to 1e-15.

**(f) Example sketch** — `examples/08_polar_coordinates.cpp`:
- Block 1: the projectile of test 6. Prints a 10-row table of `t, r, phi(deg), r_dot, phi_dot, a_r, a_phi, resid_r, resid_phi` and the max residuals over the whole run, with the O(dt²) note printed beside them. Writes `build/polar_projectile.csv` with those columns via the generic `io::write_csv` overload (the `Trajectory` overload's fixed 7-column header cannot carry polar columns).
- Block 2: **Taylor problem 1.40** as a live demo. Sweep the launch angle and find the largest θ for which `d|r|/dt` (use `dot(r, v)`, same sign, no division) never dips negative during the flight. Bisect on θ to 40 iterations. Measured with `UniformGravity`, v₀ = 20 m/s, RK4 `dt = 1e-5`: **θ_max = 70.5287794°**, against `arcsin(sqrt(8/9))` = 70.5287794° — **7.0e-13 relative**. This is the most satisfying number in the chapter: a non-obvious closed form your simulator finds by search.
- Block 3: `r̂`/`φ̂` printed at four angles, with `d r̂/dt` next to `φ̇ φ̂`, so the reader sees eq. 1.42 as numbers.

**(g)** M, 5–8 h.

**(h) Learning.** Physics: why a rotating basis produces acceleration terms that are *not* forces — the single most common conceptual error in the rest of the course; why Cartesian is the right internal representation and polar the right *reporting* representation; that `atan2` is a physics bug fix, not a C detail. C++: designing a small aggregate struct instead of overloading a vector type with a second meaning (and *why* — `norm()` of a polar pair is nonsense); `constexpr` where arithmetic allows and `inline` where `<cmath>` forbids; zero-guarding a function whose mathematical domain excludes a point, the way `normalized()` does; `std::hypot` over `sqrt(x*x+y*y)` for overflow safety; the harness gap — `NT_CLOSE_VEC3` is templated on `.x/.y/.z` (`tests/test_util.hpp:51`), so `PolarAccel` needs two `NT_CLOSE` calls, and bind an `expected` local first because the 3-argument macro splits on top-level commas.

---

## 5. `analytic::UniformCircle` — the circular-motion oracle ch. 1 needs

**(a)** A ch.-1-framed closed form for uniform circular motion, so `a = v²/R`, `T = 2πR/v` and problem 1.41's `T_string = mRω²` stop being reachable only through the §2.5 cyclotron functions.

**(b) Physics (§1.7, problems 1.10, 1.41, 1.45).**
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

**(c) Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/analytic/circular.hpp`
- create `/Users/jackery/Github/newtonian-library/examples/09_circular_motion.cpp`
- modify `/Users/jackery/Github/newtonian-library/tests/test_numeric_vs_analytic.cpp` (new ch.-1 block) and `/Users/jackery/Github/newtonian-library/tests/test_analytic.cpp` (internal consistency)
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — add **`#include <newton/analytic/circular.hpp>`** to the analytic block, after `kinematics.hpp`
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — bullet in the `newton::analytic` section

**(d) Signatures.**
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
Optional add-on, 20 minutes: problems 1.11/1.12's prescribed ellipse and elliptical helix, `struct EllipticalPath { double b, c, omega, vz; ... }`, which exercise `coords::to_polar` on a non-circular path and satisfy `x²/b² + y²/c² = 1` to roundoff.

**(e) Closed forms and assertions.**

Internal consistency (`test_analytic.cpp`, no integrator): `centripetal_accel() == speed()²/radius` and `== centripetal_accel_from_speed(speed(), radius)`; `period() == circular_period_from_speed(speed(), radius)`; `dot(position(t), velocity(t)) == 0` and `norm(position(t)) == radius` over a swept t — all 1e-12.

Numeric-vs-analytic (`test_numeric_vs_analytic.cpp`), R = 3 m, ω = 2 rad/s, m = 1.5 kg. The force is a **string**: constant magnitude `T = m R ω² = 18 N` directed at the centre, `F = -T * normalized(r)` — problem 1.41's actual force, not a spring. One full period, RK4, `dt = 1e-4`:

| Assertion | Tolerance | Measured |
|---|---|---|
| `norm(s.r) == R` every sample (the constraint holds) | 1e-11 | 3.7e-14 |
| `norm(s.v) == analytic.speed()` every sample | 1e-11 | 5.0e-14 |
| `dot(s.r, s.v) == 0` every sample (**problem 1.45**) | 1e-11 | 2.0e-13 |
| `norm(a) == norm(v)²/norm(r)` every sample | 1e-11 | 2.3e-13 |
| period: `2 * simulate_until(..., events::turning_point_x()).back().t` vs `2π/ω` = 3.14159265358979 | 1e-12 rel | **3.0e-15** rel |
| tension magnitude `== analytic.central_force(m)` = 18 N | 1e-12 | algebraic |

Note the period trick, which needs Proposal 6's tiny event addition: starting at `(R,0,0)` with `v = (0,Rω,0)`, `v.x` starts at 0 and next vanishes a **half** period later, so the measured full period is `2 × t_event`. (With `forces::sum`-style linear central force `-mω²r` instead of the string, the same run gives 4.5e-14 / 9.2e-14 / 4.5e-13 / 5.4e-13 — indistinguishably good, because for exact circular initial conditions the two force laws agree. Use the string; it is the one ch. 1 poses.)

**(f) Example sketch** — `examples/09_circular_motion.cpp`: problem 1.41's astronaut. Prints R, ω, speed, centripetal acceleration, period, and tension, then the numeric/analytic pair for each, plus the two running invariants (`max ||r|-R|` and `max |r·v|`) over one period with a line explaining that `r·v = 0` *is* problem 1.45's theorem. Writes `build/circle.csv` (the `Trajectory` overload works here — it is ordinary Cartesian motion) and a second generic CSV of `t, r, phi_deg, r_dot, phi_dot` from `coords::polar_velocity`, demonstrating that `ṙ = 0` and `φ̇ = ω` to roundoff — which is what makes the polar decomposition collapse to `a_r = -Rω²`.

**(g)** S–M, 3–4 h (depends on Proposal 4 for the polar read-out and Proposal 6's event for the period).

**(h) Learning.** Physics: that constant *speed* is not constant *velocity*; that the force holding a circle is whatever the constraint supplies, and its magnitude is fixed by kinematics; problem 1.45's `|r|` const ⟺ `r ⊥ v` as a runnable invariant rather than an exercise. C++: a struct of pure accessors as an analytic oracle (the `analytic::SHM` pattern without any of SHM's ch.-5 content); guarding a divide inside a `const` member function; `[[nodiscard]]` on accessors and why discarding `period()` is always a bug.

---

## 6. The half-pipe: `systems::CircularTrack` + `events::turning_point_x` — ch. 1's only `[Computer]` problems

**(a)** Integrate `φ̈ = -(g/R) sin φ` from rest at 20° and at 90°, plot against the small-angle curve, and read the normal force off the radial equation — Taylor problems 1.50 and 1.51, exactly as posed.

**(b) Physics (§1.7, Example 1.2).** A point mass on a rigid frictionless circular track of radius R, φ measured from the lowest point. Polar coordinates centred on the centre of curvature: `r = R` fixed, so `ṙ = r̈ = 0`, and the two components of eq. (1.48) become
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

**(c) Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/systems/circular_track.hpp`
- create `/Users/jackery/Github/newtonian-library/examples/10_half_pipe.cpp`
- create `/Users/jackery/Github/newtonian-library/tests/test_circular_track.cpp`
- modify `/Users/jackery/Github/newtonian-library/include/newton/simulate.hpp` — one new event in the existing `events` namespace (beside `ground_impact` and `apex`, lines 222–231)
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — add **`#include <newton/systems/circular_track.hpp>`** to the systems block (beside line 48's `rocket.hpp`)
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — a `newton::systems` file-by-file entry and the new `events` function in the `simulate.hpp` spec block (line 222–226)

No new namespace: `newton::systems` is already on rule 4's list, and "self-contained special-purpose simulator" is exactly what this is.

**(d) Signatures.**
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
/// store phi [rad] in r.x, phi_dot [rad/s] in v.x, the GENERALIZED force
/// (a torque, -m g R sin phi [N m]) in the returned x component, and the
/// GENERALIZED mass (m R^2 [kg m^2]) in the integrator's `mass` slot. Then
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

    /// Generalized force (torque) -m g R sin(phi) [N m], in the x slot.
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

**(e) Closed forms, method, and assertions.** Taylor's numbers: `R = 5 m`, `g = 9.8 m/s²`, so `ω = sqrt(9.8/5) = 1.4 rad/s` exactly and `T_lin = 2π/1.4 = 4.48798951 s`; `φ̈ = -1.96 sin φ`.

**The method for every turning-point quantity — this is not optional.** Chain `simulate_until` legs:

```cpp
const systems::CircularTrack tk{5.0, 9.8, 1.0};
Particle q = tk.particle(phi0, 0.0);
double t0 = 0.0;
for (int leg = 0; leg < 12; ++leg) {            // 12 half-periods = 6 periods
    SimOptions o;  o.t0 = t0;  o.t_end = t0 + 3.2;  o.dt = 1e-3;  o.method = Method::RK4;
    const Trajectory tr = simulate_until(q, tk, o, events::turning_point_x());
    const State& tp = tr.back();
    NT_CHECK(std::fabs(systems::CircularTrack::phi_dot(tp)) < 1e-10);  // guard: the event really fired
    /* ... record phi and N at tp ... */
    q  = tk.particle(systems::CircularTrack::phi(tp), systems::CircularTrack::phi_dot(tp));
    t0 = tp.t;
}
```
Three details that make this correct. (i) `simulate_until` refines the crossing by bisecting on the **step size**, re-stepping from the previous state ~40 times (`simulate.hpp:185-202`), so the returned state sits on `φ̇ = 0` to the integrator's own accuracy, not to O(dt). (ii) It returns on the **first** crossing, so each run yields exactly **one** turning point — chaining legs is how you get a series. (iii) At `t0` the state has `φ̇ = 0` exactly, so `e_prev == 0`; the strict `e_prev*e_next < 0` test (`simulate.hpp:185`, and the comment above it) correctly declines to call that a crossing, so the run does not terminate at once. The `|φ̇| < 1e-10` guard catches the one failure mode — a leg whose `t_end` was too short for the event to fire, which would silently return the `t_end` state instead.

Measured with that method, `dt = 1e-3`, RK4, 12 legs:

| Quantity | φ₀ = 20° | φ₀ = 90° | Assert |
|---|---|---|---|
| max \| \|φ_turning\| − φ₀ \| over 6 periods | **7.4e-15 rad** | **2.1e-14 rad** | 1e-9 |
| max \| N/m − g cos φ₀ \| at turning points | **2.4e-14** | **2.1e-13** | 1e-9 |
| reference `g cos φ₀` | 9.20898768 m/s² | 6.0e-16, i.e. **0 to 1e-9** (not "exactly 0") | — |
| measured full period (2 × first leg) | **4.52240801 s** | **5.29735622 s** | see below |
| ratio to `T_lin` | 1.00766903 | 1.1803406 | see below |

**Pin `dt ≤ 1e-3` in the test and say why.** Same method at `dt = 1e-2`: amplitude drift 3.2e-11 / 1.4e-10 and normal-force error 1.1e-10 / **1.4e-9** — the 90° case would fail a 1e-9 assertion. At `dt = 1e-4` nothing improves (6.4e-15 / 1.2e-13): `dt = 1e-3` is already at the roundoff floor.

**Do not scan recorded samples for a `φ̇` sign change as a route to any of these assertions.** Measured with the identical run, taking the nearest recorded sample to each crossing: amplitude error up to **6.7e-8** (20°) and **2.1e-7** (90°), normal-force error up to **6.8e-7** and **6.1e-6** — three to four orders outside 1e-9, because the nearest sample misses the turning point by O(φ̈ dt²/8) in φ and O(φ̈ dt/2) in φ̇, and `N = m(g cos φ + R φ̇²)` is first-order in that φ̇ error. If you ever need a turning-point series without chaining legs, *scan then refine*: re-step from the sample **before** the crossing with a bisected sub-step. If you keep a raw scan, the assertion must be loosened to ~1e-5 with the O(dt²) quantization floor written in the comment.

**Period assertions, ch.-1-safe.** Do not assert an amplitude-correction formula — that is not in ch. 1. Assert instead: (i) the measured period **exceeds** `T_lin` in both cases, which follows from ch. 1 alone since `sin φ < φ` makes the true restoring torque weaker than the linear one; (ii) at 20° the measured period is within 1% of `T_lin` (measured +0.767%); (iii) at 90° it is more than 15% above `T_lin` (measured +18.03%). That trio captures "the small-angle limit is good at 20° and bad at 90°" without importing ch. 5.

**Small-angle curve comparison — assert the failure as well as the success.** Over 3 `T_lin` at `dt = 1e-3`, max `|φ_num(t) − φ₀ cos(ωt)|`:
- φ₀ = 20° (= 0.3491 rad): **0.0461 rad (2.64°)**, i.e. 13% of the amplitude → assert `< 0.05`;
- φ₀ = 90° (= 1.5708 rad): **3.112 rad (178°)** — the two curves end up essentially in antiphase → assert `> 1.0`. Problem 1.51's expected answer *is* the divergence, so the test must demand it, not tolerate it.

**Integrator quality, energy-free.** Same refined-turning-point drift metric, `dt = 1e-3`, 6 periods:

| Method | drift at 20° | drift at 90° |
|---|---|---|
| `Method::Euler` | **6.97e-3 rad** (2% of φ₀) | **3.5e-2 rad** |
| `Method::SemiImplicitEuler` | 4.4e-7 rad | 1.5e-6 rad |
| `Method::RK4` | 7.4e-15 rad | 2.1e-14 rad |

Assert Euler drifts more than 1e-3 rad and RK4 less than 1e-9 rad. This is the ch.-1-licensed replacement for `examples/06_integrator_comparison.cpp`'s energy-drift chart: the turning points of a periodic motion must not wander, and you do not need the word "energy" to say so.

**(f) Example sketch** — `examples/10_half_pipe.cpp`, reproducing both `[Computer]` problems:
- Header comment: the physical situation, §1.7 / Example 1.2, the analytic expectation (none — that is the point), and what to look at.
- Problem 1.50: φ₀ = 20°, RK4, `dt = 1e-3`, `t_end = 14 s` (just over three periods, matching the window the problem asks for). Problem 1.51: φ₀ = π/2, `t_end = 21 s`.
- Prints, per amplitude: measured period, `T_lin = 4.48799 s`, the ratio, the max deviation from `φ₀ cos(ωt)`, the turning-point drift for Euler / Euler–Cromer / RK4, and `N/m` at a turning point next to `g cos φ₀`.
- Prints `N/m` at the lowest point too — measured **10.9820246 m/s² at φ₀ = 20°** and **29.4 m/s² at 90°** — labelled honestly: *you cannot predict these two numbers yet. Ch. 4's energy methods will let you; write them down now and check them then.*
- `build/half_pipe.csv` via the generic overload, columns `t, phi20, small20, phi90, small90` on a common 0–21 s grid, so one `python3 scripts/plot_csv.py build/half_pipe.csv t phi20 small20` gives problem 1.50's figure and swapping the columns gives 1.51's.

**(g)** M, 6–9 h (depends on Proposal 4 for the polar reasoning, though not on its code).

**(h) Learning.** Physics: that a constraint turns a 2-D vector problem into one scalar ODE; that the "unused" radial equation is where the constraint force comes from; that a nonlinear equation can be perfectly well-behaved and still have no closed form; exactly how far a linearisation can be trusted, measured rather than asserted. C++: reusing a general-purpose integrator for a generalized coordinate by reinterpreting a struct's slots — and the discipline that makes such a trick safe rather than a landmine (named `pack`/`phi`/`phi_dot` accessors, a doc comment that says what the slots mean, a rule that these States never reach Cartesian code); writing a lambda-returning factory in the established `events` idiom (`[[nodiscard]] inline auto name(args) { return [captures](const State&){...}; }`); chaining stateful calls where each output seeds the next input; defensive post-conditions on a function that can silently return the wrong thing.

---

## 7. `newton::frames` — Galilean boosts and the turntable as *observation* (§1.4)

**(a)** Re-express a finished trajectory in another frame, so problems 1.26, 1.27 and 1.46 become three lines of post-processing — the one ch. 1 concept missing from both the code and the roadmap.

**(b) Physics (§1.4, §1.5).**
```
Newton I defines the inertial frames: F = 0 => v = const.
A frame moving at constant V relative to an inertial frame is inertial:
   r' = r - V t,   v' = v - V,   a' = a        =>  F = m a is Galilean invariant
A uniformly accelerating frame is NOT inertial:
   r'' = r - (1/2) A t^2,   v'' = v - A t,   a'' = a - A
A rotating frame is NOT inertial. What the turntable observer SEES:
   r' = Rot(-w t) r,    v' = Rot(-w t) (v - w z_hat x r)
   |r'| = |r|, and a puck crossing at constant v0 appears to move at
   sqrt(v0^2 + rho^2 w^2) at radius rho.
```
Chapter 1 never adds a centrifugal or Coriolis **force**; it only asks what the moving observer measures. Fictitious forces are ch. 9. This module transforms *recorded output* and must never touch an equation of motion.

**(c) Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/frames/galilean.hpp`
- create `/Users/jackery/Github/newtonian-library/tests/test_frames.cpp`
- create `/Users/jackery/Github/newtonian-library/examples/11_frames.cpp`
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — new block after the coords block:
  ```cpp
  // frames — re-expressing a recorded trajectory in another frame (Taylor §1.4)
  #include <newton/frames/galilean.hpp>
  ```
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — **hard rule 4's namespace list gains `newton::frames`** (bundle this with Proposal 4's `newton::coords` amendment if you do them together), plus a file-by-file section

**(d) Signatures.**
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

**(e) Closed forms and assertions** in `tests/test_frames.cpp`:

1. **Galilean invariance, the headline.** Simulate `forces::UniformGravity{m}` from `(r₀, v₀) = ((0,10,0), (5,2,0))`, m = 2 kg, RK4 `dt = 1e-3`, 1.5 s. Simulate again from `(r₀, v₀ − V)` with `V = (3,4,-1)`. Then `observe_from(traj1, LinearFrame{{}, V, {}})` must equal `traj2` sample by sample. Tol **1e-11**; measured worst `|Δr| + |Δv|` = **4.0e-13**. This works precisely because `UniformGravity::operator()` ignores `v` — say so in the test comment.
2. **The free-fall frame (problem 1.26's third frame).** `observe_from(traj1, LinearFrame{{}, {}, g})` must turn the parabola into the straight line `r₀ + v₀t`. Tol 1e-11; measured worst deviation **6.0e-13**. A parabola becoming a line is the most vivid possible statement of "this frame is not inertial".
3. **The turntable (problems 1.27, 1.46).** Puck crossing the centre at constant velocity, `v₀ = 2 m/s` along ŷ, `ω = 1.3 rad/s`, t ∈ [-2, 2]: a rotation preserves length, so `norm(observed.r) == norm(s.r)` to **1e-14** (measured **4.4e-16**), and the observed speed equals `sqrt(v₀² + ρ²ω²)` with `ρ = |r|` to **1e-12** (measured **1.8e-15**) — the analytic result the published solutions give for 1.27. Also assert the observed path is a spiral in and back out: `norm(observed.r)` decreases monotonically to 0 then increases.
4. **`is_inertial()`**: true for `{{}, V, {}}`, false for any nonzero `accel`, false for any nonzero `omega`.

**(f) Example sketch** — `examples/11_frames.cpp`: one kicked puck, three descriptions (problem 1.26). Prints a short table of `t, x, y` in the original frame, the boosted frame, and the uniformly-accelerating frame, with a one-line verdict under each ("inertial: straight line", "inertial: straight line, different slope", "non-inertial: parabola, yet no force acts"). Then the turntable (1.46): prints `rho, phi_ground_deg, phi_rotating_deg, observed_speed, sqrt(v0^2+rho^2 w^2)` at a few times and writes `build/frames.csv` with a `series` tag column (the idiom `examples/02_projectile_drag.cpp` already uses) so all four paths plot together.

**(g)** M, 4–6 h.

**(h) Learning.** Physics: that Newton I is a *definition of a class of frames*, not a special case of Newton II; that acceleration is the Galilean invariant and therefore `F = ma` is frame-independent among inertial frames, while position and velocity are not; that a non-inertial frame shows you motion without force, which is the entire motivation for ch. 9. C++: overloading a free function on a tag-like parameter type as an extension point; `constexpr bool` predicates on an aggregate; transforming a container into a new one without mutating the input (`[[nodiscard]]`, value semantics, and why `Trajectory` being cheap to copy makes this honest); writing a doc comment whose job is to record a *limitation* so that a future you does not over-trust the function.

---

## 8. `examples/12_atwood.cpp` — one constraint, two bodies, Newton III for the tension

**(a)** The Atwood machine from Newton II per mass plus the third law plus `a₁ = -a₂`, validated against two exact formulas — with no new library code.

**(b) Physics (§1.4, §1.5).** Two masses on a light string over a frictionless pulley. Newton II on each, the third law for the string tension (the same `T` pulls up on both), and the constraint that the string's length is fixed (`a₁ = -a₂`):
```
m1 a = m1 g - T
m2 a = T - m2 g           (taking a > 0 as m1 descending)
=>  a = (m1 - m2) g / (m1 + m2)
    T = 2 m1 m2 g / (m1 + m2)

Limits worth asserting:  m2 -> 0  =>  a -> g, T -> 0     (free fall)
                         m1 == m2 =>  a = 0,  T = m g     (balanced)
```
As one degree of freedom: generalized coordinate `x` = how far `m1` has descended, generalized mass `m1 + m2`, generalized force `(m1 - m2) g` — constant, so the motion is constant-acceleration and `analytic::position_const_accel` from Proposal 2 is the exact oracle.

**Honesty note, from INPUT A §3.6:** this is *not* a Taylor ch. 1 problem. A search of all 51 ch.-1 problem statements finds no Atwood machine and no pulley; Taylor's Atwood machine is a ch. 7 Lagrangian showcase. The physics is entirely ch. 1, so build it — but do not tell yourself you are following the book here, and say so in the file header.

**(c) Files.** Create `/Users/jackery/Github/newtonian-library/examples/12_atwood.cpp`, and a small block in `/Users/jackery/Github/newtonian-library/tests/test_numeric_vs_analytic.cpp`. No header changes, no umbrella line.

**(d) Signatures.** None public. In the example:
```cpp
// One degree of freedom: x = descent of m1. Generalized mass m1 + m2,
// generalized force (m1 - m2) g — the two Newton II equations with T
// eliminated by the third law and the constraint a1 = -a2.
const auto atwood = [m1, m2, g](double, const Vec3&, const Vec3&) {
    return Vec3{(m1 - m2) * g, 0.0, 0.0};
};
const Particle dof{m1 + m2, Vec3{}, Vec3{}};
```
The same slot-reinterpretation discipline as Proposal 6 applies, though here it is milder — `x` really is metres.

**(e) Closed forms and assertions.** `m1 = 3 kg`, `m2 = 2 kg`, `g = 9.80665`:

| Assertion | Expected | Measured |
|---|---|---|
| `a == (m1-m2)g/(m1+m2)` | 1.96133 m/s² | — |
| `x(2 s)` numeric vs `analytic::position_const_accel(2, {}, {}, a)` | 3.9226599999999996 m | numeric 3.9226600000000134, abs err **1.4e-14** at RK4 `dt = 1e-2` → tol 1e-12 |
| Newton III consistency: `T` from body 1 (`m1(g-a)`) == `T` from body 2 (`m2(g+a)`) | 23.53596 N both | identical to printed precision → tol 1e-12 |
| `T == 2 m1 m2 g/(m1+m2)` | 23.53596 N | → tol 1e-12 |
| limit `m2 = 0`: `a == g`, `T == 0` | — | tol 1e-12 |
| limit `m1 == m2`: `a == 0`, `T == m g` | — | tol 1e-12 |

The exactness at `dt = 1e-2` is the point: the motion is constant-acceleration, so RK4 is exact and any error you see is a bug in your setup, not truncation.

**(f) Example sketch.** Prints the free-body reasoning as comments, then `a` and `T` three ways (from `m1`, from `m2`, from the closed form) so the third-law consistency is visible as three equal numbers; then the numeric-vs-analytic descent at 2 s; then a small table sweeping `m2/m1` from 0 to 1 showing `a → g` and `a → 0`, written to `build/atwood.csv`.

**(g)** S, 2 h.

**(h) Learning.** Physics: that a constraint is an *equation* you impose and a constraint force is what you solve for; that the third law is what lets one symbol `T` appear in two equations; that eliminating the constraint force gives you a single-DOF problem. C++: reducing two coupled bodies to one `Particle` by choosing a generalized coordinate, and when that is legitimate; reusing an analytic oracle written for a different purpose; sweeping a parameter into a `std::vector<std::vector<double>>` for the generic `io::write_csv`.

---

## 9. `systems/nbody.hpp` — a System whose pairwise forces obey Newton III by construction

**(a)** A multi-body state, a pairwise-force concept that makes `F_ij = -F_ji` structurally true, and momentum conservation to roundoff as the sharpest invariant ch. 1 offers.

**(b) Physics (§1.5).**
```
Newton III, weak form:   F_12 = -F_21
Newton III, strong form: F_12 parallel to (r_1 - r_2)      (central)
P = sum_alpha m_alpha v_alpha
dP/dt = F_external       internal forces cancel pairwise      (1.25)-(1.29)
F_external = 0  =>  P = const
R_cm = (1/M) sum m_alpha r_alpha,   M R_cm_ddot = F_external
Perfectly inelastic collision (problem 1.30):  v' = m1 v / (m1 + m2)
```
Momentum conservation is a *structural* consequence of how you apply a pair force, not a numerical accident — which is why it should hold to roundoff rather than to truncation error, and why it is a sharper test than anything that depends on step size.

**Boundary to document now:** ch. 1 also shows the third law *failing* for magnetic forces between moving charges (problem 1.32), with a relative violation of order `v₁v₂/c²`. So if you ever add a velocity-dependent pair force, the "momentum is conserved" guarantee below is no longer automatic. Write that into the header while everything is still velocity-independent.

**(c) Files.**
- create `/Users/jackery/Github/newtonian-library/include/newton/systems/nbody.hpp`
- create `/Users/jackery/Github/newtonian-library/tests/test_system.cpp`
- create `/Users/jackery/Github/newtonian-library/examples/13_two_body.cpp`
- modify `/Users/jackery/Github/newtonian-library/include/newton/analysis/momentum.hpp` — add `total_mass` and `velocity_of_cm` (its own C++ note already flags their absence)
- modify `/Users/jackery/Github/newtonian-library/include/newton/newton.hpp` — add **`#include <newton/systems/nbody.hpp>`** to the systems block
- modify `/Users/jackery/Github/newtonian-library/DESIGN.md` — **line 15's scope sentence ("single-particle Newtonian mechanics") must change**, both new concepts must be specified beside "The Force concept" (lines 41–51), and a new file-by-file section added. `notes/04-extending.md:167` already names this as the v0.2 headline item; this moves it from roadmap to spec.

**(d) Signatures — including the two new concepts.**
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

**(e) Closed forms and assertions** in `tests/test_system.cpp`. Fixture: two bodies, `m = (1, 3) kg`, `r = ((-1,0,0), (2,0.3,0))`, `v = ((0,0.5,0), (0.1,-0.2,0))`, so `M = 4 kg`, `|P₀| = 0.3162 kg m/s`; `CentralPair{2.0}`; `dt = 1e-3`; 10⁴ steps (10 s); minimum separation during the run 0.58 m (no singularity).

1. **The concepts are enforced at compile time.** A C++17 detection idiom, which is also the verifier's fix made mechanical:
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
   I compiled both `static_assert`s against the real headers: they pass.
2. **Newton III structurally.** With `NoExternalForce`, `total_forces` must sum to zero. For **N = 2** assert **exact** `== Vec3{}` — measured `(0, 0, 0)`, bit-exact, because `F[0] = f` and `F[1] = -f` and IEEE addition of a value and its negation is exact. For **N = 3** assert 1e-15 instead, not exact: three pairs cancel pairwise but the *order* of the additions differs per body, so `(f01+f02) + (-f01+f12) + (-f02-f12)` need not round to zero — measured **2.2e-16**. That contrast is one of the best floating-point lessons in the repo; put it in a comment.
3. **Strong third law (problem 1.34).** `norm(cross(ri - rj, f_ij)) == 0` for `CentralPair` — tol 1e-15, measured **1.3e-16** with `|f| = 2 N`. Geometry only: angular momentum itself is ch. 3.
4. **Isolated system, momentum exactly conserved.** `max |P(t) - P₀|` over 10⁴ steps: tol **1e-12**, measured **2.4e-14**. Then the assertion that makes it meaningful: re-run at `dt = 5e-4` for 2×10⁴ steps and assert the drift is *still* below 1e-12. Truncation error would fall by 4×; roundoff will not — and if your drift instead *grows* with step size, your force application is asymmetric, which is the bug this test exists to catch.
5. **Centre of mass moves uniformly.** `max |R_cm(t) - (R₀ + V_cm t)|` over the same run: tol **1e-12**, measured **2.6e-14**. This is the dynamical half of §1.5 that the existing `center_of_mass` could never demonstrate, because nothing evolved more than one body.
6. **External field — the verifier's assertion, restated with the concrete object.** Same pair force, `body = systems::UniformField{{0, -9.80665, 0}}` (a `BodyForce`, *not* `forces::UniformGravity`). Then
   ```
   P(t)    = P0 + M g t
   R_cm(t) = R0 + V_cm t + ½ g t²
   ```
   Tol **1e-9 absolute** on both. Measured over 10⁴ steps: `max |P - (P₀ + M g t)|` = **1.25e-12** and `max |R_cm - ...|` = **8.5e-13**, with `|P|` reaching ~392 kg m/s by t = 10 s (so 1.25e-12 is ~3e-15 relative — roundoff, as it should be, since the field is constant and Verlet is exact for constant acceleration).
7. **Inelastic collision (problem 1.30), algebra only.** Build the before state (`m1 = 2` at 7 m/s, `m2 = 3` at rest) and the after state (both at `v'`), and assert `total_momentum` is equal: `v' = m1 v/(m1+m2)` = **2.8 m/s**, tol 1e-12. No integration — collisions are not a force law, and ch. 1 treats this as a momentum-bookkeeping problem.

**(f) Example sketch** — `examples/13_two_body.cpp`: the two-body fixture under `CentralPair` alone. Prints `|P(t) - P₀|` and `|R_cm(t) - (R₀ + V_cm t)|` at ten sample times, with a line explaining that these stay at 1e-14 because the cancellation is structural, not because the integrator is accurate. Then the same run with `UniformField` added, printing `P(t) - P₀` beside `M g t`. Writes `build/two_body.csv` through the generic overload with a `body` tag column (0, 1, and 2 for the centre of mass) so all three paths plot together and the CoM's straight line / clean parabola is visible against the two wiggling bodies.

**(g)** L, 12–20 h. This is real design work, correctly scoped as v0.2 in `notes/04-extending.md`. Do the other nine first.

**(h) Learning.** Physics: that the third law is what makes momentum conservation *derivable* rather than empirical; the difference between the weak and strong forms and why the strong one is extra; that only external forces move the centre of mass; that a conservation law can be built into the structure of a program so that violating it requires bypassing the interface. C++: defining a *concept by documentation plus a detection idiom* (`std::void_t`, `std::declval`, partial specialization, `static_assert`) — the C++17 way to get a compile-time error message instead of a physics bug, and the direct answer to "how do I stop someone plugging the wrong callable in here"; two template callable parameters with distinct arities in one function; why `forces::sum`'s fold expression cannot mix them; designing a new aggregate (`SystemState`) to be compatible with existing free functions (`analysis::total_momentum(s.bodies)`) instead of requiring an adapter; floating-point cancellation that *is* exact versus cancellation that merely rounds small, and how to write a test that distinguishes roundoff from truncation by changing `dt`.

---

## 10. `notes/05-chapter-1-problems.md` — the study page that closes the loop

**(a)** A table mapping Taylor problems 1.1–1.51 to the library call, test, or example that answers each — so the library becomes a way of *doing the problem set*, not a parallel activity.

**(b) Physics.** All of ch. 1; no new physics.

**(c) Files.**
- create `/Users/jackery/Github/newtonian-library/notes/05-chapter-1-problems.md`
- modify `/Users/jackery/Github/newtonian-library/notes/02-physics-to-code.md` — its ch. 1 "Key equations" block stops at `F = dp/dt = m a`; add §1.7's polar set, §1.4's Galilean transformations, and §1.5's system momentum
- modify `/Users/jackery/Github/newtonian-library/notes/00-start-here.md` — add the new note to the reading order at line 328
- modify `/Users/jackery/Github/newtonian-library/docs/research/taylor-ch1-7-map.md` — **fix three errors in its ch. 1 section** (lines 38–122), each verified against all 51 problem statements: (i) Atwood's machine is **not** a ch. 1 problem (zero hits for "Atwood" or "pulley" across all 51); (ii) a **conical** pendulum is not in ch. 1 — problem 1.41's twirling mass is explicitly in gravity-free space, so there is no cone; (iii) cylindrical coordinates are **ch. 1** (problems 1.47–1.49), not ch. 4 as the map defers them. Also soften the map's claim that ch. 1 is "where all the integrators should appear": ch. 1 justifies exactly one, and the `[Computer]` tag appears on exactly two problems (1.50, 1.51).

**(d) Signatures.** None. Structure: one table with columns *problem · what it asks (your own words) · library surface · file*. Group by the chapter's own grouping: 1.1–1.25 → §1.2 vectors and ODE warm-ups; 1.26–1.34 → §1.4–1.5 frames, third law, momentum; 1.35–1.40 → §1.6 Cartesian Newton II; 1.41–1.51 → §1.7 polar and cylindrical. Paraphrase every problem in your own words — do not transcribe Taylor's text.

**(e) Validation.** A checklist at the top: every row naming a test must name a test that exists and passes. Run `make check` and paste the summary line.

**(f) Example.** N/A.

**(g)** S, 3–4 h (do it last; it is the index to everything above).

**(h) Learning.** Physics: the chapter as a *structure* rather than 51 unrelated exercises — you will notice that 19 of them are the same vector-identity skill and that only two require a computer. C++/tooling: documentation as a first-class deliverable; keeping a spec (`DESIGN.md`), a manifest (`include/newton/newton.hpp`) and a research doc in sync, and what it feels like when they have silently drifted apart.

---

## Needs later chapters — don't start yet

| Tempting feature | Needs | Why you cannot validate it from ch. 1 |
|---|---|---|
| Air resistance, linear `-bv` or quadratic `-c|v|v`; terminal velocity; the optimal launch angle dropping below 45° | **Ch. 2** §2.1–2.4 | Ch. 1 has no velocity-dependent force at all, and every ch. 1 projectile explicitly neglects air resistance. You would have no analytic reference. Already implemented (`forces/drag.hpp`, `analytic/drag.hpp`, `examples/02_projectile_drag.cpp`) — leave it alone. |
| The "boosting a drag simulation changes the physics" demo | **Ch. 2** | It is the perfect §1.4 lesson and the library is positioned for it, but it needs a force ch. 1 does not have. Put it in Proposal 7's doc comment as a caveat now, and build the demo when you reach ch. 2. |
| Lorentz force; charged particle in a B field; Boris pusher | **Ch. 2** §2.5, §2.7 | Ch. 1 touches E and B only in problems 1.32/1.33, and only to discuss the third law — never to integrate motion. Already implemented. |
| Rocket / variable mass | **Ch. 3** §3.2 | `F = dp/dt` with `ṁ ≠ 0` is past §1.4's constant-mass reduction. Already implemented (`systems/rocket.hpp`). |
| **Angular momentum as a diagnostic** (`analysis::angular_momentum` exists — do not build on it) | **Ch. 3** §3.4 | Problem 1.34's own statement forward-references ch. 3. Proposal 9 asserts the *geometry* of a central pair force, which is ch. 1; `L` conservation is not. |
| Work, kinetic and potential energy, **energy drift as an integrator metric** | **Ch. 4** §4.1–4.4 | The big temptation. "Symplectic integrators conserve energy" is the standard motivation for Verlet, but ch. 1 has no concept of energy, so you cannot state it yet. Use momentum conservation (Proposal 9) and turning-point drift (Proposal 6) instead — both are energy-free and both are sharper. Already implemented (`analysis/energy.hpp`, `examples/06_integrator_comparison.cpp`). |
| The closed form for `N` at the bottom of the half-pipe | **Ch. 4** | Proposal 6 prints the measured value (10.98 and 29.4 m/s²) precisely so you can check it against your own derivation in ch. 4. |
| Damped and driven oscillators, resonance, Q, amplitude-dependent period formulas | **Ch. 5** §5.1–5.6 | Ch. 1 gives undamped SHM only as a small-angle *approximation* to one specific problem. Do not generalise Proposal 6 into an `Oscillator` module, and do not assert the `1 + φ₀²/16` period correction — that is standard pendulum theory, not Taylor ch. 1. Already implemented. |
| A general holonomic constraint solver (Lagrange multipliers, `[M Jᵀ; J 0][a; λ] = [F; -J̇v]`) | **Ch. 6–7** | Proposals 3, 6 and 8 solve their constraints *by hand*, which is what ch. 1 does. A general engine is ch. 7 and is already scoped in `notes/04-extending.md`. |
| Central-force orbits, Kepler, `PointGravity` as a ch. 1 fixture | **Ch. 8** (and §4.8 for the force law) | Problem 1.11's ellipse is a *prescribed kinematic* ellipse, not a dynamical orbit. This is why Proposal 4 uses a uniform-gravity projectile as its polar fixture and Proposal 9 uses a constant-magnitude pair force. |
| **Centrifugal and Coriolis forces** in a rotating frame | **Ch. 9** | Ch. 1's `2ṙφ̇` term is *kinematic*. Problems 1.27 and 1.46 ask only what the rotating observer *sees*, which is exactly what Proposal 7 implements. Adding fictitious forces to an equation of motion now would be jumping eight chapters. |

**One design decision to make now, while it is free.** Velocity Verlet and the symplectic family assume `F = F(r)`. Every ch. 1 force satisfies that; ch. 2's drag and ch. 5's damping do not. Add a trait flag — `static constexpr bool depends_on_velocity = false;` on each force functor, or a `depends_on_velocity()` predicate — *now*, while every answer is `false`, and have the integrator choice consult it. Discovering this in ch. 2 costs far more. Proposal 9's `PairForce`/`BodyForce` concepts are the natural place to state it for systems, where a velocity-dependent pair force also voids the momentum guarantee (problem 1.32).

---

## Photos that would help

The book's text is not accessible online, and several things I relied on are marked "could not confirm" in the research report. A phone photo of each of these would let me tighten or correct the proposals above. Highest value first.

1. **Pages 34–41 — the full problem set 1.1–1.51, including the `[Computer]` tags and the difficulty stars.** Single highest value. Two independent things ride on it: (a) whether `[Computer]` really appears on only 1.50 and 1.51 — the evidence is strong (all 51 third-party transcriptions scanned) but depends on those transcriptions preserving the tag faithfully; if a third problem carries it, Proposal 6 grows. (b) Taylor's star ratings, which no transcription preserves, so the examples above cannot be graded by difficulty at all right now.
2. **Problems 1.50 and 1.51 specifically** (likely on the last page of the set). These are the *only* two problems Proposal 6 claims to reproduce. I want to confirm the plot windows (I have 0–14 s and 0–21 s second-hand), the exact `φ₀` values, and whether Taylor asks for anything beyond `φ(t)` and the overlay.
3. **Pages 26–33 — all of §1.7, including Example 1.2 and equations (1.37)–(1.57).** This is the section Proposals 4 and 6 are built from. I particularly want: the exact form and sign convention of eq. (1.48); whether Taylor measures `φ` from the bottom of the trough or from the vertical (my sign conventions for `N = m(g cos φ + Rφ̇²)` assume the former); and whether Example 1.2 computes the normal force explicitly or only mentions it.
4. **Pages 23–26 — §1.6, including Example 1.1.** I am told, but could not confirm, that Example 1.1 is a block on an incline in §1.6. Proposal 3 is built on that assumption. Does Taylor's version include friction? Does he solve for the normal force? And is eq. (1.35) written as three scalar equations, as one vector equation, or both?
5. **The chapter's table of contents and the "Principal Definitions and Equations of Chapter 1" page (p. 33).** That one page would confirm every section title and page number and every equation-number binding I am relying on second-hand, and it would settle how many numbered examples the chapter actually contains (I can confirm Example 1.2 exists; I do not know whether the chapter stops there).
6. **Pages 17–23 — §1.5, especially the N-particle momentum derivation (1.25)–(1.29) and Figure 1.8's caption.** Proposal 9's doc comments paraphrase that derivation and the magnetic-force caveat. I would like to match Taylor's notation for the internal/external force split rather than invent my own.
7. **Page 9–16 — §1.3 and §1.4.** The only section I have no fetched corroboration for at all is §1.3 (mass, force, superposition). Proposal 9's external/internal force split and Proposal 7's framing of Newton I as a definition of inertial frames would both be firmer with the actual text in view.
8. **Problems 1.47–1.49** (cylindrical coordinates). Proposal 4 adds a cylindrical section largely on the strength of problem 1.49's closed-form helix. If 1.47's notation differs from `(ρ, φ, z)`, the struct should match the book.