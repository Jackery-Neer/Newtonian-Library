# Coverage audit: Taylor ch. 1–4 topics vs the library (2026-11-24)

You asked whether 15 topics from chapters 1–4 are covered in the codebase. Each topic had its own
Opus agent. Each agent:
- read the code, tests, examples and notes;
- ran its own verification programs against hand-derived answers;
- reported claims with file:line evidence.

I then re-checked every claim myself, by opening the cited lines and re-running the agents' programs
(and some probes of my own). Where an agent was wrong or imprecise, the correction is noted below.

Topics 1–11 audit the code as it was committed before this work (commit `bcb6267`). Topics 12–16 cover
the conservative-force, potential-energy and non-conservative-energy work built on 2026-11-24/25. They are
at the bottom, along with what changed for topics 1, 4, 9 and 10.

Verdicts:
- **Covered** means everything a ch. 1–4 student would expect is there and tested.
- **Partial** means the core is there and correct, but pieces are missing.
- **Missing** means not there.

## Summary

| # | Topic | Verdict | What exists | Biggest gaps |
|---|---|---|---|---|
| 1 | Derivatives and integrals of scalars and vectors | Partial (now being filled) | Time integration of Newton II: four integrators, tested to orders 1, 1, 2, 4 | No numeric derivative, quadrature, work or impulse integral before today. The new `fields` module adds gradient/curl/divergence and line integrals, and `work_done` adds ∫F·v dt |
| 2 | Mass and forces in Cartesian, cylindrical, spherical | Partial | Newton II in Cartesian; position conversions to and from cylindrical and spherical (`coords/`) | No unit vectors (ρ̂, φ̂, r̂, θ̂), no projection of velocities or forces onto them, no polar velocity/acceleration (§1.7's F_r = m(r̈ − rφ̇²)) |
| 3 | Momentum | Partial | `momentum(Particle)`, `analysis::total_momentum` | No momentum from a `State`, no `total_mass` or CM velocity, no test of `momentum()` |
| 4 | Air resistance, range, trajectory, b → 0 | Partial, **with a bug** | Linear and quadratic drag forces; every §2.2–2.4 closed form (matches RK4 to ~1e-12); example 02 | **The linear closed forms break as b → 0** (see Bugs). No range-with-drag function, no b → 0 coverage, no f_quad/f_lin helper |
| 5 | Conservation of momentum | Partial | `total_momentum` bookkeeping (problem 1.30 works as before/after algebra) | **Can't simulate interacting bodies**: no pair forces, no N-body `System`, no collisions. This is the largest gap in ch. 1/3 |
| 6 | Rockets | Partial (solid core) | `systems::simulate_rocket` (burn + coast), Tsiolkovsky, closed-form v(t); agrees with hand-derived formulas to ~1e-15 | No altitude closed form, no burnout/max-height helpers, no input validation, Tsiolkovsky untested |
| 7 | Angular momentum, moment of inertia, torque | Partial | `angular_momentum(mass, State, about)`, correct about any point | No `torque`, no total L of several particles, no moment of inertia |
| 8 | Center of mass | Partial | `center_of_mass(vector<Particle>)`, correct; P = MṘ checks out | No continuous bodies (hemisphere 3a/8, cone h/4, …), no `total_mass` / `velocity_of_cm` |
| 9 | Kinetic energy | Partial (being filled) | `kinetic_energy` | No work or power before today. `power()` and `work_done()` are being added now |
| 10 | Energy conservation | Partial (being filled) | `total_energy`, `energy_drift`, `spring_U`, `uniform_gravity_U`; integrator energy behavior matches theory exactly | No dissipation bookkeeping before today (`energy_balance` is being added now). `uniform_gravity_U` assumes −y and is never used in a test |
| 11 | Conservation of angular momentum | Partial | ℓ is conserved by Verlet/semi-implicit Euler to ~1e-14 for central forces (measured) | No test simulates a force and checks ℓ; no Kepler's-second-law helper or test; no total L |

### Updates after the 2026-11-24/25 build

- **Topic 1 now has spatial calculus.**
  - `fields::gradient`, `jacobian`, `curl`, `divergence`: 5-point stencil, about 1e-12 accurate.
  - `fields::line_integral` and `circulation`: adaptive Gauss–Kronrod.
  - `analysis::work_done` (∫F·v dt) and `analysis::power`.
  - Still missing: polar kinematics, impulse, and velocity from samples.
- **Topic 4's b → 0 bug is fixed** (bug 1 below). New functions:
  - `analytic::linear_drag_range`, exact for any b ≥ 0;
  - `analytic::linear_drag_range_small_drag`, Taylor's R ≈ R_vac(1 − 4η/3), with the next term +14/9 η² verified;
  - `analytic::quadratic_drag_horizontal_distance`.

  The fix uses exactly the ln(1 − ε) and ln(1 + ε) series you did in class. Example 02 now prints the analytic range and Taylor's small-drag formula.
- **Topics 9 and 10:** `power`, `work_done` and `energy_balance` (the Δ(T + U) = W_nc ledger) now exist, with tests against closed forms.

## Bugs found (verified by me, not just by the agents)

1. **FIXED 2026-11-25.** The bug was that the linear-drag closed forms failed as b → 0 (`include/newton/analytic/drag.hpp:53-58, 71-76`, before the fix). For
   a 30 m/s projectile at 45° after 3 s, where the vacuum answer is y = 19.51 m:

   | b [kg/s] | `linear_drag_position(...).y` |
   |---|---|
   | 0 | **NaN** (0 × ∞) |
   | 1e-6 | 19.5098 (fine) |
   | 1e-10 | **24,405.8 m** |
   | 1e-14 | **−2.35e12 m** |

   Two things cause it. The formula multiplies τ = m/b (huge) by (1 − e^{−t/τ}) (tiny), which loses
   every significant digit. And `LinearDragParams{}` defaults to b = 0, so a default-built object returns
   NaN.

   **Fix:** rewrite with u = bt/m, φ₁(u) = (1 − e^{−u})/u (computed with `expm1`) and
   φ₂(u) = (u − 1 + e^{−u})/u² (a short series for small u). Then r = r₀ + v₀ t φ₁ + a t² φ₂, which
   becomes the vacuum formula exactly at b = 0.
2. **`simulate` can overshoot `t_end`.** With t_end = 1.0 and dt = 0.4, the run ends at t = 1.2
   (`llround(2.5) = 3` steps). This is undocumented. `dt ≤ 0` throws `length_error` instead of giving a
   clear message.
3. **`quadratic_drag_fall_distance` overflows** (it returns ∞) once g·t/v_ter is past about 710
   (`drag.hpp:127`).
4. **`analysis::energy_drift` with E₀ = 0** silently returns plain joules instead of a ratio, which is
   undocumented (HEAD `energy.hpp:65`).

## Comments and docs that are wrong
- `analytic/rocket.hpp:41-44` says the velocity formula returns −∞ at burnout of all mass. It returns
  **+∞**: ln(m₀/(m₀ − kt)) → +∞.
- `integrators/velocity_verlet.hpp:62-67` says the convergence test checks "a drag problem". It only
  uses a spring.
- `tests/test_numeric_vs_analytic.cpp:108-112` says the driven-oscillator test catches an RK4
  stage-time bug. Measured: a wrong-stage RK4 passes it (1.4e-4 error against a 5e-3 threshold).
- `analysis/oscillation.hpp:12-15` implies `angular_momentum` takes a `Trajectory`. It takes
  `(mass, State)`.
- `examples/01_free_fall.cpp:45-49` says printing `momentum(drop)` shows the mass cancelling. The drop
  starts at rest, so it prints (0, 0, 0).
- `examples/02_projectile_drag.cpp:16-17` says the range drops to "roughly half". The example's own
  output gives a ratio of 0.644.
- `forces/electromagnetic.hpp:34-35` explains when qE does work by talking about parallel vectors.
  What actually matters is whether qE·v ≠ 0.
- `notes/02-physics-to-code.md:276-279` says example 05 uses `analysis/energy.hpp`. It doesn't.
- `notes/02-physics-to-code.md:228-231` calls a parametric circle test "Kepler's second law".
- `notes/chapters/ch01-what-to-build.md:14, 112` and `docs/research/chapter-01-repo-audit.md:25` still
  say there is no coords module.
- `notes/learning-path.md:156, 278` cite lines 28–42 for `angular_momentum`. It's at 28–38.
- `docs/research/taylor-ch1-7-map.md:876` has the dissipation sign backwards. It should read
  T + U − ∫f_nc·v dt = const (line 516 has it right).
- `analysis/momentum.hpp:7-9` leaves out that internal forces must be **central** for a system's L to
  be conserved. A non-central equal-and-opposite pair keeps P but not L (measured: L changed 28×).
- `README.md` never mentions the coords module. The docs update in this build adds it.

## Tests that exist but are weak
- `test_analysis.cpp:125-135`: the "about another point" ℓ check compares a value that is exactly 0
  with a `<`, and never checks over time.
- `test_analysis.cpp:155-158`: says "empty / massless" input is tested, but only `{}` is checked.
- `test_integrators.cpp:179-184`: the energy bounds are 10³–10⁴× looser than measured.
- Untested functions:
  - `momentum()`
  - `tsiolkovsky_delta_v`
  - `linear_drag_velocity` / `linear_drag_position`
  - `DampedOscillator::v` at t > 0
  - `events::apex`
  - rocket altitude
  - `uniform_gravity_U`

## Smallest additions that would close the gaps (none implemented yet; your call)

| Topic | Addition |
|---|---|
| 2 | `coords/basis.hpp`: ρ̂, φ̂, r̂, θ̂; component projections; `polar_velocity` / `polar_acceleration`; an angle-unwrap helper |
| 3, 5, 8 | `analysis::total_mass`, `velocity_of_cm`, `momentum(mass, State)`; `systems/nbody.hpp` (pair forces with +f/−f, coupled Verlet/RK4 where every stage sees all bodies at the same instant); a two-body example and `tests/test_system.cpp` |
| 4 | The φ₁/φ₂ rewrite above; `analytic::linear_drag_range` (+ the small-drag expansion R ≈ R_vac(1 − 4v_y0/(3v_ter)), verified numerically); `quadratic_drag_horizontal_distance`; a drag-ratio helper |
| 6 | `analytic::rocket_altitude_gravity`, burnout state and max height; input validation |
| 7, 11 | `analysis::torque`, `total_angular_momentum(particles, about)`, `moment_of_inertia`, `angular_momentum_drift(traj, about)`, `swept_area`; an ℓ-conservation test on a simulated orbit; a Kepler example |
| 8 | `analytic/center_of_mass.hpp` with the continuous-body formulas |

## Where the evidence is
The full per-topic reports, verification programs and outputs are in this session's scratch folder
(`potentials-work/audit/NN-*.md` and `NN-*/`). That folder may be cleared on reboot. Everything above
was re-checked against the repo.

---

## Topics 12–16 (audited after the new modules were built)

These five audits doubled as the final review of the new `fields` / `potentials` / energy code. They found
**four real blockers**. I reproduced each one with my own program before it was fixed, and re-ran each
reproducer after the fix.

| # | Topic | Verdict after fixes | What exists | Remaining gaps |
|---|---|---|---|---|
| 12 | Conservative forces | Covered, with documented limits | `potentials::check_conservative`. Order: probe v → probe t → curl (confirmed at two step sizes and by a small Stokes loop) → closed-loop tiles → verdict. Verdicts: Conservative, TimeDependentPotential, VelocityDependent (+ NoWork/Dissipative/Feeding/Mixed/Unknown), NonzeroCurl, PathDependent, Inconclusive | Curl that exists only between time slices can be missed (raise `time_slices`). Lines that only clip a face or corner of the box can be missed (0% for lines through the inner 80% of the box) |
| 13 | Central forces | Partial | Spherically symmetric central forces classify Conservative even off-center; f(θ)r̂ forces classify NonzeroCurl; −∇U is radial; ℓ is conserved about the center | No `forces::Central{f}` helper, no radial-U helper, no spherical unit-vector API, no §4.8 example with ℓ |
| 14 | Potential energy functions | Partial | Ready-made: `uniform_gravity_U` (mgy), `spring_U` (½k\|r\|²). Any conservative force's U via `find_potential`. The U(∞) = 0 recipe for gravity-like forces | No ready-made Coulomb, Morse, LJ, Yukawa or anchored-spring U (they're one-line lambdas, or use find_potential). No §4.6 1D toolkit (equilibria, turning points, period by quadrature) |
| 15 | PE from F, and F from PE | Covered | `find_potential` / `potential_from_force` (U = U_ref − ∫F·dr, including time-dependent U(r,t)), and `force_from_potential` (F = −∇U, with `length_scale` and `center`). Round trips agree to about 1e-12 | No reference-at-infinity mode, and no analytic-gradient hook |
| 16 | Non-conservative forces, Δ(T + U) = W_nc | Covered | `work_done`, `energy_balance` (4th-order in the sample spacing, gauge-invariant scale), `power`. Checks: damped oscillator, quadratic fall, magnetic (W = 0), driven steady state W_drive = −W_drag = πbωA², and a forgotten force is detected. Force headers now say why drag, magnetic and drives are not conservative | No kinetic friction force (it needs stop-at-rest logic) |

**Blockers found and fixed during these audits** (each reproduced by me before and after):
1. **Missed holes.** The loop test (circles in three center planes) missed antiparallel wires, off-center
   wires and oblique lines: 12.7% of random lines. Now rectangular cells tile 3 planes per axis at 4 grid
   sizes. Re-check: antiparallel wires and an off-center wire are both correctly path-dependent.
2. **Silently wrong U.** `find_potential(force_from_potential(−6/r))` kept a singular reference at the
   origin and returned ΔU = 25,037 J where the truth is −3.27 J. The reference is now validated for path
   independence. Re-check: it moves the reference and returns −3.27224 J, exact.
3. **Off-origin steps.** Finite-difference steps were measured from the origin, so a potential centered
   1 km away had a 1.7e-4 force error. The fix adds a `center` to `force_from_potential`, and the
   classifier now measures its step from the box center. Re-check: 1.5e-10.
4. **A false alarm in the ledger.** The energy-ledger scale read 1.0 (total failure) for a perfect
   free-fall run starting at E = 0. The scale is now built from the energies involved. Re-check: 6e-14.

Smaller fixes from the same audits:
- Detours are clamped into the box.
- A caller's reference point outside the box is rejected.
- Time-dependent validation now includes t_min.
- A rounding-noise false "nonzero curl" is now rejected by the Stokes check.
- New velocity-power labels: Feeding for negative damping, Unknown for v/|v| at rest.
- ∂U/∂t wording.
- NaN-aware test helpers.
- Docs for the two ready-made U helpers now state their conventions (pending in the final docs pass).
