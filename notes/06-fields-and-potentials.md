# 06 — Fields and potentials: Taylor chapter 4 in code

This note explains `newton::fields` and `newton::potentials`, the modules that turn Taylor's chapter 4 into
functions, plus the energy-ledger functions added to `analysis/energy.hpp`. It is for someone who has just met
§4.1–4.5 and §4.8 and wants to see what the code does, why, and what its answers cannot tell you. Read it next to
`examples/07_potential_from_force.cpp` (force to potential energy) and `examples/08_force_from_potential.cpp`
(potential energy to motion). The API itself, with short snippets, is in `notes/05-library-reference.md`.

Every number below was measured on this repository's code with Apple clang 21 on arm64, and the snippets can be
re-run. Line references like `include/newton/fields/derivatives.hpp:176`
were checked on the day this note was written; if a header has moved on, search for the name instead.

## 1. Physics to code

| Taylor | The statement | Where it lives |
|---|---|---|
| §4.1 | Work W = ∫F·dr; power P = F·v; work–KE theorem ΔT = W | `fields::line_integral` (a path you choose), `analysis::power`, `analysis::work_done` (the simulated path) |
| §4.2 | Conservative: (i) F depends on position only, (ii) W is path-independent, i.e. ∮F·dr = 0; then U(r) = U_ref − ∫F·dr | `potentials::probe_dependence` (i), `fields::circulation` (ii), `potentials::potential_from_force` (U) |
| §4.3 | F = −∇U | `fields::gradient`, `potentials::force_from_potential` |
| §4.4 | ∇×F = 0 on a simply connected region is the local form of (ii); Stokes: ∮F·dr = ∫(∇×F)·n̂ dA | `fields::curl`, `fields::detail::relative_curl`, the loops in `check_conservative` |
| §4.5 | U(r, t): F = −∇U at each instant, but dE/dt = ∂U/∂t | `U(t, r)` accepted by `total_energy`, `energy_drift`, `energy_balance`, `force_from_potential`, `LineIntegralPotential` |
| §4.8 | Central force F = f(r) r̂ = −(dU/dr) r̂; but a force along r̂ whose size depends on direction, such as F = x r, has curl (0, −z, y) | `fields::gradient` of a function of \|r\| is purely radial; `coords::` projects onto r̂, θ̂, φ̂ |
| §4.1–4.2, §5.4 | Δ(T + U) = W_nc: drag and drives change E by exactly the work they do | `analysis::work_done`, `analysis::energy_balance` |

Two ideas run through the design. Taylor's two conditions become **tests you can run on any callable**, so no table of
"known conservative forces" is needed. And every test is an experiment that samples, so the answer is evidence.
Sections 4 and 5 are about what that evidence is worth.

## 2. Numerical derivatives

### 2.1 The 5-point stencil is Richardson extrapolation

The plain central difference D(h) = [f(x+h) − f(x−h)] / 2h has error (h²/6) f‴ + O(h⁴). The same formula with step 2h has
four times that h² error, so 4D(h) − D(2h) cancels it, and dividing by 3 leaves an error of order h⁴ (about h⁴f⁽⁵⁾/30).
Written out, that is `central_difference`
(`include/newton/fields/derivatives.hpp:176`):

```
f'(x) ≈ [ 8 (f(x+h) − f(x−h)) − (f(x+2h) − f(x−2h)) ] / 12h
```

It costs four evaluations per derivative, so `gradient` and `jacobian` cost 12. It is templated on the *value type*, so one
body serves a `double` field (a gradient component) and a `Vec3` field (a Jacobian column).

### 2.2 The step-size trade-off, measured

A smaller h cuts truncation error (like h⁴) but amplifies rounding error (like ε/h, with ε = 2.2×10⁻¹⁶, since you subtract
two nearly equal numbers and divide by a small one). The best h is where they meet. This table comes from the snippet below,
a Yukawa potential U = −e^(−r/a)/r with a = 0.5 m at |r| = 0.37 m. Errors are relative to the exact ∇U.

| h / \|r\| | 2-point | 5-point |
|---|---|---|
| 1e-2 | 6.1e-5 | 2.4e-8 |
| 1e-3 | 6.1e-7 | 2.3e-12 |
| **1e-4** | 6.1e-9 | **1.4e-12** |
| 1e-6 | 7.4e-11 | 9.1e-11 |
| 1e-8 | 9.9e-9 | 1.4e-8 |
| 1e-10 | 1.2e-6 | 1.7e-6 |

```cpp
auto U = [](const Vec3& r) { const double d = norm(r); return -std::exp(-d / 0.5) / d; };   // Yukawa, a = 0.5 m
const Vec3 r{0.3, 0.2, 0.1};
const double d = norm(r);
const Vec3 exact = std::exp(-d / 0.5) * (1.0 / (d * d) + 1.0 / (0.5 * d)) * (r / d);       // by hand
for (double c : {1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10}) {
    const double h = c * d;
    const Vec3 two_point{(U(r + Vec3{h, 0, 0}) - U(r - Vec3{h, 0, 0})) / (2.0 * h),
                         (U(r + Vec3{0, h, 0}) - U(r - Vec3{0, h, 0})) / (2.0 * h),
                         (U(r + Vec3{0, 0, h}) - U(r - Vec3{0, 0, h})) / (2.0 * h)};
    const Vec3 five_point = fields::gradient(U, r, h);
    std::printf("%.0e  %.1e  %.1e\n", c, norm(two_point - exact) / norm(exact),
                norm(five_point - exact) / norm(exact));
}
```

The 5-point error falls 10⁴ per decade (h⁴), bottoms out near 10⁻¹² around h = 10⁻³|r|…10⁻⁴|r|, then **rises** as rounding
takes over. The 2-point formula falls only 10² per decade and its best case (6e-11) needs a tiny, already noisy step. So
the 5-point stencil wins at any step, and its sweet spot is wide.

The default is h = 10⁻⁴·max(|r|, length_scale) (`default_step`, `include/newton/fields/derivatives.hpp:67`),
a decade below the optimum here, because real fields change faster than |r| (a Yukawa range shorter than |r|, gravity
at |r| < 1 m). In an earlier study of harder fields (not re-run for this note), h = 10⁻³ gave a worst-case Jacobian
error of 1e-6 and h = 10⁻⁴ gave 1e-10. Over the standard fields the worst case is about 3e-12, far below an
integrator's own error. The step grows with |r| so it also works at Earth-orbit scale (|r| ~ 10⁷ m): a fixed 10⁻⁴ m step
there would be 10⁻¹¹ of the length over which gravity changes, and rounding (about ε times that ratio) would spoil
2×10⁻⁵ of the answer.

### 2.3 Why `length_scale` exists: the atomic Morse potential

Near the origin |r| → 0 would give h → 0, so |r| is floored by `length_scale`, 1 m by default. For everyday mechanics
the floor is invisible. For a bond-sized potential it is a trap. Take U = D (1 − e^(−(r−r₀)/S))² with D = 10⁻¹⁹ J,
r₀ = 10⁻¹⁰ m, S = 3×10⁻¹¹ m, at r = 1.2×10⁻¹⁰ m:

- **Default step**, h = 10⁻⁴ m. At the four stencil points the exponent is about −3×10⁶, so the exponential underflows
  to exactly 0.0 and U = D at all four. The gradient is **exactly zero**: no warning, no NaN, a molecule that feels no
  force. In a simulation the energy drift is then exactly 0, which looks perfect, while the position is badly wrong.
- **`default_step(r, 1e-10)`**, h = 1.2×10⁻¹⁴ m. The gradient is 1.66547×10⁻⁹ N, equal to the exact value in every printed digit.

A power law behaves differently from an exponential: with the default step, a Lennard-Jones force at atomic scale comes
out as a wrong, tiny number (3e-54 N where 1.3e-12 N is right) rather than exactly 0, which is harder to notice.

`force_from_potential(U, length_scale, center)` takes the same knob
(`include/newton/potentials/force_from_potential.hpp:72`), and `check_conservative` derives it from
the box: L = half the largest side (`box_length_scale`, `include/newton/potentials/conservative.hpp`),
so a molecule-sized box gets molecule-sized steps automatically.

**The `center` argument.** The step is 10⁻⁴·max(|r − center|, length_scale), measured from `center`, which defaults to the
origin. For a potential centred at the origin that is exactly right. For a well centred 1 km away, evaluated 1 m from the
well, the step is 0.1 m, far too coarse for a structure 1 m across: measured force error 4e-4 without `center` and
2.5e-10 with `center = {1000, 0, 0}`. Rule: if your potential's structure (the nucleus, the well, the singularity) is
not near the origin, pass its position as `center`, or write U about the origin. The classifier does the same for its
box by measuring distances from the box centre.

### 2.4 Where finite differences are wrong

They cannot see structure smaller than about 2h.

- **A kink.** For U = |x| the exact slope is ±1. The stencil returns, at x/h = 0, 0.5, 1, 1.5, 2: 0, 0.5833,
  **1.1667**, 1.0833, 1.0000: a peak 7/6 of the truth (17% high), gone by 2h.
- **A singularity.** For −GMm/r the stencil straddles r = 0 within 2h and returns garbage.
- **A multi-valued U.** U = −atan2(y, x) jumps by 2π across the negative x axis. Away from the cut its gradient is the
  hole field (−y, x, 0)/(x² + y²): curl-free, but not the gradient of any single-valued U. A stencil straddling the cut
  returns a spike of order 2π/h.

## 3. Line integrals

### 3.1 Gauss–Kronrod in plain words

A Gauss rule with n points integrates every polynomial up to degree 2n − 1 exactly, if the points sit at special spots.
A Kronrod extension adds n + 1 points, keeping the old ones, so all 2n + 1 integrate degree 3n + 1 exactly. With n = 7,
15 evaluations (`kGkNodes`, `include/newton/fields/line_integral.hpp:75`)
give two answers: K15, very accurate, and G7, less so, built from 7 of the same points. Their difference |K15 − G7| is a
free error estimate, deliberately pessimistic (10²–10⁴ too big for a smooth integrand), which is the right bias for "do I
trust this?". A smooth integrand needs one piece (15 evaluations). Measured: a spring along a segment that misses its
anchor needs 15, and PointGravity from r = 1 m to r = 2 m needs 45.

### 3.2 Why the answer is a struct

`LineIntegral` (`include/newton/fields/line_integral.hpp:48`) holds `value`, `error`, `scale`,
`evaluations`, `converged`. A line integral can fail in ways a bare `double` cannot say: the path hits a singularity, the
budget runs out, `F` returns NaN. Throwing would turn a singularity on one of the 351 loops of a time slice into a crash of
the whole test, when it is really *information*; a quiet wrong number would be worse. So the detailed functions return
everything, and the plain-`double` overloads return NaN on failure, so a failure cannot pass for a result. `scale` is
S = ∫|F||dr|, "the work if F always pointed along the path", the yardstick for "is this tiny?": W = 10⁻¹² J is zero if
S = 1 J and a real effect if S = 10⁻¹⁰ J.

### 3.3 Global versus recursive adaptivity

The textbook adaptive scheme is recursive: if a piece's error is too big, split it and integrate each half to *half the
tolerance*. Near a singularity this explodes: once the per-piece tolerance falls below roundoff, every piece fails and each
failure spawns two more. The research for this library measured **24 to 44 million evaluations** for a path ending near a
point mass.

`integrate_adaptive` (`include/newton/fields/line_integral.hpp:197`) is *global*: it keeps a
list of pieces, sums their errors, and if the total is too big splits only the **piece with the largest error**. Work goes
where the error is, and `kMaxPieces` = 200 (`include/newton/fields/line_integral.hpp:93`)
caps the cost at 15 + 30 × 199 = 5985 evaluations. The stopping rule is

```
error ≤ max( rel_tol · |W|,  scale_tol · S )
```

The first term is the usual relative goal. The second is a floor: when the true answer is 0 (a conservative loop), no
routine reaches a relative goal, because positive and negative parts cancel only down to roundoff times S. A piece
narrower than 10⁻¹² of the parameter range "is a point": splitting further only chases the singularity, so the search
stops with `converged = false`.

### 3.4 The symmetric-singularity story

Integrate PointGravity along the x axis from (1, 0, 0) to (−1, 0, 0), straight through the centre. The true work is
infinite on each half but cancels by symmetry. The integrand F·dr is **odd** about the centre, and the Gauss–Kronrod nodes
are symmetric about the segment's midpoint, which is the centre. K15 and G7 then cancel to the same number, |K15 − G7| = 0,
and the integral "converges" to the principal value 0: a number from an integral that does not exist.

The fix is a second criterion (`kScaleErrorTol`, `include/newton/fields/line_integral.hpp:111`):
the integral of |F||dr| must *also* be resolved. That integrand is **even**, never negative, so it cannot cancel and
exposes the divergence by refusing to settle. The pieces keep splitting until the minimum width stops the search.
Measured: both paths through the centre, to (−1, 0, 0) and to (−2, 0, 0), return `converged = false` after 1215 evaluations.

The price: |F||dr| can have a kink on a perfectly smooth path. A spring segment through its anchor has |F| = k|x|, a V,
and needs extra splits: 45 evaluations instead of 15 (105 with a natural length of 0.5 m). The answer stays accurate; only
the cost rises. The library accepts that, because it is what makes a symmetric path through a singularity report failure
instead of a plausible number. `potentials` then recovers with detours (§5).

### 3.5 Circles

`circulation` integrates r(θ) = c + R(cos θ û + sin θ ŵ) over θ ∈ [0, 2π] with the same routine. The frame (û, ŵ, n̂) is
right-handed (`circle_frame`, `include/newton/fields/line_integral.hpp:256`), so the loop runs
counterclockwise seen from the tip of `normal`, and flipping the normal flips the sign. A zero normal or a radius ≤ 0 returns
`converged = false` without evaluating anything. (The classifier uses circles for its small Stokes check, §4.2; its big loops are
rectangles built from four segments, §4.3.)

## 4. How `check_conservative` decides

The function (`include/newton/potentials/conservative.hpp`) applies Taylor's
conditions in his order and stops at the first one that fails:

0. **A degenerate box or invalid evidence configuration** gives `Inconclusive`, with the reason in `report.note`.
   The curl, loop and time-slice counts must be positive; a callable capable of time or velocity dependence also
   needs a positive probe count, and its relevant range and tolerances must be finite and positive.
1. **Probe v, then t** (`run_probe`, `include/newton/potentials/conservative.hpp`). By default it uses 64 samples, each a random r
   in the box, two random times and two random velocities in [−v_max, v_max]³. Change only v (or only t) and compare the
   exact `Vec3` values: any difference between two finite results *proves* dependence. Dependence on v gives
   `VelocityDependent` (and `velocity_power`, below). Stop.
2. **Curl at 256 random points** per time slice. `jacobian_with_error` (13 evaluations) gives the Jacobian plus an error
   estimate, scored by `relative_curl`. A score above 10⁻⁶ is a *candidate*; it becomes `NonzeroCurl` only after two
   confirmations (below). Stop at the first confirmed one.
3. **Closed loops**: 351 rectangles per slice. A converged loop with |∮F·dr| / ∮|F||dr| above 10⁻⁸ gives `PathDependent`.
   All loops run, so `max_loop_ratio` is the true worst, and `worst_loop_center` says where to look.
4. **Otherwise** `TimeDependentPotential` if F depends on t, else `Conservative`.

### 4.1 The curl criterion and its error term

|curl F| has units of N/m, so "is it small?" needs "compared with what?". The score
(`fields::detail::relative_curl`, `include/newton/fields/derivatives.hpp:246`) is

```
relative_curl = max(0, |curl J| − error) / ( ‖J‖_F + F_max / ℓ )
```

The denominator is the field's own scale: its derivatives, plus |F|/ℓ, the variation expected from a field of size |F|
changing over ℓ = max(|r − c|, L), where c is the **box centre** and L half its largest side. (Measuring from the box centre,
not the origin, keeps the step sensible for a box that sits far from the origin: the same lesson as `center` in §2.3.)
Every term is linear in F, so the score ignores units and overall size. Measured: a spring 7e-14, uniform gravity exactly 0,
and a spring plus a swirl of relative strength δ scores about δ (δ = 10⁻⁷ gave 9.8e-8). The threshold is 10⁻⁶. The header
records, on the test battery, a worst conservative score of 2.7e-13 and a weakest non-conservative field it must catch at
4.9e-5, with the threshold between the two.

The **error term** is the 4th difference F(x+2h) − 4F(x+h) + 6F(x) − 4F(x−h) + F(x−2h)
(`jacobian_with_error`, `include/newton/fields/derivatives.hpp:212`). It is zero for any
cubic polynomial, so it measures only what is *not* the smooth field: rounding noise, or a singularity the step cannot
resolve. A curl smaller than that could have been manufactured by the method, so it is subtracted.

The order of operations matters. The stencil computes `8 * (p1 - m1) - (p2 - m2)`, differences first, so a constant field
(p1 = m1 = p2 = m2) gives exactly 0. Summing `c − 8c + 8c − c` instead leaves roundoff that looks like a curl; the
research found uniform gravity misclassified as `NonzeroCurl` before this was fixed.

### 4.2 Two confirmations: h/2 and a small Stokes loop

A candidate must pass both before it counts as `NonzeroCurl`.

(a) **Recompute with the step halved**
(`include/newton/potentials/conservative.hpp`). A **true** curl does not depend on h, so
the two agree to 10%. Truncation error shrinks 16-fold when h is halved and rounding noise changes completely, so most
artifacts disagree.

(b) **A small Stokes loop** (`include/newton/potentials/conservative.hpp`). Integrate F around a circle of radius 0.05 L centred on the
sample, with its normal along the measured curl. Stokes says the circulation must equal |curl|·πR², to 25%. Coherent
rounding noise can fool (a), because it can look the same at both step sizes, but it does not integrate to a circulation.

Candidates that fail either are counted in `unconfirmed_curl_candidates` and ignored.

### 4.3 Loops, Stokes, and holes

Why loops, when curl is local and cheap? Taylor's §4.4 says curl F = 0 implies path independence only on a **simply
connected** region, where every loop can shrink to a point. Remove a point from 3D space and it stays simply connected
(PointGravity is fine). Remove a *line* and it does not. The hole field (−y, x, 0)/(x² + y²) has curl 0 everywhere off
the z axis, yet ∮F·dr = 2π around any loop that encircles it. Measured in the default box [−1, 1]³: worst curl score 0,
worst loop ratio **0.89**, verdict `PathDependent`. Stokes' theorem ∮F·dr = ∫(∇×F)·n̂ dA is the bridge: a loop sees all the
curl inside it, including a delta function on the axis that no sampled point lands on.

The layout is deterministic. In each of the 3 coordinate directions there are 3 parallel planes (the box centre, and a
quarter of the side to either side). In each plane, a g×g grid of rectangular cells tiles the box's cross-section, for
g = 1, 2, 3 and 5. That is 3 × 3 × (1 + 4 + 9 + 25) = **351 loops per slice**, each four straight segments integrated with
`line_integral_detailed` (`detail::cell_circulation`, `include/newton/potentials/conservative.hpp`).
The grid sizes are coprime, so a singular line cannot sit on the cell edges of all of them at once.

Why rectangles in many planes, and not a few circles? A hole is seen by a loop only if the line pierces the loop's disc,
and a straight line crosses the box's cross-section somewhere, not necessarily near the centre. The earlier design, a
fixed set of circles in the three central planes, missed 12.7% of random oblique lines, and missed antiparallel wire
pairs whose circulations cancel in one big loop (the header's account). Cell loops cancel a wire pair only if both wires are
in the same cell, which the finer grids and shifted planes avoid. Measured with the new layout: a pair of antiparallel
wires at x = ±0.5 gives `PathDependent` (worst ratio 0.88, at a loop centred on one wire), and a hole displaced to
(0.37, 0.21) gives `PathDependent` (0.86).

Each loop is integrated with `scale_tol = 10⁻²·loop_tol`, so a converged loop has error below 10⁻²·loop_tol·S and a ratio
above `loop_tol` can never be quadrature noise. A loop with a segment that touches a singularity does not converge; it is
skipped and counted, never taken as evidence (the hole field in the default box skips 38 of 351; PointGravity skips 12).

### 4.4 Time slices, and why not t = 0

The induced swirl F = sin(t)(−y, x, 0) has curl 2 sin(t) ẑ, exactly zero at t = 0, π, 2π, … A test at t = 0, or at nice
round instants, calls it conservative. Measured: with the time window collapsed to [0, 0] the verdict is `Conservative`;
with the default window it is `NonzeroCurl` at t = 1.2361 s. `slice_times`
(`include/newton/potentials/conservative.hpp`) puts one golden-ratio offset in each of K = 5
equal sub-windows, t_k = t_min + (k + 0.618…)/K · (t_max − t_min): irrational, so never on a round instant, and
deterministic. A force that does not depend on t gets one slice, at t_min. Only those K instants are tested: a curl that exists
only between them is missed (raise `time_slices`). The random points are seeded
(`std::mt19937_64` plus a hand-written conversion to [0, 1), `include/newton/potentials/conservative.hpp`) because
`std::uniform_real_distribution` may differ between standard libraries.

### 4.5 Inconclusive, and velocity_power

`Inconclusive` means no verdict can be given, and `report.note` says why. It happens for a degenerate box, or if no positive
evidence was found *and* the probe skipped every sample, or more than half the curl samples of a slice were skipped, or every
loop of a slice was skipped (non-finite force values, or integrals that did not converge). Positive evidence wins: a confirmed
curl or a loop with work proves non-conservativeness even if parts of the box could not be evaluated. A force that returns NaN
everywhere gives `Inconclusive` (measured: 256 of 256 curl samples skipped), and so does a flat box.

`velocity_power` uses the same probe samples, and it labels the velocity-dependent **part** F_v = F(t, r, v) − F(t, r, 0), not the
whole force. It takes the cosine between v and F_v (`include/newton/potentials/conservative.hpp`):

| Label | Cosine at every usable sample | Example |
|---|---|---|
| `NoWork` | \|cos\| ≤ 10⁻⁹ | magnetic force q v×B, P = F·v = 0 |
| `Dissipative` | cos < −10⁻⁹ | linear or quadratic drag, kinetic friction |
| `Feeding` | cos > +10⁻⁹ | negative damping, F = +b v |
| `Mixed` | both signs genuinely occur | the made-up force (v_y, 0, 0) |
| `Unknown` | no usable sample | friction written −μ v/\|v\| is NaN at rest; or F_v too small to measure |

A sample is "usable" only if F_v is bigger than the rounding noise of the subtraction F(v) − F(0) (about 4e-7 of |F|), so a tiny
magnetic part on a huge electric field gives `Unknown` rather than a wrong label. Measured: `MagneticForce` and a `LorentzForce`
with B ≠ 0 give NoWork; `LinearDrag` and `QuadraticDrag` give Dissipative; +0.5 v gives Feeding; friction gives Unknown.
`LorentzForce` with B = 0 computes q(E + v × 0), bitwise independent of v (−0.0 equals +0.0), so it is reported as
velocity-independent: the probe judges what the force *does*, not what its signature allows. The label does not say what the whole
force does: for drag in a wind the label is `Dissipative` though the true drag power can be positive (§7).

## 5. The honest limits

A verdict is **evidence**. `NonzeroCurl`, `PathDependent` and `VelocityDependent` are strong evidence that a force is not
conservative (a curl confirmed twice, a loop that does work, a measured dependence on v). `Conservative` means only "no evidence
found".

- **Probes prove dependence, never independence.** A force that switches on at t = 100 s looks static on the default window
  [0, 10] s: measured, the verdict is `Conservative` and `depends_on_time` is false. A speed threshold is found reliably up to
  about 1.4·v_max (the corners of the v_max cube), and not beyond. Set `t_max` and `v_max` to cover the motion you care about.
- **The box matters.** "Conservative" is a statement about a region. The hole field is `PathDependent` in [−1, 1]³ and
  `Conservative` in [0.5, 3] × [−1, 1]², which avoids the axis (loop ratio 6e-17).
- **Lines that clip the box.** A singular line that only crosses a face or a corner of the box can be missed: the header's
  measurement is 4.8% of random lines over the whole box, and 0% for lines through the inner 80%. Make the box
  generous enough to contain the motion and the sources.
- **Dilution.** The loop ratio divides by ∮|F||dr|, so a huge uniform field added to a circulating one hides it. Measured with
  the hole field plus a constant push: 10⁶ N gives ratio 3.9e-6, 10⁸ N gives 3.9e-8 (still above the 10⁻⁸ tolerance), and 10⁹ N
  gives 3.9e-9: `Conservative`, though the force is not.
- **Noisy evaluators.** The Stokes confirmation rejects a false *local*
  `NonzeroCurl` caused by coherent rounding noise, but it cannot make a
  discontinuous numerical evaluator path-independent: its rectangular loops
  can still report `PathDependent`. Treat a force whose own values are noisy
  near the classifier tolerances as `Inconclusive` in your own analysis, or
  improve the evaluator.
- **Weak non-conservative parts.** A spring plus a swirl of relative size 10⁻⁷ has curl score 9.8e-8, below the curl threshold,
  and is caught by the loops (`PathDependent`). At 10⁻⁸ the loop ratio is 8.7e-9, just under its tolerance, and the verdict is
  `Conservative`. The sensitivity floor is therefore between 10⁻⁸ and 10⁻⁷ relative. A force whose own evaluation is noisy at
  the 10⁻⁹ relative level may be misjudged, because noise can survive the curl confirmation (usually the Stokes loop rejects it,
  not always) and noisy loops do not converge.
- **The reference point.** U(r) = U_ref − ∫F·dr from r_ref, and only differences of U are physical, so *where U is zero is
  your choice*. `find_potential` defaults to the origin when it is in the box (so mgy, ½k|r|² and −F·r come out in textbook
  form), otherwise the box centre, and **validates** it by checking *path independence*: from the reference to each of the 8
  box corners (at t_min and every time slice), the work along the straight path and along an independent off-axis detour must
  both converge and agree to 10⁻⁹. Convergence alone proves too little. A finite F at the reference proves nothing
  (PointGravity returns 0 exactly at its centre), and a numerical gradient (`force_from_potential`) of a singular U is *finite*
  near the singularity but is wrong, in a direction-dependent way, inside a core of a few steps. Two different paths to the
  same corner expose both. PointGravity at the origin fails, so the reference is moved, and `describe()` says so (measured:
  to (−0.5, −0.5, −0.5) m; for a numerical −1/r potential, to a point off the axes). A reference *you* choose is never moved: if
  it is outside the box or fails, `reference_ok` is false and `note` says why.
- **U(∞) = 0.** Gravity and Coulomb are conventionally zero at infinity, which no finite integral can start from. Pass any
  regular reference with its textbook value, `potential_from_force(gravity, r_ref, -GM*m / norm(r_ref - center))`, and U matches
  −GMm/r *inside the box you checked*, to about 10⁻¹². For softened gravity the textbook value is
  `-GM*m / sqrt(norm2(r_ref - center) + eps*eps)`.
- **U is certified only inside the tested box.** `work_between` (`include/newton/potentials/potential_from_force.hpp:77`)
  retries with four sideways detours when a straight path fails or is ill-conditioned (S more than 10³ times |W|), keeps
  the converged candidate with the smallest S, and returns NaN if none converge, never a silently wrong number. A detour is
  legitimate only where F is path-independent: around the hole field, detours on opposite sides of the axis differ by exactly
  2π. So `find_potential` stores the box in the potential (`LineIntegralPotential::region`) and clamps every detour into it.

## 6. Time-dependent U and gauge freedom

If F = −∇U(r, t), then along the motion dE/dt = d(T + U)/dt = F·v + ∇U·v + ∂U/∂t = ∂U/∂t, since F = −∇U cancels the
middle terms. Energy is not conserved though nothing dissipates: whoever moves the potential does work. Note that
`force_from_potential` takes only the position gradient; ∂U/∂t is not part of the force.

U is defined only up to U → U + c(t) for **any function** c(t), because −∇U does not change: a whole function's worth of
freedom, not one additive constant. `LineIntegralPotential`
(`include/newton/potentials/potential_from_force.hpp:141`) freezes t during the integral and so fixes
U(r_ref, t) = U_ref at every t, one particular c(t). E(t) = T + U then depends on the gauge, but

```
E(t) − E(0) = ∫₀ᵗ (∂U/∂t) dt′
```

holds in **every** gauge, because both sides shift by c(t) − c(0). Compare differences within one gauge, or check the
identity, which is what `energy_balance` does.

Measured on a moving spring support, F = −k(r − A cos(ωt) x̂) with k = 4, A = 0.5, ω = 3, from rest at the origin: the
line-integral gauge U_LI = 2|r|² − 2 cos(3t) x gives E in [0, 2.88] J; the natural gauge ½k|r − s(t)|² gives E in
[0.1495, 3.38] J; they differ by ½kA² cos²(ωt) to 1.4×10⁻¹⁵; and `energy_balance` is about 3×10⁻¹¹ in both. Same physics,
different E(t).

## 7. Energy of non-conservative forces

Split the net force as F = −∇U + F_nc. The work–KE theorem ΔT = W_net and W_cons = −ΔU give

```
Δ(T + U) = W_nc,      i.e.      dE/dt = F_nc · v       (plus ∂U/∂t if U depends on t)
```

so drag and drives change E by exactly the work they do (Taylor §4.1–4.2, §5.4). `power(F, state)` is F·v at one state.
`work_done(traj, F)` follows the **recorded motion**: W[i] = ∫₀^{tᵢ} F(t′, r, v)·v dt′. That is the difference from
`line_integral`, which follows a path you pick and cannot handle a v-dependent force (the work of drag depends on how fast
you traverse the path). `energy_balance(traj, U, F_nc)` returns the residual of Δ(T + U) = W_nc (+ ∫∂U/∂t) at every
sample. It is near 0 whenever −∇U + F_nc is the force you actually simulated.

**Sign convention and scale.** W_nc < 0 means F_nc removes energy (drag: −b v²). The residual is (energy change) − (booked
work), so a **negative** residual means energy was lost that you did not book (a forgotten drag), and a positive one that energy
appeared from nowhere (a forgotten drive). The scale is the largest energy that actually crosses the ledger: the maximum over
samples of (T_i, |U_i − U_0|, |W_nc,i|, |∫∂U/∂t|_i), or 1 J if all are zero. It is built only from kinetic energy and
*differences* of U, so it is gauge-invariant. An earlier scale included |E₀|, which is 0 for a body falling from rest at
U = 0: a perfect free-fall ledger then read as total failure (1.0) because the scale vanished. Measured now:
a free fall from rest gives residuals around 8e-14.

For a time-dependent U, the 5-point ∂U/∂t stencil starts at 1% of the local
recorded time gap. It is rounded upward only enough that its shifted
timestamps are representable, so moving the same samples from t = 0 to a
large absolute clock value does not turn a millisecond-scale derivative into
a seconds-scale one.

**Accuracy is O((record_every·dt)⁴).** `cumulative_integral`
(`include/newton/analysis/energy.hpp`) is composite Simpson on pairs of intervals, written for **uneven**
spacing, so `simulate_until`'s short last interval or a `record_every` that does not divide the step count are fine. At odd
samples it adds one last interval from the parabola through the last three. The power is evaluated *exactly* at each
sample, so the only error is quadrature plus the trajectory's own. Measured on a damped oscillator (m = 1, k = 4, b = 0.4,
x₀ = 1; exact W_drag(10 s) = −1.9596345276 J), RK4, dt = 0.01, largest ledger residual as `record_every` grows:

| record_every | sample spacing | ledger residual |
|---|---|---|
| 1 | 0.01 s | 1.1e-8 |
| 2 | 0.02 s | 1.7e-7 |
| 4 | 0.04 s | 2.6e-6 |
| 8 | 0.08 s | 4.2e-5 |

Each doubling multiplies the error by about 16 = 2⁴. (The period is 3.16 s, so even 0.08 s is 40 samples per period.) The
trapezoid rule at dt = 0.01 gave 3.6e-7 in W_drag(10), against 1.2e-8 for this Simpson. For bookkeeping, record every step
and resolve at least about 20 samples per period.

**The ledger finds mistakes.** Simulate the damped oscillator but give `energy_balance` the spring's U and a *zero* F_nc
(forgetting the drag): the residual grows to −0.98, negative because the drag removed energy that the ledger did not book. A
forgotten or double-counted force shows up as a growing residual.

**Splitting at rest.** `velocity_independent_part(F)` = F(t, r, 0) and `velocity_dependent_part(F)` = F(t, r, v) − F(t, r, 0)
split any force exactly. If the rest part is curl-free it has a potential U₀, and d(T + U₀)/dt = F_v·v + ∂U₀/∂t. That is exact
**bookkeeping**, not Taylor's physical classification, which sorts forces by what they are. It is clean for a magnetic force
(F_v·v = 0, so T + U₀ is conserved) and for drag in still air. It misleads for drag in a wind, F = −b(v − w): the split books
the constant push b w as a "potential" −b w·r and −b v as "dissipation". Measured with b = 0.5, w = (3, 0, 0), v = (1, 0, 0):
the true drag power is **+1.0 W** (the wind pushes), the split's "dissipation" is −0.5 W, and the heat produced, −b|v − w|²,
is −2.0 W. The ledger balances, but U₀ is not stored energy. For a wind, pass the whole drag as `F_nc` to `energy_balance`.
A drive is a similar choice: put it in U(t, r), or treat it as `F_nc` (Taylor ch. 5's natural bookkeeping). Both balance.

## 8. The C++ this module leans on

**`if constexpr` plus `std::is_invocable_r_v`: shape detection.** The classifier accepts three signatures and works out
which one it got, at compile time (`shape_v`, `include/newton/potentials/conservative.hpp`):

```cpp
template <class F> const char* shape_of(const F&) {
    if constexpr (std::is_invocable_r_v<Vec3, const F&, double, const Vec3&, const Vec3&>) return "force F(t, r, v)";
    else if constexpr (std::is_invocable_r_v<Vec3, const F&, double, const Vec3&>) return "time field F(t, r)";
    else if constexpr (std::is_invocable_r_v<Vec3, const F&, const Vec3&>) return "static field F(r)";
    else static_assert(sizeof(F) == 0, "needs F(t, r, v), F(t, r) or F(r) returning Vec3");
}

int main() {
    std::puts(shape_of([](double, const Vec3& r, const Vec3& v) { return -r - v; }));   // force F(t, r, v)
    std::puts(shape_of([](const Vec3& r) { return -r; }));                               // static field F(r)
    // shape_of([](int) { return 0; });  // compile error: needs F(t, r, v), F(t, r) or F(r) returning Vec3
}
```

`is_invocable_r_v<R, F, Args...>` asks whether `F` can be called with `Args...` and the result converts to `R`. In a false
`if constexpr` branch the code is *discarded*, not compiled, so a one-argument lambda never meets a three-argument call.
The library wraps whichever shape it found into the 3-argument Force concept (`AsForce`,
`include/newton/potentials/conservative.hpp`), and `energy.hpp` uses the same trick for `U(r)` versus `U(t, r)`
(`include/newton/analysis/energy.hpp`), so every existing call compiles unchanged. The test is on
`const F&`, so a `mutable` lambda is rejected on purpose: a force that changes between calls would break "same input, same output".

**`static_assert` with a message.** A wrong callable gets one sentence, not a page of template errors
(`include/newton/fields/derivatives.hpp:299`, and `kShapeMessage` in `conservative.hpp`).
The `jacobian` message shows how to freeze a Force with a lambda.

**`std::optional`.** `find_potential` takes `std::optional<Vec3> r_ref = std::nullopt`
(`include/newton/potentials/potential_from_force.hpp:345`) to tell "you did not choose a reference" from "you
chose the origin". A sentinel like "{0, 0, 0} means unset" would make the origin impossible to request.

**Factories that return function objects.** `velocity_independent_part(F)`
(`include/newton/potentials/conservative.hpp`) returns a lambda that captured the adapted force; its return type
is `auto` because a lambda's type has no name. `force_from_potential(U)`
(`include/newton/potentials/force_from_potential.hpp:114`) is a factory because C++17 cannot deduce the
template argument of an aggregate (`GradientForce{lambda, 1.0}` does not compile; C++20 added that), while a function
template can.

**`enum class`.** `Verdict` and `VelocityPower` (`include/newton/potentials/conservative.hpp`) are scoped:
`Verdict::Conservative` cannot clash with another enum's `Conservative` and never converts to an int by accident.

**`inline constexpr` tables.** The Gauss–Kronrod nodes and weights are namespace-scope `inline constexpr double` arrays:
`constexpr` makes them compile-time constants and `inline` gives one definition across every translation unit (see
`notes/01-cpp-for-this-library.md`). `0x1.0p-53` is a C++17 hexadecimal floating-point literal for exactly 2⁻⁵³.

**Templates, not `std::function`.** Every public function is a template over the callable's type, so a lambda is called
directly and can be inlined, with no heap allocation or indirect call. That matters when `gradient` calls the field 12
times, a loop test thousands, and `energy_drift` with a line-integral potential hundreds per sample. The cost is real: each
lambda instantiates new code, definitions must live in headers, and template errors are hard to read (hence the
`static_assert`s). Forwarding references (`ScalarField&&`) accept lambdas, function objects and function pointers without
copying.

## 9. Where to go next

Run both examples and compare `energy_drift` (large: it contains real physics) with `energy_balance` (tiny: it accounts for
it). Then break something on purpose: add `1e-8` of a swirl to a spring, collapse the time window to [0, 0] for the induced
swirl, or give `gradient` the wrong `length_scale`. What is still unbuilt (analytic gradients supplied by a potential, a
reference-at-infinity mode) is in `notes/04-extending.md`.
