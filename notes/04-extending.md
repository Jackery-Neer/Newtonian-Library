# Extending the Library: A Semester Roadmap

This library was built to grow alongside Mechanics 1, and it's built to
keep growing after. This note has three parts: a worked example of adding
a new force, the general recipe for adding a new analytic reference and
the test that validates it, and a paragraph on each v0.2+ item from
DESIGN.md's Roadmap section — what it is, and roughly what it would take.

## Adding a new force: worked example (a Yukawa potential)

Every force in this library follows the same recipe:
`forces::Spring`, `forces::PointGravity`, and every other functor in
`include/newton/forces/` all satisfy one concept —

```cpp
Vec3 operator()(double t, const Vec3& r, const Vec3& v) const;
```

— and are plain aggregates that compute a physical formula from `t`, `r`,
`v`. Adding a new one means finding the formula, writing the struct, and
hooking it into the test/analytic layers that check it.

Say you want a **Yukawa force** — a screened attraction that falls off
faster than gravity, used as a toy model for the strong nuclear force or
for a screened Coulomb interaction in a plasma. Taylor's chapter-4
potentials table (see `docs/research/taylor-ch1-7-map.md` §4.2) lists it
alongside gravity and Coulomb:

```
U(r) = -(k/r) e^(-r/a)
```

where `k` sets the interaction strength and `a` is the screening length —
at `r ≫ a` the potential is exponentially suppressed rather than falling
off as a bare `1/r`.

**Step 1 — derive the force.** `U` depends only on the radius, so this is
a central force (Taylor §4.8): `F = -dU/dr · r̂`. Differentiate:

```
dU/dr = k e^(-r/a) [1/r^2 + 1/(a r)]
F_r   = -dU/dr = -k e^(-r/a) [1/r^2 + 1/(a r)]
```

For `k > 0`, `F_r < 0` — the force points back toward the center, i.e. it's
attractive, matching the physical Yukawa interaction.

**Step 2 — write it in Cartesian form, matching the Force concept.** With
`d = r - center` and `dist = |d|`, `r̂ = d/dist`:

```cpp
struct Yukawa {
    double k{1.0};       // interaction strength [N*m], the Coulomb-like prefactor
    double a{1.0};        // screening length [m]
    Vec3   center{};      // source location [m]

    [[nodiscard]] Vec3 operator()(double /*t*/, const Vec3& r, const Vec3& /*v*/) const {
        const Vec3 d = r - center;
        const double dist = norm(d);
        if (dist <= 0.0) return {};  // zero-guard: particle exactly at the source
        const double decay = std::exp(-dist / a);
        return (-k * decay * (1.0 / (dist * dist * dist) + 1.0 / (a * dist * dist))) * d;
    }
};
```

Note the pattern this borrows from the existing forces: the zero-guard on
`dist <= 0.0` mirrors `PointGravity`'s guard on its softened denominator
(`include/newton/forces/gravity.hpp`) and `Spring`'s zero-guarded
`normalized()` call (`include/newton/forces/spring.hpp`) — a particle
sitting exactly on the source is a genuine singularity of this force, and
returning zero rather than `NaN` or an `inf`-propagated vector is the
library-wide convention (see `core/vec3.hpp`'s `normalized()` doc comment
for why this pattern exists at all). `t` and `v` are unused and commented
out, matching every other position-only force.

**Step 3 — doc-comment it** the way every force in this library is
commented: a `Physics:` paragraph citing the Taylor section and explaining
what the force represents physically, and (if there's a non-obvious C++
decision, like the zero-guard above) a `C++ note:` paragraph explaining it.
Every parameter gets its unit in square brackets, per DESIGN.md's hard
rule 5.

**Step 4 — place it and register it.** A new physical force belongs in
`include/newton/forces/`, as its own file (`yukawa.hpp`) if it's a whole
new physical idea, or added to an existing file if it's a close cousin of
something already there. Add its include line to
`include/newton/forces/forces.hpp`, the module umbrella, so
`#include <newton/forces/forces.hpp>` (or the top-level
`newton/newton.hpp`) picks it up automatically.

**Step 5 — test it.** Two independent checks, matching how the existing
forces are tested:

- **Hand-computed value** in `tests/test_forces.cpp`: pick a specific
  `(k, a, center)` and a specific `r`, compute `F` by hand (or with a
  calculator), and assert the functor returns it within a tight tolerance.
  Also test the `dist == 0` edge case returns exactly the zero vector, not
  `NaN` — the same pattern `test_forces.cpp` already uses for
  `QuadraticDrag` at `v = 0`.
- **Energy conservation**, since Yukawa is conservative and its potential
  is known in closed form. Run `simulate()` with `Yukawa{...}` as the
  force, and check `analysis::energy_drift` (`include/newton/analysis/energy.hpp`)
  stays bounded, passing a hand-written potential lambda
  `[k, a, center](const Vec3& r) { return -k / norm(r - center) *
  std::exp(-norm(r - center) / a); }` — precisely the pattern
  `uniform_gravity_U` and `spring_U` already establish for `UniformGravity`
  and `Spring`.

That five-step recipe — derive, write the struct against the Force
concept, doc-comment, register in the umbrella header, test both the raw
value and (when a potential exists) energy conservation — is exactly how
every force currently in `include/newton/forces/` was built, and it's the
recipe for anything from Taylor's problem sets you want to simulate that
isn't already here (a Coulomb force, a buoyancy restoring force, a
velocity-cubed drag law — anything expressible as `F(t, r, v)`).

## Adding a new analytic reference + test

The other half of extending this library is adding a new closed-form
oracle to `newton::analytic` — the layer that lets you trust the numbers
`simulate()` produces. The general recipe:

1. **Find the closed form.** `docs/research/taylor-ch1-7-map.md` §8 is a
   consolidated table of every closed-form solution the book supplies for
   chapters 1–5, each row citing the section it comes from; check there
   first before re-deriving something already worked out.
2. **Write it as `[[nodiscard]] inline` functions or a small struct**, in
   whichever file under `include/newton/analytic/` it belongs to (a new
   file if it's a new physical system; an added function if it's a variant
   of something already there, the way `analytic/drag.hpp` groups the
   linear and quadratic drag results together). Match the existing style:
   parameters and members carry units in `//` comments; anything using
   `std::sin`/`exp`/`sqrt`/`log`/etc. is plain `inline` (not `constexpr`,
   since none of those are `constexpr` in C++17 — see any analytic file's
   own C++ note for this exact trade-off); pure-arithmetic accessors
   (`SHM::period()`, `LinearDragParams::tau()`) are `constexpr`.
3. **Add a numeric-vs-analytic test** in `tests/test_numeric_vs_analytic.cpp`
   — the library's soul, per DESIGN.md. The pattern is always: build the
   matching `forces::` functor (or `forces::sum(...)` combination), run
   `simulate()` with a suitable `Method` and `dt`, and compare the
   resulting `Trajectory` against your new closed form using
   `tests/test_util.hpp`'s tolerance-based comparison — never `==` on a
   computed double. Choose the tolerance from the integrator's own order:
   RK4 at `dt = T/100` should match an analytic oracle to something like
   `1e-6`–`1e-8` relative error; a first-order method needs a much looser
   tolerance for the same `dt`. `docs/research/numerical-integration.md`
   §7.4 has the reasoning behind combined absolute/relative tolerances if
   you need to pick one from scratch.
4. **If no closed form exists** (Taylor's own example: full 2D quadratic
   drag), you can't do step 3 as written — fall back to the self-convergence
   check instead: run the same problem at `dt` and `dt/2`, confirm the
   difference shrinks at the rate the integrator's order predicts. This
   needs no analytic solution at all and is exactly the technique
   `tests/test_integrators.cpp`'s convergence-order tests already use (see
   `notes/03-numerical-methods.md`'s "how the test suite validates all of
   this" section).

## The v0.2+ roadmap

Everything below is deliberately **not implemented** in this v0.1 library,
with one exception: the `potentials` item has since been built and is marked
as such. DESIGN.md lists these as the roadmap, condensed here from the fuller
treatment in `docs/research/taylor-ch1-7-map.md` and
`docs/research/numerical-integration.md`. Each is a real chunk of design
work, not a quick add.

### N-particle `System` with pairwise forces (Newton III by construction)

Everything in v0.1 — `Particle`, `State`, `simulate()` — is single-particle.
A multi-body `System` would hold `std::vector<Particle>` plus a list of
*pairwise* forces (`Vec3 f(t, particle_a, particle_b)`, a genuinely
different concept from the current per-particle `Force`), and the crucial
design choice is applying `+f` to one particle and `-f` to the other
*inside the same loop iteration* — so momentum conservation
(`analysis::total_momentum`, already in `include/newton/analysis/momentum.hpp`)
becomes an exact machine-precision guarantee of the code's structure,
rather than something that could silently break from a forgotten `-=` or a
double-counted pair. This is "Newton III by construction" — the
architecture makes the bug impossible instead of merely testing for it. It
would need a new `simulate()` variant (or overload) that steps every
particle each iteration, and `System::centerOfMass()` /
`totalAngularMomentum()` convenience methods alongside the free functions
`analysis::center_of_mass` and `analysis::angular_momentum` already have.
A single particle becomes the trivial `N = 1` case of this — the research
doc's explicit recommendation is to build the N-particle version once and
never maintain two integrator code paths.

### Cubic-Hermite dense output + event detection by root-finding

`simulate_until()` (`include/newton/simulate.hpp`) already refines an
event crossing accurately, but it does so by *bisecting the step size and
re-running the real integrator* up to 40 times per crossing — correct, but
it costs extra force evaluations. Cubic-Hermite dense output is a cheaper
alternative: given the two endpoint states of a step, `(r_n, v_n)` and
`(r_{n+1}, v_{n+1})`, the four Hermite basis functions
(`h00, h10, h01, h11`, a function of `s = (t - t_n)/dt ∈ [0,1]`)
interpolate `r(t)` and `v(t)` *anywhere inside the step* with no extra
force evaluations at all — and the interpolation error is `O(dt^4)`,
strictly better than Verlet's own `O(dt^2)` global error, so it costs
nothing in accuracy. Two things become possible on top of it: recording
trajectory samples at *exact* requested output times (decoupling physics
`dt` from output `dt`, rather than only the current `record_every` stride),
and root-finding an event crossing (`g(t) = 0` for a general user function,
not just the two built into `newton::events`) directly on the interpolant
via bisection or Brent's method — cheaper than re-stepping, and reusable
for any event the way the current bisection logic is currently duplicated
conceptually inside `simulate_until` alone.

### Boris pusher for magnetic forces

`forces::MagneticForce` (Taylor §2.5/2.7) is special: because `F = qv × B`
is always perpendicular to `v`, it does zero work, and `|v|` is *exactly*
conserved by the true dynamics — one of the cleanest hard invariants in the
whole book. RK4 and velocity Verlet only *approximate* that conservation;
over many cyclotron periods you'll see a slow spiral in or out of the
gyroradius even though the physics says there shouldn't be one (this is
exactly what `analytic::magnetic`'s `cyclotron_radius`/`cyclotron_period`
in `include/newton/analytic/magnetic.hpp` exist to detect). The **Boris
pusher** is a purpose-built leapfrog-style scheme for `F = qv × B` (and its
Lorentz-force extension) that rotates the velocity by exactly the right
angle each half-step using a `tan(θ/2)` half-angle trick, and conserves
`|v|` and the gyroradius to machine precision over arbitrarily long runs —
it's the standard integrator for this exact problem in plasma physics for
that reason. Adding it means a new integrator file specific to
`MagneticForce`/`LorentzForce` (it doesn't generalize to arbitrary forces
the way RK4 and Verlet do), and it pairs naturally with
`examples/03_charged_particle.cpp` as a three-way |v|-drift comparison
against the RK4 and Verlet runs already there.

### Adaptive RK45

`simulate()`'s fixed-step contract is a deliberate choice, not an
oversight — see `notes/03-numerical-methods.md`'s discussion of why the
symplectic pair (semi-implicit Euler, velocity Verlet) needs a *constant*
`dt` to keep their bounded-energy guarantee. Adaptive stepping (embedded
Dormand–Prince 5(4) with an error estimator, FSAL reuse of the last stage
as the next step's first, and a PI step-size controller of the form
`dt_new = dt · clamp(S·(tol/err)^(1/5))`) belongs in a genuinely separate
code path restricted to the RK family — it should never be reachable
through `Method::SemiImplicitEuler` or `Method::VelocityVerlet`, and its
options struct would carry `rtol`/`atol` instead of a fixed `dt`. It earns
its complexity whenever the fastest timescale in a problem changes a lot
over the run — an eccentric orbit needing a tiny step near perihelion and
a coarse one far away is the textbook case — something none of the four
current fixed-step methods handle gracefully without hand-tuning `dt` down
for the whole run.

### A `potentials` module with numerical gradients and curl tests — implemented 2026-11-24

**Status: built.** This was the first v0.2 item to move out of the roadmap, and it
came out as two modules rather than one, because the numerical calculus is useful
on its own:

- `include/newton/fields/derivatives.hpp` — `gradient`, `jacobian`, `curl`,
  `divergence` by 5-point central differences (Richardson-extrapolated, as this
  section originally proposed), for any callable field.
- `include/newton/fields/line_integral.hpp` — adaptive Gauss–Kronrod line
  integrals and circulation, which the roadmap text did not anticipate: the
  "curl sampled at random points" test alone cannot see a force like
  `(-y, x, 0)/(x² + y²)`, whose curl is zero everywhere off an axis yet which
  does work around a loop that encircles it, so the classifier also integrates
  around closed loops.
- `include/newton/potentials/conservative.hpp` — `check_conservative`, the
  "is this force conservative" test (a `Verdict` plus the evidence behind it),
  `probe_dependence`, and the split helpers `velocity_independent_part` and
  `velocity_dependent_part`.
- `include/newton/potentials/potential_from_force.hpp` — `potential_from_force`
  and `find_potential`: U from F, as a line integral from a reference point.
- `include/newton/potentials/force_from_potential.hpp` — `force_from_potential`:
  F from U, `-∇U` by numerical gradient, as a ready-made Force. This is the
  "define a new potential in one line and get a working force for free" idea
  from the original plan.
- `include/newton/analysis/energy.hpp` gained a time-dependent `U(t, r)` form for
  `total_energy`/`energy_drift`, plus `power`, `work_done` and `energy_balance`
  for the non-conservative side of the ledger.
- Tests: `tests/test_fields.cpp`, `tests/test_line_integrals.cpp`,
  `tests/test_potentials.cpp`, and additions to `tests/test_analysis.cpp`.
  Examples: `examples/07_potential_from_force.cpp` and
  `examples/08_force_from_potential.cpp`.

`notes/06-fields-and-potentials.md` explains how all of it works and, just as
important, what its answers can and cannot prove; `notes/05-library-reference.md`
has the API with compiled snippets. Two choices differ from the original sketch,
on purpose: the potential is a plain callable (`double U(const Vec3&)`, or
`U(t, r)` for a time-dependent one) rather than a `Potential` class, so it plugs
into `energy_drift` with no adapter; and the curl test is a classifier with a
report, not a bare assertion, because a single "is the curl near zero" check
needs a tolerance, and the right tolerance depends on the field.

**Still on the roadmap** from this item:

- *Analytic gradients supplied by a potential.* `force_from_potential` always
  differentiates numerically (a few parts in 1e12 on well-scaled fields, at 12
  evaluations of U per force call, wrong within about two steps of a kink or singularity, and needing
  `length_scale` and `center` hints when the structure is tiny or far from the origin).
  A potential that can also supply its own gradient (a struct with both `U` and
  `grad`, or an overload taking a second lambda) would skip the numerics where
  the formula is known, and give a way to test the numerical path against exact
  forces on every run.
- *A reference-at-infinity mode.* Gravity, Coulomb, Yukawa and Morse potentials
  are conventionally zero at infinity, which a finite line integral cannot start
  from. Today you pass a regular `r_ref` and `U_ref = U_textbook(r_ref)` (the
  recipe in `potential_from_force`'s doc comment). A built-in mode could do that
  for the common central forces, or integrate in from a large radius with a
  tail correction.

### 1D potential toolkit (turning points, period by quadrature)

Taylor §4.6's 1D energy method turns any 1D potential `U(x)` into a
description of the motion without ever integrating in time:
`v(x) = ±sqrt(2(E - U(x))/m)`, turning points where `U(x) = E`, equilibria
where `U'(x0) = 0` classified stable or unstable by the sign of `U''(x0)`,
and — the most numerically interesting piece — the period of bounded
motion by direct quadrature, `T = 2∫ dx / sqrt(2(E-U(x))/m)` between the
two turning points. That integral has a genuine `1/sqrt(...)` singularity
at each endpoint (the particle's speed goes to zero exactly there), so a
naive Simpson's rule fails right where it matters most; the standard fix
is a Gauss–Chebyshev-style substitution (`x = midpoint + halfwidth·cos θ`)
that absorbs the singularity into the change of variables. This toolkit is
the natural generalization of `analytic::SHM`'s small-oscillation formulas
already in `include/newton/analytic/oscillator.hpp` — `SHM` is the special
case `U(x) = ½kx²`; a `smallOscillationFrequency(U, m, x0) = sqrt(U''(x0)/m)`
function would recover it from an arbitrary nonlinear `U` the same way
Taylor derives it as the leading-order behavior near any smooth potential
minimum (the same motivation `forces::Spring`'s own doc comment cites for
why Hooke's law shows up "everywhere").

### Fourier tools for §5.7–5.9

Taylor's starred sections decompose a periodic (but non-sinusoidal) drive
into harmonics, `f(t) = Σ_n [a_n cos(nωt) + b_n sin(nωt)]`, and superpose
each harmonic's steady-state response using the same
`analytic::DrivenOscillator` formulas already implemented for a single
frequency (`include/newton/analytic/oscillator.hpp`) — just evaluated at
`ω → nω` for each `n` and summed. Validating that superposition needs a
small DFT (a textbook radix-2 FFT is plenty for the book's square-wave and
sawtooth examples) and a Parseval/RMS check (`⟨x²⟩ = ½Σ Aₙ²`). This
belongs in `newton::analysis` alongside the existing measurement functions
in `include/newton/analysis/oscillation.hpp` — it would consume a
`Trajectory`'s `r.x` samples exactly the way `find_peaks` and
`measure_period` already do — and pairs naturally with a new
`forces::SquareWaveDrive`/`SawtoothDrive`, which is really just a
compile-time-known sum of `forces::SinusoidalDrive` terms and could reuse
`forces::sum(...)` from `include/newton/forces/compose.hpp` almost as-is.

### The Lagrangian hook (constraint solver, then dual-number autodiff)

Covered in more depth at the end of `notes/02-physics-to-code.md`'s
Chapter 6–7 section (condensed from
`docs/research/taylor-ch1-7-map.md` §6.4–6.5 and §7.4–7.5) — the two-step
plan is worth restating here as the roadmap item it is. **Step one**: a
Cartesian holonomic-constraint solver — `g(r) = 0` constraints resolved by
a Lagrange multiplier, solving `[M Jᵀ; J 0][a; λ] = [F_applied; -J̇v]` each
step — added to the existing force engine, so pendulums, beads on wires,
and Atwood machines become reachable *without leaving the Force paradigm
or touching a single existing integrator*. This is the recommended first
step precisely because it's genuinely Newtonian: no symbolic
differentiation, no generalized coordinates, just one more kind of force
(a constraint force) computed each step alongside gravity and drag.
**Step two**, only once step one exists and only if the generalized-
coordinates formalism itself is wanted: a true `newton::lagrangian` module
needs the mass matrix `M_ij = ∂²T/∂q̇_i∂q̇_j` and the generalized-force
vector computed from a user's `T(q, q̇, t)` and `U(q, t)`, and getting
those derivatives *exactly*, without asking a student to hand-differentiate
their own kinetic energy expression, is what forward-mode automatic
differentiation via dual numbers is for — a `Dual<double>` type propagated
through ordinary arithmetic gives exact first derivatives for free, and
nesting it once (`Dual<Dual<double>>`) gives the second derivatives `M_ij`
needs. That's the reason this is flagged as a roadmap item rather than a
quick add: it would require `Vec3` and the rest of the core to be
templated on the scalar type from the start (`Vec3<T>` instead of a
hard-coded `Vec3` of `double`), which is a straightforward decision to make
on day one and a genuine rewrite to retrofit later — exactly the trade-off
`docs/research/taylor-ch1-7-map.md`'s closing recommendation calls out.
