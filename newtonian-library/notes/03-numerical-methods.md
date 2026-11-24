# Numerical Methods: The Four Integrators

`newton` ships four one-step integrators
(`include/newton/integrators/euler.hpp`,
`semi_implicit_euler.hpp`, `velocity_verlet.hpp`, `rk4.hpp`), all with the
identical shape

```cpp
template <class Force>
State xxx_step(const State& s, double mass, double dt, Force&& force);
```

and one fixed-step driver, `simulate()` in `include/newton/simulate.hpp`,
that calls one of them `n` times and hands back a `Trajectory`. This note
explains what distinguishes the four — not just the update formulas, which
you can read straight off the header comments, but *why* they behave so
differently, how to pick a step size, how to tell when you picked a bad
one, and how the test suite proves all of this rather than just asserting
it. Everything here condenses
`docs/research/numerical-integration.md`, which has the full derivations —
go there when "why" isn't enough and you want the algebra.

## Why compare four integrators at all?

Because the deceptively simple decision "how do I turn `F = ma` into
`x(t)`" has no single right answer — every method trades accuracy, cost,
and long-run behavior against each other differently, and the trade that's
right for a two-second projectile is wrong for a million-orbit planetary
simulation. `examples/06_integrator_comparison.cpp` is built to make that
trade-off visible on one plot: SHM, all four methods, `dt = T/100`, 200
periods, final relative energy error printed per method. Keep that example
in mind through the rest of this note — it's the live demonstration of
everything below.

## The §2.0 table

Every number below is for the undamped harmonic oscillator `x'' = -ω²x`,
with `k = ω·dt` — a single dimensionless number that governs essentially
every stability and accuracy statement an integrator makes. At `dt = T/100`
(a good default), `k ≈ 0.0628`.

| Method | Order | Force evals/step | Stability limit | Energy behavior | At dt = T/100 |
|---|---|---|---|---|---|
| Explicit Euler | 1 | 1 | none — unstable for every `dt` | grows as `E_n = E_0 (1+k²)^n` | **+48% energy per period** |
| Semi-implicit Euler | 1 | 1 | `k < 2` | bounded, oscillates | wobble ≈ `k/2` ≈ **3%**, forever, no trend |
| Velocity Verlet | 2 | 1 (2 in this library's standalone step — see below) | `k < 2` | bounded, oscillates | wobble ≈ `k²/4` ≈ **1.0e-3**; phase drifts linearly |
| RK4 | 4 | 4 | `k < 2√2 ≈ 2.83` | slow monotone decay | `ΔE/E ≈ -k⁶/72` ≈ `-8.5e-10` per step |

Read the "order" column as "how many powers of `dt` the *global* error
scales with" — halving `dt` for a 1st-order method halves the error;
halving `dt` for RK4 divides the error by 16. Read "force evals/step" as
the actual cost: RK4 is 4× the per-step work of the other three, which
matters when you're integrating millions of steps.

## Explicit Euler: why it spirals outward

`euler_step` (`include/newton/integrators/euler.hpp`) is the most literal
possible discretization — evaluate the force at the *old* state, take one
linear step along it:

```
a_n     = F(t_n, r_n, v_n) / m
r_{n+1} = r_n + dt * v_n          <- uses v_n, not v_{n+1}
v_{n+1} = v_n + dt * a_n
```

There are three equivalent ways to see why it always grows, and it's worth
holding all three in your head at once because they're the same fact seen
from algebra, linear algebra, and geometry:

1. **Algebra.** On `x'' = -ω²x`, direct substitution shows `E_{n+1} = E_n ·
   (1 + ω²dt²)` exactly — every single step, the energy is multiplied by a
   number strictly greater than 1. Over `n` steps, `E_n = E_0 (1+k²)^n`,
   which is exponential growth in simulated time. Shrinking `dt` shrinks
   the *rate* but never removes the exponential — there is no stable `dt`.
2. **Linear algebra.** In scaled coordinates `(ωx, v)`, one Euler step is
   the linear map with determinant `1 + k² > 1`. The true flow of a
   Hamiltonian system exactly preserves phase-space area (Liouville's
   theorem); forward Euler inflates it by `1+k²` on every step.
3. **Geometry — the picture worth drawing once and never forgetting.** In
   the `(ωx, v)` plane, the true motion of an oscillator is a circle.
   Forward Euler steps along the **tangent line** to that circle at the
   current point. A tangent line to a circle never re-enters the circle —
   it lies entirely outside except at the point of tangency — so every
   single step lands strictly farther from the origin than the state it
   started from. There's no clever choice of `dt` that fixes this; it's
   baked into using a straight line to approximate a curve that bends away
   from you.

`euler_step` exists in this library specifically as that demonstration —
the doc comment calls it a "pedagogical foil", and `tests/test_integrators.cpp`
literally asserts that its energy grows, as a check that the energy
diagnostic itself is wired correctly. Never reach for it as a simulation
default.

## What "symplectic" means, in the terms this library needs

`semi_implicit_euler_step` (a.k.a. Euler–Cromer,
`include/newton/integrators/semi_implicit_euler.hpp`) makes exactly one
change from forward Euler — advance velocity first, then use the **new**
velocity for the position update:

```
a_n     = F(t_n, r_n, v_n) / m
v_{n+1} = v_n + dt * a_n
r_{n+1} = r_n + dt * v_{n+1}       <- v_{n+1}, not v_n
```

Same cost (one force evaluation, fully explicit — "semi-implicit" is a
slight misnomer, nothing is actually solved implicitly), same first order
per step — but that one reordering makes the per-step map on `(ωx, v)`
phase space have determinant exactly 1: it neither inflates nor deflates
phase-space area. A map with that property is called **symplectic**, and
symplectic maps have a property worth understanding rather than just
citing, because it explains everything about how these integrators behave
over long runs:

A symplectic integrator does not approximately solve your actual
Hamiltonian `H`. It solves — essentially exactly, up to exponentially small
corrections — a nearby "**shadow**" Hamiltonian `H̃ = H + dt^p · (something)
+ ...`, a slightly different physical system that happens to be very close
to the one you asked for. Energy is *exactly* conserved for the shadow
system (it's a genuine Hamiltonian flow), so the *true* energy — the one
you actually care about — is trapped forever within a band around a
constant, of relative half-width roughly `k/2` for semi-implicit Euler.
It oscillates. It never trends. That's the "wobble ≈ 3%, forever, no
drift" row in the table above: the 3% is the width of the band you're
stuck orbiting inside, not an error that keeps growing.

Two things worth knowing about that guarantee:

- **It needs a constant `dt`.** The shadow-Hamiltonian argument assumes
  the map you're iterating is the same map every step. Change `dt`
  mid-simulation and you jump to a *different* shadow Hamiltonian —
  the bounded-energy property resets, and drift can come back. This is
  exactly why `include/newton/simulate.hpp` only offers fixed-step
  integration (see its own header comment) — an adaptive step size would
  silently throw away the one guarantee that makes the symplectic methods
  worth using.
- **For a central force, semi-implicit Euler conserves angular momentum
  exactly**, to machine roundoff, not just approximately — a one-line
  algebraic proof (each step's position update uses `v_{n+1}`, and
  `r_n × v_{n+1} = r_n × (v_n + dt·a(r_n))= r_n × v_n` whenever `a` is
  parallel to `r_n`). That makes it a surprisingly good *cheap* orbit
  integrator despite being "only" first order.

## Velocity Verlet, and the subtlety this library resolves for you

`velocity_verlet_step` (`include/newton/integrators/velocity_verlet.hpp`)
is the textbook second-order symplectic workhorse, written in
kick-drift-kick form:

```
a_n     = F(t_n, r_n, v_n) / m
v_half  = v_n + 0.5*dt*a_n            (half kick)
r_{n+1} = r_n + dt*v_half             (full drift)
a_{n+1} = F(t_{n+1}, r_{n+1}, ...) / m
v_{n+1} = v_half + 0.5*dt*a_{n+1}     (half kick)
```

It's exact to roundoff for constant acceleration at *any* `dt` — the very
first thing this library's tests check (`tests/test_integrators.cpp`), and
a fast way to catch a sign error or a dropped `0.5`. Its energy behavior is
one order better than semi-implicit Euler's: bounded oscillation of
relative amplitude `≈ k²/4`, still with zero secular drift. But there's a
second trade you're always making when you pick Verlet: its numerical
frequency is slightly off (`ω_num ≈ ω(1 + k²/24)`), so while energy stays
perfect, *phase* drifts linearly in time — after enough periods your
oscillator is measurably out of sync with the true one even though its
energy still looks great. The research doc's phrase for this is "gets
structure and statistics right forever, gets long-run phase wrong" — the
mirror image of RK4's trade below.

**The velocity-dependence subtlety.** The derivation above quietly assumes
`F = F(t, r)` — no `v` in the force. The moment a force depends on
velocity (drag, damping, a magnetic field), the final half-kick becomes
*implicit*:

```
v_{n+1} = v_half + 0.5*dt * F(t_{n+1}, r_{n+1}, v_{n+1}) / m
                                                  ^^^^^^^ on both sides
```

You cannot literally write that line — `v_{n+1}` is what you're solving
for. `velocity_verlet_step` in this library resolves it with a
**predictor–corrector**: predict `v_{n+1}` with an explicit,
Euler-quality estimate, use *only* that prediction to evaluate the force,
then correct with a real half-kick:

```
v_pred  = v_half + 0.5*dt*a_n                    (predictor)
a_{n+1} = F(t_{n+1}, r_{n+1}, v_pred) / m         (the ONLY place v_pred is used)
v_{n+1} = v_half + 0.5*dt*a_{n+1}                 (corrector, using the real a_{n+1})
```

The reason this doesn't quietly downgrade the method to first order (worth
actually understanding, not just trusting): `v_pred` is only a first-order
estimate, so it carries `O(dt²)` error — but that error only ever feeds
into a *force evaluation*, and the force's contribution to `v_{n+1}` is
itself scaled by `0.5*dt` in the corrector. An `O(dt²)` error, multiplied
by `dt`, is `O(dt³)` — exactly Verlet's own local truncation error, so no
order is lost. Get this wrong — e.g. skip the predictor and reuse `v_half`
or `v_n` directly as the force's velocity argument — and you silently lose
an order of accuracy; `tests/test_integrators.cpp`'s convergence-order test
(next section) is specifically designed to catch that class of bug by
watching the fitted order for a drag problem drop from ~2 toward ~1.

One cost worth knowing: this standalone step function evaluates the force
**twice** per call (once for `a_n`, once for `a_{n+1}` using `v_pred`),
because a bare step function has nowhere to cache `a_n` from the *previous*
call the way a production integrator object would. A cached-acceleration
variant that gets back to one evaluation per step is a documented v0.2
optimization (see `notes/04-extending.md`), left out here so this file's
predictor–corrector logic stays the whole story with no state hidden
elsewhere.

Symplecticity itself is not preserved once `F` depends on `v` — but a
velocity-dependent force is generally dissipative anyway (drag removes
energy from the system), so there was no exact conservation law to lose in
the first place. What the predictor–corrector protects is the thing that
actually mattered: the accuracy order.

## RK4: the other trade

`rk4_step` (`include/newton/integrators/rk4.hpp`) is the classic 4-stage
Runge–Kutta method applied to the coupled system `y' = (v, a)`. Four force
evaluations per step (`k1v..k4v`) buy fourth-order global accuracy —
sixteen times better for every halving of `dt` than Verlet's four times.
But look at the energy row in the table: RK4 is **not** symplectic and
**not** time-reversible. Its amplification factor on an oscillator is
slightly less than 1, so it introduces a slow, monotone, secular energy
*decay* — the opposite failure mode from Verlet's bounded wobble. At sane
step sizes (`dt = T/100`) that decay is minuscule, `≈ -8.5e-10` per step —
but it never reverses, and over a genuinely long run (millions of periods
of a planetary orbit, say) it accumulates into a real, visible spiral-in
that a symplectic method would never produce. This is why the research doc
frames RK4 and Verlet as mirror images: RK4 is "phase-perfect, energy-wrong
over aeons"; Verlet is "energy-perfect, phase-wrong over aeons." Pick RK4
when pointwise accuracy — hitting a specific position at a specific time,
or handling a velocity-dependent force with full accuracy rather than
Verlet's second-order patch — matters more than long-run structure. Pick
Verlet for millions of periods of genuinely conservative dynamics, where
bounded energy matters more than exact phase.

The one implementation detail worth internalizing from
`rk4.hpp` itself: **stage order matters**. At each stage you must compute
the position-slope `k_ir` before the velocity-slope `k_iv`, because `k_iv`'s
position argument uses the `k_ir` you just computed (`k2v` needs `k1r`,
`k3v` needs `k2r`, `k4v` needs `k3r`). Swap that order — reuse `k1r` where
`k2r` belongs, say — and the method still compiles, still runs, and
silently degrades from fourth order to second. Nothing crashes. Nothing
warns. Only a convergence-order test (below) catches it.

## Choosing a step size

`dt` is never picked in absolute terms — it's picked relative to the
fastest timescale in your problem. For an oscillator, that timescale is
the period `T = 2π/ω`, and the standard rule of thumb is

```
dt = T / 100        (good accuracy — the default this library's examples use)
dt = T / 20         (rough minimum before results get visibly wrong)
dt = T / 1000       (high precision)
```

This is *not* the same thing as the stability limit. Stability tells you
where the numerics explode (`k < 2` for the symplectic pair, `k < 2√2` for
RK4); accuracy tells you where the numbers are actually *right*. At `k =
1.9` (just under Verlet's stability wall), Verlet is technically stable but
its numerical frequency is off by ~32% (exact: the per-step map gives
`cos(ω_num·dt) = 1 − k²/2`, so `ω_num/ω − 1 = 2·asin(k/2)/k − 1 ≈ 0.319`
at `k = 1.9`) and its energy wobble is ~90% — a
stable-but-useless simulation. `dt = T/100` gives `k ≈ 0.063`, deep inside
the region where both stability and accuracy hold comfortably. Pick `dt`
from the accuracy rule; treat the stability limit purely as a safety net
you should never be operating anywhere near.

For a problem that isn't a pure oscillator, estimate the fastest relevant
timescale directly: linear drag's `τ = m/b` (Taylor §2.2, exposed as
`LinearDragParams::tau()` in `include/newton/analytic/drag.hpp`), or a
spring-mass system's `T = 2π√(m/k)`. When several timescales are present at
once (e.g. a stiff spring on a light mass, or strong damping `β ≫ ω0`),
the *fastest* one sets `dt` for the whole run even if it's not the
timescale you actually care about — that's what "stiff" means, and it's
worth knowing the name even though this v0.1 library has no stiff-system
remedy beyond "pick a smaller `dt`" (implicit integrators are a v0.2+ idea,
not implemented here).

## Detecting a bad `dt`

Two cheap, general-purpose diagnostics, both of which this library gives
you the building blocks for:

1. **Energy monitor.** Run `include/newton/analysis/energy.hpp`'s
   `energy_drift(trajectory, U)` over your simulated run and look at the
   shape of the result, not just its final value:
   - **Bounded oscillation with no trend** is the expected, healthy
     signature for the symplectic methods (semi-implicit Euler, velocity
     Verlet) — this is the "sawtooth" pattern the research doc calls
     normal, not a bug.
   - **A slow, monotone slope** is RK4's expected signature — tiny at a
     sane `dt`, but a real, non-reversing trend.
   - **Unbounded exponential growth** means either explicit Euler (working
     as designed — never use it for anything that should conserve energy)
     or a `dt` past a stability limit on one of the other three methods.
   - **A sudden trend that starts small and grows without bound** on a
     method that's supposed to be bounded (symplectic or RK4) is the alarm
     that says `dt` is too big — shrink it and re-check.
2. **Halve-and-compare — the definitive test.** Run the same simulation to
   the same final time with `dt` and with `dt/2`. If the two trajectories
   differ by more than the tolerance you actually need, `dt` was too big,
   full stop — this needs no analytic solution at all (it's
   *self*-convergence), so it works even on problems this library has no
   closed form for, like full 2D quadratic drag.

## How the test suite validates all of this

`tests/test_integrators.cpp` turns every claim in this note into an
assertion, using three independent, complementary strategies — all three
have to pass, because each one catches a different class of bug:

- **Analytic comparison.** Constant force is the sharpest first check:
  Verlet and RK4 must reproduce it to roundoff (`~1e-12`) at *any* `dt`,
  and forward Euler must **not** be exact — a single test that already
  discriminates a correct implementation from a subtly broken one. SHM
  against `analytic::SHM` (`include/newton/analytic/oscillator.hpp`) is the
  workhorse case used for everything below.
- **Convergence order.** Halve `dt` several times on a problem with a known
  exact answer, compute the error at each `dt`, and fit the observed order
  `p = log2(e_j / e_{j+1})`. The suite asserts the fitted order lands
  within `±0.2` of each method's advertised order — 1 for both Euler
  variants, 2 for Verlet, 4 for RK4. This is the test that catches the
  stage-order bug in RK4 and the missing-predictor bug in Verlet described
  above: both are silent, both still "run", and both show up as nothing
  else *except* a fitted order that's wrong.
- **Energy behavior over many periods.** Run 100 periods of SHM with each
  method and check the *shape* of the energy-drift curve matches the table
  at the top of this note: Euler's relative energy growth exceeds 10%,
  semi-implicit Euler and Verlet stay bounded under 1%, RK4's drift stays
  tiny. This directly exercises `analysis::energy_drift`.

`tests/test_integrators.cpp` also checks the time-bookkeeping rule from
`include/newton/simulate.hpp`'s header comment — that `t` is recomputed as
`t0 + i*dt` every iteration rather than accumulated with `t += dt` — by
asserting the final state's `t` equals `t_end` to `1e-12`. That rule
matters because `dt` is generally not exactly representable in binary
floating point, so a naive `t += dt` loop commits a fresh rounding error
every step, and those compound; recomputing `t` from scratch each time
caps the error instead of letting it grow with the number of steps.

`examples/06_integrator_comparison.cpp` is the live version of the "energy
behavior over many periods" test above, scaled up to 200 periods and all
four methods on one CSV (`build/energy_drift.csv`, columns `t, E_euler,
E_semi, E_verlet, E_rk4`) — run it, plot it with
`python3 scripts/plot_csv.py build/energy_drift.csv`, and you'll see this
entire note as one picture: Euler's curve runs away, semi-implicit Euler's
and Verlet's oscillate in a bounded band, and RK4's decays almost
imperceptibly.
