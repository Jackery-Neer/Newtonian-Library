# Numerical Integration of Newtonian ODEs
### A design guide for a hand-rolled C++ simulation library

---

## 0. Scope and notation

We integrate

```
m x''(t) = F(t, x, v),    x(t0) = x0,  v(t0) = v0
```

where `x, v, F` live in R^d (d = 1, 2, 3 for a single particle; d = 3N for N particles).

Throughout: `dt` = step size, `n` = step index, `a = F/m`, `ω` = angular frequency of the fastest oscillation, `T = 2π/ω` its period, and `k ≡ ω·dt` the dimensionless step (this single number controls nearly every stability and accuracy statement below).

Convention: `x_n ≈ x(t_n)`, `t_n = t0 + n·dt`.

---

## 1. Recasting the second-order ODE as a first-order system

Every integrator below is really an integrator for `y' = f(t, y)`. The recast is mechanical:

```
Let y = [ x ]        (dimension 2d)
        [ v ]

Then y' = f(t, y) = [ v              ]
                    [ F(t, x, v) / m ]
```

That is: **`dx/dt = v` and `dv/dt = a(t, x, v)`**, a coupled pair. The first equation is trivial (its "force" is just `v`); this is why second-order-specific integrators (Verlet, Störmer) can be cheaper than generic ones — they exploit the fact that half the system is trivial.

### Representation choices in C++

Three viable layouts, in increasing generality:

```cpp
// (A) Explicit pair — clearest, best for a teaching/first library.
struct State {
    Vec3 x;   // position, metres
    Vec3 v;   // velocity, m/s
};

// (B) N particles, array-of-structs. Simple, cache-decent, easy to read.
struct System {
    std::vector<Vec3>  x, v, a;      // size N
    std::vector<double> invMass;     // 1/m, size N; 0.0 == immovable
};

// (C) Flat state vector — generic, plugs into generic RK code.
using Y = std::vector<double>;   // [x0..x_{d-1}, v0..v_{d-1}]
```

**Recommendation:** build around (B). A single particle is `N = 1`; you never need a second code path. Layout (C) is tempting for genericity but forces you to hand-write index arithmetic and hurts readability; do not adopt it unless you later want to plug in a generic ODE solver.

**Store `invMass`, not `mass`.** You always need `F/m`, never `m` alone (except for kinetic energy and momentum diagnostics, where you can store both). `invMass = 0.0` gives you infinitely heavy / pinned particles for free with no branches.

**Variable mass** (rockets): `m` becomes `m(t)`. Do **not** integrate `d(mv)/dt = F` naively — see §7.8. Provide an overload where the user supplies acceleration directly, `a(t, x, v)`, so they can write `a = u*mdot/m(t) - g` themselves.

---

## 2. The integrators

### 2.0 Quantitative summary (memorise this table)

All entries for the undamped harmonic oscillator `x'' = -ω²x`, with `k = ω·dt`.

| Method | Global order | Force evals / step | Stability limit | Energy behaviour | Concrete numbers at `dt = T/100` (k ≈ 0.0628) |
|---|---|---|---|---|---|
| Forward Euler | 1 | 1 | **none** — unstable for all `dt` | grows as `E_n = E_0 (1+k²)^n` | **+48% energy per period**; amplitude ×1.22 per period |
| Symplectic Euler | 1 | 1 | `k < 2` (`dt < 0.318 T`) | bounded, oscillates | relative energy wobble ≈ `k/2` ≈ **3%**, forever, no drift |
| Velocity Verlet | 2 | 1 (with cached `a`) | `k < 2` (`dt < 0.318 T`) | bounded, oscillates | relative energy wobble ≈ `k²/4` ≈ **1.0e-3**; phase error ≈ `1.0e-3` rad/period, *accumulating* |
| RK4 | 4 | 4 | `k < 2√2 ≈ 2.828` (`dt < 0.450 T`) | slow monotone **decay** | `ΔE/E ≈ -k⁶/72` per step ≈ `-8.5e-10`; ≈ `-8.5e-5` after 1000 periods |
| Backward Euler | 1 | 1 + solve | unconditional | monotone **decay** (artificial damping) | strongly damps; unusable for conservative physics |
| Implicit midpoint | 2 | 1 + solve | unconditional | bounded (symplectic) | good stiff + conservative option, but needs a nonlinear solve |

The last two are included for context (§7.1); implement the first four.

---

### 2.1 Explicit (forward) Euler

**Update equations** — everything evaluated at the *old* state:

```
a_n     = F(t_n, x_n, v_n) / m
x_{n+1} = x_n + dt * v_n
v_{n+1} = v_n + dt * a_n
t_{n+1} = t_n + dt
```

Note carefully: `x_{n+1}` uses `v_n`, **not** `v_{n+1}`. That one detail is the entire difference from §2.2.

- **Order:** local truncation error `O(dt²)`, global `O(dt)`.
- **Cost:** 1 force evaluation, no extra storage.
- **Energy:** unconditionally divergent on oscillators (below).

#### Why it spirals outward — three ways to see it

**(a) Exact algebra.** For `x'' = -ω²x` with `E = ½(v² + ω²x²)`:

```
E_{n+1} = ½(v_n - dt ω² x_n)² + ½ω²(x_n + dt v_n)²
        = ½v_n² + ½ω²x_n²  +  ½dt²ω²(ω²x_n² + v_n²)     [cross terms cancel exactly]
        = E_n · (1 + ω² dt²)
```

So `E_n = E_0 (1 + k²)^n = E_0 (1+k²)^{t/dt} ≈ E_0 · exp(ω² dt · t)`. Energy grows **exponentially**, with a rate proportional to `dt`. Halving `dt` halves the growth rate — but never eliminates it. There is no stable `dt`.

**(b) Phase-space area.** In scaled coordinates `u = ωx`, the step is the linear map `M = [[1, k], [-k, 1]]`, with `det M = 1 + k² > 1`. Liouville's theorem says the true flow preserves phase-space area; forward Euler *inflates* it by `1+k²` every step. Eigenvalues are `1 ± ik`, modulus `√(1+k²) > 1`.

**(c) Geometry.** In the `(ωx, v)` plane the exact trajectory is a circle. Forward Euler steps along the **tangent** to that circle. A tangent line always lies outside the circle, so every step lands strictly farther from the origin. This is the picture to draw for yourself once and never forget.

**When to use it:** essentially never for mechanics. Legitimate uses: (i) as a pedagogical foil in your test suite — assert that its energy *does* grow at the predicted rate, which is a surprisingly good check that your energy diagnostic is correct; (ii) purely dissipative, non-oscillatory relaxation where the exact solution decays and you only need a rough answer. Do not ship it as a default.

---

### 2.2 Semi-implicit (symplectic) Euler — Euler–Cromer

**Update equations** — velocity first, then position with the *new* velocity:

```
a_n     = F(t_n, x_n, v_n) / m
v_{n+1} = v_n + dt * a_n
x_{n+1} = x_n + dt * v_{n+1}          // <-- v_{n+1}, not v_n
```

(The adjoint variant — position first with old `v`, then velocity using `F(x_{n+1})` — is equally symplectic. Pick one, document it, be consistent.)

"Semi-implicit" is a slight misnomer: nothing is solved implicitly, it is fully explicit and costs exactly as much as forward Euler. The name reflects that `x_{n+1}` is defined by an implicit-looking formula that happens to be evaluable.

- **Order:** global `O(dt)` — same as forward Euler.
- **Cost:** 1 force evaluation. Identical to forward Euler.
- **Energy:** bounded, no secular drift (for `F = F(x)`).

#### Why energy is bounded

The map in scaled coordinates is `M = [[1-k², k], [-k, 1]]`, `det M = 1` exactly — area preserving, i.e. **symplectic**. More strongly, the quadratic form

```
Q(x, v) = ω²x² + v² - ω²·dt·x·v
```

is **exactly conserved**, to machine precision, for all `n`. Since `E = ½(v² + ω²x²) = ½Q + ½ω²·dt·x·v`, and `ω|x||v| ≤ E`, the computed energy oscillates within a band of relative half-width `≈ k/2 = ω·dt/2` around a constant, forever. It never trends.

This is a special case of **backward error analysis**: a symplectic integrator does not approximately solve your Hamiltonian `H`; it *exactly* solves (to exponentially small error) a nearby "shadow" Hamiltonian `H̃ = H + dt^p·H_p + ...`. Energy of the shadow system is conserved exactly, so the true energy stays within `O(dt^p)` of it forever. **This argument requires a constant `dt`.** Change `dt` mid-run and you jump between shadow Hamiltonians, and drift returns.

#### Other exact conservation properties

- **Linear momentum:** if you compute pairwise forces and apply `+f` to i and `-f` to j (Newton's third law by construction), then `Σ m_i v_i` changes by `dt·Σ F_i = 0` exactly — conserved to floating-point summation error.
- **Angular momentum, central forces:** `L_{n+1} = x_{n+1} × v_{n+1} = (x_n + dt·v_{n+1}) × v_{n+1} = x_n × v_{n+1} = x_n × (v_n + dt·a(x_n)) = x_n × v_n = L_n`, using `a ∥ x`. **Exact.** This makes it a great cheap orbital integrator.

**When to use it:** the default "cheap and safe" integrator. Games, real-time interactive physics, large N-body where per-step cost dominates and 1st-order accuracy is acceptable, anything where you need long-time qualitative correctness more than pointwise accuracy. It handles velocity-dependent forces without any special casing (it just isn't symplectic anymore — and neither is the physics, so nothing is lost).

---

### 2.3 Velocity Verlet (and leapfrog)

**Update equations, textbook form** (assumes `F = F(t, x)`):

```
a_n     = F(t_n, x_n) / m
x_{n+1} = x_n + dt*v_n + 0.5*dt*dt*a_n
a_{n+1} = F(t_{n+1}, x_{n+1}) / m
v_{n+1} = v_n + 0.5*dt*(a_n + a_{n+1})
```

**Update equations, kick–drift–kick form** (implement this one):

```
v_half  = v_n + 0.5*dt*a_n            // half kick
x_{n+1} = x_n + dt*v_half             // full drift
a_{n+1} = F(t_{n+1}, x_{n+1}) / m     // <-- the ONLY force evaluation
v_{n+1} = v_half + 0.5*dt*a_{n+1}     // half kick
// carry a_{n+1} into the next step as a_n
```

The two forms are algebraically identical. KDK makes the **one force evaluation per step** obvious: `a_{n+1}` becomes the next step's `a_n`, so you evaluate `F` exactly once per step provided you keep `a` in your integrator state. A common beginner bug is evaluating `F` twice per step and concluding Verlet costs 2×; it does not.

- **Order:** local `O(dt³)`, global `O(dt²)` in **both** `x` and `v`.
- **Cost:** 1 force eval/step, one extra `Vec3` array (`a`) of storage.
- **Energy:** bounded, `O(dt²)` oscillation, no drift. **Time-reversible** and symplectic.

#### Relationship to leapfrog and position Verlet

Three algorithms, one trajectory:

| Name | Update | Where velocity lives |
|---|---|---|
| **Leapfrog (drift-kick-drift / KDK)** | `v_{n+1/2} = v_{n-1/2} + dt·a_n` ; `x_{n+1} = x_n + dt·v_{n+1/2}` | half-integer steps |
| **Velocity Verlet** | as above | integer steps (synchronised) |
| **Position (Störmer) Verlet** | `x_{n+1} = 2x_n − x_{n−1} + dt²·a_n` | not stored |

They generate the **same sequence of positions** (in exact arithmetic). Velocity Verlet is the preferred implementation because:
1. `x` and `v` are synchronised at the same time, so energy `E = ½mv² + U(x)` is meaningful at every step. With raw leapfrog, using `v_{n+1/2}` in the energy gives a spurious `O(dt)` offset — a classic source of "why is my energy wrong" confusion.
2. Position Verlet needs `x_{n-1}`, needs special startup, cannot restart from `(x, v)`, and recovering velocity as `(x_{n+1} − x_{n−1})/(2dt)` suffers catastrophic cancellation for small `dt`. Avoid it.

#### Accuracy character: bounded energy, drifting phase

The exactly conserved quadratic for the harmonic oscillator is

```
Q = (1 − k²/4)·ω²x² + v²
```

so relative energy error oscillates with amplitude `≈ k²/4` and never grows. But the *frequency* is wrong: the numerical frequency satisfies `cos(ω_num·dt) = 1 − k²/2`, giving

```
ω_num = ω · (1 + k²/24 + O(k⁴))
```

Phase error accumulates **linearly in time**: `Δφ ≈ 2πN·k²/24` after N periods. At `dt = T/100` that is ~`1.0e-3` rad per period, so after ~1000 oscillations your particle is a full radian out of phase while its energy is still perfect to 0.1%. **This is the defining trade of symplectic integrators, and you must state it in your docs**: they get long-term *statistics and structure* right, not long-term *phase*. If you need to know where the pendulum is after ten thousand swings, use RK4 with a small step or a high-order symplectic method, not Verlet.

Verlet is **exact** for constant acceleration (free fall, uniform E-field). Make that your first unit test — it should match analytic to machine precision at any `dt`.

Angular momentum for central forces is conserved exactly, by the same drift/kick argument as §2.2 (each drift preserves `x × v_half`, each kick preserves `x × v` when `a ∥ x`).

#### The velocity-dependence subtlety (drag, magnetic forces, Langevin)

Classic Verlet is derived assuming `F = F(x)` (or `F(t, x)`). If `F` also depends on `v`, the final half-kick becomes **implicit**:

```
v_{n+1} = v_half + 0.5*dt*F(t_{n+1}, x_{n+1}, v_{n+1})/m     // v_{n+1} on both sides
```

Four standard resolutions, in the order you should offer them:

**(1) Closed-form solve for linear drag — do this when you can.** If `F = F_c(t, x) − b·v` (`b` a scalar or matrix):

```
v_{n+1} = ( v_half + 0.5*dt*F_c(t_{n+1}, x_{n+1})/m )  /  ( 1 + 0.5*dt*b/m )
```

Exact, one force evaluation, keeps full second order, and is unconditionally stable in the drag term. If your library special-cases linear drag at all, do it this way.

**(2) Predictor–corrector (the common approximation).** Predict the velocity with an explicit Euler-quality estimate, use it to evaluate the force, then correct:

```
v_pred  = v_half + 0.5*dt*a_n                      // or v_n + dt*a_n
a_{n+1} = F(t_{n+1}, x_{n+1}, v_pred) / m
v_{n+1} = v_half + 0.5*dt*a_{n+1}
```

The predictor has local error `O(dt²)`; it is multiplied by `½dt` in the correction, so the local error contribution is `O(dt³)` — **second order is retained globally**. Cost stays at 1 force evaluation. Symplecticity is lost, but a velocity-dependent force is non-Hamiltonian anyway, so there was nothing to preserve. This is the pragmatic default for a student library.

**(3) Fixed-point iteration.** Repeat step (2)'s force evaluation with the corrected `v` 1–3 times, or until `|Δv| < tol`. Converges when `½·dt·|∂a/∂v| < 1`; in the limit it is the implicit trapezoidal rule in `v`. Costs one force eval per iteration. Use when drag is strongly nonlinear (quadratic drag at high speed) and you refuse to shrink `dt`.

**(4) Operator splitting.** Apply the drag as an exact exponential factor in its own sub-step: `v ← v·exp(−b·dt/(2m))` before and after the conservative KDK. This is the "BAOAB"-family trick, is unconditionally stable for arbitrarily large `b`, and is the right answer for thermostatted MD. Beyond a first library, but worth a comment in the header.

**(5) Just use RK4.** When forces are genuinely velocity-dependent and dissipative, the problem is not Hamiltonian, so Verlet's structure preservation buys you nothing and its 2nd order is simply worse than RK4's 4th. **Recommend RK4 for drag-dominated projectile problems and Verlet for conservative ones.**

**API implication:** make the distinction visible at compile time. Either (a) two force concepts — `ConservativeForce: (double t, Vec3 x) -> Vec3` and `GeneralForce: (double t, Vec3 x, Vec3 v) -> Vec3` — with `velocity_verlet` accepting only the former; or (b) a `static constexpr bool depends_on_velocity` trait on the force functor, `static_assert`ed inside `velocity_verlet` unless the user opts into the predictor–corrector variant `velocity_verlet_pc`. Silently applying plain Verlet to a `v`-dependent force is a correctness bug that produces plausible-looking wrong answers; make it un-writable.

---

### 2.4 Classic RK4 for the coupled `[x, v]` system

Generic RK4 for `y' = f(t, y)`; specialised to `y = (x, v)`, `f = (v, a)` with `a(t, x, v) = F(t, x, v)/m`.

**Stages, written out for the second-order case:**

```
k1x = v_n
k1v = a( t_n,           x_n,                v_n               )

k2x = v_n + 0.5*dt*k1v
k2v = a( t_n + 0.5*dt,  x_n + 0.5*dt*k1x,   v_n + 0.5*dt*k1v  )

k3x = v_n + 0.5*dt*k2v
k3v = a( t_n + 0.5*dt,  x_n + 0.5*dt*k2x,   v_n + 0.5*dt*k2v  )

k4x = v_n + dt*k3v
k4v = a( t_n + dt,      x_n + dt*k3x,       v_n + dt*k3v      )

x_{n+1} = x_n + (dt/6)*(k1x + 2*k2x + 2*k3x + k4x)
v_{n+1} = v_n + (dt/6)*(k1v + 2*k2v + 2*k3v + k4v)
```

Order of evaluation matters: **compute `k_ix` before `k_iv` at each stage** — note that `k2v` needs `k1x`, `k3v` needs `k2x`, `k4v` needs `k3x`. Getting this wrong (e.g. using `k1x` where `k2x` belongs) silently degrades the method to 2nd order; your convergence test in §4 will catch it, which is exactly why you write that test.

**Simplification worth implementing.** Because `k_ix` are simple linear combinations of `v_n` and the `k_iv`, the position update collapses:

```
x_{n+1} = x_n + dt*v_n + (dt*dt/6)*(k1v + k2v + k3v)
```

Verify: `(dt/6)(k1x + 2k2x + 2k3x + k4x) = (dt/6)(6v_n + dt·k1v + dt·k2v + dt·k3v)`. This saves arithmetic and makes the "`x` advances by `v·dt` plus an acceleration correction" structure obvious. The velocity update keeps all four `k_iv` and the `2,2` weights.

- **Order:** local `O(dt⁵)`, global `O(dt⁴)`.
- **Cost:** **4 force evaluations per step** (no FSAL in classic RK4). Memory: 4 stage vectors per component, or reuse two temporaries.
- **Energy:** *not* symplectic, *not* time-reversible. On an oscillator the amplification factor is `|R(ik)|² = 1 − k⁶/72 + k⁸/576`, so for `k < 2√2` it is slightly **contractive** — RK4 numerically damps oscillators. The damping is `O(dt⁵)` per unit time (`ΔE/E ≈ −ω⁶dt⁵·t/72`), i.e. astonishingly small for reasonable `dt`, but **monotone and secular**. Run a Kepler orbit for 10⁷ periods with RK4 and the planet spirals in; run it with Verlet and it doesn't.

**Cost-per-accuracy.** For target error `ε`: Verlet needs `dt_V ≈ √(ε/C₂)`, RK4 needs `dt_R ≈ (ε/C₄)^{1/4}`. Work is `1/dt_V` vs `4/dt_R` force evaluations per unit time. RK4 wins decisively for `ε ≲ 1e-4` and the margin widens as `ε` shrinks. For `ε ~ 1e-2` and long horizons, Verlet wins on structure.

**When to use RK4:** general-purpose accuracy over moderate horizons; any velocity-dependent force; trajectory problems where you want the *position at a specific time* to be right (projectile with drag, orbit transfer targeting); anything where you'd otherwise fight Verlet's velocity subtlety. **When not to:** millions of periods of conservative dynamics (use Verlet), or hard per-step cost budgets (4× the force evaluations).

---

## 3. Timestep selection

### 3.1 Identify the fastest timescale first

`dt` is never chosen in absolute terms; it is chosen relative to the shortest characteristic time in the problem.

| System | Fastest timescale | Rule of thumb |
|---|---|---|
| Spring/oscillator | `T = 2π√(m/k)` | `dt = T/100` for engineering accuracy, `T/20` minimum, `T/1000` for high precision |
| Pendulum, small angle | `T = 2π√(L/g)` | same |
| Linear drag | `τ = m/b` | `dt ≲ τ/10`; hard stability limit `dt < 2τ` for explicit |
| Gravity / orbit | `T_orbit = 2π√(a³/GM)` | `dt = T/200` at perihelion; for N-body use `dt = η·min_i √(ε/|a_i|)` or `η·min_i |v_i|/|a_i|`, `η ≈ 0.01–0.03` |
| Collision / contact | contact duration `t_c ≈ π√(m/k_contact)` | `dt < t_c/10`, else tunnelling |
| Generic | `ω_max = √(λ_max(−(1/m)∂F/∂x))` | `dt < 0.1 / ω_max` |

For coupled systems, `ω_max` is the largest eigenvalue of the stiffness matrix `−M⁻¹ ∂F/∂x`. You rarely compute it; estimate it from the stiffest single spring, `ω ≈ √(k_stiffest / m_lightest)`. Note that the **lightest** mass matters, not the average — one light particle attached to a stiff spring sets `dt` for the whole simulation.

### 3.2 Stability limits (hard ceilings)

For `x'' = −ω²x`, `k = ω·dt`:

- Forward Euler: **no stable `dt`**.
- Symplectic Euler, velocity Verlet: `k < 2`, i.e. `dt < T/π ≈ 0.318·T` (`≈ 3.18` steps per period).
- RK4: `k < 2√2 ≈ 2.8284`, i.e. `dt < 0.450·T`.
- Explicit treatment of linear drag `v' = −v/τ`: `dt < 2τ` (Euler-type), else `v` oscillates in sign and diverges.

**Stability is not accuracy.** At `k = 1.9`, velocity Verlet is stable but the numerical frequency is off by ~19% and the energy wobble is ~90%. The stability limit tells you where things explode; the accuracy rule of thumb (`k ≈ 0.06`, i.e. `T/100`) tells you where things are *right*. Design your API so users pick `dt` from the accuracy rule and the stability limit is only a safety net.

### 3.3 Detecting a bad `dt`

Build these in from day one; they cost almost nothing and save hours.

1. **Energy monitor.** Compute `E(t)` every K steps. Flag if `|E − E₀|/max(|E₀|, E_scale) > tol` (say 1e-2) for a conservative system, or if a least-squares line fit to `E(t)` has a significantly nonzero slope. A *bounded oscillation* is expected for symplectic methods; a *trend* is the alarm.
2. **Non-finite guard.** `if (!std::isfinite(...)) throw / return error`. Check once per K steps, not every step — cheap enough either way, but the branch in the inner loop costs. Blow-up is unmistakable: energy goes `1e3 → 1e30 → inf` within a handful of steps.
3. **Step-size sanity ratios.** Per step, `|Δx| / L_char` and `|Δv| / (|v| + v_eps)` should both be ≲ 0.05. Cheap, dimensionally sensible, catches near-collisions and stiff spikes.
4. **Halve-and-compare (the definitive test).** Run to the same final time with `dt` and `dt/2`. If the two trajectories differ by more than your tolerance, `dt` is too big — full stop. This is self-convergence (§4.2) used as a runtime diagnostic.
5. **Frequency check.** For oscillatory problems, count zero crossings and compare the measured period to the analytic one. A period error above a few percent means `dt` is too big.
6. **Symptom dictionary:**
   - Amplitude growing on an undamped oscillator → forward Euler, or `dt` past the stability limit.
   - Amplitude decaying on an undamped oscillator → RK4 with huge `dt`, backward Euler, or a bug that added damping.
   - Energy sawtooths with constant amplitude → normal for symplectic methods. Not a bug.
   - Energy jumps at one instant then stays flat → a discrete event (collision, force discontinuity), not an integrator problem.
   - Results change when you change output/sampling frequency → your sampling is modifying `dt`. Real bug (§5.6).

---

## 4. Validation methodology

Your test suite is the actual product here. Three independent families of checks; all three must pass.

### 4.1 Compare against analytic solutions

Have each of these as a fixture with an exact `x(t), v(t)`.

**(a) Constant force / free fall.** `a = −g ĵ`. `x(t) = x0 + v0 t − ½g t² ĵ`.
Verlet and RK4 are **exact** here (to roundoff), at any `dt`. Symplectic Euler has an `O(dt)` offset (`x_n = x0 + n·dt·v0 − ½g dt² n(n+1)`, vs exact `−½g dt²n²`), forward Euler the mirror-image offset. First smoke test; catches sign errors and `½dt²` typos instantly.

**(b) Simple harmonic motion.** `m x'' = −k x`, `ω = √(k/m)`.
```
x(t) = x0·cos(ωt) + (v0/ω)·sin(ωt)
v(t) = −x0·ω·sin(ωt) + v0·cos(ωt)
E    = ½m v² + ½k x² = const
A    = √(x0² + (v0/ω)²)
```
The workhorse test. Use it for convergence order, energy behaviour, and the quantitative predictions in the §2.0 table (assert forward Euler grows by `(1+k²)^n`; assert Verlet's `Q = (1−k²/4)ω²x² + v²` is constant to ~1e-13 — a beautifully sharp test that catches almost any Verlet bug).

**(c) Damped oscillator.** `m x'' + b x' + k x = 0`; `γ = b/(2m)`, `ω₀ = √(k/m)`.
```
underdamped (γ < ω₀):   ω_d = √(ω₀² − γ²)
   x(t) = e^{−γt} [ x0 cos(ω_d t) + ((v0 + γ x0)/ω_d) sin(ω_d t) ]
critically damped (γ = ω₀):
   x(t) = e^{−γt} [ x0 + (v0 + γ x0) t ]
overdamped (γ > ω₀):    r± = −γ ± √(γ² − ω₀²)
   x(t) = C₊ e^{r₊t} + C₋ e^{r₋t},  C± = ±(v0 − r∓ x0)/(r₊ − r₋)
```
Exercises velocity-dependent forces. Test all three regimes; the critically damped case has a repeated root and catches formula bugs.

**(d) Driven damped oscillator (tests time dependence).** `F = F₀cos(Ω t)`; steady state
```
x_ss(t) = (F₀/m)/√((ω₀²−Ω²)² + 4γ²Ω²) · cos(Ω t − δ),   tan δ = 2γΩ/(ω₀²−Ω²)
```
Run past the transient (`t ≫ 1/γ`) and compare amplitude and phase. This is the test that catches "I used `t_n` where I should have used `t_n + dt/2`" in RK4 stages — a very common bug that is invisible for autonomous systems.

**(e) Projectile with linear drag + terminal velocity.** `m v' = −mg ĵ − b v`; `τ = m/b`, `v_t = gτ = mg/b`.
```
v_x(t) = v_x0 · e^{−t/τ}
x(t)   = x0 + v_x0·τ·(1 − e^{−t/τ})          → horizontal range asymptote x0 + v_x0·τ
v_y(t) = (v_y0 + v_t)·e^{−t/τ} − v_t          → −v_t as t → ∞
y(t)   = y0 + (v_y0 + v_t)·τ·(1 − e^{−t/τ}) − v_t·t
```
Terminal-velocity assertion: after `t = 10τ`, `|v_y + v_t| < 1e-3·v_t`. Also assert the horizontal-range asymptote.

**(f) Quadratic drag, vertical fall from rest (nonlinear).** `m dv/dt = −mg + c v²`, `v_t = √(mg/c)`:
```
speed(t)    = v_t · tanh(g t / v_t)
distance(t) = (v_t²/g) · ln cosh(g t / v_t)
```
The only closed-form nonlinear-drag case; use it to validate your predictor–corrector Verlet against RK4.

**(g) Rocket equation (variable mass).** Burn rate `ṁ > 0`, exhaust speed `u`, `m(t) = m₀ − ṁ t` for `t < t_burnout`:
```
no gravity:        v(t) = v0 + u·ln( m₀ / m(t) )
constant gravity:  v(t) = v0 + u·ln( m₀ / m(t) ) − g·t
```
Integrate `dv/dt = u·ṁ/m(t) − g`. **Do not** integrate `d(mv)/dt = F_ext` — see §7.8. This test exists specifically to force your API to handle time-varying mass correctly.

**(h) Kepler two-body (the long-run structure test).** `r'' = −μ r/|r|³`, `μ = GM`:
```
ε = |v|²/2 − μ/|r|        (specific energy, conserved)
a = −μ/(2ε)               (semi-major axis)
T = 2π√(a³/μ)             (period)
L = r × v                 (conserved, exactly, for Verlet/symplectic Euler)
e_vec = (v × L)/μ − r̂     (Laplace–Runge–Lenz; conserved ⇒ no apsidal precession)
```
Run an `e = 0.5` ellipse for 10³–10⁵ periods. Assert: Verlet's `a` shows no trend; RK4's `a` shrinks monotonically at the predicted rate; `|L|` is conserved to ~1e-14 by Verlet. Measure apsidal precession per orbit as a phase-error metric. This single test discriminates every integrator property that matters.

### 4.2 Convergence testing (the most important test you will write)

**Protocol.** Fix the final time `T_end` and the initial condition. For `dt_j = dt₀ / 2^j`, `j = 0..5`, compute

```
e_j = || y_numeric(T_end; dt_j) − y_exact(T_end) ||
```

Then the observed order is

```
p_j = log2( e_j / e_{j+1} )
```

Assert `p_j ∈ [p_expected − 0.15, p_expected + 0.15]` for the middle of the range. A cleaner statistic: least-squares fit `log e` vs `log dt` over the usable range and check the slope.

**Rules that make this test reliable:**

- **Fix `T_end`, not the step count.** `n_steps = llround(T_end / dt)` — and choose `dt₀` so `T_end/dt₀` is an integer, so halving stays exact.
- **Use an error norm over the whole trajectory,** not a single endpoint. `e = max_t |x_num(t) − x_exact(t)|` over the recorded samples, or the discrete L2 norm. A single endpoint can accidentally sit near a zero of the error and produce a nonsense order estimate.
- **Test `x` and `v` separately.** Some buggy schemes are 2nd order in `x` and 1st in `v`.
- **Watch the roundoff floor.** Error cannot go below `≈ ε_mach · |x| · √(n_steps)` (random-walk accumulation; worst case `n_steps·ε_mach·|x|`). With RK4, `dt = T/1000` already gives `e ~ 1e-12` and the next halving hits the floor, producing a bogus `p ≈ 1`. **Choose your `dt` range so all points are 3+ decades above the floor** — for RK4 that means fairly large `dt` (`T/10` to `T/160`), for Verlet a wider window. Plot the L-shaped curve once so you recognise it.
- **Self-convergence when no analytic solution exists.**
  ```
  p ≈ log2( ||y(dt) − y(dt/2)|| / ||y(dt/2) − y(dt/4)|| )
  ```
  Works for any problem. Use it to validate your integrators on the actual systems your users will run.
- **Richardson extrapolation** gives you a reference solution of order `p+1` for free: `y_ref ≈ (2^p·y(dt/2) − y(dt)) / (2^p − 1)`.

### 4.3 Conservation checks

| Quantity | Expression | Expected |
|---|---|---|
| Energy | `E = Σ ½m_i|v_i|² + U(x)` | symplectic: bounded oscillation, zero trend, amplitude `O(dt^p)`. RK4: monotone drift `O(dt^4)` per unit time in position error, `O(dt^5·t)` in energy. Forward Euler: exponential growth |
| Linear momentum | `P = Σ m_i v_i` | **machine precision** for internal pairwise forces, all methods, if you apply `+f/−f` |
| Angular momentum | `L = Σ m_i x_i × v_i` | **machine precision** for central forces with symplectic Euler / Verlet; `O(dt⁴)` drift with RK4 |
| Centre of mass | `X_cm` linear in `t` | machine precision for isolated systems |

**Reporting.** Use `max_t |E(t) − E(0)| / E_scale` for the bound, and the **slope of a linear fit** to `E(t)` for the drift. Reporting only the max conflates the two, and the distinction (bounded vs drifting) is exactly what you're testing.

**Momentum conservation is the single best structural unit test for an N-body force routine** — it fails loudly the moment you forget a `-=` or double-count a pair.

Also test: **time-reversibility.** Integrate forward `N` steps, negate `v`, integrate `N` more steps. Verlet and symplectic-Euler-adjoint pairs should return to the start to within roundoff (`~1e-12`); RK4 will not (`O(dt⁴)`); forward Euler will be wildly off. Cheap, needs no analytic solution, and catches asymmetries in your kick/drift ordering.

---

## 5. Practical API considerations

### 5.1 Fixed-step first — and why

**Recommend fixed-step for v1, unambiguously.** Reasons:

1. **Symplectic integrators require constant `dt`.** Adaptive stepping destroys the shadow-Hamiltonian argument and reintroduces energy drift; you would be building a feature that breaks your two best integrators.
2. **Reproducibility.** Same input → same output, bit for bit. Essential for the test suite and for a student debugging physics rather than infrastructure.
3. **Convergence testing is trivial** with fixed steps and awkward with adaptive ones (you have to compare tolerances, not step sizes).
4. **Uniform output** makes plotting, FFTs, and animation straightforward.
5. Implementation is ~30 lines instead of ~300 (embedded error estimator, PI controller, step rejection, min/max clamping, dense output).

Add adaptive stepping in v2, as a **separate** code path for RK-family methods only: Dormand–Prince 5(4) with FSAL, a PI controller (`dt_new = dt·min(f_max, max(f_min, S·(tol/err)^{1/5}))`, `S ≈ 0.9`), step rejection, and free dense output. Document that adaptive + symplectic is a contradiction.

### 5.2 The force callable

```cpp
// Concept-level contract, single particle:
//   Vec3 operator()(double t, const Vec3& x, const Vec3& v) const;
//
// N particles, fill-in style (no allocation in the inner loop):
//   void operator()(double t,
//                   std::span<const Vec3> x,
//                   std::span<const Vec3> v,
//                   std::span<Vec3>       f_out) const;   // MUST zero f_out first
```

Design points:

- **Template on the callable, don't use `std::function` in the hot loop.** `std::function` is an indirect call the compiler cannot inline; for a 20-flop force this can cost 5–10×. Template the integrator; offer a `std::function`-based wrapper at the outermost API boundary only if you need runtime polymorphism.
  ```cpp
  template <class Force>
  void velocity_verlet_step(System& s, const Force& F, double t, double dt);
  ```
- **Fill-in, not return-by-value,** for N particles. Returning `std::vector<Vec3>` allocates once per force evaluation — that is four allocations per RK4 step.
- **Zero the accumulator inside the force functor,** not at the call site, so a user-written force can never forget. Document it as part of the contract.
- **Advertise velocity-independence.** A `static constexpr bool depends_on_velocity = false;` member (or two distinct concepts) lets `velocity_verlet` `static_assert` and lets you dispatch to the exact linear-drag solve when a `LinearDrag` tag is present.
- **Composition.** Provide a `SumForce<F1, F2, ...>` variadic combinator so users write `Gravity{} + LinearDrag{b} + Spring{k, anchor}` instead of one monolithic lambda. Its `depends_on_velocity` is the OR of its parts — computed at compile time.

### 5.3 Single particle vs systems

Write the N-particle version and let `N = 1` be a special case. Two thin conveniences on top:

```cpp
// Convenience overload — internally constructs an N=1 view, no duplicated integrator logic.
Trajectory simulate(Vec3 x0, Vec3 v0, double mass, ForceFn F, double t0, double t1, double dt);
```

Do **not** write two integrators. Every bug you fix you will fix once.

For N particles, `std::vector<Vec3>` (array-of-structs) is the right first choice: readable, contiguous, fine for N up to ~10⁵. Structure-of-arrays (`vector<double> x, y, z`) only pays off once you are vectorising, and it obscures the physics. Reserve capacity up front.

### 5.4 Integrator state

Velocity Verlet must persist `a` between steps to keep 1 force eval/step. That means integrators are **objects**, not free functions:

```cpp
class VelocityVerlet {
    std::vector<Vec3> a_;   // cached a_n
    bool primed_ = false;
public:
    template <class Force>
    void step(System& s, const Force& F, double t, double dt);
    void invalidate();      // call if the user mutates state or the force out-of-band
};
```

`invalidate()` matters: if a user teleports a particle or changes a spring constant between steps, the cached `a_` is stale and the trajectory is silently wrong. Make the cache explicitly invalidatable, and document it.

### 5.5 Time stepping loop

```cpp
const long long n_steps = std::llround((t_end - t0) / dt);
for (long long i = 0; i < n_steps; ++i) {
    const double t = t0 + static_cast<double>(i) * dt;   // NOT t += dt   (see §7.2)
    integrator.step(system, force, t, dt);
    observer(t + dt, system);
}
```

Integer step count, multiplication for `t`, observer after the step. Never `while (t < t_end) { ...; t += dt; }` — the number of steps then depends on accumulated roundoff and can differ by one between builds.

### 5.6 Recording / sampling trajectory points

**Physics `dt` and output `dt` must be decoupled.** Storing every step of a 10⁷-step run is 240 MB of `Vec3` and usually useless; changing `dt` to hit output times corrupts symplecticity.

Three mechanisms, offer all three:

1. **Stride.** `every_k = max(1, llround(sample_dt / dt))`; record when `i % every_k == 0`. Simplest; output times are `k·every_k·dt`, which may not be the round numbers the user asked for.
2. **Cubic Hermite dense output** — record at *exact* requested times. You have `x, v` at both ends of the step, which is exactly the data cubic Hermite needs:
   ```
   s   = (t − t_n)/dt,  in [0,1]
   h00 =  2s³ − 3s² + 1     h10 =  s³ − 2s² + s
   h01 = −2s³ + 3s²         h11 =  s³ − s²
   x(t) = h00·x_n + h10·dt·v_n + h01·x_{n+1} + h11·dt·v_{n+1}
   v(t) = ( (6s²−6s)·x_n + (3s²−4s+1)·dt·v_n
          + (−6s²+6s)·x_{n+1} + (3s²−2s)·dt·v_{n+1} ) / dt
   ```
   Interpolation error is `O(dt⁴)` — strictly better than Verlet's `O(dt²)` global error, so it costs you nothing in accuracy. **Recommended default.**
3. **Event detection.** Root-find `g(t) = 0` (e.g. `g = y`, ground impact; or `g = v·x`, apsis passage) on the Hermite interpolant using bisection or Brent, then record/act at the exact crossing. Detect a sign change of `g` across a step, then bracket within it. This is how you get an exact projectile range instead of "somewhere within one step of the ground."

**Observer pattern for the API:**

```cpp
template <class Force, class Observer>
void integrate(System& s, const Force& F, double t0, double t1, double dt, Observer&& obs);
```

with library-provided observers: `NullObserver`, `StrideRecorder`, `DenseRecorder{times}`, `EnergyMonitor`, `TeeObserver<A,B>`. Template it — a `NullObserver` should compile to nothing.

### 5.7 Miscellaneous

- Return a `Trajectory { std::vector<double> t; std::vector<std::vector<Vec3>> x, v; }` or a column-major flat buffer; either way `reserve()` from the known sample count.
- Avoid virtual dispatch inside the step loop. If you want runtime integrator selection, put the virtual at the `run_whole_trajectory` level, not `step`.
- Sum forces in a **fixed order** (e.g. always ascending pair index). Floating-point addition is not associative; a parallel or reordered sum breaks bit-reproducibility and makes tight-tolerance tests flaky.
- Provide `dt` as an explicit required argument. Never default it — a silently defaulted `dt` is a silently wrong simulation.

---

## 6. Units

**Policy: SI everywhere, `double` everywhere, no exceptions in v1.**

Rationale:

- **Textbook constants work directly.** `G = 6.67430e-11`, `g = 9.80665`, `c = 2.99792458e8`, `k_B = 1.380649e-23`. Students copy formulas from their mechanics text and they just work.
- **No conversion layer means no conversion bugs.** Unit mixups (feet/metres, degrees/radians) are the single most common source of "the physics is wrong" in student codes, and the cure is to have exactly one unit system.
- **Doubles have the range.** IEEE-754 double spans `~1e±308` with 15–16 significant digits. Even solar-system dynamics in SI (`M_sun = 1.989e30 kg`, `AU = 1.496e11 m`) uses only a fraction of that. Dynamic range is *not* a reason to abandon SI.
- **Angles in radians, always.** State it in every doc comment. Provide `constexpr double deg = M_PI/180.0;` so users write `45*deg`.

**Concretely:**

```cpp
namespace units {  // all SI
    inline constexpr double G  = 6.67430e-11;   // m^3 kg^-1 s^-2
    inline constexpr double g0 = 9.80665;       // m s^-2
    inline constexpr double deg = 3.14159265358979323846 / 180.0;
    // convenience multipliers for INPUT ONLY:
    inline constexpr double km = 1e3, AU = 1.495978707e11, day = 86400.0, yr = 3.15576e7;
}
```
Convenience multipliers (`5*units::km`) let users write readable inputs without introducing a second internal unit system. Document every function parameter with its unit in the comment: `/// @param b linear drag coefficient [kg/s]`.

**Where natural / reduced units earn their keep** — mention in docs, don't build in v1:

- **Gravitational N-body:** `G = 1`, lengths in AU, times in years, masses in `M_sun` (then `G = 4π²`). Keeps numbers `O(1)`, improves conditioning, makes tolerances interpretable.
- **Molecular dynamics:** Lennard-Jones reduced units (`σ = ε = m = 1`).
- **Nondimensionalisation as a *validation* tool:** rescaling an SI problem to `O(1)` numbers and re-running is a good way to check that your tolerances are not accidentally absolute-scale-dependent.

If you want to support these later, do it by making constants **data**, not literals:

```cpp
struct UnitSystem { double G, g0; /* ... */ };
inline constexpr UnitSystem SI       {6.67430e-11, 9.80665};
inline constexpr UnitSystem Astro    {4.0*M_PI*M_PI, 0.0};   // AU, yr, M_sun
```

Pass a `UnitSystem` to force constructors. This is a 20-line change and avoids ever hard-coding `6.674e-11` inside `Gravity`.

**Compile-time dimensional analysis** (strong typedefs, or `mp-units` in C++23) is genuinely valuable and would catch a whole class of bugs — but it triples the API's syntactic weight and the error messages are brutal. **Not for v1.** Note it in the roadmap.

---

## 7. Common pitfalls

### 7.1 Stiff systems

**The problem.** Explicit methods have `dt < C/ω_max`. A "rigid" rod modelled as `k = 1e6 N/m` on `m = 1e-3 kg` gives `ω = 3.16e4 rad/s`, so velocity Verlet needs `dt < 6.3e-5 s` — 16,000 steps per second of simulated time — even though the interesting motion happens on a 1-second timescale. **Stiffness = a wide separation between the fastest timescale you must resolve for stability and the slowest one you care about.**

**Symptoms:** blow-up at any `dt` you consider reasonable; wildly oscillating velocities in one component; "it works at `dt=1e-5` but not `1e-4`" with a sharp cliff between them (that cliff is the stability boundary, and it's exactly at `dt = 2/ω`).

**Remedies, in order of preference for a student library:**

1. **Question the model.** A very stiff spring is usually a badly-modelled rigid constraint. Replace it with an actual constraint (distance constraint solved by RATTLE, or a reduced-coordinate formulation). This removes the stiffness rather than fighting it.
2. **Soften the stiffness or increase the mass** if physical fidelity permits.
3. **Implicit methods.**
   - *Backward Euler:* `v_{n+1} = v_n + dt·a(t_{n+1}, x_{n+1}, v_{n+1})`, `x_{n+1} = x_n + dt·v_{n+1}`. Unconditionally stable, 1st order, and **numerically dissipative** — it damps everything, including motion you wanted. Great for games (nothing ever explodes), bad for physics (your orbit spirals in). Requires solving a nonlinear system per step (Newton + Jacobian, or fixed-point).
   - *Implicit midpoint:* `y_{n+1} = y_n + dt·f(t_n + dt/2, (y_n + y_{n+1})/2)`. 2nd order, A-stable, **symplectic**, time-reversible. The right answer when you need both stiffness tolerance and conservation. Still needs a nonlinear solve.
4. **Multiple time stepping** (r-RESPA): integrate the stiff force with `m` small sub-steps inside each large step of the soft forces. Cheap when the stiff force is also cheap to evaluate.
5. **Exponential / split integrators** for stiff *linear* parts (drag, harmonic bonds): solve that part exactly (§2.3 option 4).

For v1: **detect and report**, don't solve. Estimate `ω_max` from the user's force (or let them declare it), and warn if `ω_max·dt > 1`. A clear diagnostic message beats a silent divergence.

### 7.2 Accumulating `t` by repeated addition

```cpp
double t = t0;
for (long long i = 0; i < n; ++i) { step(...); t += dt; }   // WRONG
```

`dt` is generally not exactly representable (`0.001` is not a binary fraction), so each `t += dt` commits a rounding error. Over `n` steps these accumulate — as a random walk `~ ε·|t|·√n` typically, but up to `~ ε·|t|·n` in the worst case. For `n = 10⁷` and `t = 10⁴`, the drift can reach `1e-5 s`. That is enough to visibly dephase a driven oscillator and to corrupt a comparison against an analytic `x(t)`.

```cpp
const double t = t0 + static_cast<double>(i) * dt;          // RIGHT: one rounding, always
```

One multiplication, one rounding error, error bounded by `ε·|t|` for all `i`. Also use an **integer step count** `n = llround((t1-t0)/dt)` rather than a `while (t < t1)` condition, so the number of steps is deterministic across compilers and optimisation levels.

If you need both exactness and incremental updating, use Kahan/Neumaier compensated summation — but `t = t0 + i*dt` is simpler and strictly better here.

### 7.3 `float` vs `double`

Use **`double` for all state and all arithmetic.** `float` has ~7 decimal digits (`ε ≈ 1.2e-7`); after 10⁵ steps, accumulated roundoff in position is comparable to a `dt = T/100` truncation error, so you gain nothing from a smaller `dt` — you have hit the roundoff floor and your convergence test flatlines. Worse, energy conservation tests become meaningless because roundoff drift swamps the integrator's own behaviour.

`float` is acceptable only for (a) rendering/visualisation buffers, (b) on-disk trajectory storage where you accept the loss, (c) GPU work where you have measured the trade. Never for the integrator state.

`long double` (80-bit on x86) is rarely worth it — 2–4× slower, non-portable width, and if you need it your algorithm is probably the problem.

### 7.4 Comparing floating point in tests

**Never `==`** on computed doubles. Use a combined absolute/relative tolerance:

```cpp
inline bool close(double a, double b, double rtol = 1e-9, double atol = 1e-12) {
    if (!std::isfinite(a) || !std::isfinite(b)) return false;   // NaN != NaN, be explicit
    return std::abs(a - b) <= atol + rtol * std::max(std::abs(a), std::abs(b));
}
```

- `atol` handles quantities that legitimately pass through zero (an oscillator's `x` at the turning point of `v`). Choose `atol` from the **problem scale**: `atol = 1e-12 * amplitude`, not a magic constant.
- `rtol` handles everything else.
- **Convergence-order assertions need loose tolerances.** Assert `|p_observed − p_expected| < 0.15`, not `p == 4.0`. The observed order approaches the true order only asymptotically.
- **NaN comparisons are always false**, so `assert(x == x)` is a NaN check but `assert(!(a > b))` does not mean `a <= b`. Check `isfinite` explicitly.
- Catch2's `Approx` / `WithinRel` / `WithinAbs` and doctest's `doctest::Approx` implement this; use them rather than rolling your own, but understand the default tolerances (Catch2's `Approx` default `epsilon` is `~1e-8` *relative*, which is too loose for machine-precision assertions and too tight for `dt`-truncation assertions — always set it explicitly).
- **Do not assert bit-exact results across compilers.** FMA contraction, `-ffast-math`, and different libm implementations of `exp`/`sin` change last bits. Compare against analytic solutions with physics-scaled tolerances instead.

### 7.5 Verlet-specific traps

- **Computing energy from the half-step velocity.** `v_half` is `O(dt)` away from `v_n`; using it gives a spurious energy offset that looks like a bug in your integrator. Always use the synchronised `v_{n+1}`.
- **Evaluating the force twice per step** (once for `a_n`, once for `a_{n+1}`) instead of caching. Doubles your cost for nothing.
- **Forgetting `invalidate()`** after external state mutation (§5.4).
- **Applying plain Verlet to a velocity-dependent force** without the §2.3 handling. Silently reduces to 1st order or worse.
- **Dropping the `½a·dt²` term** in the position update. Reduces Verlet to symplectic Euler; the convergence test catches it (order 2 → 1).

### 7.6 Force-routine traps

- **Zero-vector normalisation.** `F_drag = −c·|v|²·(v/|v|)` divides by zero when `v = 0`. Write `F_drag = −c·|v|·v` instead — algebraically identical, singularity-free. When you truly need a unit vector, guard: `if (n2 < eps2) return Vec3::zero();`.
- **Gravitational singularity.** `F = −Gm₁m₂ r/|r|³` blows up as `r → 0` and will produce `inf` on a close encounter. Use Plummer softening `|r|³ → (|r|² + ε²)^{3/2}` with `ε` chosen as a fraction of the typical separation, and document that this changes the physics at short range.
- **Forgetting to zero the accumulator** before summing pairwise forces. Produces forces that grow linearly with step number — a spectacular and confusing failure.
- **Breaking Newton's third law** in the pair loop (`f[i] += f_ij` but forgetting `f[j] -= f_ij`, or double-counting by looping over all `i,j` instead of `j > i`). Momentum conservation test catches this immediately (§4.3).

### 7.7 Build and numerics hygiene

- **`-ffast-math` is not free.** It permits reassociation (breaking reproducibility), assumes no NaN/Inf (breaking your `isfinite` guards, which may be optimised out entirely), and flushes denormals. Do not enable it in a library whose correctness tests depend on IEEE semantics. `-O2 -fno-fast-math` is the right default; if you need speed, profile first.
- **Denormals** near zero cause 10–100× slowdowns on some hardware. If you see mysterious slowdowns as a system relaxes to rest, that is why; flush-to-zero (`_MM_SET_FLUSH_ZERO_MODE`) is the fix, applied deliberately.
- **Catastrophic cancellation** in `x_{n+1} − x_n` for small `dt`. Prefer accumulating displacements, and never derive velocity by differencing positions when you already carry `v`.

### 7.8 Variable mass — the Newton's-second-law trap

For a rocket, `F = ma` still holds *for the instantaneous mass*, but `F_ext = d(mv)/dt` does **not** — the expelled mass carries momentum away. The correct equation of motion is

```
m(t)·dv/dt = F_ext + u·ṁ_burn        where ṁ_burn = −dm/dt > 0, u = exhaust speed
```

i.e. thrust `= u·ṁ_burn` is a *force*, and you divide by the current `m(t)`. If you integrate `d(mv)/dt = F_ext` you will get a rocket that accelerates when nothing pushes it. Encode this by making `m` a function of `t` in your `System` (or by offering an acceleration-callback overload), and validate against the closed-form rocket equation in §4.1(g).

---

## Appendix A: Recommended module layout

```
include/newton/
  vec3.hpp            // Vec3: +, -, *, /, dot, cross, norm, norm2, normalized (guarded)
  units.hpp           // SI constants, input-convenience multipliers, UnitSystem
  system.hpp          // System { x, v, a, invMass, mass }, N=1 convenience view
  forces/
    concepts.hpp      // ForceLike, depends_on_velocity trait, SumForce combinator
    uniform_gravity.hpp   linear_drag.hpp   quadratic_drag.hpp
    spring.hpp            newtonian_gravity.hpp (with softening)
    driven.hpp        // F0*cos(Omega*t) -- exercises time dependence
  integrators/
    forward_euler.hpp        // ships, documented as pedagogical only
    symplectic_euler.hpp
    velocity_verlet.hpp      // static_assert !depends_on_velocity
    velocity_verlet_pc.hpp   // predictor-corrector for general F
    rk4.hpp
  integrate.hpp       // fixed-step driver + observer dispatch
  observers.hpp       // Null, Stride, DenseHermite, EnergyMonitor, Tee
  diagnostics.hpp     // energy, momentum, angular momentum, finiteness guard
  interpolate.hpp     // cubic Hermite + event root-finding

tests/
  analytic/           // §4.1 (a)-(h)
  convergence/        // §4.2 order tables per (integrator, problem)
  conservation/       // §4.3 energy / P / L / time-reversibility
  regression/         // pinned trajectories with explicit tolerances
```

## Appendix B: Default recommendations, one line each

- **Default integrator:** velocity Verlet for conservative forces; RK4 the moment velocity dependence appears.
- **Default `dt`:** `T_fastest / 100`. Provide a helper that estimates `T_fastest` from the force where possible.
- **Default precision:** `double`, SI.
- **Default stepping:** fixed. Adaptive is v2 and RK-only.
- **Default output:** cubic Hermite dense sampling at user-specified times, decoupled from `dt`.
- **Always on:** energy monitor with a drift-slope check and a finiteness guard, once every ~1000 steps.
- **First three tests to write:** constant-acceleration exactness (Verlet, RK4); SHM convergence order table; N-body momentum conservation to machine precision.