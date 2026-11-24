# Physics to Code: A Chapter Map

This note is the bridge between the equation you just derived on paper in
Mechanics 1 and the function in this repo that implements it. For every
section of Taylor's *Classical Mechanics* chapters 1–5, you get: the key
equations in plain text, the real file and function that implements them,
the example that runs them, and the closed form this library uses to check
its own numerics are right.

Nothing here is a substitute for `docs/research/taylor-ch1-7-map.md`, which
is the full derivation-level research report this note condenses — go there
for the "why", come here for the "where is it in the code".

Conventions used below: `path/to/file.hpp` is a real path under
`include/newton/` or `examples/`; `function_name(...)` and `Struct{...}` are
real C++ identifiers you can grep for. Every equation uses `^` for powers,
`sqrt()`, and `.` for a time derivative is never implied — velocities are
just named `v`.

---

## Chapter 1 — Newton's Laws of Motion (Taylor §1.1–1.7)

### Key equations

```
r = (x, y, z)                     position
v = dr/dt,  a = dv/dt              velocity, acceleration
p = m v                            momentum
F = dp/dt = m a                    Newton's second law (constant mass)
```

Every integrator in this library works by recasting that single
second-order equation as a coupled first-order pair:

```
dr/dt = v
dv/dt = F(t, r, v) / m
```

That pair *is* the `State{r, v, t}` struct — see below.

### Where it lives

- `include/newton/core/vec3.hpp` — `Vec3`, the R^3 vector type: `+ - * /`,
  `dot`, `cross`, `norm`, `normalized`. Every position, velocity,
  acceleration, and force in the library is one of these.
- `include/newton/core/particle.hpp` — `State{r, v, t}` is exactly the
  `(r, v)` pair above; `Particle{mass, r, v}` is a particle's initial
  conditions; `momentum(particle)` returns `p = m v`.
- `include/newton/core/constants.hpp` — `constants::g_earth` and friends.
- The Force concept itself, `Vec3 operator()(double t, const Vec3& r, const
  Vec3& v) const`, is Taylor's `F(r, ṙ, t)` written as a callable — see
  `notes/03-numerical-methods.md` for how the four integrators turn a Force
  plus a `State` into the next `State`.

### Examples

- `examples/01_free_fall.cpp` — a straight drop and a 45° launch, both in
  vacuum (gravity only, no drag).

### Analytic validation

`include/newton/analytic/projectile.hpp` is Taylor §1.6's vacuum-projectile
closed form:

```
x(t) = v0 cos(theta) t
y(t) = v0 sin(theta) t - (1/2) g t^2
range(v0, theta, g)       = v0^2 sin(2 theta) / g
apex_height(v0, theta, g) = v0^2 sin^2(theta) / (2g)
time_of_flight(v0, theta, g) = 2 v0 sin(theta) / g
```

This is the library's sharpest first sanity check: a degree-2 polynomial in
`t` is something `velocity_verlet_step` and `rk4_step` must reproduce to
machine roundoff at *any* `dt` (constant acceleration is exact for both —
see `notes/03-numerical-methods.md`), while `euler_step` will show a small
but nonzero error. `examples/01_free_fall.cpp` prints exactly that
numeric-vs-analytic comparison, and also demonstrates
`simulate_until(events::ground_impact())` — the event system that stops a
run precisely when a falling particle reaches `y = 0`, rather than
somewhere inside the last time step.

---

## Chapter 2 — Projectiles and Charged Particles (Taylor §2.1–2.7)

This is the chapter that justifies the whole library. Taylor's own
observation — that quadratic drag in 2D has no closed-form solution — is
this project's founding argument for numerical integration in the first
place.

### Key equations

```
Drag magnitude:  f(v) = b v + c v^2       b = beta*D, c = gamma*D^2 (sphere, Taylor's STP coefficients)

Linear drag (vector):     F = -b v
Quadratic drag (vector):  F = -c |v| v            (NOT -c|v|^2 v_hat — see below)

Linear drag, 1D closed forms (tau = m/b, v_ter = m g/b):
   vx(t) = vx0 e^(-t/tau)
   x(t)  = vx0 tau (1 - e^(-t/tau))                    -> x_inf = vx0 tau
   vy(t) = -v_ter + (vy0 + v_ter) e^(-t/tau)
   y(x)  = ((vy0+v_ter)/vx0) x + v_ter tau ln(1 - x/(vx0 tau))   [exact 2D trajectory]

Quadratic drag, decoupled 1D cases only (v_ter = sqrt(m g / c)):
   horizontal:  v(t) = v0 / (1 + c v0 t / m)
   vertical (from rest):  v(t) = v_ter tanh(g t / v_ter)
                           y(t) = (v_ter^2/g) ln cosh(g t / v_ter)

Magnetic force:  F = q v x B         (does no work: |v| exactly constant)
cyclotron omega = q B / m
cyclotron radius = m v_perp / (|q| B)
cyclotron period = 2 pi m / (|q| B)
Lorentz force:   F = q (E + v x B)
```

### Where it lives

- `include/newton/forces/drag.hpp` — `LinearDrag{b}` (`F = -b v`) and
  `QuadraticDrag{c}` (`F = -c |v| v`), plus `linear_drag_sphere_air(D)` and
  `quadratic_drag_sphere_air(D)`, which use `constants::sphere_air_beta`
  and `constants::sphere_air_gamma` from
  `include/newton/core/constants.hpp` to turn a sphere's diameter into `b`
  or `c`. Note the algebraic form: `QuadraticDrag` is written `-c|v|v`, not
  the equivalent `-c|v|^2 v̂` — the latter needs `normalized(v)`, which has
  to special-case `v = 0`; the form actually used evaluates to exactly zero
  at `v = 0` with no division at all (see `core/vec3.hpp`'s `normalized()`
  doc comment for the same zero-guard idea applied to `Spring`).
- `include/newton/forces/gravity.hpp` — `UniformGravity{mass, g}` (`F = m
  g`) for the projectile baseline.
- `include/newton/forces/electromagnetic.hpp` — `MagneticForce{q, B}` (`F
  = q v × B`) and `LorentzForce{q, E, B}` (`F = q(E + v × B)`).
- `include/newton/forces/compose.hpp` — `forces::sum(f1, f2, ...)`, how a
  drag-and-gravity projectile becomes `sum(UniformGravity{m},
  QuadraticDrag{c})` instead of a bespoke struct.

### Examples

- `examples/02_projectile_drag.cpp` — Taylor's flagship baseball (`D = 7
  cm`, `m = 0.145 kg`, `30 m/s` at 45°): vacuum vs. linear vs. quadratic
  drag on one CSV, plus a launch-angle sweep for the quadratic case showing
  the optimal angle drops below 45°.
- `examples/03_charged_particle.cpp` — a helix in a uniform `B` field,
  verifying `|v|` stays constant and the radius/period match
  `analytic::magnetic`, plus a crossed-`E×B`-drift bonus.

### Analytic validation

`include/newton/analytic/drag.hpp` implements every closed form above:
`LinearDragParams{m, b, g}` with `.tau()` and `.v_ter()`;
`linear_drag_velocity`, `linear_drag_position`,
`linear_drag_trajectory_y` (the exact 2D `y(x)` — "the best 2D integrator
oracle in the book", per the research map); and the four quadratic-drag
scalar functions `quadratic_drag_terminal_speed`,
`quadratic_drag_fall_speed`, `quadratic_drag_fall_distance`,
`quadratic_drag_horizontal_speed`.

`include/newton/analytic/magnetic.hpp` implements `cyclotron_omega`,
`cyclotron_radius`, `cyclotron_period` — the standard oracle for checking
that an integrator preserves `|v|` on a force that does no work (see
`forces/electromagnetic.hpp`'s doc comment: this is one of the few forces
in the library where "energy conservation" is an exact invariant of the
true dynamics, not just a property of a good integrator).

There is deliberately **no** closed form for full 2D quadratic drag — that
gap *is* the chapter's point. Where the library can't validate against an
oracle (2D quadratic drag), it validates by convergence order instead (see
`notes/03-numerical-methods.md`).

---

## Chapter 3 — Momentum and Angular Momentum (Taylor §3.1–3.5)

### Key equations

```
P = sum(m_i v_i),  dP/dt = F_ext          (internal forces cancel, Newton III)

Rocket (Taylor's form):
   m dv/dt = -mdot v_ex + F_ext,  mdot = dm/dt < 0
   free space:            v - v0 = v_ex ln(m0/m)                    (Tsiolkovsky)
   uniform gravity g:     v(t) = v_ex ln(m0/m(t)) - g t
   m(t) = m0 - burn_rate*t while burning

Center of mass:  R = (1/M) sum(m_i r_i),  M = sum(m_i)
   internal forces only  ->  R(t) = R0 + V0 t  exactly

Angular momentum:  ell = r x p,  torque Gamma = r x F,  d(ell)/dt = Gamma
   central force F = f(r) r_hat  ->  Gamma = 0  ->  ell conserved
   ell conserved  ->  motion is planar; areal velocity dA/dt = ell/(2m) is
   constant (Kepler's second law)
```

### Where it lives

- `include/newton/core/particle.hpp` — `momentum(particle)` = `m v`, the
  single-particle building block.
- `include/newton/analysis/momentum.hpp` — `angular_momentum(mass, state,
  about)` = `(r - about) × (m v)`; `center_of_mass(particles)`;
  `total_momentum(particles)`, each taking a `std::vector<Particle>` since
  this v0.1 library has no multi-body `System` type yet (see
  `notes/04-extending.md` for that roadmap item).
- `include/newton/analytic/rocket.hpp` — `tsiolkovsky_delta_v(v_ex, m0, m)`
  and `rocket_velocity_gravity(v_ex, m0, burn_rate, g, t)`, the two closed
  forms above.
- `include/newton/systems/rocket.hpp` — `systems::simulate_rocket(params,
  dt, t_end)`. This is a **self-contained** simulator, not built on
  `Particle`/`simulate()`: a burning rocket's mass is a function of time,
  and folding `d(mv)/dt = F_ext` into the constant-mass machinery is
  simply the wrong equation (that extra `v dm/dt` term is a fictitious
  force with no physical source — see the file's own header comment and
  `docs/research/numerical-integration.md` §7.8). `simulate_rocket` runs
  its own RK4 on the scalar pair `(y, v)`, computing `m(t)` fresh from the
  closed form every call rather than accumulating it step by step.

### Examples

- `examples/04_rocket.cpp` — vertical launch, comparing
  `systems::simulate_rocket` against `analytic::rocket_velocity_gravity`
  during the burn, reporting burnout velocity and altitude.

### Analytic validation

The rocket functions in `analytic/rocket.hpp` are the direct oracle for
`systems::simulate_rocket`. There is no dedicated angular-momentum example
in this library's v0.1 example set, but `analysis::angular_momentum` is
exercised in the test suite on circular motion, where `ell` is exactly
constant (Kepler's second law with a trivial circular orbit).

---

## Chapter 4 — Energy (Taylor §4.1–4.10)

### Key equations

```
T = (1/2) m v^2                              kinetic energy
Delta T = W = integral(F . dr)               work-energy theorem
Conservative force:  F = -grad(U)
E = T + U = const                            (all forces conservative)
dE/dt = F_nonconservative . v                 (drag/damping: dissipation rate)

Common potentials:
   Uniform gravity:  U = m g y
   Hooke's spring:   U = (1/2) k x^2   (or (1/2) k |r - anchor|^2)
```

### Where it lives

- `include/newton/analysis/energy.hpp` — `kinetic_energy(mass, v)` = `½ m
  |v|²`; `total_energy(mass, state, U)`, a template over a potential
  callable `double U(const Vec3& r)`; `energy_drift(trajectory, U)`,
  which returns `(E(t) - E(0)) / |E(0)|` at every recorded sample — "the
  single most useful integrator diagnostic" per
  `docs/research/numerical-integration.md` §4.3 (see
  `notes/03-numerical-methods.md` for how to read it). Two ready-made
  potentials ship as small lambda factories: `uniform_gravity_U(mass, g)`
  and `spring_U(k)`.
- `include/newton/forces/spring.hpp` — `Spring{k, anchor, natural_length}`
  is the force side of Hooke's law; its potential is `spring_U(k)` above
  (isotropic case, `natural_length = 0`).

Note what's *not* here yet: Taylor §4.2–4.4's general `Potential` interface
with automatic `-∇U` gradients, the curl test for "is this force
conservative", and the 1D potential toolkit (turning points, period by
quadrature) are v0.2 roadmap items — see `notes/04-extending.md`'s
"potentials module" section. In v0.1, a potential is just a plain callable
you write by hand and pass to `total_energy`/`energy_drift`, matching the
handful of forces (uniform gravity, spring) that actually have one.

### Examples

No example is dedicated solely to Chapter 4 — instead, energy diagnostics
from `analysis/energy.hpp` are the load-bearing tool inside
`examples/05_oscillations.cpp` and especially
`examples/06_integrator_comparison.cpp`, which is entirely an energy-drift
comparison across the four integrators.

### Analytic validation

Chapter 4 doesn't get its own closed-form oracle file — instead, its
central claim (`E = T + U = const` for a conservative force) is exactly
what `energy_drift` measures, and Chapter 5's `SHM::energy(m)` (below) is
the concrete `E = ½ m ω0² A²` instance used throughout the test suite and
`examples/06_integrator_comparison.cpp`.

---

## Chapter 5 — Oscillations (Taylor §5.1–5.9)

The richest chapter in the book for closed-form validation, and the
centerpiece of `newton::analytic`.

### Key equations

```
Hooke's law:  F = -k x,  U = (1/2) k x^2

SHM:  x'' + omega0^2 x = 0,  omega0 = sqrt(k/m)
   x(t) = x0 cos(omega0 t) + (v0/omega0) sin(omega0 t)
   amplitude A = sqrt(x0^2 + (v0/omega0)^2)
   period T = 2 pi / omega0                    (independent of amplitude)
   energy E = (1/2) m omega0^2 A^2

2D oscillator:  x'' = -omega_x^2 x,  y'' = -omega_y^2 y
   isotropic (omega_x = omega_y): circle/ellipse
   anisotropic: Lissajous figure, closed iff omega_x:omega_y is rational

Damped:  x'' + 2 beta x' + omega0^2 x = 0,  beta = b/(2m)
   underdamped (beta < omega0):    x = e^(-beta t) [x0 cos(omega1 t) + ... sin(omega1 t)],  omega1 = sqrt(omega0^2 - beta^2)
   critically damped (beta = omega0):  x = (C1 + C2 t) e^(-beta t)          [fastest decay]
   overdamped (beta > omega0):     x = C+ e^(r+ t) + C- e^(r- t),  r+/- = -beta +/- sqrt(beta^2-omega0^2)
   log decrement:  ln(x_n / x_{n+1}) = 2 pi beta / omega1

Driven damped (steady state):  x'' + 2 beta x' + omega0^2 x = f0 cos(omega t)
   A(omega) = f0 / sqrt((omega0^2-omega^2)^2 + 4 beta^2 omega^2)
   tan(delta) = 2 beta omega / (omega0^2 - omega^2),  delta in (0, pi)
   resonant omega = sqrt(omega0^2 - 2 beta^2)
   Q = omega0 / (2 beta)
   FWHM (of A^2) = 2 beta
```

### Where it lives

- `include/newton/forces/spring.hpp` — `Spring` (isotropic, §5.1) and
  `AnisotropicSpring{kx, ky, kz}` (§5.3, Lissajous figures).
- `include/newton/forces/drives.hpp` — `SinusoidalDrive{F0, omega, phase,
  dir}` (§5.5), the external `F0 cos(ωt + phase)` term; `ConstantForce{F}`
  for a constant bias.
- `include/newton/forces/drag.hpp` — `LinearDrag{b}` doubles as the damping
  term: `beta = b/(2m)` links it directly to `DampedOscillator::beta`
  below. Combine `Spring + LinearDrag + SinusoidalDrive` via
  `forces::sum(...)` to build Taylor's full driven-damped oscillator.
- `include/newton/analytic/oscillator.hpp` — the centerpiece file:
  - `SHM{omega0, x0, v0}` — `.amplitude()`, `.phase()`, `.x(t)`, `.v(t)`,
    `.period()`, `.energy(mass)`.
  - `DampedOscillator{omega0, beta, x0, v0}` with `enum class Regime
    {Undamped, Underdamped, CriticallyDamped, Overdamped}`; `.regime()`
    classifies from `beta` vs `omega0`; `.x(t)`/`.v(t)` dispatch correctly
    in all four regimes (including the numerically delicate
    critically-damped case, handled by a relative-epsilon test around
    `beta == omega0` rather than an exact equality — see the file's own
    comment on why an exact match would be a genuine mathematical
    singularity of the under/overdamped formulas, not a physics one);
    `.log_decrement()`.
  - `DrivenOscillator{omega0, beta, f0, omega}` — steady-state only (no
    transient): `.amplitude()`, `.phase()`, `.resonant_omega()`,
    `.quality_factor()`, `.fwhm()`.
- `include/newton/analysis/oscillation.hpp` — measurements on a *simulated*
  `Trajectory`, so the numeric run can be checked against the analytic
  formulas above: `find_peaks(traj)` (3-point local-maxima detection, for
  log-decrement extraction from a damped ring-down); `measure_period(traj)`
  (mean spacing of linearly-interpolated upward zero crossings of `r.x`);
  `measure_amplitude(traj, t_start)` (max `|r.x|` after a transient has
  died out, for steady-state driven measurements).

### Examples

- `examples/05_oscillations.cpp` — SHM period and amplitude vs.
  `analytic::SHM`; all three damping regimes plotted on one CSV; a driven
  resonance sweep (amplitude vs. `ω` over `[0.2ω0, 2ω0]`) compared to
  `analytic::DrivenOscillator::amplitude()`.
- `examples/06_integrator_comparison.cpp` — the money demo: SHM run with
  all four integrators at `dt = T/100` over 200 periods, comparing final
  relative energy error. See `notes/03-numerical-methods.md` for the full
  discussion of what this example is teaching.

### Analytic validation

Every struct in `analytic/oscillator.hpp` doubles as its own test fixture:
`SHM::x(0) == x0` and `SHM::v(0) == v0` to the last bit (a deliberate
implementation choice — see the file's C++ note on why `x(t)` is computed
as `x0 cos(ω0 t) + (v0/ω0) sin(ω0 t)` rather than through an
amplitude/phase round trip); the damped regimes must be continuous in `t`
across the `beta ≈ omega0` boundary; the driven amplitude must peak near
`resonant_omega()`. This chapter is, by a wide margin, the richest source
of "run the simulator, compare to a formula" test cases in the whole
library (see `docs/research/taylor-ch1-7-map.md` §8, rows 23–33, for the
complete table this library draws from).

---

## Chapters 6–7 — why they're out of scope, and the future hook

Chapters 6 (calculus of variations) and 7 (Lagrange's equations) are a
change of *formalism*, not of physics — and this library is deliberately
built the other way around.

**Chapter 6 has no dynamics at all.** It's a boundary-value optimization
over *paths* (`δS[y] = 0` for a functional `S[y] = ∫ f(y, y', x) dx`),
where the shortest-path, catenary, and brachistochrone problems ask "which
curve minimizes this integral", not "where is the particle at time `t`
given where it started." There is no `State`, no `Force`, nothing for an
integrator to step forward. Trying to fold it into `newton::core` would
distort the whole design — `newton` is an initial-value time-stepper, full
stop.

**Chapter 7 doesn't add new physics — it re-derives the same equations of
motion** Chapters 1–5 already solve, via a different route: instead of
summing vector forces, you differentiate a scalar `L = T - U` with respect
to generalized coordinates. A force-based engine like this one can't use
`L` directly, because doing so needs two things this library doesn't have:
constraints handled implicitly (a bead confined to a wire, a pendulum's
fixed rod length), and derivatives that are taken symbolically or by
automatic differentiation rather than evaluated at a point. Neither has
anything in common with the Chapters 1–5 core.

There is, however, a genuinely *Newtonian* path to Chapter 7's results, and
it's the better first step if you ever want one (condensed from
`docs/research/taylor-ch1-7-map.md` §6.4–6.5 and §7.4–7.5 — read those
sections for the full argument):

- **Option A — a Cartesian constraint solver (recommended first step).**
  Add holonomic constraints `g(r) = 0` to the force engine and solve for
  the constraint force with a Lagrange multiplier, still entirely in
  `(r, v)` space: `[M Jᵀ; J 0][a; λ] = [F_applied; -J̇v]`. This gets
  pendulums, beads on wires, Atwood machines, and blocks on wedges *without
  leaving the force paradigm* — it reuses the existing integrators
  unchanged, and as a bonus it directly answers "what is the tension?" (the
  constraint force `λ Jᵀ`), which the pure-Lagrangian route in option B
  deliberately hides.
- **Option B — a true `newton::lagrangian` module.** Represent `T(q, q̇, t)`
  and `U(q, t)` as user-supplied callables, build the mass matrix `M_ij =
  ∂²T/∂q̇_i∂q̇_j` and the generalized-force vector numerically, solve the
  small dense linear system `M q̈ = b` each step, then hand the resulting
  `(q, q̇)` state to the *existing* RK4/Verlet step functions — a genuinely
  new formalism sitting on top of unchanged plumbing. The natural way to
  get the required derivatives without hand-differentiating every `T` and
  `U` a student writes is forward-mode automatic differentiation via dual
  numbers (`Dual<double>`, and `Dual<Dual<double>>` nested once for the
  second derivatives `M_ij` needs) — about 150 lines of header-only C++,
  and the reason the research map recommends templating the core scalar
  type early if this hook is ever wanted.

The single most convincing demonstration that both formalisms describe the
same physics — if this ever gets built — is solving a double pendulum both
ways (Option A's Cartesian constraints vs. Option B's generalized
coordinates) and checking they agree. See `notes/04-extending.md` for what
each option would actually take to implement, as a v0.2+ roadmap item.
