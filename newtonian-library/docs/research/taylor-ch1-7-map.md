# Taylor, *Classical Mechanics* (Ch. 1–7) → C++ Newtonian Simulation Library: Requirements Map

**Source:** John R. Taylor, *Classical Mechanics*, University Science Books, 2005.
**Audience:** design input for a C++ force-based Newtonian simulation library written alongside a Mechanics 1 course using this book.

---

## 0. Executive framing

Taylor's Chapters 1–5 are almost perfectly shaped for a force-based particle simulator. The book's pedagogical arc is:

1. **Ch. 1** — set up `m r̈ = F(r, ṙ, t)` and the coordinate systems it lives in.
2. **Ch. 2** — introduce forces that make the ODE *unsolvable analytically* (quadratic drag in 2D) — this is the book's own explicit argument for numerical integration, and therefore the library's core justification.
3. **Ch. 3–4** — conserved quantities (**p, L, E**) — these are exactly the diagnostics a simulator should track to prove it is correct.
4. **Ch. 5** — linear ODEs with complete closed-form solutions — the single richest source of validation test cases and of *frequency-domain* post-processing (resonance curves, Q).
5. **Ch. 6–7** — a change of *formalism*, not of physics. Out of scope for the force engine, but worth leaving a hook.

The design implication that runs through everything: **the library should be an ODE-integrator core + a force/potential model library + a diagnostics-and-analysis layer + an analytic-reference layer**, where the analytic-reference layer exists specifically to test the other three.

### Suggested top-level module layout (chapter-driven)

| Module (namespace) | Driven by | Contents |
|---|---|---|
| `nm::core` | Ch. 1 | `Vec2`/`Vec3`, `State`, `Particle`, `System`, units/constants |
| `nm::integrate` | Ch. 1–2 | Euler, RK2/RK4, RK45 adaptive, velocity Verlet, leapfrog, Yoshida4, Boris |
| `nm::forces` | Ch. 1–2, 5 | gravity, drag (linear/quadratic), Lorentz, springs, dampers, drives |
| `nm::potentials` | Ch. 4 | `U(r)`, `-∇U`, curl test, potential library, effective potential |
| `nm::coords` | Ch. 1 | Cartesian ↔ polar/cylindrical/spherical, polar accel decomposition |
| `nm::bodies` | Ch. 3 | variable-mass/rocket, CM computation, collisions, rigid-ish aggregates |
| `nm::events` | Ch. 2–5 | root-finding event detection (ground impact, apex, turning points, zero crossings) |
| `nm::analysis` | Ch. 3–5 | conserved-quantity monitors, range/apex, resonance sweep, FFT/Fourier, Poincaré, phase portraits |
| `nm::analytic` | Ch. 1–5 | closed-form reference solutions (the validation oracle) |
| `nm::io` | all | CSV/JSON trajectory export, gnuplot/matplotlib-friendly output |
| `nm::lagrangian` *(future)* | Ch. 6–7 | see §7.5 |

---

## Chapter 1 — Newton's Laws of Motion

### 1.1 Sections
- **1.1 Classical Mechanics** — domain of validity (v ≪ c, macroscopic).
- **1.2 Space and Time** — vectors, vector algebra, dot and cross products, differentiation of vectors.
- **1.3 Mass and Force** — operational definitions.
- **1.4 Newton's First and Second Laws; Inertial Frames**.
- **1.5 The Third Law and Conservation of Momentum** — weak and strong forms.
- **1.6 Newton's Second Law in Cartesian Coordinates** — the law as three uncoupled scalar ODEs.
- **1.7 Two-Dimensional Polar Coordinates** — unit vectors that rotate; the "extra" acceleration terms.

### 1.2 Key equations

```
r = (x, y, z)                                   position
v = dr/dt = ṙ ,   a = d²r/dt² = r̈
p = m v
F = dp/dt = m r̈                                 (Newton II, constant mass)
m ẍ = F_x(r, ṙ, t),  m ÿ = F_y,  m z̈ = F_z      (Cartesian: 3 coupled 2nd-order ODEs)
F_12 = -F_21                                    (Newton III, weak form)
P = Σ_α p_α ,  dP/dt = F_ext                    (momentum of a system)
a·b = |a||b|cosθ ,  a×b = |a||b| sinθ n̂

Two-dimensional polar coordinates:
r = r r̂
d r̂/dt =  φ̇ φ̂ ,   d φ̂/dt = -φ̇ r̂
v = ṙ r̂ + r φ̇ φ̂
a = (r̈ - r φ̇²) r̂ + (r φ̈ + 2 ṙ φ̇) φ̂
   → F_r = m(r̈ - r φ̇²)         [centripetal term]
   → F_φ = m(r φ̈ + 2 ṙ φ̇)      [Coriolis-like term]
```

Standard reduction used everywhere in the library:

```
y = (r, v) ∈ R^6 ,   dy/dt = f(y, t) = ( v , F(r, v, t)/m )
```

### 1.3 Canonical worked examples and problem archetypes
- Vector identities, components, and simple kinematics.
- **Block sliding down a frictionless / frictional incline** — decomposition along tilted axes.
- **Free fall and constant-acceleration kinematics**; vacuum projectile.
- **Atwood's machine** (returns in Ch. 7 as the Lagrangian showcase).
- **Conical pendulum / uniform circular motion** — the `r φ̇²` term.
- **The oscillating skateboard in a curved half-pipe** (§1.7) — polar coordinates → `φ̈ = -(g/R) sin φ` → small-angle SHM → forward reference to Ch. 5.
- **A puck on a rotating/frictionless surface**, motion described in polar coordinates.
- Problems asking to *show* that a given force conserves momentum for a two-body pair.

### 1.4 Simulation-library capabilities implied

**Types and core loop**
- `Vec2`, `Vec3` value types: `+ - * /`, `dot`, `cross`, `norm`, `norm2`, `normalized`, `angle`, component access. Header-only, `constexpr`-friendly, `double` by default but templated on scalar type (this pays off later for autodiff in §7).
- `struct State { Vec3 r; Vec3 v; double t; }`; a flat `std::array<double,6>` view for the integrator.
- `struct Particle { double m; double q; State s; std::string label; }`.
- `class System` holding N particles + a list of forces, exposing `step(dt)`, `run(t_end)`, `run_until(predicate)`.
- **Force interface**: `using ForceFn = std::function<Vec3(const Particle&, const System&, double t)>;` — everything downstream (drag, gravity, Lorentz, springs, drives) is one of these. Support both *per-particle body forces* and *pairwise interaction forces* (`Vec3 f(const Particle& a, const Particle& b)`), because Ch. 3–4 need pairwise forces to demonstrate Newton III and internal potential energy.
- Forces must be **summed**, not replaced — `System::netForce(i, t)`.

**Integrators (`nm::integrate`)** — this is the chapter where they should all appear, so that later chapters compare them:
- `EulerExplicit`, `EulerSemiImplicit (symplectic Euler)`, `RK2 (midpoint)`, `RK4`, `RKF45`/`DormandPrince54` with adaptive step + tolerance control, `VelocityVerlet`, `Leapfrog`, `Yoshida4`.
- Common interface: `step(state, dt, deriv_fn) -> state`, plus `order()`, `isSymplectic()`, `requiresVelocityIndependentForce()`.
- **Important design note:** velocity Verlet and the symplectic family assume `F = F(r)` only. Ch. 2's drag and Ch. 5's damping are velocity-dependent, so the library needs either (a) a documented guard/warning, (b) the velocity-Verlet variant with one fixed-point/Newton iteration on `v`, or (c) automatic routing to RK4/RK45 when any registered force declares `dependsOnVelocity() == true`. Make this an explicit force-trait flag.

**Coordinates (`nm::coords`)**
- `toPolar(Vec2) -> (r, φ)`, `fromPolar`, plus velocity/acceleration transforms.
- `polarAcceleration(r, ṙ, r̈, φ, φ̇, φ̈) -> (a_r, a_φ)` and its inverse — used to *verify numerically* that a Cartesian simulation reproduces the `r̈ - rφ̇²` decomposition. This is a great teaching-and-testing function.
- Cylindrical and spherical for Ch. 4.

**Diagnostics**
- `totalMomentum()`, per-particle `momentum()`.
- `Recorder` / `TrajectoryLog`: sample state at fixed output intervals (decoupled from the internal step, so adaptive integrators still produce evenly spaced output via dense output / interpolation).

### 1.5 Analytic references for validation (Ch. 1)
```
Constant force:      r(t) = r0 + v0 t + (1/2)(F/m) t²        exact to machine precision
Vacuum projectile:   x = v0 cosθ t ,  y = v0 sinθ t - g t²/2
   apex height  h = v0² sin²θ / (2g)
   time of flight T = 2 v0 sinθ / g
   range        R = v0² sin(2θ)/g ,  max at θ = 45°
Uniform circular motion: |a| = v²/r = ω² r, exact circle, |v| constant
Two-body internal force only: P(t) = P(0) exactly
```
A degree-2 polynomial trajectory is the sharpest first test: **RK4 and Verlet must both be exact (to roundoff) for constant force**, and explicit Euler must not be — that single test already discriminates integrator implementations.

---

## Chapter 2 — Projectiles and Charged Particles

This is the chapter that *justifies the whole library*. Taylor explicitly notes that 2D quadratic drag has no analytic solution and must be attacked numerically.

### 2.1 Sections
- **2.1 Air Resistance** — the empirical form, linear vs quadratic regimes, Reynolds number.
- **2.2 Linear Air Resistance** — separable horizontal and vertical motion; terminal velocity; time constant.
- **2.3 Trajectory and Range in a Linear Medium** — closed-form trajectory, transcendental range equation, small-drag expansion.
- **2.4 Quadratic Air Resistance** — 1D closed forms; **2D coupling → numerical solution required**.
- **2.5 Motion of a Charge in a Uniform Magnetic Field**.
- **2.6 Complex Exponentials** — mathematical interlude (also the toolkit for Ch. 5).
- **2.7 Solution for the Charge in a Magnetic Field** — helical motion.

### 2.2 Key equations

**Drag models**
```
f = -f(v) v̂ ,        f(v) = b v + c v²
b = β D ,  c = γ D²           (D = diameter of a sphere)
β ≈ 1.6e-4  N·s/m²  ,  γ ≈ 0.25 N·s²/m⁴     (sphere in air, STP)

linear drag:      f_lin  = -b v                       (vector, ∝ v)
quadratic drag:   f_quad = -c v² v̂ = -c |v| v         (vector)
ratio:            f_quad/f_lin = (γ/β) D v ≈ (1.6e3 s/m²)·D·v
```

**Linear medium (§2.2–2.3)** — with `y` up, `τ = m/b`, `v_ter = m g / b`:
```
horizontal:  m v̇x = -b vx
             vx(t) = vx0 e^{-t/τ}
             x(t)  = vx0 τ (1 - e^{-t/τ})   →   x_∞ = vx0 τ

vertical:    m v̇y = -m g - b vy
             vy(t) = -v_ter + (vy0 + v_ter) e^{-t/τ}
             y(t)  = -v_ter t + (vy0 + v_ter) τ (1 - e^{-t/τ})

trajectory:  y(x) = ((vy0 + v_ter)/vx0) x + v_ter τ ln(1 - x/(vx0 τ))
             (vertical asymptote at x = vx0 τ)

range (transcendental, solve numerically):
             R = vx0 τ (1 - e^{-T/τ}) with y(T)=0
small-drag expansion:  R ≈ R_vac [ 1 - (4/3)(vy0 / v_ter) + ... ]
```

**Quadratic medium (§2.4)**
```
horizontal only:  m v̇ = -c v²
                  v(t) = v0 / (1 + c v0 t / m)
                  x(t) = (m/c) ln(1 + c v0 t / m)     (grows without bound, logarithmically)

vertical drop from rest:  m v̇ = m g - c v²
                  v_ter = sqrt(m g / c)
                  v(t)  = v_ter tanh(g t / v_ter)
                  y(t)  = (v_ter²/g) ln cosh(g t / v_ter)

2D — NOT separable, no closed form:
                  m v̇x = -c sqrt(vx² + vy²) · vx
                  m v̇y = -m g - c sqrt(vx² + vy²) · vy
```

**Charge in a uniform magnetic field (§2.5, 2.7)**
```
m v̇ = q (v × B) ,   B = B ẑ
ω = q B / m                       (cyclotron frequency)
v̇z = 0                            → uniform drift along B
η ≡ vx + i vy ,   η̇ = -i ω η  →  η = η0 e^{-iωt}
transverse motion: circle of radius  r = m v⊥ / (|q| B) = v⊥/ω
period  T = 2π m / (|q| B)   (independent of speed)
full trajectory: a helix
|v| = const exactly (magnetic force does no work)
with E field:  m v̇ = q(E + v × B) ;  drift velocity v_d = (E × B)/B²
```

### 2.3 Canonical worked examples and problem archetypes
- **Falling object with linear drag** (oil drop, small sphere in a viscous fluid) — approach to terminal velocity; "how long to reach 90 % of `v_ter`".
- **Cart/puck coasting to rest with linear drag** — exponential decay, finite total distance.
- **Projectile range in a linear medium** — solve the transcendental range equation; verify the `(4/3)(vy0/v_ter)` correction against a numeric solve. (Taylor has an explicit *computer problem* asking exactly this.)
- **The baseball with quadratic drag** — the flagship example: `D ≈ 7 cm`, `m ≈ 0.15 kg`, launched at ~30 m/s at 45°; range collapses from the vacuum ~92 m to roughly half that. Includes the classic figure comparing vacuum / linear / quadratic trajectories on one plot.
- **Which drag dominates?** Compute the crossover speed `v = β/(γD)` for a baseball, a raindrop, a cannonball, a bacterium, a millikan oil drop.
- **Optimal launch angle with drag** — no longer 45°, must be found numerically by sweeping θ.
- **Cannonball / golf ball / cyclist coasting** variants.
- **Charged particle in uniform B** — helix radius and pitch; mass spectrometer / cyclotron sizing; `E × B` drift; particle in crossed fields tracing a cycloid.

### 2.4 Simulation-library capabilities implied

**Force models (`nm::forces`)**
```cpp
UniformGravity(Vec3 g);                       // g = {0,-9.81,0}
LinearDrag(double b);                         // F = -b v
QuadraticDrag(double c);                      // F = -c |v| v
CombinedDrag(double b, double c);
// convenience constructors from physical geometry:
static LinearDrag    fromSphereInAir(double diameter);   // b = 1.6e-4 * D
static QuadraticDrag fromSphereInAir(double diameter);   // c = 0.25 * D^2
static QuadraticDrag fromDragCoefficient(double Cd, double rho, double A); // c = Cd*rho*A/2
UniformMagneticField(Vec3 B);                 // F = q v × B
UniformElectricField(Vec3 E);
LorentzForce(Vec3 E, Vec3 B);
```
Each force should expose traits: `dependsOnVelocity()`, `isConservative()`, `potentialEnergy()` (nullable). Drag returns `false`/`nullopt`; this is how the Ch. 4 energy bookkeeping knows to account for dissipation separately.

**Event detection (`nm::events`)** — essential from here on:
- `EventFn = double(const State&)` with a sign change triggering a root solve.
- Built-ins: `GroundImpact(y = y0)`, `Apex(v_y = 0)`, `MaxSpeed`, `PlaneCrossing(n̂, d)`, `Zero crossing of any user functional`.
- Implementation: detect sign change over a step, then **Brent/bisection on the dense-output interpolant** (or a Hermite cubic from the endpoint states) so the impact time is accurate to the integrator's tolerance, not to `dt`. Getting range right requires this; naive "stop at the first negative y" gives O(dt) error in R and will destroy the comparison to the analytic linear-drag range.
- Actions on event: `Stop`, `Record`, `Reflect (bounce with restitution)`, `UserCallback`.

**Projectile analysis helpers (`nm::analysis`)**
```cpp
struct FlightResult { double range, apex, timeOfFlight, impactSpeed, impactAngle; Trajectory path; };
FlightResult launch(Projectile p, double speed, double angle, DragModel m);
double rangeVsAngle(...);                  // sweep θ, return optimal angle
double optimalLaunchAngle(...);            // golden-section / Brent maximize
Trajectory compareMedia(...);              // vacuum vs linear vs quadratic on identical ICs
double crossoverSpeed(double D);           // β/(γD)
double reynoldsNumber(double rho, double v, double D, double mu);
```

**Dedicated charged-particle integrator**
- Implement the **Boris pusher**. It conserves `|v|` to machine precision in a pure magnetic field and preserves the gyroradius over long integrations, whereas RK4 slowly spirals in/out and Euler blows up. Having Euler / RK4 / Boris side by side on a single helix makes an outstanding built-in demo *and* a rigorous regression test.
- Track: gyroradius, gyrophase, pitch angle, guiding-center position, `|v|` drift.

**Quantities to track / plot for this chapter**
`x(t), y(t), vx(t), vy(t)`, the `y` vs `x` trajectory, `|v|(t)` vs `v_ter`, `v(t)/v_ter` on a log axis (should be a straight line for linear drag), range vs launch angle, range vs drag coefficient, and the 3D helix with its projection onto the plane ⟂ B.

### 2.5 Validation targets (Ch. 2)
Every 1D closed form above is a unit test. Recommended assertions:
- Linear drag horizontal: `|x_num(t) - vx0 τ(1-e^{-t/τ})| < tol` over 10 τ.
- Linear drag 2D: compare against the **exact `y(x)` trajectory** — this is a rare case of an analytically solvable 2D drag problem and is therefore the single best 2D integrator test in the book.
- Quadratic vertical: `v(t) = v_ter tanh(gt/v_ter)`; check that `v → v_ter` monotonically and that `y(t) = (v_ter²/g) ln cosh(gt/v_ter)`.
- Magnetic: assert `|v|` constant to 1e-12, radius `= m v⊥/(qB)`, period `= 2πm/(qB)`, and closure of the orbit after one period.
- 2D quadratic drag has no oracle — validate by **convergence study** (halve `dt`, confirm error ratio → `2^p`) and by cross-checking two independent integrators (RK4 vs adaptive RK45 at tight tolerance).

---

## Chapter 3 — Momentum and Angular Momentum

### 3.1 Sections
- **3.1 Conservation of Momentum** — multiparticle systems, internal forces cancel pairwise.
- **3.2 Rockets** — variable-mass systems; the rocket equation.
- **3.3 The Center of Mass** — discrete and continuous; `M R̈ = F_ext`.
- **3.4 Angular Momentum for a Single Particle** — `ℓ = r × p`, torque, central forces, Kepler's second law.
- **3.5 Angular Momentum for Several Particles** — `dL/dt = Γ_ext`, orbital + spin decomposition, moment of inertia about a fixed axis.

### 3.2 Key equations

```
P = Σ_α m_α v_α ,   dP/dt = F_ext        (internal forces cancel by Newton III)

Rocket (Taylor's form):
    m dv/dt = -ṁ v_ex + F_ext ,   with ṁ = dm/dt < 0
    thrust  = -ṁ v_ex = |ṁ| v_ex
    F_ext = 0  →  v - v0 = v_ex ln(m0/m)            (Tsiolkovsky)
    F_ext = -mg (vertical, uniform g)
            →  v(t) = v_ex ln(m0/m(t)) - g t
    with constant burn rate m(t) = m0 - |ṁ| t
    liftoff condition: |ṁ| v_ex > m0 g

Center of mass:
    R = (1/M) Σ_α m_α r_α ,   M = Σ m_α
    R = (1/M) ∫ r dm                    (continuous)
    P = M Ṙ ,  M R̈ = F_ext
    → with only internal forces: R(t) = R0 + V0 t   exactly

Angular momentum:
    ℓ = r × p ,   Γ = r × F ,   dℓ/dt = Γ
    central force F = f(r) r̂  →  Γ = 0  →  ℓ conserved
    ℓ conserved ⇒ motion is planar (r ⟂ ℓ always)
    areal velocity  dA/dt = ℓ/(2m) = const           (Kepler's 2nd law)
    in polar form:  ℓ = m r² φ̇

Many particles:
    L = Σ_α r_α × p_α ,   dL/dt = Γ_ext
    L = R × P + Σ_α r'_α × p'_α = L_orbital + L_spin(about CM)
    fixed-axis rotation:  L_z = I ω ,  I = Σ m_α ρ_α²
```

### 3.3 Canonical worked examples and problem archetypes
- **Single-stage rocket in free space** — Tsiolkovsky; mass-ratio arguments ("why chemical rockets need staging").
- **Rocket launched vertically in uniform gravity** — `v(t) = v_ex ln(m0/m) - gt`; find burnout velocity and altitude.
- **Multistage rocket** — compare final speed to a single stage of the same total mass.
- **A rocket with linear air resistance** (problem-level) — no closed form ⇒ numerical.
- **CM of composite / continuous bodies**: solid cone, hemisphere, semicircular wire, L-shaped plate.
- **Exploding shell** — projectile bursts at apex; the CM continues on the unperturbed parabola. Superb simulator demo.
- **Elastic and inelastic collisions**; 2-body recoil; a person walking on a boat/cart.
- **Kepler's second law** — equal areas in equal times, demonstrated for any central force.
- **Ice-skater / rotating-platform** angular-momentum problems.
- **A ball on a string wrapping around a post** vs **passing through a hole in a table** (torque vs no torque).

### 3.4 Simulation-library capabilities implied

**Multi-particle system state**
```cpp
class System {
  Vec3 centerOfMass() const;
  Vec3 centerOfMassVelocity() const;
  Vec3 totalMomentum() const;                    // P
  Vec3 totalAngularMomentum(Vec3 about = {0,0,0}) const;   // L
  Vec3 orbitalAngularMomentum() const;           // R × P
  Vec3 spinAngularMomentum() const;              // about the CM
  System inCMFrame() const;                      // boost to the CM frame
  double reducedMass(int i, int j) const;        // μ = m1 m2/(m1+m2)
};
```
The **CM-frame transform** is worth building properly: it is used in Ch. 3 collisions, Ch. 4 two-particle energy, and Ch. 8 two-body/Kepler if the student continues.

**Variable-mass bodies (`nm::bodies`)**
```cpp
class RocketBody {
  double m0, m_dry;
  double v_ex;                 // exhaust speed
  std::function<double(double t)> mdot;   // burn-rate schedule, ≤ 0
  Vec3 thrustDirection(const State&) const;   // fixed, prograde, or steered
};
// EOM implemented as:  a = (thrust)/m(t) + F_ext/m(t) ,  thrust = |mdot| * v_ex * n̂
// Mass must be updated inside the integrator's derivative call, not once per step.
```
Design note: variable mass means `F = ma` is **not** `F = dp/dt` in the naive sense. The library should model this as an explicit `VariableMassBody` whose derivative function returns `(ṙ, v̇, ṁ)` — i.e., the state vector gains a mass component — rather than trying to shoehorn it into a constant-mass `Particle`. Add multi-stage support (`std::vector<Stage>` with jettison events at burnout, hooked into the event system).

**Collisions / explosions**
```cpp
elasticCollision1D(m1,v1,m2,v2) -> (v1',v2');   // closed form
collide(Particle&, Particle&, double restitution);
explode(Particle, std::vector<double> fragmentMasses, std::vector<Vec3> impulses);
```
Assert `P` conserved always, `E` conserved iff `e = 1`.

**Conservation monitors (`nm::analysis`)**
```cpp
class ConservationMonitor {   // sample every N steps
  double momentumDrift() const;        // |P(t) - P(0)| / |P(0)|
  double angularMomentumDrift() const;
  double energyDrift() const;          // Ch. 4
  double arealVelocity() const;        // ℓ/(2m) for a central-force run
};
```
These should be first-class, plottable outputs — "did the integrator conserve what the physics says it should?" is the single most useful debugging view a mechanics student can have.

**Quantities to track / plot**
`P(t)`, `L(t)`, `R_cm(t)` overlaid on individual particle tracks, swept area vs time for a central force, rocket `v(t)` and `m(t)`, thrust-to-weight vs time.

### 3.5 Validation targets (Ch. 3)
```
Tsiolkovsky:        v(m) = v0 + v_ex ln(m0/m)                  exact
Rocket + gravity:   v(t) = v_ex ln(m0/(m0 - |ṁ|t)) - g t       exact
CM theorem:         with internal forces only, R(t) = R0 + V0 t exactly
                    (test: N-body gravitating cluster — CM must move in a straight line)
Exploding projectile: CM follows the undisturbed parabola after the burst
Central force:      ℓ = m r² φ̇ = const ; dA/dt = ℓ/2m = const
1D elastic collision: v1' = ((m1-m2)v1 + 2m2v2)/(m1+m2), etc.
```

---

## Chapter 4 — Energy

### 4.1 Sections
- **4.1 Kinetic Energy and Work** — work–KE theorem, line integrals, power.
- **4.2 Potential Energy and Conservative Forces** — path independence; two conditions.
- **4.3 Force as the Gradient of Potential Energy** — `F = -∇U`.
- **4.4 The Second Condition that F be Conservative** — `∇ × F = 0`.
- **4.5 Time-Dependent Potential Energy** — `dE/dt = ∂U/∂t`.
- **4.6 Energy for Linear One-Dimensional Systems** — `E = ½mv² + U(x)`, turning points, quadrature, equilibria and stability.
- **4.7 Curvilinear One-Dimensional Systems** — motion constrained to a track, arclength coordinate `s`.
- **4.8 Central Forces** — spherical symmetry ⇔ conservative; gradient in spherical coordinates.
- **4.9 Energy of Interaction of Two Particles** — `U(r1 - r2)`, `U(|r1 - r2|)`.
- **4.10 The Energy of a Multiparticle System** — `E = T + U_ext + U_int`; rigid bodies.

### 4.2 Key equations

```
T = ½ m v²  ,  ΔT = W = ∫_1^2 F · dr  ,  P = F · v

Conservative force (two conditions):
    (i)  F = F(r) only     (ii) W path-independent  ⇔  ∇ × F = 0
U(r) = -∫_{r0}^{r} F · dr'     (path-independent)
F = -∇U = -(∂U/∂x, ∂U/∂y, ∂U/∂z)
E = T + U = const     (all forces conservative)
dE/dt = ∂U/∂t         (explicitly time-dependent U)
dE/dt = F_nonconservative · v   (drag, damping → dissipation rate = -f·v)

1D energy method:
    E = ½ m ẋ² + U(x)
    v(x) = ± sqrt( 2(E - U(x)) / m )
    t(x) - t0 = ∫_{x0}^{x} dx' / v(x')                  (quadrature)
    turning points: U(x) = E
    equilibrium:    U'(x0) = 0
    stable iff U''(x0) > 0
    small oscillations: ω = sqrt( U''(x0) / m )          → bridge to Ch. 5
    period of bounded motion:
        T = 2 ∫_{x1}^{x2} dx / sqrt(2(E - U(x))/m)       (integrable singularity at endpoints)

Central force: F = f(r) r̂ ,  U = -∫ f(r) dr ,  ∇U = (∂U/∂r) r̂
Two particles: U = U(r1 - r2) ;  F_12 = -∇_1 U ,  F_21 = -∇_2 U = +∇_1 U
Multiparticle:  E = Σ ½ m_α v_α² + Σ_α U_α^ext + Σ_{α<β} U_{αβ}
```

**Potentials that appear in the chapter and problems**
```
Uniform gravity:      U = m g y
Hooke:                U = ½ k x²           (or ½ k (r - r0)²)
Newtonian gravity:    U = -G m1 m2 / r
Coulomb:              U = k q1 q2 / r
Morse (diatomic):     U(r) = A[ e^{-2(r-r0)/S} - 2 e^{-(r-r0)/S} ]
Lennard-Jones:        U(r) = U0[ (r0/r)^12 - 2 (r0/r)^6 ]
Yukawa:               U(r) = -(k/r) e^{-r/a}
```

### 4.3 Canonical worked examples and problem archetypes
- **Work done by gravity / by a spring** along different paths; show path independence.
- **Curl test problems**: given `F = (ay, bx, 0)` etc., determine for which parameters it is conservative, and find `U`.
- **Energy diagrams**: sketch `U(x)`, mark turning points for a given `E`, classify bound vs unbound motion, identify stable/unstable equilibria.
- **The diatomic molecule / Morse potential** — find `r0`, show `U'' = 2A/S²`, get the vibrational frequency `ω = sqrt(2A/(mS²))`, compare small-oscillation vs anharmonic behavior.
- **Lennard-Jones** — same treatment; equilibrium separation and well depth.
- **Bead on a wire / cycloidal track** (§4.7) — curvilinear 1D; foreshadows the brachistochrone in Ch. 6.
- **Escape velocity** from `U = -GMm/r`; radial infall.
- **A block on a spring with friction** — energy lost per cycle.
- **Time-dependent potential** (e.g. a spring whose support is being moved) — energy not conserved, `dE/dt = ∂U/∂t`.
- **Two-particle interaction energy** and the CM/relative-coordinate split.

### 4.4 Simulation-library capabilities implied

**Potential abstraction (`nm::potentials`)** — the central new concept of this chapter:
```cpp
class Potential {
public:
  virtual double  U(const Vec3& r, double t = 0) const = 0;
  virtual Vec3    gradU(const Vec3& r, double t = 0) const;   // default: central differences
  Vec3 force(const Vec3& r, double t = 0) const { return -gradU(r, t); }
  virtual bool isTimeDependent() const { return false; }
  virtual double dUdt(const Vec3& r, double t) const { return 0.0; }
};
```
- Provide **analytic gradients** for the library potentials, and a **numerical gradient fallback** (Richardson-extrapolated central differences) so a student can define a new `U` with one line and get a force for free. Ship a unit test that checks analytic gradients against numerical ones for every built-in potential.
- `CentralPotential : Potential` taking `U(r)` of the scalar radius, with `∇U = U'(r) r̂` handled automatically.
- `PairPotential` for `U(|r_i - r_j|)`, generating equal-and-opposite forces (guarantees Newton III by construction — a nice architectural nod to Ch. 3).

**Conservativeness checking**
```cpp
bool isConservative(const ForceFn& F, const Box& region, double tol);   // numerical ∇×F ≈ 0
Vec3 curl(const ForceFn& F, const Vec3& r, double h);
double workAlongPath(const ForceFn& F, const Path& p);                  // line-integral quadrature
double pathDependence(const ForceFn& F, const Path& a, const Path& b);  // W_a - W_b
```
This directly implements §4.2/4.4 and makes an excellent interactive demo: pick a force field, pick two paths between the same endpoints, print both works.

**Energy bookkeeping (`nm::analysis`)**
```cpp
double kineticEnergy() const;
double potentialEnergy() const;         // sums Potential contributions
double totalEnergy() const;
double dissipatedEnergy() const;        // running ∫ f_nc · v dt for drag/damping
double energyDrift() const;             // (E(t) - E(0)) / |E(0)|
```
Key requirement: because Ch. 2's drag and Ch. 5's damping are non-conservative, the energy accounting must track **T + U + W_dissipated = const** so that a damped run still has an exactly conserved bookkeeping quantity to test against.

**1D potential toolkit** — probably the most valuable single module for Ch. 4 homework:
```cpp
namespace nm::u1d {
  std::vector<double> findEquilibria(Fn U, double a, double b);   // roots of U'
  Stability classify(Fn U, double x0);                            // sign of U''
  double smallOscillationFrequency(Fn U, double m, double x0);    // sqrt(U''/m)
  std::pair<double,double> turningPoints(Fn U, double E, double xGuess);
  double periodByQuadrature(Fn U, double m, double E);            // handle sqrt singularity
  Curve energyDiagram(Fn U, double a, double b, double E);        // U(x) with E line + turning pts
  Curve phasePortrait(Fn U, double m, std::vector<double> energies); // v vs x level curves
}
```
For `periodByQuadrature`, use the substitution `x = ½(x1+x2) + ½(x2-x1) cos θ` (Gauss–Chebyshev) to remove the inverse-square-root endpoint singularity — a naive Simpson rule will fail here and the failure is instructive but should not be the default.

**Effective potential (forward hook to Ch. 8)**
```cpp
double U_eff(double r) { return U(r) + l*l/(2*mu*r*r); }
```
Cheap to add now, and it makes Kepler orbits nearly free later.

**Integrator implications**
- This is where **symplectic integrators earn their place**. Ship a demo comparing energy drift over 10⁴ periods: RK4 shows secular linear drift; velocity Verlet / Yoshida4 show bounded oscillation. Also implement a time-reversibility test (integrate forward `N` steps, negate `v`, integrate back; Verlet returns to the start, RK4 does not).
- Offer an optional **energy-projection** post-step (rescale `v` to restore `E`) as a comparison, and be explicit in the docs that it is a numerical hack, not physics.

### 4.5 Validation targets (Ch. 4)
```
Free fall:              ½mv² + mgy = const
Spring:                 ½mv² + ½kx² = ½kA² ; amplitude exactly recovered
Kepler/Coulomb:         E = ½μv² - GMm/r = -GMm/(2a) for a bound orbit
Escape:                 v_esc = sqrt(2GM/R)  — launch at exactly v_esc, r → ∞ with v → 0
Morse:                  r0 = minimum, U''(r0) = 2A/S², ω = sqrt(2A/(mS²));
                        small-amplitude numeric period must match 2π/ω to <1%
Lennard-Jones:          equilibrium at r = r0, U(r0) = -U0, U''(r0) = 72 U0/r0²
Turning points:         numeric x_max must satisfy U(x_max) = E to tolerance
Anharmonic period:      quadrature period must match direct time-domain period
Curl test:              ∇×F = 0 numerically for every built-in conservative force
Dissipation:            for linear drag, T + U + ∫ b v² dt = const
```

---

## Chapter 5 — Oscillations

The richest chapter for closed-form validation, and the one that most demands *analysis* features (frequency sweeps, Fourier), not just integration.

### 5.1 Sections
- **5.1 Hooke's Law** — `F = -kx` as the leading term of any smooth potential near a minimum.
- **5.2 Simple Harmonic Motion** — four equivalent solution forms; energy.
- **5.3 Two-Dimensional Oscillators** — isotropic vs anisotropic; Lissajous figures.
- **5.4 Damped Oscillations** — under/critical/over damping.
- **5.5 Driven Damped Oscillations** — particular + transient solutions.
- **5.6 Resonance** — amplitude and phase vs driving frequency; width; **Q**.
- **5.7 Fourier Series** *(starred)* — periodic driving forces.
- **5.8 Fourier Series Solution for the Driven Oscillator** *(starred)* — superposition.
- **5.9 The RMS Displacement; Parseval's Theorem** *(starred)*.

### 5.2 Key equations

```
Hooke:   F = -k x ,  U = ½ k x²
General: U(x) ≈ U(x0) + ½ U''(x0)(x-x0)²   →  k_eff = U''(x0)

SHM:     ẍ = -ω0² x ,   ω0 = sqrt(k/m)
   x(t) = A cos(ω0 t - δ) = C1 cos ω0 t + C2 sin ω0 t = Re[ B e^{i ω0 t} ]
   A = sqrt(C1² + C2²) ,  tan δ = C2/C1
   with x(0)=x0, v(0)=v0:  A = sqrt(x0² + (v0/ω0)²) , tan δ = v0/(ω0 x0)
   T = 2π/ω0   (independent of amplitude)
   E = ½ k A² = ½ m ω0² A²   (constant; T and U each average to E/2)

2D oscillator:
   ẍ = -ωx² x ,  ÿ = -ωy² y
   isotropic (ωx = ωy): ellipse (or line/circle)
   anisotropic: Lissajous figure; closed iff ωx:ωy is rational

Damped:  m ẍ = -k x - b ẋ
   ẍ + 2β ẋ + ω0² x = 0 ,   β = b/(2m)
   roots r± = -β ± sqrt(β² - ω0²)
   underdamped (β < ω0):  x = A e^{-βt} cos(ω1 t - δ) ,  ω1 = sqrt(ω0² - β²)
   critically damped (β = ω0): x = (C1 + C2 t) e^{-βt}         (fastest decay)
   overdamped (β > ω0):  x = C1 e^{-(β - sqrt(β²-ω0²))t} + C2 e^{-(β + sqrt(β²-ω0²))t}
   energy decay: E(t) ≈ E0 e^{-2βt}  (weak damping)
   logarithmic decrement:  ln(x_n / x_{n+1}) = 2π β / ω1
   decay time τ = 1/β

Driven damped:
   ẍ + 2β ẋ + ω0² x = f(t) ,  f = F/m ;  f(t) = f0 cos(ω t)
   x(t) = x_p(t) + x_tr(t) ,   x_tr → 0 like e^{-βt}
   steady state: x_p = A cos(ω t - δ)
       A(ω) = f0 / sqrt( (ω0² - ω²)² + 4 β² ω² )
       tan δ = 2βω / (ω0² - ω²) ,   δ ∈ (0, π)
       δ → 0 (ω ≪ ω0) , δ = π/2 at ω = ω0 , δ → π (ω ≫ ω0)

Resonance:
   amplitude peak at  ω_res = sqrt(ω0² - 2β²)
   A²(ω) = f0² / ((ω0² - ω²)² + 4β²ω²)     ≈ Lorentzian near ω0 for β ≪ ω0
   FWHM of A²(ω) = 2β
   Q = ω0 / (2β) ≈ 2π (energy stored) / (energy lost per cycle)
   A(ω0) = f0/(2βω0) = Q f0/ω0²  →  amplitude gain of Q over the DC response f0/ω0²
   undamped exact resonance (β=0, ω=ω0): x(t) = (f0/(2ω0)) t sin(ω0 t)   (secular growth)
   undamped near resonance: beats at |ω - ω0|/2

Fourier (starred):
   f(t) = Σ_{n} [a_n cos(nωt) + b_n sin(nωt)] ,  ω = 2π/τ
   response x(t) = Σ_n A_n cos(nωt - δ_n)  with A_n, δ_n from the formulas above at ω → nω
   Parseval:  <x²> = ½ Σ A_n²  ,  x_rms = sqrt(<x²>)
```

### 5.3 Canonical worked examples and problem archetypes
- **Mass on a spring**, vertical spring (shifted equilibrium), two springs in series/parallel.
- **Simple pendulum small-angle** `ω0 = sqrt(g/L)`; physical pendulum `ω0 = sqrt(mgd/I)`.
- **A cube/block floating in liquid** — buoyancy restoring force gives SHM with `ω = sqrt(ρ_liq g A / m)`.
- **LRC circuit** as the exact electrical analogue (`L q̈ + R q̇ + q/C = V0 cos ωt`) — the library should note the mapping so the same solver serves both.
- **Q-factor estimation examples**: a bicycle wheel spinning on its axle, a tuning fork, a plucked guitar string, an atom (very high Q), a car suspension (very low Q, near critical).
- **Shock absorbers / door closers** — why critical damping is the design target; compare settling time across the three regimes.
- **Ring-down measurement**: given a decaying trace, extract `β`, `ω1`, and `Q`.
- **Resonance disasters** — Tacoma Narrows-style narrative; sweeping ω through ω0.
- **Transient + steady-state superposition** — plot the full `x(t)` showing the transient dying out.
- **Lissajous figures** for `ωx:ωy = 1:1, 1:2, 2:3, √2:1` (the last never closes).
- **Square-wave / sawtooth driving force** (§5.8) — Fourier-decompose, superpose responses, note that the `n`th harmonic nearest `ω0` dominates.

### 5.4 Simulation-library capabilities implied

**Force / drive models**
```cpp
HookeSpring(double k, Vec3 anchor, double naturalLength = 0);
AnisotropicSpring(Vec3 k);                       // kx, ky, kz  → Lissajous
LinearDamper(double b);                          // F = -b v   (β = b/2m)
SinusoidalDrive(Vec3 F0, double omega, double phase = 0);
ArbitraryDrive(std::function<Vec3(double t)>);
FourierDrive(std::vector<std::pair<double,double>> ab, double omega);  // §5.7-5.8
SquareWaveDrive / SawtoothDrive / ImpulseTrain;
```

**Analytic reference class (`nm::analytic`)** — a first-class deliverable, not just a test fixture:
```cpp
class DampedDrivenOscillator {
public:
  DampedDrivenOscillator(double omega0, double beta, double f0, double omega);
  Regime regime() const;                 // Undamped | Under | Critical | Over
  double omega1() const;                 // sqrt(ω0² - β²)
  double x(double t) const;              // full solution incl. transient, from ICs
  double v(double t) const;
  double steadyStateAmplitude() const;   // A(ω)
  double steadyStatePhase() const;       // δ(ω)
  double resonantFrequency() const;      // sqrt(ω0² - 2β²)
  double Q() const;                      // ω0/(2β)
  double fwhm() const;                   // 2β
  double logDecrement() const;           // 2πβ/ω1
  double energy(double t) const;
};
```
Having this alongside the numeric integrator lets *every* Ch. 5 experiment be a self-checking test.

**Analysis features (`nm::analysis`)**
```cpp
struct ResonanceCurve { std::vector<double> omega, amplitude, phase; };
ResonanceCurve frequencySweep(Oscillator o, double wMin, double wMax, int n,
                              double settleTime /* ≫ 1/β */);
double measureAmplitude(const Trajectory&, double tStart);  // steady-state peak-to-peak/2
double measurePhaseLag(const Trajectory&, double omega);    // cross-correlation with the drive
double fitQ_fromFWHM(const ResonanceCurve&);
double fitQ_fromRingdown(const Trajectory&);                // exponential fit to peak envelope
struct PeakList { std::vector<double> t, x; };
PeakList findPeaks(const Trajectory&);                      // for log-decrement extraction
double measurePeriod(const Trajectory&);                    // zero-crossing based
std::vector<Complex> dft(const std::vector<double>& x);     // §5.7-5.9; a small radix-2 FFT suffices
double rms(const std::vector<double>& x);                   // Parseval check
Curve phasePortrait(const Trajectory&);                     // v vs x — spiral for damped, ellipse for SHM
PoincareSection poincare(const Trajectory&, double drivePeriod);
```
The **frequency sweep** deserves emphasis: it takes the simulator from "solve one ODE" to "reproduce Figure 5.16 of the textbook." Run N independent simulations at different ω, discard `~5/β` of transient, measure the steady-state amplitude and phase, and plot against the analytic Lorentzian. That is a complete, self-validating lab.

**Numerical caveats to build in**
- Strong damping (`β ≫ ω0`) makes the system **stiff**: the fast mode `e^{-2βt}` forces tiny steps in explicit methods. Either document a `dt < 1/(2β)` stability bound, add a check that warns when `β·dt > 0.1`, or ship an implicit/exponential integrator for the linear oscillator.
- Velocity-dependent damping breaks plain velocity Verlet's symplecticity — route damped runs to RK4/RK45 by default, or use the standard modified-Verlet formulation for a linear damper.
- Resonance amplitudes grow as `Q`; check for overflow / provide amplitude limiting in demos with `β → 0`.

**Quantities to track / plot**
`x(t)` with the `±A e^{-βt}` envelope overlaid, the three damping regimes on one axis, `E(t)` vs `E0 e^{-2βt}`, `A(ω)` and `δ(ω)` (the two-panel resonance figure), `A²(ω)` with FWHM marked, phase portraits, Lissajous curves, transient-plus-steady-state decomposition, the Fourier spectrum of the response to a non-sinusoidal drive.

### 5.5 Validation targets (Ch. 5)
Effectively everything in this chapter is a closed form — see §8 for the consolidated table. Highest-value tests:
- SHM period independent of amplitude to 1e-10; energy `½kA²` constant.
- Underdamped trace bounded by `A e^{-βt}` with the peak sequence matching the log decrement.
- Critical damping reaches within ε of zero faster than any nearby over/under-damped case (a genuine assertion, not just a plot).
- Measured steady-state `A(ω)` and `δ(ω)` match the analytic formulas to <0.1 % across a sweep.
- Measured FWHM = `2β` and `Q = ω0/(2β)` from two independent methods (sweep and ring-down) agree.
- **Nonlinear cross-check:** the full pendulum `θ̈ = -(g/L) sin θ` has exact period `T = 4 sqrt(L/g) K(sin(θ0/2))`; the numeric period must match, and must approach `T0[1 + θ0²/16 + 11θ0⁴/3072 + …]` for small `θ0`. This is the best test that the integrator handles genuine nonlinearity correctly.

---

## Chapter 6 — Calculus of Variations *(out of scope; hook only)*

### 6.1 Sections
- **6.1 Two Examples** — the shortest path between two points in a plane; the brachistochrone.
- **6.2 The Euler–Lagrange Equation**.
- **6.3 Applications of the Euler–Lagrange Equation** — shortest path, minimum surface of revolution, brachistochrone, Fermat's principle.
- **6.4 More than Two Variables** — several dependent variables → a set of EL equations (the direct bridge to Ch. 7).

### 6.2 Key equations
```
S[y] = ∫_{x1}^{x2} f(y, y', x) dx        (a functional)
stationary S  ⇔  ∂f/∂y - d/dx (∂f/∂y') = 0        (Euler–Lagrange)
several variables:  ∂f/∂y_i - d/dx (∂f/∂y'_i) = 0 ,  i = 1..n

If f has no explicit x:  f - y' ∂f/∂y' = const     (Beltrami first integral)
If f has no explicit y:  ∂f/∂y' = const

Shortest path:      f = sqrt(1 + y'²)   →  y' = const  →  straight line
Minimum surface of revolution: f = y sqrt(1+y'²)  →  y = y0 cosh((x - x0)/y0)   (catenary)
Brachistochrone:    f = sqrt((1+y'²)/y)  →  cycloid
                    x = a(θ - sin θ) ,  y = a(1 - cos θ)
                    descent time  t = sqrt(a/g) · θ_max   (tautochrone property)
Fermat:             δ∫ n ds = 0   →  Snell's law  n1 sin θ1 = n2 sin θ2
Geodesic on a sphere: great circle
```

### 6.3 Canonical examples and problem archetypes
Shortest path in a plane and on a cylinder/cone/sphere; brachistochrone; catenary / soap-film between two rings (including the two-solution / Goldschmidt no-solution regime); Fermat's principle → reflection and refraction; the shortest path in polar coordinates; Dido's isoperimetric problem (as a constrained/starred extension).

### 6.4 Why this is out of scope for a Newtonian engine
Chapter 6 contains **no dynamics**. It is a *boundary-value* optimization over paths, whereas the entire library architecture built for Ch. 1–5 is an *initial-value* time-stepper. There is no force, no state vector evolving in `t`, and nothing for an integrator to do. Attempting to fold it into `nm::core` would distort the design.

### 6.5 What a future hook could look like
If the student wants Ch. 6 support later, it belongs in a separate, decoupled module that shares only `Vec` and the quadrature routines:

```cpp
namespace nm::variational {

  // A functional S[y] = ∫ f(y, y', x) dx
  using Integrand = std::function<double(double y, double yp, double x)>;

  double action(Integrand f, const Path& y, int quadOrder = 5);   // Simpson/Gauss-Legendre

  // (a) Direct method: discretize the path into N free interior nodes and minimize S
  Path minimizePath(Integrand f, BoundaryConditions bc, int nNodes,
                    Optimizer opt = Optimizer::LBFGS);

  // (b) Indirect method: form the EL equation and solve the resulting 2nd-order BVP
  //     y'' = G(y, y', x)   via shooting + Brent root-find on the initial slope,
  //     reusing the existing IVP integrators — this is the one place Ch.6 can
  //     genuinely reuse nm::integrate.
  Path solveEulerLagrangeBVP(Integrand f, BoundaryConditions bc);

  double eulerLagrangeResidual(Integrand f, const Path& y);  // ≈0 for a true extremal
}
```
Validation is unusually clean because all three showcase problems have exact answers: `minimizePath` on `sqrt(1+y'²)` must return a straight line; on `y sqrt(1+y'²)` a catenary; on the brachistochrone integrand a cycloid with the analytically known descent time. Derivatives `∂f/∂y`, `∂f/∂y'` are obtained by finite differences or, better, by the same forward-mode dual-number autodiff recommended in §7.5 — which is the main reason to template the core scalar type early.

A second, more physical hook: the **brachistochrone as a Newtonian simulation** — build the cycloid track, put a bead on it, integrate the constrained motion with the force engine, and confirm the descent time matches `sqrt(a/g)·θ_max` and the tautochrone property (release from any height, same arrival time). That test exercises Ch. 4's curvilinear-1D machinery and needs no variational code at all.

---

## Chapter 7 — Lagrange's Equations *(out of scope; hook only)*

### 7.1 Sections
- **7.1 Lagrange's Equations for Unconstrained Motion**.
- **7.2 Constrained Systems; an Example**.
- **7.3 Constrained Systems in General** — holonomic constraints, generalized coordinates, degrees of freedom.
- **7.4 Proof of Lagrange's Equations with Constraints** — D'Alembert / virtual work.
- **7.5 Examples of Lagrange's Equations**.
- **7.6 Generalized Momenta and Ignorable Coordinates**.
- **7.7 Conclusion**.
- **7.8 More about Conservation Laws** *(starred)* — symmetry → conservation.
- **7.9 Lagrange's Equations for Magnetic Forces** *(starred)*.
- **7.10 Lagrange Multipliers and Constraint Forces** *(starred)*.

### 7.2 Key equations
```
L(q, q̇, t) = T - U
d/dt ( ∂L/∂q̇_i ) = ∂L/∂q_i        i = 1..n        (Lagrange's equations)

Degrees of freedom: n = 3N - (number of holonomic constraints)
Holonomic constraint: f(r1,...,rN, t) = 0
Generalized coordinates q1..qn ; r_α = r_α(q1..qn, t)

Generalized momentum:   p_i = ∂L/∂q̇_i
Generalized force:      Q_i = Σ_α F_α · ∂r_α/∂q_i
Ignorable (cyclic) coordinate: ∂L/∂q_i = 0  ⇒  p_i = const
Hamiltonian:  H = Σ_i p_i q̇_i - L ;  dH/dt = -∂L/∂t  ⇒  H conserved if L has no explicit t
              H = T + U when constraints are time-independent and T is quadratic in q̇

Constraint forces via multipliers (§7.10):
  d/dt(∂L/∂q̇_i) - ∂L/∂q_i = λ ∂f/∂q_i

Magnetic (§7.9):  L = ½ m v² - q(V - v·A)

Kinetic energy in generalized coordinates (the key computational object):
  T = ½ Σ_{ij} M_ij(q, t) q̇_i q̇_j + Σ_i a_i(q,t) q̇_i + a0(q,t)
  M_ij = ∂²T/∂q̇_i ∂q̇_j            (the mass matrix)
EOM in matrix form:
  M(q) q̈ = Q(q, q̇, t) - C(q, q̇) q̇ - ∂U/∂q       →  solve a small dense linear system each step
```

### 7.3 Canonical examples and problem archetypes
- **Single particle in 2D polar coordinates** — recovers the `r̈ - rφ̇²` result of §1.7 in three lines.
- **Atwood's machine** (and the double / "swinging" Atwood machine in problems) — `a = (m1-m2)g/(m1+m2)`.
- **Simple and double pendulum**; pendulum with a **moving support** (support driven horizontally, or attached to a sliding cart).
- **A block sliding on a movable wedge** — two coupled coordinates.
- **A bead on a rotating wire / spinning hoop** — the classic: effective potential, equilibria at `cos θ = g/(ω²R)`, and a pitchfork bifurcation at the critical spin `ω_c = sqrt(g/R)`.
- **Particle confined to a cone, paraboloid, cylinder, or sphere** (spherical pendulum).
- **Two masses connected by a spring** — normal modes (forward reference to Ch. 11).
- **Cyclic coordinates**: central force (`φ` ignorable → `ℓ` conserved), free particle (`x` ignorable → `p_x` conserved).

### 7.4 Why this is out of scope for a Newtonian engine
Chapter 7 does not add physics — it re-derives the *same* equations of motion the force engine already solves, using a different route that trades vector force bookkeeping for symbolic differentiation of a scalar. A force-based engine cannot use `L` directly. The two genuine deltas are (i) **constraints handled implicitly** and (ii) **derivatives that are symbolic**, and both require machinery (symbolic algebra or autodiff, and a linear solver for the mass matrix) that has nothing in common with the Ch. 1–5 core.

There is, however, a *Newtonian* way to reach Chapter 7's results, and it is the better first step (see below).

### 7.5 What a future hook could look like

**Option A (recommended, and squarely Newtonian): constraint solver in the force engine.**
Add holonomic constraints to `nm::core` and solve for the constraint forces with Lagrange multipliers *in Cartesian coordinates*. This gets pendulums, beads on wires, Atwood machines, blocks on wedges, and the double pendulum without leaving the force paradigm:
```cpp
namespace nm::constraints {
  struct Holonomic {
    std::function<double(const SystemState&)> g;      // g(r) = 0
    std::function<VecN(const SystemState&)>  gradG;   // ∂g/∂r  (autodiff or FD)
  };
  // Solve  [ M   Jᵀ ] [ a ]   [ F_applied ]
  //        [ J   0  ] [ λ ] = [ -J̇ v - stabilization ]
  // with Baumgarte stabilization (α, β) or post-step projection to control drift.
  Vec constraintForces(...);   // λ Jᵀ  — directly answers "what is the tension/normal force?"
}
```
This reproduces §7.10's constraint forces, which the pure-Lagrangian route deliberately hides — a genuine pedagogical win, and it reuses the existing integrators.

**Option B: a true `nm::lagrangian` module.**
```cpp
namespace nm::lagrangian {

  using Scalar = double;                    // templated: Dual<double> enables autodiff
  template <int N>
  struct System {
    std::function<Scalar(const VecN<N>& q, const VecN<N>& qd, Scalar t)> T;  // kinetic
    std::function<Scalar(const VecN<N>& q, Scalar t)> U;                     // potential
    Scalar L(q, qd, t) const { return T(q,qd,t) - U(q,t); }
  };

  // Build the EOM numerically — no symbolic algebra needed:
  //   M_ij = ∂²L/∂q̇_i∂q̇_j
  //   b_i  = ∂L/∂q_i - ∂²L/∂q̇_i∂q_j q̇_j - ∂²L/∂q̇_i∂t
  //   solve M q̈ = b   (small dense LDLᵀ / LU; N is 1–4 for every textbook problem)
  template <int N> VecN<N> accelerations(const System<N>&, q, qd, t);

  // Then reuse nm::integrate on the 2N-dimensional state (q, q̇).
  template <int N> Trajectory integrate(const System<N>&, State2N ic, ...);

  VecN<N> generalizedMomenta(...);         // p_i = ∂L/∂q̇_i
  Scalar  hamiltonian(...);                // H = Σ p_i q̇_i - L
  std::vector<int> detectCyclicCoordinates(...);  // numerically: |∂L/∂q_i| ≈ 0
}
```
The enabling technology is **forward-mode automatic differentiation via dual numbers** — roughly 150 lines of header-only C++, with second derivatives from nested duals (`Dual<Dual<double>>`). This is why the Ch. 1 core should be templated on the scalar type from day one: it is nearly free now and expensive to retrofit. (Finite differences also work and are simpler, but lose 4–6 digits on the second derivatives that form `M`.)

**Validation targets for a future Lagrangian module** (all with Newtonian or analytic counterparts):
```
2D polar particle:      must reproduce  m(r̈ - rφ̇²) = F_r  identically to the Cartesian engine
Atwood machine:         a = (m1 - m2) g / (m1 + m2)
Simple pendulum:        θ̈ = -(g/L) sinθ ; exact period 4 sqrt(L/g) K(sin(θ0/2))
Bead on spinning hoop:  equilibria at θ=0, π and cosθ = g/(ω²R) for ω > ω_c = sqrt(g/R)
Block on movable wedge: a_block, a_wedge closed forms
Central force:          φ ignorable ⇒ p_φ = m r² φ̇ = const to machine precision
Double pendulum:        no closed form → test H conservation and cross-check against
                        the Cartesian constrained solver of Option A (two independent
                        formalisms agreeing is the strongest available check)
Any autonomous L:       H(t) = H(0)
```
The **double pendulum solved two ways** — Option A's Cartesian constraint solver vs Option B's generalized coordinates — is the single most convincing demonstration a student can build that both formalisms are the same physics.

---

## 8. The analytic validation suite (Ch. 1–5)

These are the closed-form solutions that should back the library's regression tests. Ordered roughly by increasing difficulty. Every one belongs in `nm::analytic` as a callable reference, not just as a hard-coded expected value.

| # | Ch. | System | Closed form | What it tests |
|---|---|---|---|---|
| 1 | 1 | Constant force / free fall | `r = r0 + v0 t + ½(F/m)t²` | RK4 & Verlet must be **exact** (roundoff only); Euler must not |
| 2 | 1 | Vacuum projectile | `R = v0² sin2θ/g`, `h = v0² sin²θ/2g`, `T = 2v0 sinθ/g` | Event detection (ground impact) accuracy; θ_max = 45° |
| 3 | 1 | Uniform circular motion | `a = ω²r`, orbit closes after `T = 2π/ω` | Long-run geometric fidelity; radius drift |
| 4 | 2 | Linear drag, horizontal | `v = v0 e^{-t/τ}`, `x = v0τ(1-e^{-t/τ})`, `x_∞ = v0τ` | Exponential decay over many τ; finite total distance |
| 5 | 2 | Linear drag, vertical | `v_y = -v_ter + (v_y0+v_ter)e^{-t/τ}`; `v_ter = mg/b` | Terminal-velocity approach; asymptotic accuracy |
| 6 | 2 | **Linear drag, full 2D trajectory** | `y(x) = ((v_y0+v_ter)/v_x0)x + v_ter τ ln(1 - x/(v_x0τ))` | The best 2D drag oracle in the book |
| 7 | 2 | Linear-drag range | transcendental; `R ≈ R_vac[1 - (4/3)(v_y0/v_ter)]` for small drag | Root-finding on `y(T)=0`; small-parameter expansion |
| 8 | 2 | Quadratic drag, horizontal | `v = v0/(1 + c v0 t/m)`, `x = (m/c)ln(1 + c v0 t/m)` | Power-law (non-exponential) decay |
| 9 | 2 | Quadratic drag, vertical from rest | `v = v_ter tanh(gt/v_ter)`, `y = (v_ter²/g)ln cosh(gt/v_ter)`, `v_ter = sqrt(mg/c)` | Nonlinear 1D; hyperbolic asymptote |
| 10 | 2 | Charge in uniform **B** | helix; `r = mv⊥/qB`, `ω = qB/m`, `T = 2πm/qB`, **`\|v\|` exactly constant** | Boris vs RK4 vs Euler; speed drift; orbit closure |
| 11 | 2 | Crossed **E**,**B** | drift `v_d = (E×B)/B²`; cycloidal path | Superposed velocity-dependent forces |
| 12 | 3 | Rocket, free space | `v = v0 + v_ex ln(m0/m)` | Variable-mass EOM; mass state variable |
| 13 | 3 | Rocket in uniform gravity | `v(t) = v_ex ln(m0/(m0-\|ṁ\|t)) - gt` | Combined thrust + external force |
| 14 | 3 | Isolated N-body system | `R_cm(t) = R0 + V0 t` **exactly**; `P` constant | Newton III in the pair-force implementation |
| 15 | 3 | Exploding projectile | CM continues on the unperturbed parabola | Impulse/collision handling preserves `P` |
| 16 | 3 | 1D elastic collision | `v1' = ((m1-m2)v1 + 2m2v2)/(m1+m2)` | Collision resolver; `e = 1` conserves `E` |
| 17 | 3 | Any central force | `ℓ = m r²φ̇ = const`; `dA/dt = ℓ/2m = const` | Kepler's 2nd law; angular-momentum drift |
| 18 | 4 | Any conservative force | `T + U = const`; `∇×F = 0` numerically | Energy drift; symplectic vs RK4 over 10⁴ periods |
| 19 | 4 | Linear/quadratic drag | `T + U + ∫f_nc·v dt = const` | Dissipation bookkeeping |
| 20 | 4 | Escape from `-GMm/r` | `v_esc = sqrt(2GM/R)`; `E = 0` ⇒ `v→0` as `r→∞` | Long-range, slowly varying force; adaptive stepping |
| 21 | 4 | Morse / Lennard-Jones | `ω = sqrt(U''(r0)/m)`; LJ: `U''(r0) = 72U0/r0²` | Small-oscillation limit of a nonlinear potential |
| 22 | 4 | Bounded 1D motion | quadrature period `T = 2∫dx/sqrt(2(E-U)/m)` | Time-domain period must match the quadrature |
| 23 | 5 | **SHM** | `x = A cos(ω0 t - δ)`, `A = sqrt(x0² + (v0/ω0)²)`, `T = 2π/ω0` amplitude-independent, `E = ½kA²` | The canonical integrator benchmark |
| 24 | 5 | Underdamped | `x = A e^{-βt}cos(ω1 t - δ)`, `ω1 = sqrt(ω0²-β²)`; envelope `±Ae^{-βt}`; `E ≈ E0 e^{-2βt}`; log decrement `2πβ/ω1` | Velocity-dependent force handling |
| 25 | 5 | Critically damped | `x = (C1 + C2 t)e^{-βt}`; **fastest** return to equilibrium | Repeated-root regime; assert it beats neighbors |
| 26 | 5 | Overdamped | `x = C1 e^{-(β-s)t} + C2 e^{-(β+s)t}`, `s = sqrt(β²-ω0²)` | Stiffness; step-size stability bound |
| 27 | 5 | Driven steady state | `A(ω) = f0/sqrt((ω0²-ω²)² + 4β²ω²)`, `tanδ = 2βω/(ω0²-ω²)` | Frequency sweep vs analytic Lorentzian |
| 28 | 5 | Resonance | peak at `ω = sqrt(ω0²-2β²)`; FWHM `= 2β`; `Q = ω0/2β`; `A(ω0) = f0/(2βω0)` | Q extracted from sweep **and** ring-down must agree |
| 29 | 5 | Undamped exact resonance | `x = (f0/2ω0) t sin(ω0 t)` (linear growth) | Secular/degenerate forcing case |
| 30 | 5 | Undamped near resonance | beats at `\|ω-ω0\|/2`; `x = (2f0/(ω0²-ω²))sin(...)sin(...)` | Two-timescale accuracy |
| 31 | 5 | Anisotropic 2D oscillator | Lissajous; closed iff `ωx:ωy ∈ ℚ` | Multi-DOF decoupled integration; long-run phase drift |
| 32 | 5 | Fourier-driven oscillator | `x = Σ_n A_n cos(nωt - δ_n)`; Parseval `<x²> = ½ΣA_n²` | Superposition + spectral post-processing |
| 33 | 5 | **Full nonlinear pendulum** | `T = 4 sqrt(L/g) K(sin(θ0/2))`, `≈ T0[1 + θ0²/16 + 11θ0⁴/3072]` | The best genuine-nonlinearity test (needs an elliptic-K implementation, ~20 lines by AGM) |

### 8.1 Numerical-quality tests that need no analytic solution
These complement the table and matter just as much:
- **Order-of-convergence study.** Halve `dt` repeatedly on a problem with a known solution; fit `log(error)` vs `log(dt)`; assert the slope equals each integrator's advertised order (1 for Euler, 2 for Verlet/RK2, 4 for RK4/Yoshida4). This is the definitive integrator-correctness test and catches subtle coefficient typos that a single-`dt` accuracy test will not.
- **Time-reversibility.** Integrate forward `N` steps, negate all velocities, integrate `N` more; velocity Verlet returns to the initial state to roundoff, RK4 does not.
- **Long-time energy behavior.** Over 10⁵ periods of an SHO: Verlet's energy error stays bounded and oscillatory; RK4's drifts secularly. Ship this as a figure.
- **Adaptive-stepper tolerance honoring.** Tighten `rtol` by 10×, confirm the true error drops accordingly.
- **Cross-integrator agreement** for problems with no oracle (2D quadratic drag, double pendulum): RK45 at `rtol=1e-12` serves as the reference solution.
- **Invariant drift as a first-class output** — every `System` run should be able to emit `ΔP/P`, `ΔL/L`, `ΔE/E` vs time without extra user code.

---

## 9. Suggested build order, aligned to the course

| Phase | Weeks (typical Mech-1 pace) | Deliverable | Backed by |
|---|---|---|---|
| 1 | Ch. 1 | `Vec3`, `State`, `Particle`, `System`, Euler + RK4, uniform gravity, CSV out | Tests 1–3 |
| 2 | Ch. 2 | Drag models, Lorentz force, RK45 adaptive, **event detection with root-finding**, Boris pusher, range/apex analysis | Tests 4–11 — this is where the library becomes genuinely useful |
| 3 | Ch. 3 | Multi-particle systems, CM/`P`/`L` monitors, pair forces, variable-mass rocket, collisions | Tests 12–17 |
| 4 | Ch. 4 | `Potential` abstraction with numerical gradients, curl test, energy accounting, **velocity Verlet / Yoshida4**, 1D potential toolkit, phase portraits | Tests 18–22, plus §8.1 |
| 5 | Ch. 5 | Springs/dampers/drives, `nm::analytic::DampedDrivenOscillator`, frequency sweep, peak finding, Q extraction, small FFT | Tests 23–33 |
| 6 | Ch. 6–7 *(optional, after the course)* | Cartesian holonomic constraint solver (Option A) → dual-number autodiff → `nm::lagrangian` (Option B) | §7.5 targets |

Two cross-cutting recommendations:
1. **Template the scalar type from the very first commit** (`template<class T> struct Vec3T`). It costs almost nothing in Phase 1 and is what makes dual-number autodiff — and therefore any future Lagrangian or variational module — possible without a rewrite.
2. **Build `nm::analytic` in parallel with `nm::integrate`, not after.** Every chapter above supplies closed forms; a library whose test suite *is* the textbook's analytic results is both more correct and a better study aid than one tested against hard-coded golden numbers.

---

**Sources consulted for table-of-contents verification:**
- [Classical Mechanics, John R. Taylor — chapter listing (JHU)](https://henry.pha.jhu.edu/chapters.html)
- [Topics covered in John R. Taylor, Classical Mechanics — Physics Forums](https://www.physicsforums.com/threads/topics-covered-in-john-r-taylor-classical-mechanics.1015957/)
- [Taylor Classical Mechanics Ch. 2 problem solutions (StemJock)](https://stemjock.com/STEM%20Books/Taylor%20Classical%20Mechanics/Chapter%202/TaylorCMCh2p24.pdf)