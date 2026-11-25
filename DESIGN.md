# newton — Design Specification (v0.1)

This is the authoritative API spec for the `newton` library. Implementation
must match the signatures, file paths, and semantics here exactly — the
umbrella header `include/newton/newton.hpp` already lists every public header
and is the manifest. The research reports in `docs/research/` are the source
of truth for formulas and numerical methods:

- `docs/research/taylor-ch1-7-map.md` — physics: every formula, chapter by chapter
- `docs/research/numerical-integration.md` — integrator equations and pitfalls
- `docs/research/cpp-library-design.md` — C++ conventions, build tooling

## What this library is

A **header-only C++17** library for simulating single-particle Newtonian
mechanics, built alongside Taylor's *Classical Mechanics* ch. 1–5 by a student
taking Mechanics 1. It is a *learning library*: correctness and readability
beat generality and performance. Every public function carries a doc comment
with two labelled sections — `Physics:` and `C++ note:` — as demonstrated in
the core headers (`include/newton/core/*.hpp`). Match that style and density.

## Hard rules (apply to every file)

1. C++17. Compiles clean with `clang++ -std=c++17 -Wall -Wextra -Wpedantic`.
2. Header-only: every namespace-scope function is `inline`, a template, or
   `constexpr`; every namespace-scope variable is `inline constexpr`.
3. `#pragma once`, `.hpp` extension, angle-bracket internal includes
   (`#include <newton/core/vec3.hpp>`), each header includes what it uses.
4. Namespace per module: `newton` (core types), `newton::constants`,
   `newton::forces`, `newton::integrators`, `newton::analytic`,
   `newton::analysis`, `newton::systems`, `newton::io`, `newton::coords`,
   `newton::fields`, `newton::potentials`.
   Never `using namespace` at namespace scope in a header.
5. SI units, radians, `double`. Every function parameter documents its unit:
   `/// @param b linear drag coefficient [kg/s]`.
6. `[[nodiscard]]` on every pure function. No `noexcept` on anything calling
   user callables.
7. No dependencies beyond the C++ standard library. No `std::function` in any
   hot loop — force parameters are template callables.
8. Close namespaces with `}  // namespace newton::xyz`.

## The Force concept

A force is any callable with this shape:

```cpp
Vec3 operator()(double t, const Vec3& r, const Vec3& v) const;
```

Arguments are time [s], position [m], velocity [m/s]; the return is the force
[N]. All integrators and drivers take `Force&&` as a template parameter.
Lambdas, the functors in `newton::forces`, and `std::function` all satisfy it.

---

## File-by-file specification

### Already written (do not modify): core

- `include/newton/core/vec3.hpp` — `Vec3` aggregate; `+ - * /`, `dot`,
  `cross`, `norm`, `norm2`, `normalized` (zero-guarded).
- `include/newton/core/particle.hpp` — `State {Vec3 r; Vec3 v; double t;}`,
  `Particle {double mass; Vec3 r; Vec3 v;}` with `initial_state(t0)`,
  free `momentum(particle)`.
- `include/newton/core/constants.hpp` — `pi`, `deg`, `g_earth`, `G`,
  `sphere_air_beta`, `sphere_air_gamma` in `newton::constants`.
- `include/newton/core/trajectory.hpp` — `Trajectory {double mass;
  std::vector<State> states;}` with `size/empty/front/back`.
- `include/newton/newton.hpp` — umbrella/manifest.

### Agent A1 — forces (`newton::forces`)

`include/newton/forces/forces.hpp` is a tiny module umbrella that only
includes the six files below.

`include/newton/forces/gravity.hpp`
```cpp
/// F = m g, uniform field. Default g points in -y (2D convention: y is up).
struct UniformGravity {
    double mass{1.0};                                  // [kg]
    Vec3   g{0.0, -constants::g_earth, 0.0};           // [m/s^2]
    Vec3 operator()(double t, const Vec3& r, const Vec3& v) const;  // = mass * g
};

/// F = -G M m r̂ / |r - center|^2, Newtonian point-mass gravity (Taylor §4.8).
/// Optional Plummer softening eps [m] tames the r -> 0 singularity:
/// denominator (|d|^2 + eps^2)^(3/2). Softening changes the physics at short
/// range — say so in the doc comment.
struct PointGravity {
    double GM{1.0};        // G * (central mass) [m^3/s^2]
    double mass{1.0};      // orbiting particle's mass [kg]
    Vec3   center{};       // position of the attracting body [m]
    double softening{0.0}; // Plummer eps [m]
};
```

Every non-conservative force's doc comment (drag, magnetic, drives) says which of Taylor's §4.2 conditions it fails.

`include/newton/forces/drag.hpp` (Taylor §2.1–2.4)
```cpp
/// F = -b v (linear / viscous drag).
struct LinearDrag  { double b{0.0}; };   // [kg/s]

/// F = -c |v| v (quadratic drag; singularity-free algebraic form, NOT -c|v|² v̂).
struct QuadraticDrag { double c{0.0}; }; // [kg/m]

/// Factory helpers using Taylor's sphere-in-air coefficients:
/// linear_drag_sphere_air(D)    -> LinearDrag{ sphere_air_beta * D }
/// quadratic_drag_sphere_air(D) -> QuadraticDrag{ sphere_air_gamma * D * D }
```

`include/newton/forces/spring.hpp` (Taylor §5.1)
```cpp
/// Hooke's law toward an anchor: F = -k (|d| - natural_length) d̂ with
/// d = r - anchor. With natural_length = 0 this reduces to F = -k d, the
/// isotropic oscillator. Zero-guard d̂ via normalized().
struct Spring {
    double k{1.0};              // [N/m]
    Vec3   anchor{};            // [m]
    double natural_length{0.0}; // [m]
};

/// Anisotropic spring for Lissajous figures (Taylor §5.3):
/// F = (-kx*x, -ky*y, -kz*z) about the origin.
struct AnisotropicSpring { double kx{1.0}, ky{1.0}, kz{1.0}; };
```

`include/newton/forces/drives.hpp` (Taylor §5.5)
```cpp
/// Sinusoidal driving force F(t) = F0 cos(omega t + phase), direction dir
/// (unit vector, default {1,0,0}).
struct SinusoidalDrive {
    double F0{1.0};      // amplitude [N]
    double omega{1.0};   // [rad/s]
    double phase{0.0};   // [rad]
    Vec3   dir{1.0, 0.0, 0.0};
};

/// Constant force, e.g. a constant thrust or a uniform E-field seen by a charge.
struct ConstantForce { Vec3 F{}; };
```

`include/newton/forces/electromagnetic.hpp` (Taylor §2.5, §2.7)
```cpp
/// Magnetic force F = q v × B. Does no work: |v| is exactly constant.
struct MagneticForce { double q{1.0}; Vec3 B{0.0, 0.0, 1.0}; };

/// Full Lorentz force F = q (E + v × B).
struct LorentzForce  { double q{1.0}; Vec3 E{}; Vec3 B{}; };
```

`include/newton/forces/compose.hpp`
```cpp
/// sum(f1, f2, ...) -> a callable returning the vector sum of all forces.
/// Variadic template over callables; C++17 fold or recursion. This is how a
/// projectile gets gravity + drag: sum(UniformGravity{m}, QuadraticDrag{c}).
template <class... Fs> [[nodiscard]] auto sum(Fs... fs);
```

### Agent A2 — integrators + driver

Each step function has the identical shape (t for force evaluation comes from
`s.t`; the returned State carries `t = s.t + dt`):

```cpp
template <class Force>
[[nodiscard]] State xxx_step(const State& s, double mass, double dt, Force&& force);
```

All update equations must follow `docs/research/numerical-integration.md` §2
exactly. Files:

- `include/newton/integrators/euler.hpp` — explicit Euler. Document loudly:
  first-order, energy grows as (1+ω²dt²)ⁿ on oscillators, ships as a
  pedagogical foil, never as a default.
- `include/newton/integrators/semi_implicit_euler.hpp` — Euler–Cromer:
  v first, then x with the NEW v. First-order but symplectic for F(x).
- `include/newton/integrators/velocity_verlet.hpp` — kick-drift-kick with the
  predictor–corrector velocity handling (research §2.3 resolution (2)) so
  velocity-dependent forces (drag, damping) stay second-order:
  `v_pred = v_half + 0.5*dt*a_n` evaluated force at (t+dt, x_{n+1}, v_pred),
  then `v_{n+1} = v_half + 0.5*dt*a_{n+1}`. Document: this standalone step
  evaluates the force twice per step for clarity; a cached-acceleration
  variant is a v0.2 optimization.
- `include/newton/integrators/rk4.hpp` — classic RK4 on the coupled (x, v)
  system with the k1..k4 stages written out exactly as research §2.4 (compute
  k_ix before k_iv; the position update may use the collapsed form
  `x + dt*v + dt²/6*(k1v+k2v+k3v)` with a comment deriving it).

`include/newton/simulate.hpp` — the driver:

```cpp
namespace newton {

enum class Method { Euler, SemiImplicitEuler, VelocityVerlet, RK4 };

struct SimOptions {
    double t0{0.0};                  // [s]
    double t_end{1.0};               // [s]
    double dt{1e-3};                 // [s] — required, never defaulted sensibly; doc says so
    Method method{Method::RK4};
    int    record_every{1};          // record every Nth step (first & last always)
};

/// Fixed-step time loop. CRITICAL numerics detail (research §7.2): time is
/// computed as t = t0 + i*dt each iteration — never accumulated t += dt —
/// and the step count is n = llround((t_end - t0)/dt). After each step the
/// state's t is overwritten with the exact t0 + (i+1)*dt.
template <class Force>
[[nodiscard]] Trajectory simulate(const Particle& p, Force&& force,
                                  const SimOptions& opt);

/// Event-terminated run. `event` is a scalar function of State whose sign
/// change stops the run: double event(const State&). On sign change between
/// steps, the crossing is refined by bisection ON THE STEP SIZE dt' ∈ (0, dt]
/// re-stepping from the previous state (~40 iterations), so the final state
/// lands on event ≈ 0 to integrator accuracy, not O(dt). (Naive "stop at
/// first sign change" gives O(dt) range errors — research taylor map §2.4.)
/// The refined final state is appended and the trajectory returned.
/// If the event never fires, runs to t_end like simulate().
template <class Force, class EventFn>
[[nodiscard]] Trajectory simulate_until(const Particle& p, Force&& force,
                                        const SimOptions& opt, EventFn&& event);

namespace events {
/// event functions for simulate_until
[[nodiscard]] inline auto ground_impact(double y_level = 0.0);  // returns s.r.y - y_level
[[nodiscard]] inline auto apex();                               // returns s.v.y
}  // namespace events

}  // namespace newton
```

Method dispatch inside the loop is a `switch` calling the step functions.

### Agent A3 — analytic (`newton::analytic`)

Closed-form reference solutions; formulas from the taylor map (§ per chapter,
consolidated table §8). Every function `[[nodiscard]] inline` (or constexpr
where possible), tolerating t in seconds. Small structs group results.

`include/newton/analytic/projectile.hpp` — vacuum projectile (Taylor §1.6):
`range(v0, theta, g)`, `apex_height(...)`, `time_of_flight(...)`,
`position(t, v0, theta, g) -> Vec3` (x horizontal, y vertical).

`include/newton/analytic/drag.hpp` — Taylor §2.2–2.4:
```cpp
/// Linear drag. tau = m/b, v_ter = m g / b.
struct LinearDragParams { double m{1.0}; double b{0.0}; double g{constants::g_earth};
                          double tau() const; double v_ter() const; };   // b = 0 is the vacuum limit
Vec3 linear_drag_velocity(const LinearDragParams&, Vec3 v0, double t);
Vec3 linear_drag_position(const LinearDragParams&, Vec3 r0, Vec3 v0, double t);
double linear_drag_trajectory_y(const LinearDragParams&, double vx0, double vy0, double x); // the exact y(x)
double linear_drag_range(const LinearDragParams&, double vx0, double vy0);             // level-ground range, y(R) = 0, by bisection
double linear_drag_range_small_drag(const LinearDragParams&, double vx0, double vy0);  // Taylor's R_vac (1 - 4η/3); next term +14/9 η²
/// Quadratic drag, 1D closed forms:
double quadratic_drag_terminal_speed(double m, double c, double g);          // sqrt(mg/c)
double quadratic_drag_fall_speed(double v_ter, double g, double t);          // v_ter tanh(gt/v_ter)
double quadratic_drag_fall_distance(double v_ter, double g, double t);       // (v_ter²/g) ln cosh(gt/v_ter)
double quadratic_drag_horizontal_speed(double m, double c, double v0, double t); // v0/(1 + c v0 t/m)
double quadratic_drag_horizontal_distance(double m, double c, double v0, double t); // (m/c) ln(1 + c v0 t/m)
```
Numerics convention: the textbook forms multiply τ = m/b (huge for small drag) by 1 − e^(−t/τ) (tiny), a
catastrophic cancellation (b = 1e-10 once put a projectile 24 km high). Every function is therefore written with
the small dimensionless ratios u = b t/m, ε = b x/(m vx0), z = g t/v_ter and Taylor's own series in them
(`detail::phi1`, `phi2`, `psi`, `ln_cosh`, `lambda`, `tanh_over_z`, `ln1p_over`), so each is exact at b = 0 (c = 0)
and none divides by b; the quadratic fall functions do not overflow for large t.

`include/newton/analytic/oscillator.hpp` — Taylor ch. 5. The centerpiece:
```cpp
/// x(t) = A cos(w0 t - delta) from ICs (Taylor §5.2).
struct SHM {
    double omega0;                       // sqrt(k/m) [rad/s]
    double x0, v0;
    double amplitude() const;            // sqrt(x0² + (v0/w0)²)
    double phase() const;                // atan2(v0/w0, x0)
    double x(double t) const;  double v(double t) const;
    double period() const;  double energy(double m) const;   // ½ m w0² A²
};

/// Damped (+ optionally driven) oscillator x'' + 2β x' + ω0² x = f0 cos(ω t)
/// (Taylor §5.4–5.6). Regime from β vs ω0; x(t) includes the transient with
/// coefficients fitted to (x0, v0); steady-state amplitude/phase closed-form.
enum class Regime { Undamped, Underdamped, CriticallyDamped, Overdamped };
struct DampedOscillator {
    double omega0, beta;                 // [rad/s]
    double x0{1.0}, v0{0.0};
    Regime regime() const;
    double omega1() const;               // sqrt(ω0² - β²) (underdamped)
    double x(double t) const;  double v(double t) const;   // all three regimes
    double log_decrement() const;        // 2πβ/ω1
};
struct DrivenOscillator {                // steady state only
    double omega0, beta, f0, omega;      // f0 = F0/m [m/s²]
    double amplitude() const;            // f0 / sqrt((ω0²-ω²)² + 4β²ω²)
    double phase() const;                // atan2(2βω, ω0²-ω²) ∈ (0, π)
    double resonant_omega() const;       // sqrt(ω0² - 2β²)
    double quality_factor() const;       // ω0/(2β)
    double fwhm() const;                 // 2β (of A²(ω))
};
```

`include/newton/analytic/magnetic.hpp` — Taylor §2.5/2.7: `cyclotron_omega(q,
B, m)`, `cyclotron_radius(m, v_perp, q, B)`, `cyclotron_period(m, q, B)`.

`include/newton/analytic/rocket.hpp` — Taylor §3.2:
`tsiolkovsky_delta_v(v_ex, m0, m)`, `rocket_velocity_gravity(v_ex, m0,
burn_rate, g, t)` = v_ex ln(m0/(m0 - burn_rate·t)) - g t.

### Agent A4 — analysis, io, systems

`include/newton/analysis/energy.hpp`:
```cpp
double kinetic_energy(double mass, const Vec3& v);                 // ½ m |v|²
template <class PotentialFn>   // double U(const Vec3& r)
double total_energy(double mass, const State& s, PotentialFn&& U);
/// Relative energy drift over a trajectory: (E(t) - E(0)) / |E(0)|, returned
/// as vector<double> aligned with trajectory.states. The single most useful
/// integrator diagnostic (research §4.3): bounded oscillation = symplectic
/// health; a monotone trend = drift; growth = wrong dt or explicit Euler.
template <class PotentialFn>
std::vector<double> energy_drift(const Trajectory&, PotentialFn&& U);
```
Common potentials as small inline helpers: `uniform_gravity_U(m, g) -> lambda
m*g*y`, `spring_U(k) -> lambda ½k|r|²`. They match ONE setup each: `uniform_gravity_U` assumes g
along −y, and `spring_U` assumes anchor = origin and natural_length = 0. For any other
setup (gravity along z, an anchored or natural-length spring) write the U, or use `potentials::find_potential`,
otherwise `energy_drift` reads order 1 and the blame lands on `dt`.

`include/newton/analysis/momentum.hpp`:
```cpp
Vec3 angular_momentum(double mass, const State& s, const Vec3& about = {}); // (r-about) × m v
Vec3 center_of_mass(const std::vector<Particle>&);
Vec3 total_momentum(const std::vector<Particle>&);
```

`include/newton/analysis/oscillation.hpp` — measurements on a Trajectory
using a component accessor (default: x): `find_peaks(traj)` (local maxima of
r.x by 3-point comparison) returning times+values; `measure_period(traj)`
(mean spacing of upward zero crossings of r.x, linearly interpolated);
`measure_amplitude(traj, t_start)` (max |r.x| after t_start — for
steady-state measurement after transients die).

`include/newton/io/csv.hpp`:
```cpp
/// Write "t,x,y,z,vx,vy,vz\n" rows. Returns false on failure to open.
/// Precision: %.12g. Plot with scripts/plot_csv.py or any tool.
bool write_csv(const std::string& path, const Trajectory&);
/// Generic columns variant: write_csv(path, {"omega","amplitude"}, rows)
/// with rows = vector<vector<double>>.
```

`include/newton/systems/rocket.hpp` (Taylor §3.2) — self-contained 1D
vertical rocket (variable mass CANNOT reuse the constant-mass Particle path;
research §7.8):
```cpp
namespace newton::systems {
struct RocketParams {
    double m0;                 // initial (wet) mass [kg]
    double m_dry;              // mass at burnout [kg]
    double v_ex;               // exhaust speed [m/s]
    double burn_rate;          // |dm/dt| while burning [kg/s]
    double g{constants::g_earth};
};
struct RocketSample { double t, m, v, y; };
/// Integrates dv/dt = v_ex*burn_rate/m(t) - g, dy/dt = v with RK4 on the
/// scalar pair (y, v); m(t) = m0 - burn_rate*t until m_dry, then thrust = 0.
/// Runs from t=0 to t_end. Includes exact burnout instant as a sample.
[[nodiscard]] inline std::vector<RocketSample>
simulate_rocket(const RocketParams&, double dt, double t_end);
}  // namespace newton::systems
```

### Agent A5 — examples (each: own `main()`, reads like a lesson)

Every example states, in a top-of-file comment: the physical situation, the
Taylor section it comes from, the analytic expectation, and what to look at
in the output. Print comparisons numeric-vs-analytic to stdout; write CSVs to
`build/` (create nothing — assume `build/` exists since binaries live there).

1. `examples/01_free_fall.cpp` — drop + vacuum projectile at 45°; compare
   range/apex/flight time against `analytic::projectile`; demonstrate
   `simulate_until(events::ground_impact())`.
2. `examples/02_projectile_drag.cpp` — Taylor's baseball (D=7cm, m=0.145kg,
   30 m/s at 45°): vacuum vs linear vs quadratic drag on one CSV
   (`build/projectile_drag.csv`); sweep launch angle for the quadratic case,
   report optimal angle (< 45°).
3. `examples/03_charged_particle.cpp` — helix in uniform B; verify |v| const,
   radius & period vs `analytic::magnetic`; crossed E×B drift bonus.
   CSV `build/helix.csv`.
4. `examples/04_rocket.cpp` — vertical launch: `systems::simulate_rocket` vs
   `analytic::rocket_velocity_gravity` during burn; report burnout v & alt.
5. `examples/05_oscillations.cpp` — SHM period & amplitude vs `analytic::SHM`;
   all three damping regimes on one CSV (`build/damping_regimes.csv`);
   driven resonance sweep: amplitude vs ω over [0.2ω0, 2ω0], compare to
   `analytic::DrivenOscillator::amplitude()`, CSV `build/resonance.csv`.
6. `examples/06_integrator_comparison.cpp` — THE money demo: SHM with all
   four methods at dt = T/100 over 200 periods; print final relative energy
   error per method; CSV `build/energy_drift.csv` (t, E_euler, E_semi,
   E_verlet, E_rk4). Expect: Euler grows, semi-implicit & Verlet bounded,
   RK4 tiny monotone decay.

`scripts/plot_csv.py` — small stdlib+matplotlib script:
`python3 scripts/plot_csv.py build/foo.csv [xcol ycol...]`; degrade
gracefully with a message if matplotlib is missing.

### Agent A6 — tests + build tooling

`tests/test_util.hpp` — copy the verified harness from
`docs/research/cpp-library-design.md` §7 (NT_CHECK / NT_CLOSE / report, with
its comment block explaining tolerance-based float comparison).

Test files (each with `main()`, returning `newton_test::report(...)`):
- `tests/test_vec3.cpp` — arithmetic, dot/cross identities (a×b ⊥ a, b),
  norm of {1,2,2} = 3, normalized zero-guard.
- `tests/test_forces.cpp` — each functor returns the hand-computed vector;
  `sum` adds; quadratic drag at v=0 is exactly zero (no NaN).
- `tests/test_integrators.cpp` —
  (a) constant force: Verlet & RK4 exact to ~1e-12 at dt=0.1; Euler NOT exact;
  (b) SHM convergence orders by halving dt (research §4.2: fit observed order,
  assert within ±0.2 of 1, 1, 2, 4);
  (c) energy behavior over 100 periods: Euler grows >10%, semi & Verlet
  bounded < 1%, RK4 |drift| tiny;
  (d) time computed as t0+i*dt: final state t equals t_end to 1e-12.
- `tests/test_analytic.cpp` — internal consistency: SHM x(0)=x0, v(0)=v0;
  damped regimes continuous in t; driven amplitude peaks near resonant_omega;
  vacuum projectile range at 45° = v0²/g.
- `tests/test_numeric_vs_analytic.cpp` — the library's soul:
  linear-drag 2D trajectory vs exact y(x) (rtol 1e-6 with RK4 dt=1e-4);
  quadratic fall vs tanh profile; terminal velocity approach; underdamped
  oscillator vs closed form; driven steady-state amplitude within 0.5%;
  rocket velocity vs Tsiolkovsky+gravity; magnetic |v| conservation < 1e-12
  and cyclotron period; `simulate_until` range vs analytic vacuum range to
  1e-9 (event refinement works).
- `tests/test_analysis.cpp` — energy_drift on exact SHM data ≈ 0; peaks and
  period on synthetic cos data; angular momentum of circular motion constant.

Build files (base them on the verified versions in
`docs/research/cpp-library-design.md` §6):
- `Makefile` — targets: all (examples), tests, check (build+run tests),
  clean, debug (sanitizers), install/uninstall (PREFIX=/usr/local, cp -R);
  auto-discovery via wildcard; `HEADERS := $(shell find include -name '*.hpp')`.
- `build.sh` — the verified bash script (shebang, set -euo pipefail, cd to
  script dir, compile examples + tests, run tests), heavily commented for a
  first-time .sh reader.
- `.gitignore` — `build/`, `.DS_Store`, `*.csv`.

### Agents A7/A8 — notes and README (see agent prompts)

Notes live in `notes/`, written for a student with good basic C++ and decent
advanced C++, currently taking Mechanics 1. Reference actual file paths and
real function names from this spec. README.md follows the structure in
`docs/research/cpp-library-design.md` §8.

### coords (`newton::coords`), added 2026-11-24

Cylindrical and spherical polar coordinates, physics convention: spherical
`theta` is the polar angle down from +z in [0, pi]; `phi` is the azimuth from
+x in [-pi, pi] and is the same angle in both systems. Each coordinate header
includes only `core/vec3.hpp`; `cartesian.hpp` includes both and is the only
coords header allowed to. Functions are named after the type they return.

- `include/newton/coords/cylindrical_polar.hpp` — `CylindricalPolar {double
  rho, phi, z;}`; `to_cylindrical(const Vec3&)`.
- `include/newton/coords/spherical_polar.hpp` — `SphericalPolar {double r,
  theta, phi;}`; `to_spherical(const Vec3&)`, with theta from
  `atan2(rho, z)` rather than `acos(z / r)`.
- `include/newton/coords/cartesian.hpp` — `to_cartesian` for both types, plus
  the direct `to_cylindrical(const SphericalPolar&)` and
  `to_spherical(const CylindricalPolar&)`.
- `tests/test_coords.cpp` — hand-computed points, the atan2 quadrant check
  (problem 1.42), round trips in both directions, direct vs via-Cartesian
  agreement, and the acos precision loss near the poles.

### fields (`newton::fields`), added 2026-11-24

Numerical calculus on fields given as plain callables (Taylor §4.1–4.4). A
*scalar field* is any callable `double U(const Vec3& r)`; a *vector field* is
any callable `Vec3 F(const Vec3& r)`. Fields depend on position only: a
library Force `F(t, r, v)` becomes a field by freezing `t` and `v` in a lambda,
`[&](const Vec3& r) { return force(t, r, Vec3{}); }`. Both headers include
only `core/vec3.hpp` (plus `core/constants.hpp` for `pi`) and are independent
of each other. Public templates take the callable as a forwarding reference
and `static_assert` its shape with `std::is_invocable_r_v`, so a wrong
callable gets one readable sentence instead of a page of template errors.

- `include/newton/fields/derivatives.hpp` — 5-point (4th-order) central
  differences, documented as Richardson extrapolation of the 2-point formula.
  `default_step(r, length_scale = 1.0)` = `1e-4 * max(|r|, length_scale)` [m];
  `Jacobian {Vec3 dx, dy, dz}` stores the *columns* ∂F/∂x, ∂F/∂y, ∂F/∂z, so
  `J.dx.y` is ∂Fy/∂x; `constexpr` `curl(const Jacobian&)`,
  `divergence(const Jacobian&)`, and `frobenius_norm`; and the numerical
  `gradient`, `jacobian`, `curl`, `divergence`, each with an explicit-step
  overload `(…, r, h)` and a default-step overload `(…, r)`. Overloads rather
  than an "`h <= 0` means automatic" sentinel. `detail::` holds the stencil,
  and the classifier's per-point sample: `jacobian_with_error` (13 evaluations)
  and `relative_curl`. Convention: the step follows `|r|` but is floored by
  `length_scale`, which defaults to 1 m. At atomic scale that floor silently
  makes the gradient exactly zero, so pass `default_step(r, 1e-10)`.
- `include/newton/fields/line_integral.hpp` — global adaptive Gauss–Kronrod
  (G7–K15, QUADPACK's QAG strategy). `LineIntegral {value, error, scale,
  evaluations, converged}`; `line_integral_detailed(F, a, b, rel_tol = 1e-12,
  scale_tol = 1e-14)`; `line_integral(F, a, b)` and `line_integral(F, path)`
  (a `std::vector<Vec3>` polyline); `circulation_detailed(F, center, normal,
  radius, ...)` and `circulation(F, center, normal, radius)`. The plain-`double`
  overloads return NaN if the integral did not converge, so a failure cannot
  pass for a number. Convention: a circle runs counterclockwise seen from the
  tip of `normal` (right-hand rule). Limits: at most 200 pieces
  (5985 evaluations), and any non-finite `F` on the path means
  `converged = false`.
- `tests/test_fields.cpp` and `tests/test_line_integrals.cpp` — gradient,
  Jacobian, curl and divergence against closed forms (including atomic-scale
  fields and the F = x r counterexample of §4.8); segments, polylines and
  circulation against exact work, the hole field's 2π, and the failure flags.

### potentials (`newton::potentials`), added 2026-11-24

"Is this force conservative, and what is its potential energy?" and the
reverse, "what force does this potential energy give?" (Taylor §4.2–4.5).
Everything here accepts a force in any of three shapes, detected at compile
time with `if constexpr` + `std::is_invocable_r_v`, checked in this order:
`Vec3 F(double t, const Vec3& r, const Vec3& v)` (the library's Force),
`Vec3 F(double t, const Vec3& r)` (a time-dependent field), and
`Vec3 F(const Vec3& r)` (a static field). Anything else is a `static_assert`.
Depends on `newton::fields`.

- `include/newton/potentials/conservative.hpp` — `Box {lo, hi}` (the region
  tested; it must have positive size in all three directions); `CheckOptions`
  (region, time window, `v_max`, sample counts, `loop_grid_levels`, tolerances,
  a fixed RNG seed); `Verdict` (`Conservative`, `TimeDependentPotential`,
  `VelocityDependent`, `NonzeroCurl`, `PathDependent`, `Inconclusive`);
  `VelocityPower` (`NotApplicable`, `NoWork`, `Dissipative`, `Feeding`, `Mixed`,
  `Unknown`); `ConservativeReport` (the verdict plus the evidence: worst curl
  and where/when, worst loop ratio and `worst_loop_center`/`worst_loop_time`,
  sample and loop counts, and `note`, a string literal saying why a verdict is
  `Inconclusive`), with `has_potential()` and `is_conservative()`;
  `probe_dependence(F, opt)`; `check_conservative(F, opt)`; `to_string(Verdict)`,
  `to_string(VelocityPower)`, `describe(report)`; and the split helpers
  `velocity_independent_part(F)` and `velocity_dependent_part(F)`.
  Check order: (0) a degenerate box is `Inconclusive` with a note; (1) probe `v`,
  then `t`; (2) curl at 256 sampled points per frozen time slice, with the
  finite-difference step measured from the BOX CENTRE
  (`h = 1e-4 max(|r − centre|, L)`); a candidate counts only if it survives
  recomputation at `h/2` AND a small Stokes loop (circulation around a circle
  of radius 0.05 L, normal along the curl, must equal |curl|·πR² to 25%);
  (3) closed loops: in each of 3 parallel planes per axis (the box centre and
  centre ± a quarter of the side), g×g grids of rectangular cell loops
  with g = 1, 2, 3, 5 (pairwise coprime) tile the cross-section, 3 × 3 × 39 =
  351 loops per slice, each four straight-segment integrals; (4) the verdict.
  The verdict is evidence, not proof: it applies to the box and time window
  given, can show a force is *not* conservative, and can only gather evidence
  that it is. Convention: time slices use golden-ratio offsets, never just
  `t = 0` (the induced swirl `sin(t)(-y, x, 0)` has zero curl there).
  `VelocityPower` describes the velocity-dependent *part* F(v) − F(0): `NoWork`
  (magnetic), `Dissipative` (drag, friction), `Feeding` (negative damping),
  `Mixed` (both signs genuinely occur), `Unknown` (no usable sample, e.g.
  friction written −μ v/|v| is NaN at rest). The split helpers are exact
  bookkeeping, `F = F(t, r, 0) + [F − F(t, r, 0)]`, not Taylor's physical
  classification (drag in a wind is the caveat: the label is `Dissipative`
  even though the true drag power can be positive). Remaining documented limit:
  a singular line that only clips the box near a face or corner can be missed
  (measured 4.8% of random lines over the whole box, 0% for lines through the
  inner 80%).
- `include/newton/potentials/potential_from_force.hpp` —
  `LineIntegralPotential<Force> {force, r_ref, U_ref, region}`, a function
  object with `U(t, r)` and `U(r)` (t = 0), computing `U_ref − ∫ F(t, r', 0)·dr'`
  from `r_ref` with `t` frozen; the optional `region` (a `std::optional<Box>`,
  set by `find_potential`) clamps the integration detours into the tested box;
  `potential_from_force(F, r_ref, U_ref)`;
  `PotentialResult<Force> {report, U, r_ref_used, reference_ok, reference_moved,
  note}`; `find_potential(F, opt, r_ref, U_ref)`, which classifies and, if a
  potential exists, validates the reference point; and `describe(PotentialResult)`.
  Validation checks PATH INDEPENDENCE: from the reference to each of the 8 box
  corners, at `t_min` and every time slice, the work along the recovered straight
  path and along an independent generic detour must both converge and agree to
  1e-9. Conventions: the default reference is the origin when it is in the box
  (so `mgy`, `½k|r|²` and `-F·r` come out in textbook form), otherwise the box
  centre; a default reference that fails validation is moved; a caller's `r_ref`
  is never moved, and one outside the box gives `reference_ok = false` with a
  `note`; an integral that cannot be done returns NaN, never a wrong number. The
  U(∞) = 0 recipe: pass a regular `r_ref` and `U_ref = U_textbook(r_ref)`, which
  makes U match the textbook form *inside the checked box* (for softened gravity,
  `-GMm/sqrt(|r_ref − c|² + ε²)`). For a time-dependent force
  `U(r_ref, t) = U_ref` at every `t`, one gauge among many
  (see `notes/06-fields-and-potentials.md`).
- `include/newton/potentials/force_from_potential.hpp` —
  `GradientForce<PotentialFn> {U, length_scale, center}`, a Force-concept function
  object `F(t, r, v) = -∇U(r, t)`, and the factory
  `force_from_potential(U, length_scale = 1.0, center = {})` (a factory because
  C++17 cannot deduce the template argument of an aggregate). `U` may be `U(r)` or
  `U(t, r)`; the time-dependent form takes priority. The finite-difference step is
  `default_step(r − center, length_scale)`: it follows the distance from `center`,
  so set `center` for a potential centred away from the origin (measured: a
  well 1 km away, force 1 m from it, error 4e-4 without `center`, 3e-10 with).
- `tests/test_potentials.cpp` — the verdict table for every library force
  and the homework fields, time slices, regions (the hole field in two boxes,
  off-centre holes, antiparallel wires), `Inconclusive`, velocity-power labels,
  potentials against textbook forms, reference validation, time-dependent U with
  the gauge identity, round trips F→U→F and U→F→U, and the split helpers.

### analysis/energy.hpp additions, added 2026-11-24

Backward compatible: every existing call compiles and behaves as before.

- `total_energy` and `energy_drift` accept a time-dependent potential
  `double U(double t, const Vec3& r)` as well as `U(r)`, chosen with
  `if constexpr`. With a time-dependent `U` the drift includes real physics,
  `E(t) - E(0) = ∫ ∂U/∂t dt` (Taylor §4.5), not just integrator error.
- `power(F, state)` = `F(t, r, v)·v` [W] (§4.1).
- `work_done(traj, F)` — cumulative work `W[i] = ∫ F·v dt`, aligned with
  `traj.states` (`W[0] = 0`). The power is evaluated at each sample and
  integrated with non-uniform composite Simpson (`detail::cumulative_integral`),
  so a short last interval from `simulate_until` and a `record_every` that does
  not divide the step count both work. Accuracy O((record_every·dt)⁴): record
  every step for bookkeeping.
- `energy_balance(traj, U, F_nc)` — the residual of `Δ(T + U) = W_nc`
  (plus `∫ ∂U/∂t dt` for a time-dependent `U`) at every sample, divided by the
  largest energy that crosses the ledger: the maximum over samples of
  `(T_i, |U_i − U_0|, |W_nc,i|, |∫ ∂U/∂t dt|_i)`, or 1 J only if all vanish. It
  is gauge-invariant (only differences of `U` enter) and contains no `|E₀|`
  (a correct free fall has `E₀ ≈ 0` and still reads ~0, where an `E₀`-based scale
  read 1.0). `∂U/∂t` comes from a 5-point difference in `t` whose step follows
  the sample spacing. Sign convention: `W_nc < 0` means `F_nc` removes energy
  (drag); the residual is (energy change) − (booked work), so a NEGATIVE residual
  means energy was lost that `F_nc` does not book (a forgotten drag) and a
  positive one that energy appeared from nowhere (a forgotten drive). Requirement:
  `-∇U + F_nc` must equal the force that was simulated. Accuracy: quadrature
  O((record_every·dt)⁴) plus the integrator's own energy error (RK4 O(dt⁴),
  velocity Verlet O(dt²), Euler O(dt)).
- `tests/test_analysis.cpp` gained checks for each of these.

---

## Roadmap (v0.2+, documented in notes/04-extending.md, not implemented)

N-particle `System` with pairwise forces (Newton III by construction);
cubic-Hermite dense output + event detection by root-finding; Boris pusher
for magnetic forces; adaptive RK45; 1D potential toolkit (turning points,
period by quadrature); analytic gradients supplied by a potential, and a
reference-at-infinity mode for `potential_from_force`; Fourier tools for §5.7–5.9; the Lagrangian hook (constraint
solver, then dual-number autodiff) per taylor map §7.5.
