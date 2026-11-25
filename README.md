# newton

A small, header-only C++17 library for single-particle Newtonian mechanics —
vectors, forces, four numerical integrators, the closed-form textbook
solutions to check them against, and tools for the energy chapter (is this
force conservative? what is its potential energy?) — written alongside
Taylor's *Classical Mechanics*, chapters 1–5.

## Quick start

```cpp
#include <newton/newton.hpp>
#include <cstdio>

int main() {
    using namespace newton;

    // Taylor's baseball: D = 7 cm, m = 0.145 kg, launched at 30 m/s, 45 deg.
    Particle ball{0.145, {0.0, 0.0, 0.0}, {21.213, 21.213, 0.0}};

    forces::UniformGravity gravity{ball.mass};
    forces::QuadraticDrag  drag{forces::quadratic_drag_sphere_air(0.07)};

    SimOptions opt;
    opt.dt     = 1e-4;
    opt.t_end  = 5.0;
    opt.method = Method::RK4;

    // Integrate gravity + drag together (Newton's second law: forces add),
    // stopping exactly when the ball returns to y = 0.
    Trajectory traj = simulate_until(ball, forces::sum(gravity, drag), opt,
                                      events::ground_impact());

    std::printf("range = %.3f m, flight time = %.3f s\n", traj.back().r.x, traj.back().t);
    return 0;
}
```

```
clang++ -std=c++17 -Wall -Wextra -Wpedantic -I include quickstart.cpp -o quickstart
```

```
range = 59.136 m, flight time = 3.796 s
```

(Compare to the vacuum range at the same launch speed/angle,
`analytic::range(30.0, constants::pi / 4.0)` ≈ 91.8 m — the ~35% shorter
range is drag, exactly the effect Taylor §2.4 introduces this problem to
illustrate.)

## Requirements

- **C++17** — no later standard needed, none of the code reaches for C++20.
- **No dependencies** beyond the C++ standard library — header-only, nothing
  to link.
- Verified with **Apple clang 21.0.0** (arm64, Darwin 27.2.0), and earlier with Apple clang 15.0.0. Any
  reasonably modern clang or gcc with `-std=c++17` should work.
- **Python 3 + matplotlib**, optional — only needed for
  `scripts/plot_csv.py`; it degrades to a message and exits cleanly if
  matplotlib isn't installed.

## Installation

The library is one directory, `include/newton/`, with no build step. Three
ways to make `#include <newton/newton.hpp>` resolve, in the order you should
reach for them — full detail and rationale in `notes/00-start-here.md`:

1. **`-I include`** on the compile command (what this repo's own `Makefile`
   and `build.sh` do). Simplest, most explicit, no system changes — use this
   for a long time, including from a shell alias for scratch files outside
   the repo.
2. **Symlink into `/usr/local/include`** once the API you're using has
   stabilized: `sudo ln -s "$PWD/include/newton" /usr/local/include/newton`
   — edits in this repo are then live everywhere, with no reinstall step.
   `sudo make install` / `sudo make uninstall` do the same job with a `cp
   -R` instead, for a frozen rather than live copy.
3. **`CPLUS_INCLUDE_PATH`** — an environment variable clang/gcc both honor
   as an extra search path. Works, but it's invisible to anyone (including a
   future you in a new terminal) who doesn't have it set — reach for it only
   for quick, temporary experiments.

## Modules

Every header is reachable individually, or all at once via the umbrella
`#include <newton/newton.hpp>`.

**Core** — vocabulary types used everywhere.

| Header | Namespace | Contents |
|---|---|---|
| `newton/core/vec3.hpp` | `newton` | `Vec3`; `+ - * /`; `dot`, `cross`, `norm`, `norm2`, `normalized` |
| `newton/core/particle.hpp` | `newton` | `State {r, v, t}`, `Particle {mass, r, v}`, `momentum(particle)` |
| `newton/core/constants.hpp` | `newton::constants` | `pi`, `deg`, `g_earth`, `G`, `sphere_air_beta`, `sphere_air_gamma` |
| `newton/core/trajectory.hpp` | `newton` | `Trajectory {mass, states}` — `size/empty/front/back` |

**Coords** — cylindrical and spherical polar coordinates (Taylor §1.7, §4.8).
Simulate in Cartesian, convert at the edges.

| Header | Namespace | Contents |
|---|---|---|
| `newton/coords/cylindrical_polar.hpp` | `newton::coords` | `CylindricalPolar {rho, phi, z}`, `to_cylindrical(Vec3)` |
| `newton/coords/spherical_polar.hpp` | `newton::coords` | `SphericalPolar {r, theta, phi}` (θ down from +z), `to_spherical(Vec3)` |
| `newton/coords/cartesian.hpp` | `newton::coords` | `to_cartesian` for both types, plus direct cylindrical ↔ spherical conversions |

**Forces** — `Vec3 operator()(double t, const Vec3& r, const Vec3& v) const`
functors, composable with `forces::sum(...)`.

| Header | Namespace | Contents |
|---|---|---|
| `newton/forces/forces.hpp` | `newton::forces` | module umbrella — includes the six below |
| `newton/forces/gravity.hpp` | `newton::forces` | `UniformGravity`, `PointGravity` (with Plummer softening) |
| `newton/forces/drag.hpp` | `newton::forces` | `LinearDrag`, `QuadraticDrag`, `linear_drag_sphere_air`, `quadratic_drag_sphere_air` |
| `newton/forces/spring.hpp` | `newton::forces` | `Spring`, `AnisotropicSpring` |
| `newton/forces/drives.hpp` | `newton::forces` | `SinusoidalDrive`, `ConstantForce` |
| `newton/forces/electromagnetic.hpp` | `newton::forces` | `MagneticForce`, `LorentzForce` |
| `newton/forces/compose.hpp` | `newton::forces` | `sum(f1, f2, ...)` — variadic superposition of forces |

**Integrators and driver** — one-step methods `State -> State`, and the
fixed-step time loop that drives them.

| Header | Namespace | Contents |
|---|---|---|
| `newton/integrators/euler.hpp` | `newton::integrators` | `euler_step` — explicit Euler, pedagogical foil only |
| `newton/integrators/semi_implicit_euler.hpp` | `newton::integrators` | `semi_implicit_euler_step` — symplectic Euler–Cromer, the cheap default |
| `newton/integrators/velocity_verlet.hpp` | `newton::integrators` | `velocity_verlet_step` — 2nd-order symplectic, predictor–corrector for v-dependent forces |
| `newton/integrators/rk4.hpp` | `newton::integrators` | `rk4_step` — classic 4th-order Runge–Kutta, the accuracy workhorse |
| `newton/simulate.hpp` | `newton` | `Method`, `SimOptions`, `simulate(...)`, `simulate_until(...)`, `events::ground_impact()`, `events::apex()` |

**Analytic** — closed-form textbook solutions, the validation oracle for
everything above.

| Header | Namespace | Contents |
|---|---|---|
| `newton/analytic/projectile.hpp` | `newton::analytic` | `range`, `apex_height`, `time_of_flight`, `position` (vacuum) |
| `newton/analytic/drag.hpp` | `newton::analytic` | `LinearDragParams`; linear-drag velocity/position/`y(x)`; `linear_drag_range` and Taylor's small-drag `linear_drag_range_small_drag`; quadratic-drag 1D closed forms (speed, fall distance, `quadratic_drag_horizontal_speed` and `_distance`). All are exact at zero drag |
| `newton/analytic/oscillator.hpp` | `newton::analytic` | `SHM`, `Regime`, `DampedOscillator`, `DrivenOscillator` |
| `newton/analytic/magnetic.hpp` | `newton::analytic` | `cyclotron_omega`, `cyclotron_radius`, `cyclotron_period` |
| `newton/analytic/rocket.hpp` | `newton::analytic` | `tsiolkovsky_delta_v`, `rocket_velocity_gravity` |

**Fields** — numerical calculus on any scalar or vector field you can write as a
function: `double U(const Vec3& r)` or `Vec3 F(const Vec3& r)` (Taylor §4.1–4.4).

| Header | Namespace | Contents |
|---|---|---|
| `newton/fields/derivatives.hpp` | `newton::fields` | `default_step`, `Jacobian`, `gradient`, `jacobian`, `curl`, `divergence` (5-point central differences; explicit-step and default-step overloads), and `curl`/`divergence`/`frobenius_norm` of a `Jacobian` |
| `newton/fields/line_integral.hpp` | `newton::fields` | `LineIntegral {value, error, scale, evaluations, converged}`, `line_integral` (segment or polyline), `circulation` (circle), and the `_detailed` versions that return the whole struct — adaptive Gauss–Kronrod |

**Potentials** — "is this force conservative, and what is its U?", and the reverse
(Taylor §4.2–4.5). Accepts forces `F(t, r, v)`, `F(t, r)` or `F(r)`.

| Header | Namespace | Contents |
|---|---|---|
| `newton/potentials/conservative.hpp` | `newton::potentials` | `check_conservative` → `ConservativeReport` (a `Verdict`, a `VelocityPower` label for the velocity-dependent part, the evidence, and a `note` when inconclusive); `probe_dependence`; `describe`, `to_string`; `Box`, `CheckOptions`; `velocity_independent_part`, `velocity_dependent_part` |
| `newton/potentials/potential_from_force.hpp` | `newton::potentials` | `potential_from_force` → `LineIntegralPotential`; `find_potential` → `PotentialResult` (classify, validate the reference point by checking path independence, then hand back U, certified inside the tested box) |
| `newton/potentials/force_from_potential.hpp` | `newton::potentials` | `force_from_potential(U, length_scale, center)` → `GradientForce`, F = −∇U as a Force you can pass to `simulate` (set `length_scale` for atomic sizes, `center` for a potential centred away from the origin) |

**Analysis, systems, io**

| Header | Namespace | Contents |
|---|---|---|
| `newton/analysis/energy.hpp` | `newton::analysis` | `kinetic_energy`, `total_energy`, `energy_drift` (U may depend on time), `power`, `work_done`, `energy_balance` (the ledger for drag and drives), `uniform_gravity_U` (g along −y only), `spring_U` (anchor at the origin, natural length 0 only; for any other setup write the U or use `find_potential`) |
| `newton/analysis/momentum.hpp` | `newton::analysis` | `angular_momentum`, `center_of_mass`, `total_momentum` |
| `newton/analysis/oscillation.hpp` | `newton::analysis` | `PeakList`, `find_peaks`, `measure_period`, `measure_amplitude` |
| `newton/systems/rocket.hpp` | `newton::systems` | `RocketParams`, `RocketSample`, `simulate_rocket` — self-contained variable-mass flight |
| `newton/io/csv.hpp` | `newton::io` | `write_csv` — trajectory and generic-column CSV export |

## Building examples and tests

```
make            # compile every examples/*.cpp into build/
make check      # compile + run every tests/*.cpp
make clean      # remove build/
make debug      # rebuild with sanitizers + -O0 -g
./build.sh      # the same job as a plain shell script instead of a Makefile
```

Examples live in `examples/`, one physical situation per file, each printing
a numeric-vs-expected comparison to stdout (and, for most, writing a CSV to
`build/` for plotting with `python3 scripts/plot_csv.py
build/whatever.csv`). `notes/00-start-here.md` walks through the Makefile and
the `clang++` invocation flag by flag if any of that is unfamiliar.

| Example | Taylor | What it shows |
|---|---|---|
| `01_free_fall.cpp` | §1.6 | A drop and a 45° projectile against the closed forms; `simulate_until` |
| `02_projectile_drag.cpp` | §2.1–2.4 | Taylor's baseball under vacuum, linear, and quadratic drag; the best launch angle; the exact linear-drag range against the small-drag formula |
| `03_charged_particle.cpp` | §2.5–2.7 | The helix in a magnetic field, and the E×B drift |
| `04_rocket.cpp` | §3.2 | The rocket against the Tsiolkovsky equation with gravity |
| `05_oscillations.cpp` | §5.1–5.6 | Simple harmonic motion, the three damping regimes, a resonance sweep |
| `06_integrator_comparison.cpp` | numerics (energy is ch. 4) | All four integrators on one oscillator for 200 periods |
| `07_potential_from_force.cpp` | §4.1–4.5, §5.4 | "Is this force conservative, and what is its U?": the same endpoints along two paths, the curl test, a verdict table, `find_potential` against textbook U, the hole field, a time-dependent U, and a closed energy ledger for gravity plus drag (no CSV) |
| `08_force_from_potential.cpp` | §4.3–4.5, §5.1, §5.4 | "From U to F to motion": a Morse potential through `force_from_potential` into `simulate`, small oscillations against ω = √(U″/m), a moving-support U(t, r), and the ledger with drag added (CSV: `build/morse_damped.csv`) |

## Learning more

- **`notes/00-start-here.md`** — a guided tour of this repository: every
  directory and file type, the build tooling, and how to make the include
  path work from outside this repo. Start here; it also gives a suggested
  order for the notes below.
- **`notes/01-cpp-for-this-library.md`** — every non-basic C++ feature the
  library uses (template callables, fold expressions, `inline constexpr`,
  and more), each grounded in a real snippet from `include/newton/`.
- **`notes/02-physics-to-code.md`** — the bridge from a Taylor section to
  the function that implements it, chapter by chapter.
- **`notes/03-numerical-methods.md`** — the four integrators in depth: what
  each buys and costs, and how to pick.
- **`notes/04-extending.md`** — the v0.2 roadmap below, expanded into a
  semester-length plan for growing past v0.1.
- **`notes/05-library-reference.md`** — every function in the library with a
  compiled usage example, plus how the repo and the tests fit together.
- **`notes/06-fields-and-potentials.md`** — the physics, numerics, and C++
  behind `fields` and `potentials`: how the derivatives and line integrals
  work, how `check_conservative` decides, and what its verdicts can and
  cannot prove.
- **`DESIGN.md`** — the authoritative API specification: exact signatures,
  file paths, and semantics for every header above.
- **`docs/research/`** — the physics (`taylor-ch1-7-map.md`), the numerics
  (`numerical-integration.md`), and the C++/tooling conventions
  (`cpp-library-design.md`) this library was designed from.

## Status and roadmap

**v0.1 (current)** covers Taylor chapters 1–5 for a **single particle**:
vacuum and drag-affected projectiles, charged-particle motion in uniform E/B
fields, a self-contained variable-mass rocket, and the full oscillator
family (undamped, damped in all four regimes, driven steady-state) —
against four hand-written integrators (explicit Euler, semi-implicit
Euler, velocity Verlet, RK4) and a fixed-step driver with bisection-refined
event stopping. Since the first cut it has also gained polar coordinates
(`coords`), numerical derivatives and line integrals (`fields`), a
conservative-force classifier with potential-from-force and force-from-
potential (`potentials`), and an energy ledger for non-conservative forces
(`analysis::work_done`, `energy_balance`). Deliberately out of scope for
v0.1: multi-particle systems, adaptive step sizes, and anything beyond
Taylor ch. 5.

**v0.2 roadmap** (design only — not implemented; see `DESIGN.md`'s roadmap
section and `notes/04-extending.md` for the expanded plan; the `potentials`
item from the original list has been built):

- An N-particle `System` with pairwise forces (Newton's third law by
  construction)
- Cubic-Hermite dense output plus event detection by root-finding
- A Boris pusher for magnetic forces
- Adaptive RK45
- A 1D potential toolkit (turning points, period by quadrature)
- Analytic gradients supplied by a potential, and a reference-at-infinity
  mode for `potential_from_force` (the rest of the `potentials` plan is
  done: see `notes/04-extending.md`)
- Fourier tools for Taylor §5.7–5.9
- A Lagrangian hook (constraint solver, then dual-number autodiff), per the
  Taylor map's §7.5
