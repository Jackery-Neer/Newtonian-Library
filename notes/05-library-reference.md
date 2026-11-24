# 05 — Library reference

Every function that `make install` copies, what it's for, and how to call it. After that:
how this library fits the way C++ libraries have always been built, what everything else in
the repo is for, and how the tests work.

*Checked against the code on 2026-11-24: 30 headers, 53 functions (counting overloads
separately), 26 types. Every usage snippet below was compiled with strict warnings
(`-Wall -Wextra -Wpedantic -Wconversion -Wshadow`) and run. If this file and a header ever
disagree, the header wins; tell the next update which one changed.*

---

## Part 1 — The API: what gets installed

`make install` copies exactly one folder, `include/newton/`, into your system include
directory. Nothing else goes with it. That folder *is* the library; everything else in the repo
exists to build, test, explain, or demonstrate it.

### The shortest complete program

```cpp
#include <newton/newton.hpp>
#include <cstdio>

int main() {
    using namespace newton;                       // fine in your own .cpp, never in a header
    Particle ball{0.145, {}, {20.0, 15.0, 0.0}};  // mass [kg], start [m], launch velocity [m/s]
    auto force = forces::sum(forces::UniformGravity{ball.mass},
                             forces::QuadraticDrag{forces::quadratic_drag_sphere_air(0.07)});
    SimOptions opt;
    opt.t_end = 10.0;
    opt.dt = 1e-3;
    Trajectory flight = simulate_until(ball, force, opt, events::ground_impact());
    std::printf("range = %.2f m after %.2f s\n", flight.back().r.x, flight.back().t);
}
```

It prints `range = 45.16 m after 2.81 s`. Almost every program follows the same five steps:

1. **Describe the body:** a `Particle` (mass and initial conditions).
2. **Describe the forces:** force structs, added together with `forces::sum`.
3. **Choose how to step:** a `SimOptions` (time span, step size, method).
4. **Run it:** `simulate` or `simulate_until`.
5. **Read the results:** the `Trajectory`, the `analysis` functions, or `io::write_csv`.

**Rules that apply everywhere:**
- **SI units, angles in radians.** Multiply degrees by `constants::deg`.
- **y is up.** `UniformGravity`'s default `g` points along −y, and `events::ground_impact()` watches `r.y`.
- **Folder = namespace.** `include/newton/forces/` holds `newton::forces`, `include/newton/analytic/` holds `newton::analytic`, and so on.
- **Every pure function is `[[nodiscard]]`.** Throwing away a result is a compiler warning, because it's always a mistake.

### core — `newton` and `newton::constants`

| Name | What it is |
|---|---|
| `constants::pi`, `constants::deg` | π, and π/180: `30.0 * constants::deg` is 30° in radians |
| `constants::g_earth` | 9.80665 m/s² |
| `constants::G` | 6.67430×10⁻¹¹ m³ kg⁻¹ s⁻² |
| `constants::sphere_air_beta`, `sphere_air_gamma` | Taylor's drag coefficients for a sphere in air (§2.1) |
| `Vec3 {x, y, z}` | A 3D vector, (0, 0, 0) by default. Supports `+ − * /`, unary `−`, `+= −= *= /=`, `== !=` |
| `dot(a, b)` | a·b |
| `cross(a, b)` | a×b (order matters: `cross(b, a)` is the negative) |
| `norm(v)`, `norm2(v)` | \|v\| and \|v\|² (`norm2` skips the square root) |
| `normalized(v)` | Unit vector in v's direction; returns the zero vector for zero input instead of NaN |
| `State {r, v, t}` | One instant of motion: position, velocity, time |
| `Particle {mass, r, v}` | A body's mass and initial conditions. `initial_state(t0)` turns it into a `State` |
| `momentum(particle)` | p = m v |
| `Trajectory {mass, states}` | A simulation's output: every recorded `State`, in time order. Has `size()`, `empty()`, `front()`, `back()` |

```cpp
Vec3 a{1, 2, 3}, b{0, 1, 0};
Vec3 c = 2.0 * a + b;                           // (2, 5, 6)
double ab = dot(a, b);                          // 2
Vec3 n = cross(a, b);                           // (-3, 0, 1)
double angle = 30.0 * constants::deg;           // 0.5236 rad
Particle ball{0.145, {0, 0, 0}, {20, 15, 0}};
Vec3 p = momentum(ball);                        // (2.9, 2.175, 0)
```

### coords — `newton::coords`

| Name | What it is |
|---|---|
| `CylindricalPolar {rho, phi, z}` | ρ = distance from the z axis, φ = angle around z measured from +x, z unchanged |
| `SphericalPolar {r, theta, phi}` | Physics convention: θ measured *down from +z*, in [0, π]; φ around z, as above |
| `to_cylindrical(Vec3)`, `to_spherical(Vec3)` | Cartesian → curvilinear |
| `to_cartesian(CylindricalPolar)`, `to_cartesian(SphericalPolar)` | Back to a `Vec3`; the argument's type decides which one runs |
| `to_cylindrical(SphericalPolar)`, `to_spherical(CylindricalPolar)` | Directly between the two |

```cpp
auto cyl = coords::to_cylindrical(Vec3{3, 4, 5});   // {5, 0.9273, 5}
auto sph = coords::to_spherical(Vec3{0, 0, 2});     // {2, 0, 0}: on the +z axis
Vec3 back = coords::to_cartesian(sph);              // (0, 0, 2)
Vec3 q = coords::to_cartesian(coords::SphericalPolar{2.0, constants::pi / 2, 0.0});  // (2, 0, 0)
```

These types have no arithmetic on purpose: simulate in Cartesian and convert at the edges. With
literal numbers, write the type out. `to_cartesian({2.0, 1.0, 0.0})` doesn't compile, because
the same three numbers name different points in the two systems.

### forces — `newton::forces`

A force is any callable with the signature `Vec3 operator()(double t, const Vec3& r, const Vec3& v) const`
that returns newtons. Each library force is a small struct: you set its parameters with braces,
and then it acts like a function.

| Force | Law | Parameters |
|---|---|---|
| `UniformGravity` | F = m g | `mass` [kg]; `g` [m/s²], default (0, −9.80665, 0) |
| `PointGravity` | F = −GM m d/\|d\|³, with d = r − center | `GM` [m³/s²], `mass`, `center`, `softening` [m] |
| `LinearDrag` | F = −b v | `b` [kg/s] |
| `QuadraticDrag` | F = −c \|v\| v | `c` [kg/m] |
| `Spring` | F = −k(\|d\| − L₀) d̂, with d = r − anchor | `k` [N/m], `anchor`, `natural_length` |
| `AnisotropicSpring` | F = (−kx·x, −ky·y, −kz·z) | `kx`, `ky`, `kz` |
| `SinusoidalDrive` | F = F0 cos(ωt + phase) · dir | `F0` [N], `omega`, `phase`, `dir` |
| `ConstantForce` | F = constant | `F` [N] |
| `MagneticForce` | F = q v×B | `q` [C]; `B` [T], default ẑ |
| `LorentzForce` | F = q(E + v×B) | `q`, `E`, `B` |

Helpers: `linear_drag_sphere_air(D)` and `quadratic_drag_sphere_air(D)` turn a sphere's
diameter [m] into `b` or `c`. `sum(f1, f2, ...)` combines any number of forces into one.

```cpp
forces::UniformGravity gravity{ball.mass};
forces::QuadraticDrag air{forces::quadratic_drag_sphere_air(0.07)};
auto net = forces::sum(gravity, air);
Vec3 F = net(0.0, ball.r, ball.v);       // forces are callable: F(t, r, v)

// Your own force: any lambda with the same signature works wherever a library force does.
auto wind = [](double, const Vec3&, const Vec3&) { return Vec3{0.5, 0.0, 0.0}; };
auto net2 = forces::sum(gravity, air, wind);
```

That lambda is the main way to extend the library: you never have to edit it to try a new force law.

### integrators — `newton::integrators`

| Function | Accuracy | Use it for |
|---|---|---|
| `euler_step` | 1st order; energy grows on oscillators | Demonstrating why it's bad |
| `semi_implicit_euler_step` | 1st order; energy stays bounded | Cheap, stable runs |
| `velocity_verlet_step` | 2nd order; energy stays bounded | Long runs of conservative systems |
| `rk4_step` | 4th order; 4 force evaluations per step | The default; anything with drag |

All four have the signature `State xxx_step(const State& s, double mass, double dt, Force&& force)` and return the
state one `dt` later. You rarely call them directly, because `simulate` calls them for you.
Use them when you want to write your own time loop:

```cpp
State next = integrators::rk4_step(ball.initial_state(), ball.mass, 1e-3, net);
```

### simulate — `newton`

| Name | What it does |
|---|---|
| `Method` | `Euler`, `SemiImplicitEuler`, `VelocityVerlet`, `RK4` |
| `SimOptions` | `t0` (default 0), `t_end` (1), `dt` (1e-3), `method` (RK4), `record_every` (1) |
| `simulate(particle, force, opt)` | Steps from `t0` to `t_end` and returns a `Trajectory` |
| `simulate_until(particle, force, opt, event)` | The same, but stops the moment `event(state)` changes sign. The stopping point is refined to the integrator's accuracy, not just to the nearest step |
| `events::ground_impact(y_level = 0)` | Stops when the particle comes down through `y_level` |
| `events::apex()` | Stops at the top of the arc (where v.y = 0) |

An event is any function from a `State` to a `double`, so you can write your own:
`[](const State& s) { return s.r.x - 10.0; }` stops when x crosses 10 m.

```cpp
SimOptions opt;
opt.t_end = 10.0;
opt.dt = 1e-3;
opt.method = Method::RK4;
opt.record_every = 10;                                                     // keep every 10th step
Trajectory traj = simulate(ball, net, opt);                                // 1001 samples
Trajectory flight = simulate_until(ball, net, opt, events::ground_impact());
double range = flight.back().r.x;                                          // 45.16 m
```

Pick `dt` from the fastest thing in the problem: about a hundredth of a period for an
oscillator, and a tenth of τ = m/b for drag. `notes/03-numerical-methods.md` explains why.

### analytic — `newton::analytic` (closed-form answers)

These functions don't simulate anything. They are the textbook's solutions, written as
functions, so simulations have answer keys to check against.

| Function or type | Gives | Taylor |
|---|---|---|
| `range(v0, θ, g)`, `apex_height(...)`, `time_of_flight(...)`, `position(t, v0, θ, g)` | Vacuum projectile results (`g` defaults to `g_earth`) | §1.6 |
| `LinearDragParams {m, b, g}` | Holds the parameters; `tau()` = m/b, `v_ter()` = mg/b | §2.2 |
| `linear_drag_velocity(p, v0, t)`, `linear_drag_position(p, r0, v0, t)` | Exact v(t) and r(t) under gravity plus linear drag | §2.2–2.3 |
| `linear_drag_trajectory_y(p, vx0, vy0, x)` | The exact path y(x); only valid for x < vx0·τ | §2.3 |
| `quadratic_drag_terminal_speed(m, c, g)` | √(mg/c) | §2.4 |
| `quadratic_drag_fall_speed(v_ter, g, t)`, `quadratic_drag_fall_distance(v_ter, g, t)` | The tanh and ln cosh results for a fall from rest | §2.4 |
| `quadratic_drag_horizontal_speed(m, c, v0, t)` | v0 / (1 + c·v0·t/m) | §2.4 |
| `cyclotron_omega(q, B, m)`, `cyclotron_radius(m, v_perp, q, B)`, `cyclotron_period(m, q, B)` | Circular motion in a magnetic field | §2.5 |
| `tsiolkovsky_delta_v(v_ex, m0, m)`, `rocket_velocity_gravity(v_ex, m0, burn_rate, g, t)` | Rocket equation, without and with gravity | §3.2 |
| `SHM {omega0, x0, v0}` | `x(t)`, `v(t)`, `amplitude()`, `phase()`, `period()`, `energy(m)` | §5.2 |
| `DampedOscillator {omega0, beta, x0, v0}` | `regime()`, `omega1()`, `x(t)`, `v(t)`, `log_decrement()`, correct in all four regimes | §5.4 |
| `Regime` | `Undamped`, `Underdamped`, `CriticallyDamped`, `Overdamped` | §5.4 |
| `DrivenOscillator {omega0, beta, f0, omega}` | Steady state: `amplitude()`, `phase()`, `resonant_omega()`, `quality_factor()`, `fwhm()` | §5.5–5.6 |

```cpp
double R = analytic::range(20.0, 45.0 * constants::deg);        // 40.79 m
analytic::LinearDragParams ldp{1.0, 0.5};                        // m, b (g defaults)
Vec3 v2 = analytic::linear_drag_velocity(ldp, {10, 5, 0}, 2.0);  // velocity at t = 2 s
analytic::SHM shm{2.0, 1.0, 0.0};                                // ω0, x0, v0
double T = shm.period();                                         // π s
analytic::DampedOscillator d{2.0, 0.5, 1.0, 0.0};                // ω0, β, x0, v0: set all four
analytic::DrivenOscillator drv{2.0, 0.1, 1.0, 1.9};              // ω0, β, f0, ω
double Q = drv.quality_factor();                                 // 10
```

### analysis — `newton::analysis` (measuring a simulation)

| Function | What it returns |
|---|---|
| `kinetic_energy(m, v)` | ½ m v² |
| `total_energy(m, state, U)` | T + U(r), for any potential-energy function `U` |
| `energy_drift(traj, U)` | (E(t) − E(0)) / \|E(0)\| for every sample. The best "is my `dt` small enough?" check |
| `uniform_gravity_U(m, g)`, `spring_U(k)` | Ready-made `U` functions to pass into the two above |
| `angular_momentum(m, state, about)` | (r − about) × m v; `about` defaults to the origin |
| `center_of_mass(particles)` | Σ mr / Σ m over a `std::vector<Particle>` |
| `total_momentum(particles)` | Σ m v |
| `find_peaks(traj)` | The local maxima of x, as a `PeakList {t, x}` |
| `measure_period(traj)` | The average time between upward zero crossings of x |
| `measure_amplitude(traj, t_start)` | The largest \|x\| after `t_start` (for a steady state) |

The last three only look at the x component (`r.x`).

```cpp
auto U = analysis::uniform_gravity_U(ball.mass, constants::g_earth);
Trajectory vac = simulate(ball, gravity, opt);                   // no drag, so energy is conserved
std::vector<double> drift = analysis::energy_drift(vac, U);      // about 1e-12 with RK4
Vec3 L = analysis::angular_momentum(ball.mass, flight.back());   // about the origin
Particle other{1.0, {4.0, 0.0, 0.0}, {}};
Vec3 R_cm = analysis::center_of_mass({ball, other});
```

### systems — `newton::systems`

| Name | What it is |
|---|---|
| `RocketParams {m0, m_dry, v_ex, burn_rate, g}` | Wet mass, dry mass, exhaust speed, burn rate [kg/s]; `g` defaults to `g_earth` |
| `RocketSample {t, m, v, y}` | One instant: time, mass, vertical velocity, altitude |
| `simulate_rocket(params, dt, t_end)` | Vertical flight, burn then coast. The burnout instant is always one of the samples |

The rocket has its own simulator because its mass changes, and `Particle` assumes a fixed mass.

```cpp
systems::RocketParams rp{1e5, 4e4, 3000.0, 1000.0};       // m0, m_dry, v_ex, burn_rate
auto samples = systems::simulate_rocket(rp, 0.01, 90.0);  // dt, t_end
```

### io — `newton::io`

| Function | Writes |
|---|---|
| `write_csv(path, traj)` | The columns `t,x,y,z,vx,vy,vz` |
| `write_csv(path, names, rows)` | Any columns you choose |

Both return `false` if the file can't be written. They're `[[nodiscard]]`, so check the result:

```cpp
if (!io::write_csv("build/flight.csv", flight)) std::puts("write failed");
std::vector<std::vector<double>> rows{{30.0, 50.1}, {45.0, 59.1}};
if (!io::write_csv("build/sweep.csv", {"angle", "range"}, rows)) std::puts("write failed");
```

Then plot it with `python3 scripts/plot_csv.py build/flight.csv`.

### Three things to watch for

1. **A few structs have no default values.** `DampedOscillator` (`omega0`, `beta`), `RocketParams`
   (`m0`, `m_dry`, `v_ex`, `burn_rate`) and `RocketSample` leave those members uninitialized if you
   don't set them. `analytic::DampedOscillator d;` gives you garbage values. Always fill them in with
   braces. Every other struct in the library gives each member a default, so these are worth fixing.
2. **A `[[nodiscard]]` warning means you threw away a result.** Use the result, or decide you really don't need it.
3. **Nothing simulates two or more interacting bodies yet.** `center_of_mass` and `total_momentum`
   measure particles you set up yourself, but nothing integrates bodies pushing on each other.
   That's the biggest gap for Chapter 3.

---

## Part 2 — Philosophy: is a library just its functions?

You write the functions, and you also show how to use them. The functions are the product, but
every C++ library that people actually use ships three layers:

| Layer | Question it answers | In newton |
|---|---|---|
| **Interface** | What can I call? | The headers in `include/newton/`, plus their `///` doc comments |
| **Demonstration** | How do I use it? | The README quick start, `examples/`, `notes/` |
| **Evidence** | Can I trust it? | `tests/` |

How this came about:

- **C, in the 1970s.** The Unix C library documented every function with a manual page: a SYNOPSIS
  (the exact declaration from the header), a DESCRIPTION, and often an example. They still exist;
  try `man 3 printf`. *The C Programming Language* (Kernighan and Ritchie) taught the language almost
  entirely through complete example programs.
- **The STL, in 1994.** Alexander Stepanov's Standard Template Library introduced *generic programming*:
  write an algorithm once, separately from the data it works on, and let the user plug in the
  specifics. `std::sort(begin, end, compare)` takes *your* comparison function. newton follows this
  directly: `simulate(particle, force, options)` takes *your* force. The library supplies the algorithm,
  and you supply the problem. That's also the answer to the earlier "a function per problem, or general
  machinery?" question.
- **Boost, from 1998 on.** Boost is a peer-reviewed collection of C++ libraries, and its review expects
  documentation and tests along with the code. Many Boost libraries later joined the standard library
  (smart pointers, `std::function`, regex, filesystem). This is where "a library ships examples and
  tests" became the normal expectation in C++.
- **Modern header-only libraries.** Libraries like nlohmann/json, which newton's folder layout is
  modeled on, have READMEs made almost entirely of usage examples, because people decide whether to use
  a library from its first example.

Two principles to take from this:

1. **Your examples are your API's first user.** If an example is awkward to write, the function is
   badly designed, so fix the function, not the example. Writing the example is how you find out.
2. **Everything public is a promise.** Once someone, including future you, writes code against a
   function, renaming it or changing its parameters breaks their code. That's why helpers live in a
   `detail` namespace: it means "internal, may change". Changing things is still cheap now; once a
   library has users, it uses version numbers to announce breaking changes.

**Where the demonstration lives.** It isn't in the installed headers. `make install` copies only
`include/`, so someone using the installed library gets the doc comments and nothing else. The doc
comments are the only documentation that travels with the code, which is why every function in
newton has one, with units. The README, examples, and tests are read in the repo (or on GitHub, once
it's published).

**Your coords module covers two of the three layers.** The headers have doc comments, `test_coords.cpp`
is the evidence, and `DESIGN.md` is updated. The missing layer is demonstration: no example uses
coords yet, and the README's module table doesn't mention it.

---

## Part 3 — Everything else in the repo

| Path | What it is | Installed? |
|---|---|---|
| `include/newton/` | The library (Part 1) | **Yes, the only thing installed** |
| `README.md` | The front page: quick start, module table, install instructions. Its module table predates coords | No |
| `DESIGN.md` | The spec: the hard rules, the Force concept, every file's API | No |
| `Makefile`, `build.sh` | Build the examples and tests, run the tests, install | No |
| `examples/` | Six lesson programs (table below) | No |
| `tests/` | Seven test programs plus the harness, `test_util.hpp` (Part 4) | No |
| `scripts/plot_csv.py` | Plots any CSV the library writes | No |
| `docs/research/` | Six research reports: the Taylor ch. 1–7 map, numerical methods, C++ library design, and three Chapter 1 reports | No |
| `notes/` | Your study notes: 00–05, `learning-path.md`, `chapters/`, `progress/` | No |
| `.vscode/` | Editor settings for the C/C++ Runner extension and IntelliSense | No |
| `.gitignore` | Tells git to ignore `build/`, `.DS_Store`, and `*.csv`. Does nothing until you run `git init` | No |
| `build/` | Generated programs and CSVs. Safe to delete (`make clean`) | No |

| Example | Taylor | What it shows |
|---|---|---|
| `01_free_fall.cpp` | §1.6 | A drop and a 45° projectile against the closed forms; `simulate_until` |
| `02_projectile_drag.cpp` | §2.1–2.4 | Taylor's baseball under vacuum, linear, and quadratic drag; the best launch angle |
| `03_charged_particle.cpp` | §2.5–2.7 | The helix in a magnetic field, and the E×B drift |
| `04_rocket.cpp` | §3.2 | The rocket against the Tsiolkovsky equation with gravity |
| `05_oscillations.cpp` | §5.1–5.6 | Simple harmonic motion, the three damping regimes, a resonance sweep |
| `06_integrator_comparison.cpp` | numerics (energy is ch. 4) | All four integrators on one oscillator for 200 periods |

**About `.vscode/`.** The C/C++ Runner extension compiles single files its own way, separately from the
Makefile: its debug setting points at `examples/build/Debug/outDebug`, not at `build/`. Its warning flags
are stricter than the Makefile's (`-Wconversion`, `-Wshadow`, and others), and the library compiles cleanly
under those too (checked 2026-11-24). If the extension and `make` ever disagree, trust `make`, because
that's what the tests run under.

---

## Part 4 — Tests

### What a test is

Each file in `tests/` is an ordinary C++ program with its own `main()`. It calls library functions
and checks the results with three macros from `tests/test_util.hpp`:

| Macro | Passes when |
|---|---|
| `NT_CHECK(condition)` | The condition is true |
| `NT_CLOSE(a, b, tol)` | `a` and `b` agree within `tol` |
| `NT_CLOSE_VEC3(a, b, tol)` | The same, for vectors, component by component |

The last line is always `return newton_test::report("name");`. It prints a summary and makes the
program exit with 0 if everything passed, or 1 if anything failed. That exit code is the whole
mechanism: there's no test framework, just programs that exit with 0 or 1, and `make` stops when one
exits with 1.

A complete test file:

```cpp
// tests/test_example.cpp
#include "test_util.hpp"        // quotes: this file sits right next to it
#include <newton/newton.hpp>    // angle brackets: found through -I include

using namespace newton;

int main() {
    Particle p{2.0, {0.0, 100.0, 0.0}, {}};      // dropped from rest at 100 m
    SimOptions opt;
    opt.t_end = 3.0;
    opt.dt = 0.01;
    const Trajectory traj = simulate(p, forces::UniformGravity{p.mass}, opt);

    // RK4 is exact for a constant force, so this holds to roundoff at any dt.
    NT_CLOSE(traj.back().r.y, 100.0 - 0.5 * constants::g_earth * 9.0, 1e-10);
    NT_CHECK(traj.back().t == 3.0);

    return newton_test::report("example");
}
```

When a check fails, it prints the file, the line, and the exact expression that failed. Adding
`NT_CLOSE(norm(Vec3{1.0, 2.0, 2.0}), 3.5, 1e-12);` (the norm is really 3) prints:

```
  FAIL tests/test_example.cpp:19  norm(Vec3{1.0, 2.0, 2.0}) ~= 3.5
example: 3 checks, 1 failed
```

### Choosing a tolerance

Computer arithmetic rounds, so results are almost never *exactly* equal.

- **Exact `==`** only for values that were copied or assigned, never computed. `test_coords.cpp` checks
  `c345.z == 5.0` exactly because `z` is copied straight through.
- **About 1e-12** for a handful of arithmetic operations.
- **For simulations, match the method.** RK4's error shrinks like dt⁴ and Euler's like dt. A tolerance
  tighter than the method's real error fails for no reason; a much looser one hides bugs.

### Four kinds of checks in this repo

| Kind | Example | Where |
|---|---|---|
| Exact identity | a×b is perpendicular to both a and b | `test_vec3.cpp` |
| Simulation vs closed form | A simulated quadratic-drag fall matches v_ter·tanh(gt/v_ter) | `test_numeric_vs_analytic.cpp` |
| Method contract | Halving dt shrinks RK4's error 16×, so it really is 4th order | `test_integrators.cpp` |
| Round trip or invariant | Cartesian → spherical → Cartesian returns the same point | `test_coords.cpp` |

### Why they're worth it

- **A safety net.** Every program depends on every header, so any edit could break anything. After
  each change, `make check` tells you within seconds whether you broke something anywhere.
- **They check the physics, not just the code.** A simulation that agrees with Taylor's closed form
  to 1e-9 means your derivation, the force, and the integrator are all right at once.
- **They catch real bugs.** The rocket's burnout bug, where thrust cut off one RK4 stage too early,
  showed up as a mismatch against the closed-form rocket equation, which is exactly the comparison
  the tests make. `test_numeric_vs_analytic.cpp` now checks that instant, so the bug can't come back.
- **They're precise documentation.** `test_coords.cpp` shows exactly what `to_spherical` returns at
  the origin, and why `atan2` beats `acos`, more precisely than a comment could.

### The current suites (2026-11-24)

| Suite | Checks | Covers |
|---|---|---|
| `analysis` | 529 | Energy, momentum, and peak/period/amplitude measurement on made-up data |
| `analytic` | 139 | The closed forms agree with themselves (x(0) = x0, each regime, the resonance peak) |
| `coords` | 2,779 | Your new module: hand-computed points, round trips, `atan2` vs `acos` |
| `forces` | 29 | Each force law against a hand-computed vector |
| `integrators` | 18 | Exactness for constant force, convergence order, energy behavior, time bookkeeping |
| `numeric_vs_analytic` | 15 | Simulations against closed forms: drag, oscillators, rocket, magnetic field, stopping events |
| `vec3` | 33 | Vector algebra identities |

Total: **3,542 checks, all passing.**

### How to run them

| Command | What it does |
|---|---|
| `make check` | Builds any tests that changed, then runs all of them. Stops at the first suite that fails; otherwise ends with "All tests passed." |
| `make tests` | Builds the tests without running them |
| `./build.sh` | Builds everything and runs **every** suite, even after one fails, so you see all failures at once |
| `make build/tests/test_coords && ./build/tests/test_coords` | Runs one suite. The quickest loop while you work on one module |
| `make debug` | Deletes `build/`, rebuilds with debug info and sanitizers (memory-error and undefined-behavior detectors), then runs everything. Use it when something crashes or acts strangely |
| `echo $?` | Run right after any of the above: 0 means success |

### Adding a test

Create `tests/test_<name>.cpp` using the skeleton above. The Makefile and `build.sh` find it
automatically. Write the test right after the function, while you still know what the right answer
is, then run `make check`.

---

## Loose ends found while writing this

1. **No example uses coords yet.** An `examples/07_coordinates.cpp` would complete the demonstration
   layer: for instance, report a projectile's path in spherical coordinates, or show that a uniform
   circular orbit has constant ρ.
2. **The README's module table doesn't list `coords/`.** It's a one-section edit.
3. **Give `DampedOscillator`, `RocketParams` and `RocketSample` default member values**, so they match
   the rest of the library (see "Three things to watch for").
