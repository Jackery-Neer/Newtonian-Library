# 05 — Library reference

Every function that `make install` copies, what it's for, and how to call it. After that:
how this library fits the way C++ libraries have always been built, what everything else in
the repo is for, and how the tests work.

*Checked against the code on 2026-11-25: 35 headers, 87 functions (counting overloads
separately), 37 types. Every usage snippet below was compiled with strict warnings
(`-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion`; `-Wno-unused-variable` for the
snippets that compute a value only to show it) and run. If this file and a header ever
disagree, the header wins; tell the next update which one changed.*

*How the counts were made.* Headers: `find include -name '*.hpp' | wc -l`, which includes
`newton.hpp` and `forces.hpp`. Functions and types: from clang's AST of a one-line program
that includes `<newton/newton.hpp>` (`clang++ -fsyntax-only -Xclang -ast-dump=json
-Xclang -ast-dump-filter=newton`). A function is any namespace-scope function or function
template outside a `detail` namespace, each overload counted separately, with member functions
and operators left out (so `Vec3`'s `+` and `Spring::operator()` are not counted). A type is
any struct, class, class template or enum at namespace scope outside `detail`. The same rule
gives the previous 53 functions and 26 types for the 30 headers that existed before `fields`
and `potentials`; the new work adds 5 headers, 34 functions (17 in `fields`, 11 in
`potentials`, 3 in `analysis/energy.hpp`, 3 in `analytic/drag.hpp`) and 11 types.

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
| `SimOptions` | `t0` (default 0), `t_end` (1), `dt` (1e-3), `method` (RK4), `record_every` (1). Times must be finite with `t_end >= t0`; `dt` must be finite and positive |
| `simulate(particle, force, opt)` | Takes full `dt` steps and at most one shorter final step, ending exactly at `t_end`; invalid timing inputs throw `std::invalid_argument` |
| `simulate_until(particle, force, opt, event)` | The same, but stops at a later event zero or sign change. Interior crossings are refined to integrator accuracy; an initial zero is ignored so a projectile can launch from ground |
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
| `linear_drag_range(p, vx0, vy0)` | The exact level-ground range under gravity plus linear drag (y(R) = 0, solved by bisection; no closed form exists) | §2.4 |
| `linear_drag_range_small_drag(p, vx0, vy0)` | Taylor's first-order result R_vac (1 − 4η/3), η = vy0·b/(m g). The next term is +(14/9) η² | §2.4 |
| `quadratic_drag_horizontal_speed(m, c, v0, t)` | v0 / (1 + c·v0·t/m) | §2.4 |
| `quadratic_drag_horizontal_distance(m, c, v0, t)` | (m/c) ln(1 + c·v0·t/m), the coasting distance | §2.4 |
| `cyclotron_omega(q, B, m)`, `cyclotron_radius(m, v_perp, q, B)`, `cyclotron_period(m, q, B)` | Circular motion in a magnetic field | §2.5 |
| `tsiolkovsky_delta_v(v_ex, m0, m)`, `rocket_velocity_gravity(v_ex, m0, burn_rate, g, t)` | Rocket equation, without and with gravity | §3.2 |
| `SHM {omega0, x0, v0}` | `x(t)`, `v(t)`, `amplitude()`, `phase()`, `period()`, `energy(m)` | §5.2 |
| `DampedOscillator {omega0, beta, x0, v0}` | `regime()`, `omega1()`, `x(t)`, `v(t)`, `log_decrement()`, correct in all four regimes | §5.4 |
| `Regime` | `Undamped`, `Underdamped`, `CriticallyDamped`, `Overdamped` | §5.4 |
| `DrivenOscillator {omega0, beta, f0, omega}` | Steady state: `amplitude()`, `phase()`, `resonant_omega()`, `quality_factor()`, `fwhm()` | §5.5–5.6 |

```cpp
double R = analytic::range(20.0, 45.0 * constants::deg);        // 40.79 m
analytic::LinearDragParams ldp{1.0, 0.5};                        // m, b (g defaults; b = 0 is exact vacuum)
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
| `total_energy(m, state, U)` | T + U for any potential-energy function: `U(r)`, or `U(t, r)` if the potential changes with time |
| `energy_drift(traj, U)` | (E(t) − E(0)) / \|E(0)\| for every sample. The best "is my `dt` small enough?" check. With a `U(t, r)` it also contains real physics, ∫∂U/∂t dt |
| `power(F, state)` | P = F·v in watts, for one force at one state |
| `work_done(traj, F)` | The cumulative work W = ∫F·v dt of one force along the run, one value per sample (W[0] = 0), from Simpson's rule on the samples. `F` may depend on v or t |
| `energy_balance(traj, U, F_nc)` | The energy ledger: [(T+U)(t) − (T+U)(0) − W_nc(t) − ∫∂U/∂t dt] / scale, for every sample. Near 0 whenever `-∇U + F_nc` is the force you simulated, even with drag or a drive. The 5-point ∂U/∂t stencil follows the local sample gap, not the absolute timestamp. Pass a zero force for a purely conservative run |
| `uniform_gravity_U(m, g)`, `spring_U(k)` | Ready-made `U` functions to pass into the two above. Each matches one setup only: `uniform_gravity_U` assumes g along −y, and `spring_U` assumes anchor = origin and natural length 0. For anything else write the U (or use `potentials::find_potential`), or `energy_drift` reads order 1 and `dt` gets the blame |
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

Energy with a non-conservative force (drag), and with a potential that changes in time:

```cpp
Particle stone{1.0, {0.0, 0.0, 0.0}, {10.0, 10.0, 0.0}};
forces::UniformGravity weight{stone.mass};
forces::LinearDrag drag{0.5};
SimOptions run;
run.t_end = 2.0;
run.dt = 1e-3;                                                  // record_every stays 1: the ledger wants every step
Trajectory fall = simulate(stone, forces::sum(weight, drag), run);
auto Ug = analysis::uniform_gravity_U(stone.mass, constants::g_earth);

std::vector<double> lost = analysis::energy_drift(fall, Ug);       // about −0.73 at the end: drag took 73% of E(0)
std::vector<double> W = analysis::work_done(fall, drag);           // W.back() ≈ −72.76 J, exactly the energy lost
std::vector<double> ledger = analysis::energy_balance(fall, Ug, drag);   // ≈ 1e-13 everywhere: the books balance
double P = analysis::power(drag, fall.front());                    // −b v² = −100 W at the start

// A spring whose stiffness pulses: U depends on t, so E is not constant even with no drag.
auto pulse = [](double t, const Vec3& r, const Vec3&) { return -4.0 * (1.0 + 0.5 * std::sin(t)) * r; };
auto Ut = [](double t, const Vec3& r) { return 2.0 * (1.0 + 0.5 * std::sin(t)) * norm2(r); };
auto no_force = [](double, const Vec3&, const Vec3&) { return Vec3{}; };
SimOptions pulsing;
pulsing.t_end = 20.0;
pulsing.dt = 1e-3;
Trajectory pulsed = simulate(Particle{1.0, {1.0, 0.0, 0.0}, {}}, pulse, pulsing);
std::vector<double> d = analysis::energy_drift(pulsed, Ut);              // up to 28% of E(0): real, it is ∫∂U/∂t
std::vector<double> b = analysis::energy_balance(pulsed, Ut, no_force);  // about 1e-11: the ∂U/∂t term accounts for it
```

Record every step (`record_every = 1`) when you use `work_done` or `energy_balance`: their accuracy
is O((record_every·dt)⁴), and they are only as good as the samples they are given.

### fields — `newton::fields` (calculus on a field you can write as a function)

A *scalar field* is any callable `double U(const Vec3& r)`, such as a potential energy. A *vector
field* is any callable `Vec3 F(const Vec3& r)`, such as a force that depends on position only. These
functions differentiate and integrate such fields numerically, so you can test an idea without doing
the algebra by hand. A library force also takes `t` and `v`; to use one here, freeze them in a lambda
(`[&](const Vec3& r) { return force(t, r, Vec3{}); }`), as the snippet below does. Passing a callable
of the wrong shape is a compile error with a message that says what shape is needed.

| Name | What it does |
|---|---|
| `default_step(r, length_scale = 1)` | The finite-difference step the functions below use: 1e-4·max(\|r\|, `length_scale`) [m] |
| `gradient(U, r)`, `gradient(U, r, h)` | ∇U of a scalar field, by 5-point central differences |
| `jacobian(F, r)`, `jacobian(F, r, h)` | All nine ∂Fᵢ/∂xⱼ, as `Jacobian {dx, dy, dz}`. The members are columns: `J.dx.y` is ∂Fy/∂x |
| `curl(F, r)`, `divergence(F, r)` | ∇×F and ∇·F at a point, each with an explicit-`h` overload |
| `curl(J)`, `divergence(J)`, `frobenius_norm(J)` | The same quantities read off a `Jacobian` you already have (`curl(J)` and `divergence(J)` are `constexpr`) |
| `line_integral(F, a, b)` | ∫F·dr along the straight segment a → b. Returns NaN if it could not converge |
| `line_integral(F, path)` | The same along a polyline: a `std::vector<Vec3>` of corner points |
| `circulation(F, center, normal, radius)` | ∮F·dr once around a circle, counterclockwise seen from the tip of `normal` (right-hand rule) |
| `line_integral_detailed(F, a, b)`, `circulation_detailed(F, c, n, R)` | The same two integrals, returning the whole `LineIntegral` |
| `LineIntegral {value, error, scale, evaluations, converged}` | The answer; its estimated absolute error; ∫\|F\|\|dr\| (how big the answer could have been); how many times `F` was called; and whether to trust it |

```cpp
auto U = [](const Vec3& r) { return 0.5 * 4.0 * norm2(r); };       // a spring, k = 4 N/m
Vec3 r{1.0, 2.0, 0.5};
Vec3 g = fields::gradient(U, r);                                   // ∇U = k r = (4, 8, 2)

auto swirl = [](const Vec3& p) { return Vec3{-p.y, p.x, 0.0}; };
Vec3 c = fields::curl(swirl, r);                                   // (0, 0, 2)
double d = fields::divergence(swirl, r);                           // 0

// A library force becomes a field once t and v are frozen:
forces::Spring spring{4.0};
auto frozen = [&](const Vec3& x) { return spring(0.0, x, Vec3{}); };
fields::Jacobian J = fields::jacobian(frozen, r);                  // −4 times the identity matrix
double dFx_dx = J.dx.x;                                            // −4
double trace = fields::divergence(J);                              // −12 = −3k

// Atomic scale: the default step assumes lengths of about 1 m and finds nothing in a bond-sized well.
auto morse = [](const Vec3& p) {
    const double s = 1.0 - std::exp(-(norm(p) - 1e-10) / 3e-11);
    return 1e-19 * s * s;                                          // joules
};
Vec3 at{1.2e-10, 0.0, 0.0};
Vec3 wrong = fields::gradient(morse, at);                          // exactly (0, 0, 0): no error, no warning
Vec3 right = fields::gradient(morse, at, fields::default_step(at, 1e-10));   // (1.665e-9, 0, 0) N
```

```cpp
auto swirl = [](const Vec3& p) { return Vec3{-p.y, p.x, 0.0}; };   // curl = 2 ẑ, so work depends on the path
Vec3 a{1.0, 0.0, 0.0}, b{0.0, 1.0, 0.0};
double direct = fields::line_integral(swirl, a, b);                // 1 J
std::vector<Vec3> corner{a, {1.0, 1.0, 0.0}, b};
double around = fields::line_integral(swirl, corner);              // 2 J: same endpoints, different work
double loop = fields::circulation(swirl, Vec3{}, Vec3{0.0, 0.0, 1.0}, 2.0);   // 2πR² = 8π = 25.13 J

fields::LineIntegral I = fields::line_integral_detailed(swirl, a, b);
std::printf("%.6f J, error %.1e, %d evaluations, converged: %s\n", I.value, I.error, I.evaluations,
            I.converged ? "yes" : "no");
```

Two cautions. First, a finite-difference derivative is wrong within about two steps of a kink (`|x|`
at 0) or a singularity (`1/r` at 0), and the default step is 10⁻⁴ of `max(|r|, length_scale)`. Second,
a line integral that passes through a singularity comes back with `converged == false` (or NaN from the
plain-`double` overloads), never with a made-up number, so check it. `notes/06-fields-and-potentials.md`
explains both.

### potentials — `newton::potentials` (is it conservative? and what is U?)

These functions answer Taylor chapter 4's two questions about a force: does it have a potential
energy, and if so, what is it? And the reverse: given a potential energy, what force is that? A
force can be any callable of the shapes `F(t, r, v)` (a library force), `F(t, r)` (a time-dependent
field) or `F(r)` (a static field); the code works out which at compile time.

| Name | What it does |
|---|---|
| `check_conservative(F, opt)` | Tests a force and returns a `ConservativeReport`. It checks velocity and time dependence, then the curl at sampled points (confirmed at a second step size and by a small Stokes loop), then the work around 351 rectangular closed loops per time slice |
| `Verdict` | `Conservative`, `TimeDependentPotential`, `VelocityDependent`, `NonzeroCurl`, `PathDependent` (zero curl but a loop does work: a hole), `Inconclusive` |
| `VelocityPower` | For a `VelocityDependent` force, what its velocity-dependent *part* F(v) − F(0) does: `NoWork` (magnetic), `Dissipative` (drag, friction), `Feeding` (negative damping), `Mixed` (both signs occur), `Unknown` (no usable sample), or `NotApplicable` |
| `ConservativeReport` | `verdict`, `velocity_power`, `depends_on_time`, `depends_on_velocity`, the evidence (`max_relative_curl`, `worst_curl_point`, `worst_curl_time`, `max_loop_ratio`, `worst_loop_center`, `worst_loop_time`, sample and loop counts, and `note`, the reason for an `Inconclusive`), and `has_potential()` / `is_conservative()` |
| `describe(report)`, `to_string(verdict)` | A sentence or two explaining the verdict, or its short name |
| `CheckOptions`, `Box` | Where to test (`region`, a box [m] with positive size in all three directions; a flat box gives `Inconclusive`), the time window (`t_min`, `t_max` [s]), the largest velocity probed (`v_max` [m/s]), and sample counts and tolerances (the defaults are sound) |
| `probe_dependence(F, opt)` | Only the first step: a `Dependence {on_time, on_velocity}` |
| `find_potential(F, opt, r_ref, U_ref)` | Classifies the force and, if it has a potential, returns a `PotentialResult {report, U, r_ref_used, reference_ok, reference_moved, note}`. The reference point is validated by checking that two different paths to each box corner give the same work. A reference you choose that is outside the box, or fails, gives `reference_ok = false`. Check `has_potential()` before using `U`; `describe(result)` explains the reference point |
| `potential_from_force(F, r_ref, U_ref)` | The potential U(r) = U_ref − ∫F·dr from `r_ref`, as a function object: call `U(r)`, or `U(t, r)` for a time-dependent force. It does **not** check that the force is conservative |
| `force_from_potential(U, length_scale, center)` | The force F = −∇U as an ordinary force (`F(t, r, v)`), ready for `simulate`. `U` is `U(r)` or `U(t, r)`. `length_scale` [m] is the same atomic-scale knob as in `fields`; `center` [m] is where U's structure sits, because the step grows with the distance from it (default: the origin) |
| `velocity_independent_part(F)`, `velocity_dependent_part(F)` | Split a force into F(t, r, 0) and F(t, r, v) − F(t, r, 0), each a force. Bookkeeping, not Taylor's classification |

```cpp
potentials::CheckOptions opt;                                      // the box to test: default [−1, 1]³ m
opt.region = {{-2.0, -2.0, -2.0}, {2.0, 2.0, 2.0}};

auto swirl = [](const Vec3& p) { return Vec3{-p.y, p.x, 0.0}; };
potentials::ConservativeReport rep = potentials::check_conservative(swirl, opt);
bool curly = rep.verdict == potentials::Verdict::NonzeroCurl;      // true
std::puts(potentials::describe(rep).c_str());                      // "not conservative: curl F is nonzero at ..."

bool spring_ok = potentials::check_conservative(forces::Spring{4.0}, opt).is_conservative();   // true
auto drag_rep = potentials::check_conservative(forces::LinearDrag{0.5}, opt);
// drag_rep.verdict == Verdict::VelocityDependent, drag_rep.velocity_power == VelocityPower::Dissipative
std::printf("%s: %s\n", potentials::to_string(drag_rep.verdict), potentials::to_string(drag_rep.velocity_power));
```

```cpp
potentials::CheckOptions opt;
opt.region = {{-2.0, -2.0, -2.0}, {2.0, 2.0, 2.0}};

// U from F: classify first, then use the potential.
auto found = potentials::find_potential(forces::Spring{4.0}, opt);
if (found.has_potential()) {
    double U1 = found.U(Vec3{1.0, 0.0, 0.0});                      // 2 J = ½ k r²: the reference is the origin
    double U2 = found.U(Vec3{0.0, 2.0, 0.0});                      // 8 J
}

// Gravity's U is zero at infinity, which a finite integral cannot start from. Anchor it at any regular
// point instead, giving that point its textbook value:
forces::PointGravity gravity{10.0, 1.0};                           // GM = 10 m³/s², m = 1 kg
const Vec3 r_ref{1.0, 0.0, 0.0};
auto Ug = potentials::potential_from_force(gravity, r_ref, -10.0 / norm(r_ref));   // U_ref = −GMm/|r_ref|
double U3 = Ug(Vec3{2.0, 0.0, 0.0});                               // −5 J = −GMm/r

// F from U: a one-line potential, and the force it produces.
auto F = potentials::force_from_potential([](const Vec3& r) { return 2.0 * norm2(r); });   // U = 2 r²
Vec3 f = F(0.0, Vec3{1.0, 0.0, 0.0}, Vec3{});                      // (−4, 0, 0): F = −4 r
// A potential centred away from the origin: say where, so the step follows the distance from the centre.
const Vec3 nucleus{1000.0, 0.0, 0.0};
auto well = [nucleus](const Vec3& r) { return -1.0 / norm(r - nucleus); };           // U = −1/d, units J
auto Fwell = potentials::force_from_potential(well, 1.0, nucleus);                   // 1 m from the well: 2.5e-10 error; without `nucleus`, 4e-4
Vec3 pull = Fwell(0.0, Vec3{1001.0, 0.0, 0.0}, Vec3{});                             // (−1, 0, 0)
SimOptions sim;
sim.t_end = 3.0;
Trajectory run = simulate(Particle{1.0, {1.0, 0.0, 0.0}, {}}, F, sim);   // a spring with k = 4
bool derived_ok = potentials::check_conservative(F, opt).is_conservative();   // true: −∇U is always conservative

// Split a force at rest: gravity stays, drag is the part that exists only because the particle moves.
auto both = forces::sum(forces::UniformGravity{1.0}, forces::LinearDrag{0.5});
auto at_rest = potentials::velocity_independent_part(both);        // F(t, r, 0) = (0, −9.81, 0)
auto moving = potentials::velocity_dependent_part(both);           // F(t, r, v) − F(t, r, 0) = −0.5 v
Vec3 drag_only = moving(0.0, Vec3{}, Vec3{3.0, 0.0, 0.0});         // (−1.5, 0, 0)
```

**What the verdicts mean.** `check_conservative` samples; it does not prove. A `NonzeroCurl`,
`PathDependent` or `VelocityDependent` verdict is a proof that the force is not conservative. A
`Conservative` verdict means "no evidence against it" inside the box and the time window you gave, and
nothing outside them. `describe(report)` quotes the evidence, and `notes/06-fields-and-potentials.md`
explains exactly what was tried and what could be missed. Two more things to know: a potential from
`find_potential` or `potential_from_force` is **NaN** at a point where its integral cannot be
computed (never a silently wrong number), and `U` is only certified inside the box that was tested
(`find_potential` hands the box to `U`, which keeps its integration detours inside it). One limit worth knowing:
a singular line (a "hole") that only clips the box near a face or corner can be missed; over random lines
the measured miss rate is 4.8% across the whole box and 0% for lines through the inner 80%.

### systems — `newton::systems`

| Name | What it is |
|---|---|
| `RocketParams {m0, m_dry, v_ex, burn_rate, g}` | Wet mass, dry mass, exhaust speed, burn rate [kg/s]; `g` defaults to `g_earth`. Require finite `0 < m_dry < m0`, finite `v_ex >= 0`, finite positive `burn_rate`, and finite `g` |
| `RocketSample {t, m, v, y}` | One instant: time, mass, vertical velocity, altitude |
| `simulate_rocket(params, dt, t_end)` | Vertical flight, burn then coast. Burnout and `t_end` are exact samples. `dt` must be finite and positive; `t_end` finite and non-negative; invalid inputs (or an unrepresentable step count) throw `std::invalid_argument` |

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

**Your coords module covers all three layers on paper, but one is thin.** The headers have doc comments,
`test_coords.cpp` is the evidence, and the README's module table now lists it. The demonstration layer is
still missing: no example uses coords yet. `fields` and `potentials` are the model to copy: headers with
doc comments, three test files, two examples (`07` and `08`), and a note (`06-fields-and-potentials.md`).

---

## Part 3 — Everything else in the repo

| Path | What it is | Installed? |
|---|---|---|
| `include/newton/` | The library (Part 1) | **Yes, the only thing installed** |
| `README.md` | The front page: quick start, module tables (one per module, including coords, fields and potentials), install instructions | No |
| `DESIGN.md` | The spec: the hard rules, the Force concept, every file's API | No |
| `Makefile`, `build.sh` | Build the examples and tests, run the tests, install | No |
| `examples/` | Eight lesson programs (table below) | No |
| `tests/` | Ten test programs plus the harness, `test_util.hpp` (Part 4) | No |
| `scripts/plot_csv.py` | Plots any CSV the library writes | No |
| `docs/research/` | Six research reports: the Taylor ch. 1–7 map, numerical methods, C++ library design, and three Chapter 1 reports | No |
| `notes/` | Your study notes: 00–06, `learning-path.md`, `chapters/`, `progress/` | No |
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
| `07_potential_from_force.cpp` | §4.1–4.5, §5.4 | Is this force conservative, and what is its U? Two paths, the curl test, a verdict table, `find_potential`, the hole field, a time-dependent U, and a closed energy ledger for gravity plus drag |
| `08_force_from_potential.cpp` | §4.3–4.5, §5.1, §5.4 | From U to F to motion: a Morse potential through `force_from_potential`, small oscillations against ω = √(U″/m), a moving-support U(t, r), and the ledger with drag |

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

### Five kinds of checks in this repo

| Kind | Example | Where |
|---|---|---|
| Exact identity | a×b is perpendicular to both a and b | `test_vec3.cpp` |
| Simulation vs closed form | A simulated quadratic-drag fall matches v_ter·tanh(gt/v_ter) | `test_numeric_vs_analytic.cpp` |
| Method contract | Halving dt shrinks RK4's error 16×, so it really is 4th order | `test_integrators.cpp` |
| Round trip or invariant | Cartesian → spherical → Cartesian returns the same point | `test_coords.cpp` |
| Identity or exact value from physics | curl(∇U) = 0; the swirl's circulation is 2πR²; E(t) − E(0) = ∫∂U/∂t | `test_fields.cpp`, `test_line_integrals.cpp`, `test_potentials.cpp` |

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

### The current suites (2026-11-25)

| Suite | Checks | Covers |
|---|---|---|
| `analysis` | 10,666 | Energy (including `U(t, r)`), `power`, `work_done` and `energy_balance` against closed forms; momentum; peak/period/amplitude measurement on made-up data |
| `analytic` | 376 | The closed forms agree with themselves (x(0) = x0, each regime, the resonance peak) |
| `coords` | 2,779 | Hand-computed points, round trips, `atan2` vs `acos` |
| `fields` | 171 | Gradient, Jacobian, curl and divergence against closed forms, atomic-scale fields, the identities curl(∇U) = 0 and div(∇×A) = 0 |
| `forces` | 29 | Each force law against a hand-computed vector |
| `integrators` | 29 | Exactness for constant force, convergence order, energy behavior, time bookkeeping, invalid timing inputs and endpoint events |
| `line_integrals` | 102 | Segments, polylines and circulation against exact work; the hole field's 2π; the failure flags |
| `numeric_vs_analytic` | 68 | Simulations against closed forms: drag, oscillators, rocket, magnetic field, stopping events and rocket input validation |
| `potentials` | 469 | The verdict table for every library force, velocity-power labels, time slices, regions (off-centre holes, antiparallel wires), reference validation, potentials against textbook forms, round trips and invalid classifier configurations |
| `vec3` | 33 | Vector algebra identities |

Total: **14,722 checks, all passing** (from `make check`, 2026-11-26). The counts depend on how many cases a
test loops over, so they are a rough measure of size, not of quality.

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

1. **No example uses coords yet.** An `examples/09_coordinates.cpp` would complete coords' demonstration
   layer: for instance, report a projectile's path in spherical coordinates, or show that a uniform
   circular orbit has constant ρ. (The README's module table now lists `coords/`, and examples `07` and
   `08` now exist for fields and potentials.)
2. **Give `DampedOscillator`, `RocketParams` and `RocketSample` default member values**, so they match
   the rest of the library (see "Three things to watch for").
3. **`notes/02-physics-to-code.md` still calls the potentials work "v0.2 roadmap"** (in its Chapter 4
   "Where it lives" and "Examples" sections, which also say no example covers Chapter 4). It is built
   now, so those passages should point to `notes/06-fields-and-potentials.md` and examples `07` and `08`.
