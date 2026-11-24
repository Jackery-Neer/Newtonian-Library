# Learning Path — physics ↔ code

Your running map of the library, in the order you're learning the physics. Each **stage** is
one topic from Taylor: what you learned, then the exact files and line ranges that
*implement* it (forces), *solve* it (integrators, `simulate`), and *check* it (analytic
solutions, tests, examples). Tooling — the test harness, the Makefile, `build.sh` — appears at
the stage where you first need it, so by the end of the course you'll have read the whole
library, not only the physics.

**Status:** ✅ covered in class and the code walked through · 📖 covered in class, read the code now ·
🔜 next in the book · 🔒 later chapter (listed so you can see the whole map)

**Line numbers** were last verified on **2026-10-03**. If a file has changed since, search for the
named function or struct; the name is the stable anchor, the line number is a convenience.

**Updated every Friday at 3 pm** by a scheduled check-in that asks what you covered that week.

---

## Progress log

| Week of | What you covered | What changed in this file |
|---|---|---|
| 2026-09-09 | Finished Taylor ch. 1 (§1.1–1.7). Worked problems 1.18 (triangle area by cross products), 1.30 (perfectly inelastic collision), 1.37 (puck up an incline). | Stages 0–6 |
| 2026-10-02 | Taylor §2.1–2.4: linear drag in 1D and 2D, quadratic drag in 1D, and solving those ODEs. | Stages 7–8 added as 📖; 9–10 queued |

## Your contributions

Code you've written yourself goes here: what it is, where it lives, and which stage it belongs to.

| Date | What you built | Files | Stage |
|---|---|---|---|
| — | *(nothing yet; the "Your turn" items below are good first ones)* | | |

---

## Stage 0 — Orientation: how the repo fits together ✅

**What it is.** A header-only C++17 library: every `.hpp` under `include/newton/` is pasted into
whatever `.cpp` includes it, so there's no library binary to build. Examples and tests are each
one `.cpp` with its own `main()`.

| File | Lines | Look at |
|---|---|---|
| [notes/00-start-here.md](00-start-here.md) | 10–216 | Guided tour: every directory and file type |
| [include/newton/newton.hpp](../include/newton/newton.hpp) | 1–51 | The umbrella header. This is the **manifest**: every public header, grouped by module |
| [DESIGN.md](../DESIGN.md) | 22–54 | The 8 hard rules and the Force concept every force must satisfy |

---

## Stage 1 — Vectors (Taylor §1.2) ✅

**What you learned.** Vector algebra: dot product (projection, angle), cross product (area,
right-hand rule, anticommutative), magnitude, unit vectors.

| File | Lines | Look at |
|---|---|---|
| [core/vec3.hpp](../include/newton/core/vec3.hpp) | 19–30 | `Vec3` as an aggregate: default member initializers, compound operators |
| | 32–46 | `+ - * /` as free functions, so `2.0 * v` and `v * 2.0` both work |
| | 48–79 | `dot`, `cross`, `norm2`, `norm`, `normalized` (zero-guard at 74–79) |
| [tests/test_vec3.cpp](../tests/test_vec3.cpp) | 15–93 | Every check is a §1.2 identity; 57–69 is the cross product |
| [notes/01-cpp-for-this-library.md](01-cpp-for-this-library.md) | 133–218 | Aggregates, default member initializers, `[[nodiscard]]` |
| | 364–402 | `constexpr` functions: why `dot` can be one and `norm` can't |

**Problem tie-in.** 1.18 (triangle area from cross products): the build plan's Proposal 1 turns it
into a unit test.
**Your turn.** [ch01-what-to-build.md](chapters/ch01-what-to-build.md) Proposal 1: add `angle`,
`scalar_triple`, `vector_triple` to `vec3.hpp`, plus their tests.

---

## Stage 2 — Newton's second law becomes a data structure (Taylor §1.4–1.6) ✅

**What you learned.** `F = dp/dt = m r̈`, a second-order ODE. Carry velocity as its own unknown
and it becomes two first-order equations: `ṙ = v`, `v̇ = F/m`.

| File | Lines | Look at |
|---|---|---|
| [core/particle.hpp](../include/newton/core/particle.hpp) | 1–12 | The header derives `ṙ = v`, `v̇ = F/m` from `m r̈ = F` |
| | 19–23 | `State {r, v, t}`: that first-order pair, stamped with a time |
| | 27–36 | `Particle {mass, r, v}`: initial conditions |
| | 38–39 | `momentum()`: `p = m v` |
| [core/trajectory.hpp](../include/newton/core/trajectory.hpp) | 16–24 | A simulation's output: mass + a vector of `State`s |
| [core/constants.hpp](../include/newton/core/constants.hpp) | 1–25 | `pi`, `deg`, `g_earth`, `G`; the `inline constexpr` note at 8–12 |
| [notes/01-cpp-for-this-library.md](01-cpp-for-this-library.md) | 103–132 | Why constants are `inline constexpr` |

---

## Stage 3 — Forces are functions, and forces add (Taylor §1.3) ✅

**What you learned.** The net force is the vector **sum** of every force acting; it never
replaces one force with another.

| File | Lines | Look at |
|---|---|---|
| [forces/gravity.hpp](../include/newton/forces/gravity.hpp) | 1–16 | Why unused parameters are commented out (`/*t*/`) instead of deleted |
| | 30–37 | `UniformGravity`: the simplest force, `F = m g` |
| [forces/drives.hpp](../include/newton/forces/drives.hpp) | 43–50 | `ConstantForce` |
| [forces/compose.hpp](../include/newton/forces/compose.hpp) | 1–38 | `forces::sum(f1, f2, ...)`: superposition as a C++17 fold expression |
| [tests/test_forces.cpp](../tests/test_forces.cpp) | 22–35, 143–149, 171–193 | `UniformGravity`, `ConstantForce`, and `sum()` tests |
| [DESIGN.md](../DESIGN.md) | 41–54 | The Force concept: `Vec3 operator()(double t, const Vec3& r, const Vec3& v) const` |
| [notes/01-cpp-for-this-library.md](01-cpp-for-this-library.md) | 219–327 | Template callables vs `std::function` (the central design decision), and fold expressions |

---

## Stage 4 — Solving the ODE numerically (Taylor §1.6 + numerics) ✅

**What you learned.** On paper you integrate `v̇ = F/m` twice. In code an **integrator** takes
tiny steps instead: you write the force, and the integrator does the solving.

| File | Lines | Look at |
|---|---|---|
| [integrators/euler.hpp](../include/newton/integrators/euler.hpp) | 1–12, 57–62 | The update equations, then the code: the most literal possible step |
| [integrators/semi_implicit_euler.hpp](../include/newton/integrators/semi_implicit_euler.hpp) | 56–63 | Compare with Euler: **one line** differs (old `v` vs new `v`) |
| [integrators/rk4.hpp](../include/newton/integrators/rk4.hpp) | 81–125 | Four force evaluations per step, plus the collapsed position update |
| [integrators/velocity_verlet.hpp](../include/newton/integrators/velocity_verlet.hpp) | 100–114 | Kick, drift, kick |
| [simulate.hpp](../include/newton/simulate.hpp) | 1–19 | What the driver is for |
| | 33–54 | `Method` and `SimOptions`: the knobs you set |
| | 83–132 | `simulate()`, and why time is `t0 + i*dt`, never `t += dt` |
| | 134–214 | `simulate_until()`: stop on an event; bisection refinement at 185–200 |
| | 216–233 | `events::ground_impact`, `events::apex` |
| [tests/test_integrators.cpp](../tests/test_integrators.cpp) | 1–12, 61–107 | Overview, then (a): RK4 and Verlet are exact for constant force, Euler isn't |

**Skip for now.** The rest of each integrator's header comment justifies the method with energy
and "symplectic" arguments. That's Chapter 4 vocabulary; read the code now and the justification then.

---

## Stage 5 — Constant force, checked against the closed form (Taylor §1.6) ✅

**What you learned.** Constant force gives a degree-2 polynomial in `t`, and the vacuum projectile
has range `v₀² sin 2θ / g`.

| File | Lines | Look at |
|---|---|---|
| [analytic/projectile.hpp](../include/newton/analytic/projectile.hpp) | 1–62 | The library's first **analytic oracle**: your textbook formulas as functions |
| [examples/01_free_fall.cpp](../examples/01_free_fall.cpp) | 1–23 | How every example opens: situation, expectation, what to look at |
| | 39–70 | Drop from rest, stopped at the ground by `simulate_until` |
| | 71–116 | 45° projectile, numeric vs analytic (the comment at 109–111 uses Chapter 4 energy language; defer it) |
| [tests/test_analytic.cpp](../tests/test_analytic.cpp) | 124–150 | Projectile formulas agree with each other |
| [tests/test_numeric_vs_analytic.cpp](../tests/test_numeric_vs_analytic.cpp) | 200–228 | Simulated range matches the closed form to 1e-9 |

**Run it:** `./build/01_free_fall`
**Problem tie-in.** 1.37 (puck up a frictionless incline): build plan Proposal 3 turns it into an example.

---

## Stage 6 — Momentum of a system (Taylor §1.5) ✅ (only partly in code)

**What you learned.** `P = Σ mᵢvᵢ`. Internal forces cancel in pairs (Newton III), so with no
external force `P` is constant.

| File | Lines | Look at |
|---|---|---|
| [analysis/momentum.hpp](../include/newton/analysis/momentum.hpp) | 1–21 | Notes that v0.1 can't simulate several bodies yet |
| | 44–64 | `center_of_mass`, `total_momentum` (skip 28–42, angular momentum is Chapter 3) |
| [tests/test_analysis.cpp](../tests/test_analysis.cpp) | 138–162 | Two-particle sanity checks |

**Problem tie-in.** 1.30 (perfectly inelastic collision): put the "before" and "after" particles
into `total_momentum` and check they agree.
**Gap.** Nothing simulates two interacting bodies yet. That's build plan Proposal 9, the biggest one.

---

## Stage 7 — Linear drag in 1D and 2D (Taylor §2.1–2.3) 📖 now

**What you learned.**
```
f = −b v                     τ = m/b         v_ter = m g / b
horizontal:  vx(t) = vx0 e^(−t/τ)            x → vx0 τ   (a finite distance)
vertical:    vy(t) = −v_ter + (vy0 + v_ter) e^(−t/τ)
the x and y equations decouple, so 2D linear drag has a closed form, including y(x) with a log term
```

| File | Lines | Look at |
|---|---|---|
| [core/constants.hpp](../include/newton/core/constants.hpp) | 27–36 | Taylor's sphere-in-air coefficients β and γ |
| [forces/drag.hpp](../include/newton/forces/drag.hpp) | 1–15 | The `f(v) = bv + cv²` model |
| | 22–34 | `LinearDrag`: the whole force is `return -b * v;` |
| | 55–61 | `linear_drag_sphere_air(D)`: `b = βD` |
| [analytic/drag.hpp](../include/newton/analytic/drag.hpp) | 1–21 | Why linear drag has closed forms and quadratic drag mostly doesn't |
| | 25–40 | `LinearDragParams`: `tau()`, `v_ter()` |
| | 42–58 | `linear_drag_velocity`: your horizontal and vertical results as **one vector formula** |
| | 60–76 | `linear_drag_position` |
| | 78–100 | `linear_drag_trajectory_y`: the exact `y(x)`, with the asymptote at `x = vx0 τ` |
| [tests/test_forces.cpp](../tests/test_forces.cpp) | 59–67, 82–88 | `LinearDrag` and the sphere-coefficient tests |
| [tests/test_numeric_vs_analytic.cpp](../tests/test_numeric_vs_analytic.cpp) | 16–44 | Your hand solution used as the **answer key** for a simulation |
| [integrators/velocity_verlet.hpp](../include/newton/integrators/velocity_verlet.hpp) | 27–70 | Why a velocity-dependent force breaks plain Verlet, and how this library handles it |
| [examples/02_projectile_drag.cpp](../examples/02_projectile_drag.cpp) | 1–25 | Taylor's baseball under vacuum, linear, and quadratic drag |
| | 41–55, 57–98 | The `launch()` helper, then the three force models built with `forces::sum` |
| [io/csv.hpp](../include/newton/io/csv.hpp) | 22–75 | Two `write_csv` overloads: one for a `Trajectory`, one for generic columns (what example 02 uses) |
| [notes/02-physics-to-code.md](02-physics-to-code.md) | 87–174 | Chapter 2 map (skip the electromagnetic parts for now) |

**Tooling unlocked: building your own code.** You're about to write code, so read how it gets built:

| File | Lines | Look at |
|---|---|---|
| [notes/00-start-here.md](00-start-here.md) | 76–121, 235–275 | What `make` does; the compile command flag by flag |
| [Makefile](../Makefile) | 1–56 | Configuration, wildcard file discovery, `.PHONY` |
| | 57–79 | Pattern rules: how each `.cpp` becomes its own binary |
| | 86–91 | `check` and `clean` |

**Run it:** `./build/02_projectile_drag`, then `python3 scripts/plot_csv.py build/projectile_drag.csv`

**Your turn.**
1. In a `sandbox/` folder (the Makefile ignores it), write `linear_drag_velocity` yourself from your
   notes *without looking*, then print it beside `simulate()` and the library's version.
2. Nothing tests `linear_drag_velocity` or `linear_drag_position` against a simulation yet. Write
   that test, using lines 16–44 as the template.

---

## Stage 8 — Quadratic drag in 1D (Taylor §2.4) 📖 now

**What you learned.**
```
f = −c |v| v                 v_ter = sqrt(m g / c)
horizontal only:   v(t) = v0 / (1 + c v0 t / m)       power-law decay, x grows like ln t
falling from rest: v(t) = v_ter tanh(g t / v_ter)     y(t) = (v_ter²/g) ln cosh(g t / v_ter)
```

| File | Lines | Look at |
|---|---|---|
| [forces/drag.hpp](../include/newton/forces/drag.hpp) | 36–53 | `QuadraticDrag`, and why it's written `-c|v|v` rather than `-c|v|² v̂` |
| | 62–67 | `quadratic_drag_sphere_air(D)`: `c = γD²` |
| [analytic/drag.hpp](../include/newton/analytic/drag.hpp) | 101–108 | Terminal speed |
| | 110–128 | The `tanh` speed and `ln cosh` distance for a fall from rest |
| | 130–142 | Horizontal speed; the doc comment writes out an `x(t)` that **no function implements** |
| [tests/test_forces.cpp](../tests/test_forces.cpp) | 68–81 | Quadratic drag is exactly zero, not NaN, at `v = 0` |
| [tests/test_numeric_vs_analytic.cpp](../tests/test_numeric_vs_analytic.cpp) | 45–79 | Fall from rest vs the `tanh` solution, plus the approach to terminal velocity |

**Tooling unlocked: the test harness.** You'll write tests here, so read the framework, which is
about 100 lines:

| File | Lines | Look at |
|---|---|---|
| [tests/test_util.hpp](../tests/test_util.hpp) | 1–16 | What a "test framework" really needs |
| | 17–72 | Counters, tolerance-based `close()`, `record()`, `report()` |
| | 74–103 | Why `NT_CHECK` and `NT_CLOSE` have to be macros |

**Your turn.** Add `quadratic_drag_horizontal_distance` to `analytic/drag.hpp`, then write one test
block checking both the horizontal speed and the distance against `simulate()`. Before you start:
which forces go in the `sum`, and which sign does your velocity have?

---

## Stage 9 — Quadratic drag in 2D (Taylor §2.4) 🔜 next

**Why it matters.** The x and y equations are coupled through `|v|`, so there's no closed form. This
is the book's own argument for numerical methods. A simulation is the only way to get the answer.

| File | Lines | Look at |
|---|---|---|
| [examples/02_projectile_drag.cpp](../examples/02_projectile_drag.cpp) | 12–19 | What *can* still be checked without a closed form |
| | 99–124 | Brute-force sweep of launch angle: the best angle comes out below 45° |

---

## Stage 10 — Charged particle in a magnetic field (Taylor §2.5–2.7) 🔒 next in the book

*Unknown whether your class covers §2.5–2.7 or skips it. Tell the Friday check-in, and this stage
will be promoted or marked skipped.*

| File | Lines | Look at |
|---|---|---|
| [forces/electromagnetic.hpp](../include/newton/forces/electromagnetic.hpp) | 1–49 | `MagneticForce` (`q v × B`), `LorentzForce` |
| [analytic/magnetic.hpp](../include/newton/analytic/magnetic.hpp) | 1–52 | Cyclotron frequency, radius, period |
| [examples/03_charged_particle.cpp](../examples/03_charged_particle.cpp) | all | The helix, and the E×B drift |
| [tests/test_forces.cpp](../tests/test_forces.cpp) | 150–170 | The magnetic force is always perpendicular to `v` |
| [tests/test_numeric_vs_analytic.cpp](../tests/test_numeric_vs_analytic.cpp) | 161–199 | Speed conserved to 1e-12, and the cyclotron period |

---

## Later chapters 🔒: the rest of the map

| Ch. | Topic | Where it lives |
|---|---|---|
| 3 | Rockets (variable mass), angular momentum | `systems/rocket.hpp`, `analytic/rocket.hpp`, `examples/04_rocket.cpp`, `analysis/momentum.hpp` 28–42, `tests/test_numeric_vs_analytic.cpp` 137–160 |
| 4 | Energy | `analysis/energy.hpp`, `examples/06_integrator_comparison.cpp`, `tests/test_integrators.cpp` 108–210, `notes/03-numerical-methods.md`, the energy discussion in each integrator's header |
| 5 | Oscillations | `forces/spring.hpp`, `forces/drives.hpp` 17–41, `analytic/oscillator.hpp`, `analysis/oscillation.hpp`, `examples/05_oscillations.cpp`, `tests/test_analytic.cpp` 13–123, `tests/test_numeric_vs_analytic.cpp` 80–136 |
| 7 | Lagrange's equations (your course skips ch. 6) | **Not built yet. This is your original Lagrangian idea.** The plan is in `notes/02-physics-to-code.md` 386–442 and the Lagrangian-hook section of `notes/04-extending.md`. Stages for it get added when you reach it, and they'll mostly be things for you to *build* |

**Tooling still to come.**
- [build.sh](../build.sh) (all 100 lines): the bash version of `make check`. Read it whenever you
  want to learn shell scripting.
- [Makefile](../Makefile) 93–103 (`make debug`, sanitizers): the first time your code crashes or
  behaves strangely.
- [Makefile](../Makefile) 105–116 and [notes/00-start-here.md](00-start-here.md) 276–327: installing,
  so `#include <newton/newton.hpp>` works from anywhere.
- [DESIGN.md](../DESIGN.md) 55 to the end: the file-by-file spec. Read it the first time you add a whole new module.

---

## Maintenance rules (for the weekly update)

1. Add one Progress log row per week, in your own words (no text copied from Taylor).
2. Move stages forward: 🔜 or 🔒 becomes 📖 when it's covered in class, and 📖 becomes ✅ once the
   code has been read or used.
3. A new stage needs **What you learned**, a **file + lines + what to look at** table, and **Your
   turn** items that only use physics already covered.
4. Re-check line numbers with `grep -n` for any file that changed since the last verified date, and
   update that date.
5. Log anything you've built under **Your contributions**.
6. Never put physics from a chapter you haven't reached into a 📖 stage. Homework problems get
   guidance, not answers.
