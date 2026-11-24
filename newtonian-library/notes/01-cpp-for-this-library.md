# 01 — The C++ this library actually uses

You know loops, classes, and basic templates already, so none of that is
here. This note covers every feature *beyond* that baseline which the
library's headers lean on — each with a real snippet copied from this
codebase (file path included) rather than a made-up example, so you can go
find the whole file afterward and recognize it immediately.

Read `notes/00-start-here.md` first if you haven't — it explains the
directory layout these paths live in.

## `#pragma once`

Every header in `include/newton/` opens with this line, e.g.
`include/newton/core/vec3.hpp:1`:

```cpp
#pragma once
```

Problem it solves: if `a.hpp` and `b.hpp` both `#include <newton/core/
vec3.hpp>`, and some `.cpp` includes both `a.hpp` and `b.hpp`, then without
protection `vec3.hpp`'s contents get pasted into that `.cpp` *twice* —
`struct Vec3` gets defined twice in one translation unit, a hard compile
error. `#pragma once`, placed at the top of a header, tells the compiler
"if you've already opened this exact file once in this translation unit,
skip it silently the second time." It is a compiler extension rather than
standard C++ — but every compiler you're likely to ever use (clang, gcc,
MSVC) supports it, which is why this library uses it instead of the
standard-but-wordier alternative, manual include guards:

```cpp
#ifndef NEWTON_CORE_VEC3_HPP
#define NEWTON_CORE_VEC3_HPP
// ... contents ...
#endif
```

Same effect, three lines instead of one, and a real failure mode of its own
(copy a header, forget to rename the macro, and the copy silently compiles
to nothing).

## Header-only, the ODR, and why so much is `inline`

This library ships **zero** `.cpp` files under `include/` — every function
body lives directly in a header. That's only legal because of a C++ rule
called the **One Definition Rule (ODR)**: a normal (non-template,
non-`inline`) function or variable must be defined in **exactly one**
place across the whole linked program, or the linker fails with "duplicate
symbol." If `include/newton/io/csv.hpp`'s `write_csv` were declared without
`inline` and two `.cpp` files both included that header, you'd get exactly
that error the moment you tried to link both into one program.

`inline` changes the rule for that one entity: it becomes legal to define
it *identically* in every translation unit that includes it, and the linker
merges the duplicates into one instead of complaining. Concretely, from
`include/newton/io/csv.hpp`:

```cpp
[[nodiscard]] inline bool write_csv(const std::string& path, const Trajectory& traj) {
```

Contrast that with `include/newton/integrators/euler.hpp`'s `euler_step`,
which has **no** `inline` keyword at all:

```cpp
template <class Force>
[[nodiscard]] State euler_step(const State& s, double mass, double dt, Force&& force) {
```

That's not an inconsistency — **templates already get the "may be defined
in multiple TUs" exemption for free**, without needing `inline` written out.
A template isn't compiled into code until it's instantiated with concrete
types, and the compiler emits each instantiation as a mergeable symbol on
its own. So the rule to internalize for this codebase is: *every
namespace-scope function in a header is either a template (no `inline`
needed) or a plain function (`inline` required)*. `[[nodiscard]] inline
auto` free functions like `write_csv` fall in the second bucket;
`euler_step`, `rk4_step`, `total_energy`, and everything else parameterized
on `Force`/`PotentialFn` fall in the first.

One more `inline`-adjacent case worth knowing: a member function defined
**inside** a class body is implicitly `inline`, no keyword needed — that's
why `Vec3`'s compound-assignment operators in `include/newton/core/
vec3.hpp` don't say `inline` even though they're plain (non-template)
functions:

```cpp
struct Vec3 {
    double x{0.0};
    double y{0.0};
    double z{0.0};

    constexpr Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    ...
};
```

(They're also `constexpr`, which is its own, separate reason they'd be
implicitly `inline` even if they weren't member functions — see the
`constexpr` section below.)

## `inline constexpr` variables (C++17)

`include/newton/core/constants.hpp`:

```cpp
namespace newton::constants {

inline constexpr double pi = 3.14159265358979323846;
inline constexpr double deg = pi / 180.0;
inline constexpr double g_earth = 9.80665;
inline constexpr double G = 6.67430e-11;
```

Before C++17, a shared constant at namespace scope in a header was
genuinely awkward: a `const`- (or `constexpr`-) qualified namespace-scope
variable has **internal linkage** by default in C++, meaning every
translation unit that includes the header silently gets its **own private
copy** of `g_earth`, at its own address. That compiles and links, but it's
a trap the moment you rely on identity — `&newton::constants::g_earth`
would differ between two `.cpp` files — and it's real (if tiny) duplication
in every object file.

`inline` on a *variable* is a C++17 feature specifically added to fix this:
it's the same relaxation of the ODR that lets `inline` functions be defined
in multiple TUs, extended to variables. `inline constexpr double g_earth`
means: one object, one address, for the entire program, no matter how many
`.cpp` files include `constants.hpp`. This is the modern, typed,
namespaced, debugger-visible replacement for the old C habit of `#define
G_EARTH 9.80665`.

## Aggregates and default member initializers

`Vec3` (`include/newton/core/vec3.hpp`) has no constructor at all:

```cpp
struct Vec3 {
    double x{0.0};
    double y{0.0};
    double z{0.0};
    ...
};
```

A type with all-public data members, no user-declared constructors, no
private/protected members, and no virtual functions is an **aggregate**,
and aggregates get a form of initialization constructors don't: brace-init
against the member list in order, `Vec3{1.0, 2.0, 3.0}`, or `Vec3{}` /
plain `Vec3 v;` to fall back to each member's own default. The `{0.0}`
after each member name — a **default member initializer** — is what makes
`Vec3 v;` and `Vec3{}` produce `{0, 0, 0}` instead of three uninitialized
(and, for `double`, potentially garbage/NaN-producing) values. That single
detail is a whole bug class closed off for free: nothing in this library
can accidentally carry an uninitialized position or velocity into a force
evaluation.

`Particle` (`include/newton/core/particle.hpp`) is the same pattern, and
shows *why* the specific default matters, not just that one exists:

```cpp
struct Particle {
    double mass{1.0};  ///< mass [kg], must be > 0
    Vec3   r{};        ///< initial position [m]
    Vec3   v{};        ///< initial velocity [m/s]
    ...
};
```

`mass` defaults to `1.0`, not `0.0` — a `0.0` default would make `a = F/m`
evaluate to infinity (or NaN) the instant you forgot to set a mass, which is
a much nastier failure than "everything moves as if `mass` were 1," the
`1.0`-default failure mode. Every force functor in `newton::forces`
(`include/newton/forces/*.hpp`) follows the identical shape — e.g.
`UniformGravity` in `include/newton/forces/gravity.hpp`:

```cpp
struct UniformGravity {
    double mass{1.0};                                  // [kg]
    Vec3   g{0.0, -constants::g_earth, 0.0};           // [m/s^2]
    ...
};
```

— which is what lets you write `UniformGravity{}` (Earth gravity, unit
mass) or `UniformGravity{ball.mass}` (override just the mass, keep the
default field) with no constructor ever written by hand.

## `[[nodiscard]]`

`DESIGN.md` mandates this on every pure function in the library, and it's
on essentially every declaration you'll open. Two examples worth contrasting
because they catch two *different* kinds of bug:

`include/newton/core/vec3.hpp`:
```cpp
[[nodiscard]] constexpr Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
```
Ignoring the result of `a + b` is almost certainly a typo for `a += b`;
`[[nodiscard]]` turns that typo into a compiler warning instead of a
silent no-op.

`include/newton/io/csv.hpp`:
```cpp
[[nodiscard]] inline bool write_csv(const std::string& path, const Trajectory& traj) {
```
This one is the more instructive case: `write_csv`'s return value isn't
incidental math you might reasonably ignore, it's a **success/failure
flag** — `false` means the file never got written. Call `write_csv(...)`
as a bare statement and ignore the `bool`, and you've silently swallowed an
I/O failure; `[[nodiscard]]` makes the compiler flag that at the call site
instead of you discovering it later as "why is this CSV empty."

The mechanical rule this library follows: `[[nodiscard]]` goes on every
function whose only job is to *compute and return* something (as opposed to
a function whose job is a side effect, where the return value really is
optional). Almost every function in `newton` is the former.

## Template callables vs. `std::function` — the central design decision

This is the one choice that shapes almost every function signature in the
library. A "force" is a value the *caller* supplies to an integrator, and
C++ gives you two fundamentally different ways to accept "a value the
caller supplies that you can call like a function":

**Type erasure (`std::function`)** — store *any* matching callable behind
one concrete type, `std::function<Vec3(double, Vec3, Vec3)>`, at the cost of
an indirect (virtual-dispatch-like) call and a possible heap allocation if
the callable's captured state doesn't fit in the small-buffer optimization.

**Template callables** — make the integrator itself a template on the
callable's *type*, so the compiler knows exactly what it's calling at
every call site and can inline it. This is what every integrator step
function in this library does; `include/newton/integrators/euler.hpp`:

```cpp
template <class Force>
[[nodiscard]] State euler_step(const State& s, double mass, double dt, Force&& force) {
    const Vec3 a = force(s.t, s.r, s.v) / mass;  // evaluated once, at the OLD (r, v)
    return State{s.r + dt * s.v, s.v + dt * a, s.t + dt};
}
```

`newton` picked template callables everywhere — condensed from the
tradeoff table in `docs/research/cpp-library-design.md` §5:

| | `std::function` | Template callable (`Force&&`) |
|---|---|---|
| Call in the hot loop | Indirect call; the user's lambda body typically can't be inlined | Exact type known at compile time; usually inlined to a few instructions |
| Allocation | May heap-allocate if captures exceed the small-buffer size | Never allocates |
| Measured cost | Often 2–10× slower for a tiny force called millions of times | Full speed |
| Storable in a `std::vector` of *different* force types | Yes, directly | No — needs `std::variant`, a base class, or `std::function` after all |
| Header cost | Needs `<functional>` | Nothing extra |
| Can a caller still pass a `std::function` in? | — | Yes: `std::function` itself satisfies the `Force` concept, so it's a strict subset, not a competing choice |

The last row is why this isn't a closed door: a `Force` template parameter
accepts *anything* callable as `Vec3(double, const Vec3&, const Vec3&)` —
a lambda, a `newton::forces` functor, or a `std::function` someone
constructed for their own reasons — so a caller who genuinely needs runtime
polymorphism (swap the force at runtime, store heterogeneous forces in one
container) can still hand in a `std::function`; they just pay its cost only
when they've chosen to.

One C++ detail worth being precise about: `Force&&` here is a **forwarding
reference**, not an rvalue-reference overload, *because* it's a deduced
template parameter (`template <class Force> ... Force&& force`) — see
`euler.hpp`'s own comment on this:

> `Force&&` is a forwarding reference on a template parameter (not an
> rvalue-reference overload) — it binds to lambdas, functors, and
> `std::function` alike, and because the concrete callable type is known at
> compile time the compiler can fully inline `force(...)` into the step
> instead of going through indirect/virtual dispatch.

If `Force` weren't a template parameter — say, this were a plain function
taking `Vec3(*)(...)&&` — `&&` would mean "rvalue reference," binding only
to temporaries. Because `Force` is *deduced* at the call site, `Force&&`
instead means "bind to anything, and deduce whether it was an lvalue or
rvalue" — that's the forwarding-reference rule (sometimes called a
"universal reference"), and it's what lets the same `euler_step` accept a
named `LinearDrag{0.4}` variable, a temporary `LinearDrag{0.4}`, and a
lambda, with no overloads to write.

## Variadic templates and fold expressions — `forces::sum`

`include/newton/forces/compose.hpp` is the sharpest example of a C++17
fold expression in the whole library:

```cpp
template <class... Fs>
[[nodiscard]] auto sum(Fs... fs) {
    return [fs...](double t, const Vec3& r, const Vec3& v) -> Vec3 {
        return (Vec3{} + ... + fs(t, r, v));
    };
}
```

`class... Fs` is a **template parameter pack** — `Fs` isn't one type, it's
a placeholder for however many types the caller supplies, so `sum(gravity,
drag)` instantiates `Fs = {UniformGravity, QuadraticDrag}` and `sum(gravity,
drag, wind)` instantiates a different, three-type `Fs`. `Fs... fs` is the
matching **function parameter pack**: `fs` is that many actual parameters,
one per type in `Fs`.

The body, `(Vec3{} + ... + fs(t, r, v))`, is a **fold expression** — new in
C++17, purpose-built for exactly this "combine a whole parameter pack with
one operator" pattern. Read left to right, this specific form is `(init op
... op pack-expression)`, a **binary left fold**, which the C++ standard
defines as expanding to:

```
((((Vec3{} + f1(t,r,v)) + f2(t,r,v)) + f3(t,r,v)) + ... )
```

i.e. it left-associates `+` over every force's result, seeded with
`Vec3{}` (the zero vector) so that `sum()` called with *zero* forces is
still well-formed, returning the zero force rather than being a special
case to guard against.

Two more things worth noticing in that snippet, both covered in their own
sections below: `[fs...]` is a **lambda capturing a parameter pack by
value** — the returned closure owns its own copies of every force, so it
stays valid after whatever temporaries were passed into `sum(...)` go out
of scope — and the `[[nodiscard]] auto` return type is unavoidable here
because a lambda's closure type has **no name you could write down** even
if you wanted to; see the `auto` returns section at the end of this note.

## Structured bindings

Not used inside `include/newton/` itself, but directly enabled by it: every
aggregate in the library — `State`, in particular — decomposes cleanly,
because structured bindings work on any aggregate (also on `std::pair`,
`std::tuple`, and C arrays, but plain aggregates are the common case here).
`State` is defined in `include/newton/core/particle.hpp` as:

```cpp
struct State {
    Vec3   r{};      ///< position [m]
    Vec3   v{};      ///< velocity [m/s]
    double t{0.0};   ///< time [s]
};
```

so, given a `Trajectory` (`include/newton/core/trajectory.hpp`, whose
`back()` returns `const State&`), you can write:

```cpp
#include <newton/newton.hpp>

const auto& [r, v, t] = traj.back();  // decomposes State's three members, in declaration order
```

instead of `traj.back().r`, `.v`, `.t` at every use. The names `r`, `v`,
`t` are yours to choose — they don't have to match `State`'s actual member
names — but the *order* is fixed to declaration order, so `auto [v, r, t] =
...` would silently bind `v` to the position and `r` to the velocity: a
real footgun worth remembering the one time it bites you.

One caveat if you go looking to use this inside a lambda: **structured
bindings cannot be captured by a lambda in C++17** (that's fixed in C++20).
If you hit "variable cannot be implicitly captured", copy the structured
binding into an ordinarily-declared variable first, then capture that.

## `constexpr` functions

A `constexpr` function can be evaluated **at compile time** when every
argument it's given is itself a compile-time constant (and still runs as an
ordinary function at runtime otherwise — it's not an either/or). It's also
implicitly `inline`, which is a nice bonus for a header-only library but not
the main reason to reach for it here: `constexpr` is a statement about
*what the function is allowed to depend on*.

`include/newton/core/vec3.hpp` draws the line precisely, in its own header
comment:

> Everything here is `constexpr` where the math allows (`std::sqrt` is not
> `constexpr` in C++17, so `norm()` and friends are plain `inline`).

Concretely:

```cpp
[[nodiscard]] constexpr double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
[[nodiscard]] constexpr double norm2(const Vec3& v) { return dot(v, v); }

[[nodiscard]] inline double norm(const Vec3& v) { return std::sqrt(norm2(v)); }
```

`dot` and `norm2` are pure arithmetic on `double`s, so they're `constexpr`.
`norm` calls `std::sqrt`, which C++17 does not guarantee is usable in a
constant expression (`<cmath>` functions only start becoming `constexpr` in
C++26, and this library targets C++17), so `norm` drops back to plain
`inline`. The same split
shows up again in `include/newton/analytic/oscillator.hpp`'s `SHM`: `period()`
and `energy()` are `constexpr` (pure arithmetic — no trig, no `sqrt`),
while `x(t)` and `v(t)`, which call `std::cos`/`std::sin`, are plain
`inline` member functions. This isn't inconsistency across the file; it's
the same rule — *"constexpr where the underlying math allows it, plain
inline the moment a `<cmath>` transcendental function is involved"* —
applied uniformly.

## Namespaces, and why `using namespace` never appears in a header

Every module gets its own namespace, matching its directory:
`newton` (core), `newton::forces`, `newton::integrators`, `newton::analytic`,
`newton::analysis`, `newton::systems`, `newton::io`. C++17 lets you write a
nested namespace on one line — every header in this library opens its
namespace this way, e.g. `include/newton/forces/gravity.hpp`:

```cpp
namespace newton::forces {
```

instead of the pre-C++17 `namespace newton { namespace forces { ... } }`.
Every header also **closes** its namespace with a comment repeating the
full name — `DESIGN.md` rule 8, followed everywhere:

```cpp
}  // namespace newton::forces
```

This is pure readability discipline (the compiler doesn't care), but it
pays for itself the moment a namespace block runs longer than a screen —
without the comment, scrolling to a lone `}` gives you no way to tell which
of several nested scopes just closed.

The rule with real teeth, though, is: **no `using namespace` at namespace
scope in any header, ever** (`DESIGN.md` rule 4) — and it's not a stylistic
nicety followed loosely; running `grep -rn "using namespace" include/` from
the repo root turns up **zero** matches across the entire library. The
reason is that `#include` is a textual paste (see `notes/00-start-here.md`
if that's new to you) — a `using namespace std;` sitting at namespace scope
in *any* header gets pasted into *every* translation unit that transitively
includes that header, with no way for the includer to opt out. Concretely,
that can silently change what an *unrelated* line of a user's own code
means: they write a variable or function named `count`, `distance`, or
`swap` (all names that live in `std`), and now there's an ambiguity — or
worse, a silent overload-resolution change — that traces back to a header
they never opened and have no reason to suspect.

What the library does instead:

- **Fully qualifies everything in headers** — `std::sqrt`, `std::vector`,
  `constants::g_earth`, throughout `include/newton/`.
- Reserves `using namespace newton;` for **leaf translation units** —
  `.cpp` files that nothing else `#include`s, i.e. exactly the `examples/`
  and `tests/` files, where it's genuinely harmless (nothing downstream can
  be affected) and, as a bonus, shows the API reading cleanly unqualified.
  The quick-start snippet in `README.md` does exactly this.

## `auto` return types — for lambda factories

`auto` as a return type is used sparingly and deliberately in this library:
only where the returned type is **unnameable**, never on the main
documented API. Every occurrence returns a closure — a lambda — and a
lambda's closure type has no spellable name in C++17 (each lambda
expression generates its own unique, compiler-invented class), so `auto` is
not a style choice there, it's the only option. Four real examples, all
"factory functions that build and return a force- or event-shaped
callable":

`include/newton/forces/compose.hpp`:
```cpp
template <class... Fs>
[[nodiscard]] auto sum(Fs... fs) {
    return [fs...](double t, const Vec3& r, const Vec3& v) -> Vec3 { ... };
}
```

`include/newton/simulate.hpp`, the two ready-made stop conditions for
`simulate_until`:
```cpp
[[nodiscard]] inline auto ground_impact(double y_level = 0.0) {
    return [y_level](const State& s) { return s.r.y - y_level; };
}
[[nodiscard]] inline auto apex() {
    return [](const State& s) { return s.v.y; };
}
```

`include/newton/analysis/energy.hpp`, the two ready-made potential-energy
functions for `total_energy`/`energy_drift`:
```cpp
[[nodiscard]] inline auto uniform_gravity_U(double mass, double g) {
    return [mass, g](const Vec3& r) { return mass * g * r.y; };
}
```

Notice every one of these still has `[[nodiscard]]` — `auto` and
`[[nodiscard]]` are independent choices; deducing the return type doesn't
excuse discarding it. The rule this library follows for when `auto` is and
isn't appropriate, from `docs/research/cpp-library-design.md` §5: reach for
it on factory functions returning lambdas and on template-heavy internals,
where there's genuinely no alternative or the type is obvious from a
one-line body with a visible `return`; avoid it on the ordinary public
API — a function like `analytic::range(v0, theta, g)` is written to return
a plain `double`, spelled out, precisely so the signature documents itself
without forcing a reader into the function body.
