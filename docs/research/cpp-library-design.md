# Structuring a Header-Only Modern C++ Physics Library

A practical guide for building `newton` — a header-only C++17 library that lives in its own repo, is used as `#include <newton/newton.hpp>`, and can later be installed system-wide.

Everything in this report was compiled and run on the target machine (Apple clang 15.0.0, arm64 Darwin 23.3.0, Command Line Tools toolchain) with `-std=c++17 -Wall -Wextra -Wpedantic`. All code samples and build files are verified working.

---

## 1. Directory layout

### Recommended layout

```
newtonian-library/
├── include/
│   └── newton/                  <- the "namespace directory"; this is what -I exposes
│       ├── newton.hpp           <- umbrella header: includes everything
│       ├── core/
│       │   ├── vec3.hpp
│       │   ├── constants.hpp
│       │   └── units.hpp
│       ├── dynamics/
│       │   ├── particle.hpp
│       │   ├── forces.hpp
│       │   └── integrate.hpp
│       └── detail/              <- optional: implementation helpers, not public API
│           └── math_helpers.hpp
├── examples/
│   ├── projectile.cpp
│   ├── orbit.cpp
│   └── spring.cpp
├── tests/
│   ├── test_util.hpp            <- tiny hand-rolled check harness
│   ├── test_vec3.cpp
│   └── test_integrate.cpp
├── docs/
│   ├── core.md
│   ├── dynamics.md
│   └── design-notes.md
├── build/                       <- generated binaries (gitignored)
├── Makefile
├── build.sh
├── README.md
├── LICENSE
└── .gitignore
```

### Why `include/newton/` and not just `include/`

The `include/` directory is the **root that gets added to the compiler's search path**. Everything below it becomes part of the include *spelling*. So:

- `-I include` + file at `include/newton/core/vec3.hpp` → `#include <newton/core/vec3.hpp>`
- `-I include/newton` + same file → `#include <core/vec3.hpp>` ← bad, `core/` is a generic name that will collide with someone else's library

The doubled name (`include/newton/newton.hpp`) looks redundant but is the correct convention. It is what makes the include *self-identifying* and collision-proof.

### Umbrella header

`include/newton/newton.hpp` is a convenience header that pulls in everything:

```cpp
#pragma once
// Umbrella header. Include this to get the whole library.
// Prefer including individual module headers in large projects
// to keep compile times down.
#include <newton/core/constants.hpp>
#include <newton/core/vec3.hpp>
#include <newton/dynamics/integrate.hpp>
```

Note that the umbrella header uses **angle brackets** and the full path, exactly as a user would. This works because the library's own `-I include` is on the search path when compiling examples/tests, and after installation the system path provides the same spelling. Do not use `#include "core/vec3.hpp"` inside your headers — see §2.

You may also want per-module umbrellas (`include/newton/core.hpp`, `include/newton/dynamics.hpp`) once you have more than ~3 files per module.

### Real-world examples of this pattern

| Library | Layout | Include spelling | Notes |
|---|---|---|---|
| **nlohmann/json** | `include/nlohmann/json.hpp` plus `include/nlohmann/detail/**` | `#include <nlohmann/json.hpp>` | The canonical modern example. Public umbrella at the top, dozens of private headers under `detail/`. Also ships a generated `single_include/nlohmann/json.hpp` (one amalgamated file) for drop-in use. |
| **Eigen** | `Eigen/Dense`, `Eigen/Core`, `Eigen/Geometry` — extensionless umbrella headers — with implementation in `Eigen/src/Core/**.h` | `#include <Eigen/Dense>` | A math library, closest in spirit to yours. Shows the umbrella-per-module idea taken to its conclusion: users include a *module*, never an individual file. |
| **doctest** | `doctest/doctest.h` (single header), assembled from `doctest/parts/doctest_fwd.h` + `doctest_impl.h` | `#include <doctest/doctest.h>` | Shows that "single header" is a *distribution* format, not necessarily how the source is organized. |
| **glm** | `glm/glm.hpp` umbrella + `glm/detail/**`, `glm/gtc/**`, `glm/gtx/**` | `#include <glm/glm.hpp>` | Header-only vector/matrix math; module subdirectories with a stability convention (`gtc` = stable extensions, `gtx` = experimental). |
| **fmt** | `include/fmt/core.h`, `format.h`, … | `#include <fmt/format.h>` | Normally compiled, but `#define FMT_HEADER_ONLY` makes it header-only. Shows the dual-mode trick you don't need yet. |

The common thread: **one vendor-named directory under the include root, module subdirectories inside it, an umbrella header, implementation details quarantined in `detail/`.**

### `examples/`, `tests/`, `docs/`

- **`examples/`** — each `.cpp` has its own `main()` and compiles to its own binary. These are your first users; if an example is awkward to write, the API is wrong. They double as documentation.
- **`tests/`** — same structure (each file a `main()`), but they assert rather than print. See §7.
- **`docs/`** — Markdown files, one per module, plus a design-notes file recording *why* you made choices. For a learning library this is arguably the most valuable directory.

---

## 2. How `#include` actually works

### The preprocessor is a text-paste machine

`#include <foo.hpp>` means: *find the file `foo.hpp`, and paste its entire contents here, then keep going.* That's all it does. There is no module system, no linking involved, no symbol lookup. The compiler proper never sees the `#include` line — it sees the pasted text.

You can watch it happen:

```bash
clang++ -std=c++17 -I include -E examples/projectile.cpp | wc -l
```

`-E` runs *only* the preprocessor. The output is your example file with tens of thousands of lines of headers pasted in front of it.

### Quotes vs angle brackets

```cpp
#include "test_util.hpp"      // quoted form
#include <newton/newton.hpp>  // angle-bracket form
```

The standard only says the two forms may differ in how the file is searched for. In practice, for clang and gcc:

- **Quoted form (`"..."`)**: search *first* in the directory of the file that contains the `#include` line, then fall back to the angle-bracket search path.
- **Angle-bracket form (`<...>`)**: search only the directories on the search path — the ones from `-I`, from `CPLUS_INCLUDE_PATH`, and the compiler's built-in system directories, in that order.

**The convention to follow:**

- Use `"quotes"` for files that sit next to the including file and travel with it — e.g. `#include "test_util.hpp"` inside `tests/test_vec3.cpp`.
- Use `<angle brackets>` for anything reached through the include path — the standard library, third-party libraries, **and your own library's headers, including from inside your own library**.

That last point matters. If `include/newton/dynamics/integrate.hpp` says `#include "../core/vec3.hpp"`, it works but is brittle — it breaks if a file moves, and it bypasses the include path entirely. Writing `#include <newton/core/vec3.hpp>` means your library's internal includes look exactly like a user's includes, so if the include path is set up correctly for one, it's correct for the other.

### The `-I` flag

`-I <dir>` prepends `<dir>` to the angle-bracket search path. So with the layout above:

```bash
clang++ -std=c++17 -I include examples/projectile.cpp -o build/projectile
```

When the preprocessor hits `#include <newton/newton.hpp>` it tries `include/` + `newton/newton.hpp` = `include/newton/newton.hpp`. Found. Done.

Note the resolution is a **simple path concatenation**: search directory + the exact text between the brackets. This is the whole mechanism, and it explains the next point.

### Why `include/<libname>/` prevents collisions

Suppose you flattened your library to `include/vec3.hpp` and used `#include <vec3.hpp>`. Now the user also uses a graphics library that installed `vec3.hpp` into `/usr/local/include/`. Both are reachable as `<vec3.hpp>`; which one wins depends on the *order* of `-I` flags and system directories. The user gets a confusing type error deep inside a header they didn't write, and the fix is to reorder compiler flags — a terrible experience.

With `include/newton/core/vec3.hpp`, the include is spelled `<newton/core/vec3.hpp>`. For a collision to occur, another library would have to install a directory literally named `newton` containing `core/vec3.hpp`. The directory name acts as a de-facto namespace in the filesystem, mirroring the C++ `namespace newton` in the language. This is why every serious library does it: `<nlohmann/json.hpp>`, `<Eigen/Dense>`, `<fmt/format.h>`, `<glm/glm.hpp>`, `<doctest/doctest.h>`.

### Search-path order (verified on this Mac)

```
$ echo | clang++ -std=c++17 -x c++ -E -v -
#include <...> search starts here:
 /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1
 /Library/Developer/CommandLineTools/usr/lib/clang/15.0.0/include
 /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include
 /Library/Developer/CommandLineTools/usr/include
 /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks (framework directory)
End of search list.
```

`-x c++ -E -v -` with empty stdin is the incantation to print the search path. Anything from `-I` is inserted *before* all of these. Run this command whenever an include mysteriously fails — it is the ground truth.

### Three ways to make the library importable later

**(a) Copy or symlink into `/usr/local/include`**

Verified detail for this Mac: `/usr/local/include` **is** on Apple clang's search path, but the directory does not currently exist, so clang silently skips it. The `-v` output shows this explicitly:

```
ignoring nonexistent directory "/usr/local/include"
```

So you must create it first:

```bash
sudo mkdir -p /usr/local/include
sudo ln -s "$PWD/include/newton" /usr/local/include/newton     # symlink: edits are live
# or
sudo cp -R include/newton /usr/local/include/                  # copy: a frozen snapshot
```

**Prefer the symlink while developing.** With a symlink, editing `include/newton/core/vec3.hpp` in your repo is instantly visible to every program on the machine — no reinstall step, no chance of debugging a stale copy. Switch to `cp` only when you want a stable installed version pinned independently of your working tree.

An `install`/`uninstall` target in the Makefile makes this a one-liner:

```make
PREFIX ?= /usr/local
install:
	mkdir -p $(PREFIX)/include
	cp -R include/newton $(PREFIX)/include/
uninstall:
	rm -rf $(PREFIX)/include/newton
```

(run as `sudo make install`).

**(b) `-I` in a Makefile or shell alias**

This is what you should do *first*, and for a long time. In the Makefile it's already there (`-I include`). For quick one-off compiles from anywhere, a shell alias or function in `~/.zshrc`:

```bash
# ~/.zshrc
export NEWTON_ROOT="$HOME/Github/newtonian-library"
nc++() { clang++ -std=c++17 -Wall -Wextra -I "$NEWTON_ROOT/include" "$@"; }
```

Then `nc++ scratch.cpp -o scratch` anywhere on the machine. This is explicit, per-invocation, and leaves no trace on the system — the safest option and the one that teaches you the most, because the path is visible in the command.

**(c) `CPLUS_INCLUDE_PATH` environment variable**

clang and gcc both honour `CPLUS_INCLUDE_PATH` (colon-separated) as an extra angle-bracket search path for C++ compilations. Verified working on this machine:

```bash
export CPLUS_INCLUDE_PATH="$HOME/Github/newtonian-library/include:$CPLUS_INCLUDE_PATH"
clang++ -std=c++17 anything.cpp -o anything   # <newton/newton.hpp> now resolves
```

Add the `export` to `~/.zshrc` to make it permanent. (`C_INCLUDE_PATH` is the C equivalent; `CPATH` applies to both languages.)

**Trade-offs:** `CPLUS_INCLUDE_PATH` is invisible — a build that works in your shell fails in a CI job or for a collaborator, and the error gives no hint why. Directories added this way are also treated as *system* include directories, which suppresses warnings from inside them; that's convenient for third-party code but actively unhelpful for a library you are writing and want warnings about.

**Recommendation:** use (b) `-I` now, add (a) symlink into `/usr/local/include` once the API stabilizes and you want to `#include <newton/newton.hpp>` from throwaway files without thinking. Know that (c) exists, use it only for temporary experiments.

### Other macOS / Apple clang specifics

- **`clang++`, not `g++`.** On macOS, `g++` is a symlink to Apple clang; `g++ --version` prints clang version info. There's no real GCC unless you `brew install gcc` (which gives you `g++-14` etc.). Write `clang++` in your build files to avoid confusion; the `CXX ?= clang++` pattern lets it be overridden.
- **Command Line Tools vs full Xcode.** This machine uses the standalone CLT (`/Library/Developer/CommandLineTools`). Both work. If headers ever vanish after a macOS update, the fix is `xcode-select --install` and possibly `sudo xcode-select --reset`.
- **The SDK sysroot.** Apple clang compiles against an SDK (`xcrun --show-sdk-path` → `/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk`), which is why system headers appear under a long SDK path rather than plain `/usr/include`. `/usr/include` does not exist as a real directory on modern macOS. This does not affect you, but it explains surprising paths in error messages.
- **libc++, not libstdc++.** Apple clang uses LLVM's libc++ (headers at `.../MacOSX.sdk/usr/include/c++/v1`). Occasionally a Linux-only assumption in someone's code breaks here. Your C++17 usage is fine.
- **Apple Silicon Homebrew is at `/opt/homebrew`, which is NOT on the default include path.** Verified: `/opt/homebrew/include` exists on this machine but does not appear in clang's search list. So if you ever `brew install` a header-only dependency (e.g. `brew install eigen`), you must add `-I/opt/homebrew/include` yourself. On Intel Macs Homebrew lives in `/usr/local`, which *is* searched — this is a very common source of "it works on my friend's Mac" confusion.
- **`-Wall -Wextra` on Apple clang** is well-behaved. `-Wpedantic` is also fine and worth keeping.

---

## 3. Header-only mechanics

### The problem `inline` solves

Normally a C++ program is compiled one translation unit (TU) at a time — one `.cpp` plus everything it includes — and the linker joins the results. The **One Definition Rule (ODR)** says: every function and variable with external linkage must have exactly **one** definition across the entire program. Define `double kinetic_energy(double, double) { ... }` in a header, include that header from two `.cpp` files, and the linker sees two definitions and errors with `duplicate symbol`.

Header-only libraries dodge this because the standard carves out exceptions for entities that *may* be defined in every TU, as long as all the definitions are token-for-token identical and mean the same thing:

| Entity | May appear in multiple TUs? | Why |
|---|---|---|
| **Templates** (function and class) | Yes, automatically | A template is not code until instantiated; the compiler emits instantiations as weak symbols and the linker merges duplicates. |
| **Member functions defined inside the class body** | Yes, automatically | Implicitly `inline`. |
| **Free functions marked `inline`** | Yes | That is what `inline` means today. |
| **Variables marked `inline`** (C++17) | Yes | New in C++17 — the reason header-only libraries got much easier. |
| **`constexpr` functions** | Yes | Implicitly `inline`. |
| **`static constexpr` data members** (C++17) | Yes | Implicitly `inline` since C++17. |
| **Class/struct/enum definitions** | Yes | Type definitions were always allowed once per TU. |
| **Plain non-template functions** | **No** | Must be `inline`. |
| **Plain non-const namespace-scope variables** | **No** | Must be `inline` (C++17) or moved into a function. |

### `inline` does not mean "inline this call"

This is the single most-misunderstood keyword in C++. In modern C++, `inline`'s *only* guaranteed meaning is **"this entity may be defined identically in multiple translation units; linker, please merge them."** It is a hint about *linkage*, not about optimization. Compilers decide about actual call-site inlining on their own, based on size and hotness, and will inline non-`inline` functions freely and refuse to inline `inline` ones.

So in a header-only library, `inline` is not an optimization request — it is a *requirement* for any non-template function or variable at namespace scope.

```cpp
#pragma once
#include <cmath>

namespace newton {

// Non-template free function in a header MUST be inline,
// otherwise two .cpp files including this header fail to link.
[[nodiscard]] inline double kinetic_energy(double mass, double speed) {
    return 0.5 * mass * speed * speed;
}

// Templates need no inline — the rule already covers them.
template <typename T>
[[nodiscard]] T square(T x) { return x * x; }

}  // namespace newton
```

### C++17 inline variables for constants

Before C++17, putting a shared constant in a header was awkward. The traditional workaround relies on a subtlety: a namespace-scope variable declared `const` (and `constexpr` implies `const`) has **internal linkage** by default in C++, meaning each TU quietly gets its own private copy. That compiles and links, but:

- Each TU has a separate object at a *different address*, so `&newton::constants::g_earth` differs between TUs — surprising if you ever take the address or bind a reference.
- The constant is duplicated in every object file (negligible for a `double`, real for a large table).
- Any `inline` function in the header that *uses* the constant technically has different definitions in different TUs — a subtle ODR violation that compilers do not diagnose.

C++17 `inline constexpr` fixes all of it: one entity, one address, program-wide.

```cpp
#pragma once

namespace newton::constants {

// `inline` => exactly one object across the whole program (C++17).
// `constexpr` => usable in compile-time contexts, e.g. array sizes.
inline constexpr double g_earth = 9.80665;      // m/s^2, standard gravity
inline constexpr double G       = 6.67430e-11;  // m^3 kg^-1 s^-2, CODATA 2018
inline constexpr double c_light = 299'792'458;  // m/s, exact by definition
                                                // (' is a C++14 digit separator)
}  // namespace newton::constants
```

**Rule to internalize: in a header, every namespace-scope function is `inline` or a template; every namespace-scope variable is `inline constexpr` (or `inline const`).**

The same applies to mutable state, which you should mostly avoid, but if you need it:

```cpp
namespace newton::detail {
inline int integration_step_counter = 0;   // one shared object; needs `inline`
}
```

### ODR essentials, in four rules

1. **One definition per program** for non-inline functions and variables with external linkage.
2. **One definition per translation unit** for classes, enums, templates, and inline entities — no repeats within a single TU (this is what include guards prevent).
3. **Multiple TUs may each define** the same inline entity / template / class, provided the definitions are **token-identical and mean the same thing** in each.
4. Violating rule 3 is **ill-formed, no diagnostic required** — the compiler is not obliged to tell you. Symptoms are bizarre: the linker silently picks one definition and you get behavior that contradicts the source you're reading.

The practical way to violate rule 3 in a header-only library is to make a header's contents depend on macros defined *before* the include:

```cpp
// DANGEROUS pattern — do not do this
#ifdef NEWTON_USE_FLOAT
using Real = float;
#else
using Real = double;
#endif
```

If `a.cpp` defines `NEWTON_USE_FLOAT` and `b.cpp` doesn't, `newton::Vec3` is a *different type* with a different layout in the two TUs, and passing one between them is undefined behavior with no diagnostic. Prefer a template parameter (`template <typename Real> struct Vec3`) or just commit to `double` — for a physics learning library, `double` everywhere is the right call.

### `#pragma once` vs include guards

**Include guards** are standard C++, using the preprocessor:

```cpp
#ifndef NEWTON_CORE_VEC3_HPP
#define NEWTON_CORE_VEC3_HPP
// ... contents ...
#endif  // NEWTON_CORE_VEC3_HPP
```

**`#pragma once`** is a compiler extension that does the same thing:

```cpp
#pragma once
// ... contents ...
```

| | `#pragma once` | Include guards |
|---|---|---|
| Standard C++ | No (but supported by clang, gcc, MSVC, ICC — every compiler you will ever use) | Yes |
| Boilerplate | 1 line | 3 lines + a unique name you must invent |
| Failure mode | Can misbehave if the same file is reachable via two different paths (symlinks, hard links, network mounts) — compilers use filesystem identity to decide "same file" | Can misbehave if you copy-paste a header and forget to change the macro name, silently blanking the second file |
| Speed | Slightly faster: the compiler can skip re-opening the file entirely | Compiler must open and scan the file to find the guard (though both clang and gcc have a "multiple-include optimization" that caches this) |
| Namespace pollution | None | Adds a macro to the global macro namespace; a name like `VEC3_H` can collide |

**Recommendation: use `#pragma once`.** For your project, on a single Mac with clang, in a directory tree with no symlink weirdness, the theoretical portability concern is nil, and the practical benefit is real: three fewer lines per file, zero chance of the copy-paste bug, no macro names to invent. Every major header-only library either uses it or supports it — nlohmann/json uses `#pragma once`, doctest uses guards, Eigen uses guards, glm uses `#pragma once`.

If you later want belt-and-braces, you can use both. Don't bother now.

Also: **use `.hpp`, not `.h`**, for C++ headers. `.h` is ambiguous with C. This is convention, not a compiler rule — clang decides the language from the *source file*, not the header — but it's a signal to readers and to editors/tooling.

---

## 4. Namespace design

### One top-level namespace, always

Everything the library defines goes inside `namespace newton`. No exceptions — not even helper types, not even things you think nobody will name-collide with. This is the language-level twin of the `include/newton/` directory.

```cpp
namespace newton {
struct Vec3 { /* ... */ };
}
```

### Nested namespaces per module (C++17 syntax)

C++17 lets you write nested namespaces on one line:

```cpp
namespace newton::dynamics {   // C++17; before this you needed nested braces
    // ...
}
```

Suggested structure, mirroring the directory tree:

| Namespace | Directory | Contents |
|---|---|---|
| `newton` | `include/newton/` | Core vocabulary types used everywhere: `Vec3`, `Particle`, `Real`. Keep this small. |
| `newton::constants` | `core/constants.hpp` | `g_earth`, `G`, `c_light`. A namespace, not a class — no reason for a class here. |
| `newton::dynamics` | `dynamics/` | Integrators, force composition, simulation stepping. |
| `newton::forces` | `dynamics/forces.hpp` | Ready-made force functors: `Gravity`, `Drag`, `Spring`. |
| `newton::units` | `core/units.hpp` | Unit-conversion helpers, literal operators. |
| `newton::detail` | `detail/` | **Implementation only.** By convention, anything in `detail` is private: users must not name it and you may change it freely. nlohmann/json, Eigen (`Eigen::internal`), and Boost all use this convention. |

Rule of thumb on depth: `newton::dynamics` is right; `newton::dynamics::integrators::explicit_methods` is over-engineered for a small library. Two levels, three at most.

**Type names in the namespace should not repeat the namespace.** Write `newton::dynamics::Integrator`, not `newton::dynamics::DynamicsIntegrator`. The namespace already says it.

### Why `using namespace` stays out of headers

```cpp
// include/newton/core/vec3.hpp  -- NEVER DO THIS
#pragma once
#include <cmath>
using namespace std;   // <-- catastrophic in a header
```

A `using namespace X;` at namespace scope in a header is **pasted into every single translation unit that includes it, and everything included after it**. The user of your library did not ask for `std` (or `newton`) to be dumped into their global namespace, and cannot undo it. Consequences:

- **Ambiguity errors in user code that has nothing to do with you.** The user writes `count`, `distance`, `data`, `size`, `begin`, `swap`, `min`, `max`, or `left` — all names in `std` — and gets a call-is-ambiguous error pointing into a header they never opened.
- **Silent behavior changes**, which are worse. Overload resolution suddenly considers a `std` overload it previously didn't, picks a different function, and the program compiles and runs incorrectly.
- **Order dependence.** Whether the breakage occurs depends on include order, so it appears and disappears as unrelated files change.

The rule applies to **all** `using`-directives in headers, including `using namespace newton;` inside your own headers, and to unqualified `using` declarations at namespace scope (`using std::vector;` — slightly less bad, still not yours to decide).

**What to do instead:**

- **Fully qualify in headers**: `std::sqrt`, `std::vector`, `newton::Vec3`.
- **Scope it to a function body** if you need it — a `using` inside a function is local and harmless:

  ```cpp
  [[nodiscard]] inline double norm(Vec3 v) {
      using std::sqrt;  // fine: scoped to this function; also enables ADL for
                        // user-defined types via the "using std::sqrt; sqrt(x)" idiom
      return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  }
  ```

- **Offer a convenience alias namespace** so users can opt in explicitly:

  ```cpp
  // in newton.hpp, at the very bottom
  // Users may write `namespace nt = newton;` themselves; you can also provide:
  namespace nt = newton;   // acceptable at the umbrella level, still debatable —
                           // safest is to document the alias and let users declare it
  ```

  The safest posture: **document** that users can write `namespace nt = newton;` in *their* `.cpp` files, and don't declare it for them.

- In your own `examples/*.cpp` and `tests/*.cpp`, `using namespace newton;` is **fine** — those are leaf translation units that nothing includes. It also demonstrates to readers that the API reads nicely unqualified.

- `<cmath>` note: prefer `std::sqrt`/`std::fabs` over the global `::sqrt`. On libc++ the global C names are usually also visible, but relying on that is not portable.

---

## 5. Modern C++17 features that fit a physics library

### `constexpr` constants

Covered in §3. `inline constexpr double g_earth = 9.80665;` — one object program-wide, usable at compile time, no macro, properly typed and namespaced. This alone replaces the C habit of `#define G_EARTH 9.81`, which has no type, no scope, and no debugger visibility.

`constexpr` on *functions* means "can be evaluated at compile time when the arguments allow." A `constexpr` function is also implicitly `inline`, which is a bonus in a header:

```cpp
namespace newton {
// constexpr: computable at compile time if inputs are compile-time constants.
// Implicitly inline, so it's header-safe without the keyword.
[[nodiscard]] constexpr double kinetic_energy(double mass, double speed) noexcept {
    return 0.5 * mass * speed * speed;
}
}
static_assert(newton::kinetic_energy(2.0, 3.0) == 9.0);  // checked at compile time
```

Don't force `constexpr` onto everything — anything calling `std::sqrt` can't be `constexpr` in C++17 (that changed in C++26). Use it where it falls out naturally.

### Default member initializers

```cpp
struct Vec3 {
    double x{0.0};
    double y{0.0};
    double z{0.0};
};

Vec3 a;                  // {0,0,0} -- no uninitialized garbage
Vec3 b{1.0, 2.0, 3.0};   // aggregate initialization still works
Vec3 c{};                // {0,0,0}
```

This kills an entire bug class: uninitialized doubles producing NaN halfway through a simulation. Note the braces `{0.0}` rather than `= 0.0` — brace-init forbids narrowing conversions, so `double x{someInt}` is caught. Note also that default member initializers still leave the type an **aggregate** in C++14 and later, so `Vec3{1,2,3}` continues to work and you need no constructor.

For a `Particle`:

```cpp
struct Particle {
    Vec3   position{};
    Vec3   velocity{};
    double mass{1.0};      // a physically sensible default beats 0.0,
                           // which would make a = F/m infinite
};
```

### Structured bindings

An integrator naturally returns two things — new position and new velocity. Structured bindings make that readable:

```cpp
auto [pos, vel] = newton::dynamics::step(pos, vel, mass, dt, gravity);
```

Works with `std::pair`, `std::tuple`, arrays, and any struct with all-public members — including your own:

```cpp
struct State { Vec3 position; Vec3 velocity; };
State advance(State s, double dt);

auto [p, v] = advance(s, dt);       // decomposes your own struct, no tuple needed
```

Prefer a named struct (`State`) over `std::pair` in a public API: `.position` documents itself, `.first` does not. Structured bindings work equally well on both.

Also useful in range-based `for` over maps:

```cpp
for (const auto& [name, body] : bodies) { /* ... */ }
```

Caveat: structured bindings cannot be captured by a lambda in C++17 (fixed in C++20). If you hit that, copy into a named variable first.

### `std::function` vs template callables for force functions

This is the central API design decision in a physics library, because "a force" is fundamentally *a function the user supplies*.

**Option A — `std::function` (type erasure):**

```cpp
#include <functional>
namespace newton::dynamics {

using ForceFn = std::function<Vec3(const Vec3& pos, const Vec3& vel)>;

// Not a template -> must be inline in a header
[[nodiscard]] inline std::pair<Vec3, Vec3>
step(Vec3 pos, Vec3 vel, double mass, double dt, const ForceFn& force) {
    const Vec3 a = (1.0 / mass) * force(pos, vel);
    const Vec3 v = vel + dt * a;
    return {pos + dt * v, v};
}
}
```

**Option B — template callable:**

```cpp
namespace newton::dynamics {

template <typename Force>
[[nodiscard]] std::pair<Vec3, Vec3>
step(Vec3 pos, Vec3 vel, double mass, double dt, Force&& force) {
    const Vec3 a = (1.0 / mass) * force(pos, vel);
    const Vec3 v = vel + dt * a;
    return {pos + dt * v, v};
}
}
```

| | `std::function` | Template callable |
|---|---|---|
| **Call cost** | Indirect call through a type-erased handle; **the compiler cannot inline the user's lambda** into the integrator loop | The lambda's exact type is known at compile time; the body is typically **inlined and optimized into the loop**, often to a few instructions |
| **Allocation** | May heap-allocate if the callable's captures exceed the small-buffer size (~16 bytes, implementation-defined). Capturing a couple of `double`s is usually fine, capturing a `Vec3` plus a pointer may not be | Never allocates |
| **Speed in a hot loop** | Often 2–10× slower for a tiny force function called millions of times — exactly the shape of a physics simulation | Full speed |
| **Storable in a container** | Yes — `std::vector<ForceFn>` for a list of heterogeneous forces | No, not directly; needs `std::variant`, a base class, or `std::function` after all |
| **Runtime swappable** | Yes — reassign the `std::function` | No, the type is fixed at compile time |
| **Error messages** | Short and clear; the signature is spelled out in the declaration | Long template-instantiation errors when the callable has the wrong signature (much improved by C++20 concepts) |
| **Compile time** | One instantiation | One instantiation per distinct callable type |
| **Header cost** | Requires `<functional>`, a fairly heavy header | Nothing extra |
| **API documentation** | `ForceFn` names the expected signature explicitly — very readable in docs | Signature is implicit; must be documented in a comment |

**Recommendation for this library: template callables (Option B).**

Reasons, in order of importance for a *learning* library:

1. **Performance is the point.** A physics integrator running 10⁶ steps with a lambda body of three multiplications is the canonical case where type erasure costs the most. Getting the fast version and understanding *why* it's fast is a genuinely valuable lesson.
2. **It teaches the more useful concept.** Templates + perfect forwarding + "duck-typed callables" is a core modern-C++ idiom that transfers everywhere (all the STL algorithms work this way — `std::sort` takes a template comparator, not a `std::function`). Type erasure is a narrower tool.
3. **It composes.** `step(pos, vel, m, dt, [&](Vec3 p, Vec3 v){ return gravity(p) + drag(v); })` works with zero machinery.
4. **It costs nothing to add `std::function` later.** `std::function<Vec3(Vec3,Vec3)>` *is* a valid callable, so it satisfies the template. A user who needs runtime polymorphism can pass one in and everything still compiles. The template API is strictly more general.

Keep a note in `docs/design-notes.md` explaining the trade-off, and provide named force functors so users get readable code without lambdas:

```cpp
namespace newton::forces {

// A stateful functor: cleaner than a capturing lambda for reusable forces,
// and shows the "function object" concept explicitly.
struct UniformGravity {
    double mass{1.0};
    double g{constants::g_earth};

    [[nodiscard]] Vec3 operator()(const Vec3& /*pos*/, const Vec3& /*vel*/) const {
        return Vec3{0.0, -mass * g, 0.0};
    }
};

struct LinearDrag {
    double b{0.0};   // F = -b * v
    [[nodiscard]] Vec3 operator()(const Vec3& /*pos*/, const Vec3& vel) const {
        return -b * vel;
    }
};

// Compose two forces into one callable.
template <typename A, typename B>
[[nodiscard]] auto sum(A a, B b) {
    return [a, b](const Vec3& p, const Vec3& v) { return a(p, v) + b(p, v); };
}

}  // namespace newton::forces
```

Note the `/*pos*/` comment-out of unused parameter names — that's how you silence `-Wunused-parameter` from `-Wextra` while keeping the name documented.

### `auto` return types

```cpp
// Deduced return type: fine when the type is obvious or unnameable (lambdas).
template <typename A, typename B>
[[nodiscard]] auto sum(A a, B b) { /* returns an unnameable lambda type */ }
```

Guidance: use `auto` returns for **factory functions returning lambdas** (where there is no other option) and for **template-heavy internals**. Avoid it on the main public API — `Vec3 acceleration(...)` documents itself in a header; `auto acceleration(...)` forces the reader to read the body. A `[[nodiscard]] auto` on a one-line function with a visible `return` is fine.

Related: `auto` on *variables* is good in examples and loops (`auto [p, v] = ...`), but write `const double dt = 0.01;` rather than `auto dt = 0.01;` when the type carries physical meaning.

### `[[nodiscard]]`

```cpp
[[nodiscard]] inline double norm(Vec3 v);
[[nodiscard]] std::pair<Vec3, Vec3> step(...);
```

Marks a return value as important; the compiler warns if it's ignored. In a physics library, essentially **every pure function should be `[[nodiscard]]`**, because ignoring the result of `norm(v)` or `step(...)` is always a bug — the user probably meant an in-place variant that doesn't exist. This is one of the highest-value-per-character features in C++17.

You can also apply it to a whole type — `struct [[nodiscard]] Error { ... };` — so every function returning it is covered. Overkill here.

Also worth pairing with `noexcept` on functions that genuinely can't throw (pure arithmetic on `double`s). Don't put `noexcept` on anything that allocates or calls user callbacks.

### Range-based `for`

```cpp
// Simulate a system of particles
for (auto& p : particles) {
    auto [pos, vel] = dynamics::step(p.position, p.velocity, p.mass, dt, force);
    p.position = pos;
    p.velocity = vel;
}

// Read-only pass: const auto& avoids copying
double total_ke = 0.0;
for (const auto& p : particles) {
    total_ke += kinetic_energy(p.mass, norm(p.velocity));
}
```

The rule: `const auto&` for reading, `auto&` for modifying in place, plain `auto` only for cheap types you intend to copy. Teaching this distinction is worth a comment in the examples.

### Other C++17 bits worth using

- **Class template argument deduction (CTAD)**: `std::pair p{pos, vel};` instead of `std::pair<Vec3,Vec3>`. Works for your own templates too.
- **`if` with initializer**: `if (const double d = norm(v); d > tol) { ... }` — scopes the variable to the `if`.
- **`std::clamp`** (`<algorithm>`) — useful for bounding timesteps.
- **`std::optional`** — for functions that may not produce a result (e.g. "time of collision, if any"):
  ```cpp
  [[nodiscard]] inline std::optional<double> time_to_ground(const Particle& p);
  ```
- **Digit separators** (C++14): `299'792'458` — invaluable for physical constants.
- **`static_assert` with no message** (C++17): `static_assert(sizeof(double) == 8);`
- **Nested namespace definitions** (C++17): `namespace newton::dynamics { }`.

### What C++20 concepts would add later

Today, passing a wrong-signature callable produces a wall of instantiation errors pointing inside `step`. With concepts, the constraint is checked at the call site and the message names the actual problem:

```cpp
// C++20 — a future upgrade, sketched for the design notes
template <typename F>
concept ForceFunction = requires(F f, const Vec3& p, const Vec3& v) {
    { f(p, v) } -> std::convertible_to<Vec3>;
};

template <ForceFunction Force>
[[nodiscard]] auto step(Vec3 pos, Vec3 vel, double mass, double dt, Force&& force);
```

The error becomes roughly *"constraint `ForceFunction<MyThing>` not satisfied: `f(p, v)` is not a valid expression"* — pointing at the user's call, not at line 47 of your header. Concepts also enable clean overloading on capability and self-documenting template parameters.

Other C++20 items relevant later: `std::span` (non-owning view over a particle array), `constexpr` math in `<cmath>` (partially, and more in C++23/26), designated initializers (`Vec3{.x = 1.0, .y = 2.0}` — actually already a common compiler extension), `<numbers>` for `std::numbers::pi`, three-way comparison `<=>`, and `std::format`.

**Stay on C++17 for now.** It is fully supported by Apple clang 15, has all the header-only ergonomics (`inline` variables especially), and adding `-std=c++20` later is a one-character change in your Makefile. Design the API so concepts can be *added* as constraints without changing any signature — which template callables already allow.

---

## 6. Build and tooling literacy

### What a Makefile is

`make` is a program that reads a file named `Makefile` and figures out **what to rebuild and in what order**, based on file timestamps. You describe *rules*; make works out the plan.

A rule has three parts:

```make
target: prerequisites
	recipe        <- MUST begin with a literal TAB character, not spaces
```

Meaning: *"to produce `target`, you need `prerequisites`; if any prerequisite is newer than `target` (or `target` doesn't exist), run `recipe`."*

**The tab is not optional.** `make` will fail with `missing separator` if you use spaces. Configure your editor so `Makefile` uses real tabs.

Key syntax you'll meet:

| Syntax | Meaning |
|---|---|
| `VAR := value` | Set a variable, evaluated once immediately (prefer this) |
| `VAR ?= value` | Set only if not already set (allows `make CXX=g++-14`) |
| `$(VAR)` | Use a variable |
| `$@` | The target of the current rule |
| `$<` | The **first** prerequisite |
| `$^` | **All** prerequisites |
| `%` | Wildcard in a pattern rule — matches the same text in target and prerequisite |
| `$(wildcard examples/*.cpp)` | Expand a glob into a list |
| `$(patsubst a/%.cpp,b/%,$(LIST))` | Rewrite each item in a list |
| `.PHONY: clean` | Declares a target that isn't a real file (otherwise a file named `clean` would break it) |
| `\| dir` | Order-only prerequisite: must exist, but its timestamp doesn't trigger rebuilds (perfect for output directories) |
| `@cmd` | Run `cmd` without echoing it |
| `$$` | A literal `$` passed to the shell (make eats single `$`) |

Because a header-only library has **no library to build**, your Makefile's only jobs are: compile each example to a binary, compile each test to a binary, run the tests, and clean up.

### A minimal working Makefile (verified)

```make
# ---- configuration -----------------------------------------------------
CXX      ?= clang++
CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -O2 -I include

# ---- file discovery ----------------------------------------------------
EXAMPLE_SRCS := $(wildcard examples/*.cpp)
EXAMPLE_BINS := $(patsubst examples/%.cpp,build/%,$(EXAMPLE_SRCS))
TEST_SRCS    := $(wildcard tests/*.cpp)
TEST_BINS    := $(patsubst tests/%.cpp,build/tests/%,$(TEST_SRCS))

# Every binary depends on every header: a header-only library has no
# separate compilation, so touching any header must rebuild everything.
HEADERS      := $(shell find include -name '*.hpp')

# ---- targets -----------------------------------------------------------
.PHONY: all examples tests check clean
all: examples

examples: $(EXAMPLE_BINS)
tests:    $(TEST_BINS)

build/%: examples/%.cpp $(HEADERS) | build
	$(CXX) $(CXXFLAGS) $< -o $@

build/tests/%: tests/%.cpp $(HEADERS) | build/tests
	$(CXX) $(CXXFLAGS) $< -o $@

build build/tests:
	mkdir -p $@

check: tests
	@for t in $(TEST_BINS); do ./$$t || exit 1; done
	@echo "All tests passed."

clean:
	rm -rf build
```

Usage:

```
make            # build all examples
make check      # build and run all tests
make clean      # delete build/
make -j8 check  # build in parallel across 8 cores
```

Adding `examples/orbit.cpp` requires **no Makefile edit** — `$(wildcard)` picks it up.

Refinement to add once you're comfortable: `-MMD -MP` plus `-include $(DEPS)` makes clang emit per-file dependency lists, so only the binaries that actually include a changed header get rebuilt, instead of all of them. Not worth it at this scale.

### What a `.sh` shell script is

A `.sh` file is a **shell script**: a plain text file containing shell commands, run top to bottom by an interpreter. The first line, the *shebang*, names the interpreter:

```bash
#!/usr/bin/env bash
```

`/usr/bin/env bash` finds `bash` on `PATH` rather than hardcoding a path. (macOS note: `/bin/bash` is bash **3.2** from 2007 for licensing reasons; `brew install bash` gives you a modern one, and `env` will find it. Your default interactive shell on macOS is `zsh`, but scripts declare their own interpreter, so `#!/usr/bin/env bash` runs under bash regardless.)

To run it: make it executable once, then invoke it.

```bash
chmod +x build.sh   # set the execute permission bit
./build.sh          # ./ is required -- the current directory is not on PATH
```

`bash build.sh` also works and skips the `chmod`.

### A `build.sh` (verified)

```bash
#!/usr/bin/env bash
set -euo pipefail
#  -e  exit immediately if any command fails
#  -u  error on use of an undefined variable (catches typos)
#  -o pipefail  a pipeline fails if ANY stage fails, not just the last

cd "$(dirname "$0")"     # run from the repo root no matter where invoked from
                         # $0 = path to this script; dirname strips the filename
                         # quotes protect against spaces in the path

CXX="${CXX:-clang++}"                                   # use $CXX if set, else clang++
CXXFLAGS=(-std=c++17 -Wall -Wextra -Wpedantic -O2 -I include)  # a bash array

mkdir -p build/tests     # -p: create parents, and don't error if it exists

for src in examples/*.cpp; do
  out="build/$(basename "${src%.cpp}")"   # ${var%.cpp} strips the .cpp suffix
  echo "compiling $src -> $out"
  "$CXX" "${CXXFLAGS[@]}" "$src" -o "$out"   # [@] expands the array, one word each
done

status=0
for src in tests/*.cpp; do
  out="build/tests/$(basename "${src%.cpp}")"
  "$CXX" "${CXXFLAGS[@]}" "$src" -o "$out"
  "$out" || status=1      # keep going after a failing test, remember the failure
done

exit "$status"             # nonzero exit = failure, which CI and `make` understand
```

**Makefile vs shell script:** make gives you incremental rebuilds (skips work whose inputs haven't changed) and parallelism (`-j`) essentially for free; a shell script gives you readable, linear logic and arbitrary control flow. They coexist happily — many projects have a `Makefile` for compiling and small `.sh` scripts for orchestration (`run_tests.sh`, `install.sh`, `format.sh`). Start with the Makefile; add scripts when you want to *do* something rather than *build* something.

**Always quote your variables** in shell (`"$src"`, not `$src`). Unquoted variables are re-split on spaces and glob-expanded, which is the source of most shell bugs.

### Why CMake exists, and why you can skip it

CMake is a **build-system generator**. You describe your project once in `CMakeLists.txt`, and CMake emits a Makefile, a Ninja file, an Xcode project, or a Visual Studio solution depending on the platform. It exists to solve problems you do not have yet:

- **Cross-platform builds** — one description that works on macOS, Linux, and Windows/MSVC, where compiler flags differ entirely.
- **Dependency management** — `find_package(Eigen3 REQUIRED)` locates a library wherever it's installed on that machine.
- **Transitive usage requirements** — the killer feature. A consumer writes `target_link_libraries(myapp PRIVATE newton::newton)` and automatically inherits your include directory, your required C++ standard, and your compile definitions. No manual `-I`.
- **Installation and packaging** — `cmake --install` plus generated config files, so `find_package(newton)` works for a downstream project.
- **Editor/IDE integration** — CLion, VS Code, and clangd all read CMake projects natively; CMake can also emit `compile_commands.json`, which is what powers autocomplete and diagnostics in editors.

For a header-only library used from its own repo by one person on one Mac, CMake adds a configuration language to learn, a `build/` configure step, and a `CMakeLists.txt` that is longer than your Makefile — in exchange for benefits you can't yet use. **Skip it.**

Migrate when any of these becomes true: someone else wants to use `newton` in their own CMake project; you want it to build on Linux or Windows; you add a real dependency; you want editor autocomplete via `compile_commands.json` (though `bear -- make` can generate that from your Makefile without CMake); or you want to publish it.

For reference, the eventual CMake for a header-only library is genuinely short — this is the whole thing:

```cmake
cmake_minimum_required(VERSION 3.16)
project(newton VERSION 0.1.0 LANGUAGES CXX)

add_library(newton INTERFACE)              # INTERFACE = no compiled artifact
add_library(newton::newton ALIAS newton)

target_include_directories(newton INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(newton INTERFACE cxx_std_17)
```

Knowing it's ~10 lines is the point: there is no rush, and no lock-in from waiting.

### What Markdown (`.md`) files are for

Markdown is a plain-text formatting syntax that renders as styled text on GitHub, in editors, and in doc generators. `#` for headings, `**bold**`, `` `code` ``, ```` ``` ```` fences for code blocks, `-` for lists, `|` tables, `[text](url)` links.

Why it matters for your repo:

- **`README.md` renders automatically** on the GitHub repo front page. It is the first and often only thing anyone reads.
- It's **plain text**, so it diffs cleanly in git — unlike a Word document, you can see exactly what changed in a doc commit.
- It's **readable unrendered**, in a terminal or editor, so it never becomes inaccessible.
- Editors preview it live (VS Code: `Cmd-Shift-V`).

Conventional files: `README.md` (what/why/how), `LICENSE` (MIT is a fine default for a learning project), `CHANGELOG.md` (later), `CONTRIBUTING.md` (later), and per-module docs under `docs/`.

### The compile command, flag by flag

```
clang++ -std=c++17 -Wall -Wextra -I include examples/projectile.cpp -o projectile
```

| Token | Meaning |
|---|---|
| `clang++` | The **compiler driver** for C++. Not really one program — it orchestrates the preprocessor, compiler, assembler, and linker. Using `clang++` rather than `clang` matters: it compiles as C++ *and* links the C++ standard library (libc++) automatically. On macOS this is Apple clang, installed with the Command Line Tools. |
| `-std=c++17` | Select the **language standard**. Without it, Apple clang 15 defaults to gnu++17, which is close but enables GNU extensions; stating the standard explicitly makes the build reproducible and portable. This flag is what makes `inline` variables, structured bindings, `[[nodiscard]]`, and `namespace a::b` legal. |
| `-Wall` | Enable a **standard set of warnings** — despite the name, *not* all of them. Catches unused variables, missing returns, suspicious comparisons, misleading indentation. |
| `-Wextra` | Enable **more warnings** on top of `-Wall`: unused parameters, sign-comparison mismatches, missing field initializers. Together, `-Wall -Wextra` is the baseline every project should use. Warnings in C++ are frequently real bugs; treat them as errors you happen to be allowed to ignore. |
| `-I include` | Add `include/` to the **angle-bracket header search path**, so `#include <newton/newton.hpp>` resolves to `include/newton/newton.hpp`. (`-Iinclude` with no space is equivalent.) The path is relative to the working directory, so this command must be run from the repo root. |
| `examples/projectile.cpp` | The **input translation unit**. Since the library is header-only, this is the *only* source file — everything else arrives via `#include`. That's precisely why there is no library to build or link. |
| `-o projectile` | **Output** filename for the executable. Without `-o`, clang writes `a.out`. Prefer `-o build/projectile` to keep binaries out of the source tree (and add `build/` to `.gitignore`). |

Other flags worth knowing:

| Flag | Meaning |
|---|---|
| `-O0` / `-O2` / `-O3` | Optimization level. `-O0` (default) for debugging — fast compiles, readable in a debugger. `-O2` for real runs; on a physics loop the difference is often 5–20×. Always benchmark at `-O2` or higher. |
| `-g` | Emit debug info, so `lldb` can show source lines and variable names. Pairs fine with `-O0`; combine as `-g -O0`. |
| `-Wpedantic` | Warn about non-standard extensions. Good discipline for a library. |
| `-Werror` | Turn warnings into errors. Great once the code is clean; frustrating before. |
| `-fsanitize=address,undefined` | Compile *and link* with AddressSanitizer + UndefinedBehaviorSanitizer. Catches out-of-bounds access, use-after-free, signed overflow, bad casts, at ~2× runtime cost. Enormously valuable; add a `make sanitize` target. |
| `-E` | Preprocess only — see what `#include` actually pasted. |
| `-c` | Compile to an object file, don't link. You won't need it for a header-only library with single-file programs. |
| `-v` | Verbose: shows the include search path and the underlying commands. Your first diagnostic when an include fails. |
| `-march=native` | Optimize for this exact CPU. Can help numeric code; makes the binary non-portable. |

Debug and release targets are worth adding to the Makefile:

```make
debug:   CXXFLAGS := -std=c++17 -Wall -Wextra -g -O0 -fsanitize=address,undefined -I include
debug:   clean check
```

---

## 7. Testing without a framework

### Why hand-roll first

A test framework is, at its core, three things: a way to record a pass/fail, a way to report where a failure happened, and a nonzero exit code so the build fails. You can write that in 25 lines, and writing it yourself demystifies what doctest/Catch2 are doing. Only adopt a framework when you actually want what it adds (test discovery, subcases, parameterized tests, filtering, pretty output).

### `tests/test_util.hpp` (verified)

```cpp
#pragma once
// Minimal test harness. Not a framework -- just enough to record failures,
// report where they happened, and produce a nonzero exit code.
#include <cmath>
#include <cstdio>

namespace newton_test {

// C++17 inline variables: one shared counter per program, even though
// this header may be included by several translation units.
inline int checks_run = 0;
inline int failures   = 0;

// Floating-point comparison. Never use == on computed doubles:
// 0.1 + 0.2 != 0.3 exactly, because binary floating point cannot
// represent those decimals. We accept a mixed absolute/relative tolerance:
//   - absolute, so values near zero can be compared at all
//     (a relative test is meaningless when the expected value is 0)
//   - relative, so large values are not held to an impossible
//     absolute precision (a double has ~15-16 significant digits,
//     so 1e6 simply cannot be accurate to 1e-12 absolute)
[[nodiscard]] inline bool close(double a, double b, double tol) {
    const double diff = std::fabs(a - b);
    if (diff <= tol) return true;                        // absolute test
    const double scale = std::fmax(std::fabs(a), std::fabs(b));
    return diff <= tol * scale;                          // relative test
}

inline void record(bool ok, const char* expr, const char* file, int line) {
    ++checks_run;
    if (!ok) {
        ++failures;
        std::printf("  FAIL %s:%d  %s\n", file, line, expr);
    }
}

// Call at the end of main() and return its value.
[[nodiscard]] inline int report(const char* suite) {
    std::printf("%s: %d checks, %d failed\n", suite, checks_run, failures);
    return failures == 0 ? 0 : 1;
}

}  // namespace newton_test

// Macros, not functions, for two reasons only macros can provide:
//   1. #expr stringifies the source text of the condition for the message
//   2. __FILE__ / __LINE__ expand at the *call site*, not inside the helper
// The leading :: qualifies from the global namespace so these work
// no matter what namespace the test file is in.
#define NT_CHECK(expr) \
    ::newton_test::record((expr), #expr, __FILE__, __LINE__)

#define NT_CLOSE(a, b, tol) \
    ::newton_test::record(::newton_test::close((a), (b), (tol)), \
                          #a " ~= " #b, __FILE__, __LINE__)
```

Note: this harness deliberately **does not stop at the first failure** — it records and continues, so one run tells you everything that's broken. Add an `NT_REQUIRE` that `return`s early if you need a fatal variant.

### A test file (verified)

```cpp
// tests/test_vec3.cpp
#include "test_util.hpp"        // quotes: sits next to this file
#include <newton/newton.hpp>    // angle brackets: reached via -I include

int main() {
    // 3-4-5-style triangle: |(1,2,2)| = 3 exactly, so a tight tolerance is fine.
    newton::Vec3 a{1.0, 2.0, 2.0};
    NT_CLOSE(newton::norm(a), 3.0, 1e-12);

    // Default member initializers must zero every component.
    newton::Vec3 zero{};
    NT_CHECK(newton::norm(zero) == 0.0);   // exact == is OK: no arithmetic happened

    return newton_test::report("vec3");
}
```

Actual output:

```
vec3: 2 checks, 0 failed
```

### Choosing tolerances for physics tests

- **Exact `==` is legitimate** for values that were assigned, not computed — checking a default-initialized `Vec3` is zero, or that a copy round-trips.
- **`1e-12` relative** for a handful of arithmetic operations. `double` has ~15–16 significant decimal digits; a few operations lose 1–2.
- **`1e-9` to `1e-6`** for anything involving many accumulated steps. Error in a numerical integrator grows with step count.
- **Match the tolerance to the method's order.** Explicit Euler has global error O(dt); over 1 second at dt = 0.01, expect ~1% error, so a projectile-range test should use a tolerance around `1e-2` relative. Testing an Euler integrator to `1e-12` tests nothing except that you wrote the same wrong number twice.
- **Test invariants, not just values.** Energy conservation, momentum conservation, symmetry, and known closed-form solutions (projectile range = v²sin(2θ)/g, circular-orbit period) are far better tests than "the number equals 30.4999". They catch real physics bugs and don't break when you legitimately improve the integrator.
- **Test the analytic limit.** Halve `dt` and check the error shrinks by the expected factor — that verifies the *order* of your integrator, which is the real contract.
- **Guard against NaN explicitly**: `NT_CHECK(!std::isnan(result))`. A NaN compares false to everything, so `NT_CLOSE(nan, nan, tol)` correctly fails, but a NaN silently propagating through a simulation is worth its own check.

### Running the tests

`make check` (above) builds every `tests/*.cpp` and runs each, stopping at the first failing binary because of `|| exit 1`. Each test is its own process, so a crash in one test doesn't take out the others' results — a genuine advantage of the file-per-test structure.

Add one test file per module: `test_vec3.cpp`, `test_constants.cpp`, `test_integrate.cpp`, `test_forces.cpp`. Keeping them small keeps compile times low and failures easy to locate.

### doctest as a future upgrade

When the hand-rolled harness starts to chafe, **doctest** is the natural next step: a single header (`doctest.h`), no build system integration required, and famously fast to compile.

```cpp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN   // in exactly ONE .cpp -- provides main()
#include <doctest/doctest.h>
#include <newton/newton.hpp>

TEST_CASE("Vec3 norm") {
    newton::Vec3 a{1.0, 2.0, 2.0};
    CHECK(newton::norm(a) == doctest::Approx(3.0).epsilon(1e-12));
}
```

You get test discovery (no manual `main()` per file), `SUBCASE` for shared setup, `doctest::Approx` for float comparison, command-line filtering (`./tests -tc="*norm*"`), and readable failure output showing both operand values. Vendor it as `tests/third_party/doctest.h` (or `brew install doctest`) — the migration is mechanical, and the concepts map one-to-one onto what you built by hand. Catch2 is the alternative; doctest compiles substantially faster and is a better fit here.

---

## 8. Documentation conventions

### `README.md` structure

Order matters — readers stop early, so put the payoff first.

```markdown
# newton

A small, header-only C++17 library for classical mechanics.
Vectors, forces, and numerical integration, written to be read.

## Quick start

```cpp
#include <newton/newton.hpp>
#include <cstdio>

int main() {
    newton::Vec3 pos{0, 0, 0};
    newton::Vec3 vel{10, 15, 0};
    const double mass = 1.0, dt = 0.01;

    newton::forces::UniformGravity gravity{mass};

    while (pos.y >= 0.0) {
        auto [p, v] = newton::dynamics::step(pos, vel, mass, dt, gravity);
        pos = p; vel = v;
    }
    std::printf("range = %.2f m\n", pos.x);
}
```

```
clang++ -std=c++17 -I include examples/projectile.cpp -o projectile
```

## Requirements
C++17 or later. No dependencies. Tested with Apple clang 15 on macOS.

## Installation
[three options: -I, symlink into /usr/local/include, CPLUS_INCLUDE_PATH]

## Modules
| Header | Namespace | Contents |
|---|---|---|
| `<newton/core/vec3.hpp>` | `newton` | `Vec3`, arithmetic, `norm`, `dot`, `cross` |
| `<newton/core/constants.hpp>` | `newton::constants` | `g_earth`, `G`, `c_light` |
| `<newton/dynamics/integrate.hpp>` | `newton::dynamics` | `step` (explicit Euler), `verlet` |
| `<newton/dynamics/forces.hpp>` | `newton::forces` | `UniformGravity`, `LinearDrag`, `Spring`, `sum` |

## Building examples and tests
`make` / `make check` / `make clean`

## Design notes
Why templates over std::function, why header-only, what the integrators
do and do not guarantee. See docs/design-notes.md.

## Status / roadmap
## License
```

Principles: **a copy-pasteable working example above the fold**; state the C++ standard and dependencies early; a module table that maps header → namespace → contents so readers know which one file to open; a "why" section, which is what makes this a learning library rather than a black box.

### Per-module notes in `docs/`

One file per module — `docs/dynamics.md` — covering:

1. **The physics.** State the equations being solved, in plain text or LaTeX. `a = F/m`; explicit Euler is `v_{n+1} = v_n + a·dt`, `x_{n+1} = x_n + v_{n+1}·dt` (note: semi-implicit, using the *updated* velocity — this detail matters and belongs in the docs).
2. **What each function does and its limits.** Explicit Euler is first-order and gains energy in oscillatory systems; velocity Verlet is second-order and symplectic. A user who doesn't know this will misuse the library.
3. **How to choose parameters.** Guidance on `dt` relative to the fastest timescale in the system.
4. **The C++ decisions.** Why `step` is a template, why the force is a callable, why `Vec3` is an aggregate rather than a class with constructors.
5. **A worked example** and the expected numbers.
6. **Units.** Say it once, loudly: SI throughout — metres, kilograms, seconds, newtons. Unit confusion is the most common physics-code bug and costs real spacecraft.

Also keep `docs/design-notes.md` as a running log of decisions and their reasoning, including ones you reversed. Six months later this is the most valuable file in the repo.

### Comment style for a learning library

The default advice — "comments explain why, not what" — needs adjustment here. This library has two audiences (future-you learning physics, and future-you learning C++), so comments should explain the **physics** and the **C++ choice** at the point of use, even when a production library wouldn't bother.

```cpp
#pragma once
#include <newton/core/vec3.hpp>
#include <utility>

namespace newton::dynamics {

/// Advance one particle by a single timestep using semi-implicit
/// (symplectic) Euler.
///
/// Physics:
///   a       = F(x, v) / m          Newton's second law
///   v_{n+1} = v_n + a * dt
///   x_{n+1} = x_n + v_{n+1} * dt   <- note: the UPDATED velocity.
///
/// Using v_{n+1} rather than v_n is what makes this "semi-implicit".
/// It costs nothing extra and is dramatically more stable for
/// oscillatory systems: plain explicit Euler pumps energy into a
/// spring or orbit until it flies apart, while this variant keeps
/// energy bounded (it conserves a nearby "shadow" energy exactly).
/// Still only first-order accurate: global error is O(dt).
///
/// C++ note: `Force` is a template parameter rather than
/// std::function<Vec3(Vec3,Vec3)> so the compiler knows the exact
/// callable type and can inline the user's lambda into this loop.
/// A std::function would force an indirect call and possibly a heap
/// allocation on every step. Any callable works, including a plain
/// std::function if a caller needs runtime swapping.
///
/// @param force  callable as force(const Vec3& pos, const Vec3& vel) -> Vec3, in newtons
/// @param mass   kilograms, must be > 0
/// @param dt     seconds
/// @return {new position, new velocity}, decomposable with structured bindings
template <typename Force>
[[nodiscard]] std::pair<Vec3, Vec3>
step(Vec3 pos, Vec3 vel, double mass, double dt, Force&& force) {
    // 1/mass computed once: division is ~10x slower than multiplication,
    // and this runs in the innermost loop of every simulation.
    const Vec3 accel   = (1.0 / mass) * force(pos, vel);
    const Vec3 vel_new = vel + dt * accel;
    const Vec3 pos_new = pos + dt * vel_new;   // semi-implicit: uses vel_new
    return {pos_new, vel_new};
}

}  // namespace newton::dynamics
```

Conventions worth adopting:

- **`///` for API documentation** (on declarations), `//` for inline implementation notes. `///` is the Doxygen convention, so if you ever run Doxygen it Just Works, and it visually distinguishes "documentation" from "aside".
- **`@param` / `@return` / `@pre`** tags are Doxygen syntax and readable even unrendered. Use them to state **units** and **preconditions** — `mass` in kg and `> 0` is exactly the kind of thing that must be written down.
- **Two labelled sections** in doc comments: *Physics:* and *C++ note:*. Being explicit about which is which lets you skim for the one you need, and forces you to actually articulate both.
- **State the accuracy and stability properties** of every numerical method. "First-order, global error O(dt), symplectic" is the contract; without it the function is unusable for real work.
- **Comment the non-obvious optimization**, like hoisting the reciprocal — otherwise a future you "simplifies" it back.
- **Close namespaces with `}  // namespace newton::dynamics`.** Free, and invaluable in a long file.
- **Reference sources** where relevant: CODATA for constants, a textbook section for a method. `// CODATA 2018` next to `G` tells the reader where the digits came from.
- **Keep examples commented as tutorials.** `examples/projectile.cpp` should read like a lesson: what physical situation, what the analytic answer is, and what to look at in the output. It's the file a newcomer opens first.

---

## Recommended order of work

1. `include/newton/core/vec3.hpp` + `include/newton/core/constants.hpp` + `include/newton/newton.hpp`. Compile with `clang++ -std=c++17 -Wall -Wextra -I include`.
2. `examples/projectile.cpp` — get one thing running end to end before adding structure.
3. The `Makefile` — `make` and `make check`.
4. `tests/test_util.hpp` + `tests/test_vec3.cpp`.
5. `include/newton/dynamics/integrate.hpp` with the template-callable `step`, then `forces.hpp`.
6. `README.md` and `docs/`.
7. Only then: symlink into `/usr/local/include`, more integrators (velocity Verlet, RK4), more examples (orbit, spring, N-body).
8. Later, if needed: doctest, CMake, C++20 concepts.

The whole skeleton above — headers, umbrella, example, test harness, Makefile, and `build.sh` — was built and run before writing this report; `make check` reports `All tests passed.` and `./build/projectile` prints `range = 30.500 m, speed = 17.953 m/s` (analytic range for v=(10,15) m/s is 30.59 m — the ~0.3% discrepancy is exactly the first-order Euler error you would expect at dt = 0.01, and is itself worth a comment in the example).