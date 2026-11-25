# 00 — Start here

This note is a map of the repository: what every file and directory is for,
how to build and run things, and how to make `#include <newton/newton.hpp>`
work from outside this directory once you want that. It assumes you know C++
already; it does not assume you've built a *library* (as opposed to a single
program) before, so it explains the parts that are specific to that: the
`Makefile`, `build.sh`, the include-path mechanics, and so on.

## The guided tour

```
newtonian-library/
├── DESIGN.md              the spec: every header's exact signatures
├── README.md               the front door — what the library is, quick start
├── Makefile                 build recipes for `make`
├── build.sh                  the same job, as a plain shell script
├── .gitignore                 tells git which files not to track
├── include/newton/             the library itself — this is all there is
├── examples/                    small standalone programs, one per topic
├── tests/                        assertion-based checks, one file per module
├── scripts/                       plot_csv.py, a plotting helper
├── docs/research/                   the design research this was built from
├── notes/                            you are here
└── build/                             generated binaries (never committed)
```

### `include/newton/` — the library

This is the actual product. Everything else in the repo exists to build,
test, document, or demonstrate what's in here. It's organized by module,
mirroring the namespaces:

| Directory | Namespace | What's in it |
|---|---|---|
| `core/` | `newton` | `Vec3`, `Particle`, `State`, `Trajectory`, `constants` |
| `coords/` | `newton::coords` | Cylindrical and spherical polar coordinates, and conversions |
| `fields/` | `newton::fields` | Numerical gradient/curl/divergence and line integrals of fields |
| `potentials/` | `newton::potentials` | Is this force conservative? U from F, and F from U |
| `forces/` | `newton::forces` | Force-law functors (gravity, drag, springs, ...) |
| `integrators/` | `newton::integrators` | The four `xxx_step` one-step integrators |
| (top-level) | `newton` | `simulate.hpp` — the time-stepping driver |
| `analytic/` | `newton::analytic` | Closed-form textbook solutions |
| `analysis/` | `newton::analysis` | Energy/momentum diagnostics on a `Trajectory` |
| `systems/` | `newton::systems` | Self-contained special-purpose simulators (rocket) |
| `io/` | `newton::io` | CSV export |

Every file in here is a **header** (`.hpp`), and every header ends with
`#pragma once` at the top. There are no `.cpp` files in `include/` — see the
next section for why that's possible at all.

### `.hpp` vs `.cpp` — why this library only has the former

A normal C++ program splits code into two file kinds:

- **`.cpp` (a translation unit)** — compiled *once*, independently, into an
  object file (`.o`). The compiler only ever sees one `.cpp` file at a time
  (plus whatever it `#include`s); a project with ten `.cpp` files runs the
  compiler ten times and a linker glues the results together.
- **`.hpp`** — never compiled by itself. It's text that gets pasted into
  whatever `.cpp` file `#include`s it, by the preprocessor, before the real
  compiler even starts. A header *declares* things; traditionally a matching
  `.cpp` *defines* them once, so the linker has exactly one definition to
  find.

`newton` breaks that second half of the pattern on purpose: **every**
definition — every function body, every struct — lives in a `.hpp`, and
there is no matching `.cpp` anywhere in `include/`. That's what "header-only"
means, and it's why using the library needs nothing but `-I include`: no
separate library file to compile, archive, or link against. `notes/01-cpp-
for-this-library.md` explains the C++ mechanism (`inline`, templates, ODR)
that makes this legal instead of a pile of linker errors.

The `.cpp` files that *do* exist in this repo — everything under
`examples/` and `tests/` — are ordinary translation units. Each one
`#include`s the headers it needs and defines exactly one `main()`; each
compiles to its own standalone binary.

### `Makefile` — what `make` actually does

`make` is a program that reads a file literally named `Makefile` in the
current directory and uses it to decide **what needs rebuilding**, by
comparing file timestamps. You don't tell it *how* to build the whole
project step by step; you describe a set of rules — "to produce this file,
you need these other files, and here's the shell command that makes it" —
and `make` works out an order and skips anything that's already up to date.

A rule looks like:

```make
build/01_free_fall: examples/01_free_fall.cpp $(HEADERS) | build
	$(CXX) $(CXXFLAGS) $< -o $@
```

Read as: *"to produce `build/01_free_fall`, you need
`examples/01_free_fall.cpp` and every header; if any of those is newer than
`build/01_free_fall` (or it doesn't exist yet), run the recipe line."* The
recipe line **must** start with a literal tab character — that's an old and
famously sharp edge of `make`; a Makefile indented with spaces fails with
`missing separator`.

Because this is a header-only library, there's no compiled library to
produce — the Makefile's only jobs are: compile each `examples/*.cpp` and
`tests/*.cpp` to its own binary in `build/`, run the tests, and clean up.
Targets you'll actually type:

```
make            # build every example into build/
make check      # build every test, run each one, fail if any test fails
make clean      # rm -rf build — start over from nothing
make debug      # rebuild with -O0 -g and the sanitizers on (see below)
sudo make install    # copy include/newton into /usr/local/include
sudo make uninstall  # remove it again
```

`HEADERS := $(shell find include -name '*.hpp')` inside the Makefile is
worth understanding: because the library is header-only, there's no
"has the library changed" question separate from "has any header changed" —
so *every* example and test binary depends on *every* header, and touching
any one of them triggers a full rebuild on the next `make`. `$(wildcard
examples/*.cpp)` is the other piece of the trick: it means adding
`examples/09_something.cpp` needs no Makefile edit at all — the pattern rule
picks it up automatically the next time you run `make`.

### `build.sh` — what a `.sh` file is, and why this repo has one alongside a Makefile

A `.sh` file is a **shell script**: a plain text file of shell commands, run
top to bottom by an interpreter, the same way you'd type them one at a time
at a terminal prompt. The first line (the *shebang*) names which
interpreter:

```bash
#!/usr/bin/env bash
```

To run it, either `chmod +x build.sh` once (this sets the file's "you may
execute this" permission bit) and then invoke it as `./build.sh` — the
leading `./` is required because, unlike Windows, the current directory is
not on your shell's search path by default — or skip the `chmod` and run
`bash build.sh` directly.

`build.sh` does the same job as `make`/`make check`, just written as a plain
linear script instead of a dependency graph: `set -euo pipefail` at the top
(stop on the first error, treat unset variables as errors, and make a
pipeline fail if any stage of it fails — three safety nets that are *not*
bash's default behavior), then a loop that compiles every file in
`examples/*.cpp` and `tests/*.cpp` and runs the tests. Compare it side by
side with the Makefile rule shown above if you want to see the same build
expressed both ways.

Why keep both? `make` gives you incremental builds (only recompiles what
changed) and easy parallelism (`make -j8`) essentially for free; `build.sh`
gives you a script you can read top to bottom with no `make`-specific syntax
to learn, which is a gentler on-ramp and easier to modify with arbitrary
logic (e.g. "also copy the CSVs somewhere" is a two-line addition to a
script, and a slightly more awkward one in a Makefile). Many real projects
keep both for exactly this reason — a Makefile for the day-to-day build,
small shell scripts for one-off orchestration tasks.

### `.md` files

Markdown — plain text with a light formatting syntax (`#` headings,
`` `code` ``, ```` ``` ```` fenced code blocks, `-` lists, `|` tables). It
renders as styled text on GitHub and in most editors, but — unlike a Word
document — it's still fully readable as plain text if nothing renders it,
and it diffs cleanly in `git` (you can see exactly which sentence changed in
a commit, instead of an opaque binary diff). Every prose file in this repo
(`README.md`, `DESIGN.md`, everything in `docs/` and `notes/`) is Markdown
for that reason.

### `.gitignore`

A list of filename patterns `git` should never track or offer to add, even
if they exist on disk. This repo's covers generated/local-machine content
that has no business in version control: `build/` (binaries, regenerated by
`make` any time), `.DS_Store` (a macOS Finder metadata file, not something
you ever want to commit), and `*.csv` (the example programs' output data,
also regenerated on demand). Without it, `git status` would constantly show
you dozens of rebuildable files as "untracked."

### `DESIGN.md`

The authoritative specification for the whole library: exact signatures,
exact file paths, exact semantics, cited section-by-section against
`docs/research/`. If code and `DESIGN.md` ever disagree, `DESIGN.md` is
right (or the code has a bug) — it is not just background reading, it is the
contract every header in `include/newton/` was written to satisfy.

### `docs/research/`

Three long reference reports, produced *before* any code was written, that
this library's implementation is built from:

- **`taylor-ch1-7-map.md`** — every physics formula used anywhere in the
  library, chapter by chapter through Taylor's *Classical Mechanics*, with
  the exact equation each header implements.
- **`numerical-integration.md`** — the integrator update equations
  (`euler_step`, `semi_implicit_euler_step`, `velocity_verlet_step`,
  `rk4_step` all follow this document's §2 line for line), plus the numerical
  pitfalls that motivated specific design choices you'll see called out in
  the headers (e.g. `simulate.hpp`'s "never accumulate `t += dt`" rule is
  its §7.2).
- **`cpp-library-design.md`** — the C++ and build-tooling conventions this
  whole repo follows: directory layout, `#pragma once` vs. include guards,
  `inline`/ODR mechanics, the template-callable-vs-`std::function` decision,
  the Makefile and `build.sh` this repo's actual build files are modeled on,
  and the hand-rolled test harness in `tests/test_util.hpp`.

If you want to know *why* a header does something a particular way rather
than *what* it does, one of these three documents is almost always the
answer — `notes/01-cpp-for-this-library.md` leans especially heavily on the
third one.

### `notes/` (this directory)

Explanations aimed at you, the library's author, rather than at a
hypothetical external user — the gap `README.md` and `DESIGN.md`
deliberately don't fill. See "Reading order" below.

## Building and running

Everything below assumes your terminal's current directory is the repo
root (`newtonian-library/`) — both `make` and `build.sh` look for files
relative to that directory (`build.sh` explicitly `cd`s to its own location
first, so it's forgiving about where you launch it from; plain `make`
is not, unless you pass `-C`).

```
make            # compile every examples/*.cpp into build/
make check      # compile + run every tests/*.cpp; nonzero exit if any fails
./build.sh      # the same two jobs, as a shell script instead of a Makefile
```

The very first time, `chmod +x build.sh` before `./build.sh` will be
necessary if the execute bit isn't already set (see the `build.sh` section
above for why).

### The bare compile command, flag by flag

Both `make` and `build.sh` are just wrappers that eventually run a command
that looks like this, once per `.cpp` file:

```
clang++ -std=c++17 -Wall -Wextra -Wpedantic -O2 -I include examples/01_free_fall.cpp -o build/01_free_fall
```

| Flag | Meaning |
|---|---|
| `clang++` | The C++ compiler driver (on macOS: Apple clang, from the Command Line Tools). It runs the preprocessor, compiler, assembler, and linker in sequence and hides all of that behind one command. `clang` alone would compile as C, so always `clang++` for this project. |
| `-std=c++17` | Selects the C++17 language standard explicitly. Without it, Apple clang defaults to `gnu++17` (C++17 plus GNU extensions) — stating it removes that ambiguity and is what makes `inline` variables, `namespace newton::forces { }`, and `[[nodiscard]]` legal. |
| `-Wall` | Turns on a standard set of warnings (unused variables, missing `return`s, and more) — despite the name, not literally *all* warnings. |
| `-Wextra` | More warnings on top of `-Wall`: unused function parameters, sign-comparison mismatches, and more. `-Wall -Wextra` together is this project's floor. |
| `-Wpedantic` | Warns about non-standard compiler extensions creeping into the code. |
| `-O2` | Optimization level 2 — real-speed code. (`-O0`, the default with no flag, compiles faster and is easier to step through in a debugger; always benchmark at `-O2` or higher.) |
| `-I include` | Adds `include/` to the header search path used for `#include <...>` (angle-bracket) lookups — this is the *entire* mechanism that makes `#include <newton/newton.hpp>` resolve to `include/newton/newton.hpp`. Must be run from the repo root, since the path is relative. |
| `examples/01_free_fall.cpp` | The one translation unit being compiled. Because the library is header-only, this `.cpp` plus everything it `#include`s *is* the whole program — there's no separate library object to link in. |
| `-o build/01_free_fall` | Output binary path. Without `-o`, clang would write `a.out` into the current directory instead. |

Two more flags worth knowing even though they're not in the default build:

- **`-g`** — emit debug info so a debugger (`lldb`) can show source lines and
  variable names. Pairs with `-O0`, not `-O2` (optimized code reorders and
  eliminates variables, which makes a debugger's view confusing).
- **`-fsanitize=address,undefined`** — compiles in AddressSanitizer and
  UndefinedBehaviorSanitizer, which catch out-of-bounds access and
  undefined-behavior bugs (signed overflow, bad casts, etc.) at roughly 2×
  runtime cost. This is what `make debug` turns on.

If an `#include` ever fails to resolve and you can't see why, the ground-
truth diagnostic is:

```
echo | clang++ -std=c++17 -x c++ -E -v -
```

`-v` prints the exact search path clang used, in order, with anything from
`-I` listed first.

## Three ways to make `#include <newton/newton.hpp>` work from anywhere

Inside this repo, `-I include` (baked into the Makefile and `build.sh`) is
all you need. Once you want to `#include <newton/newton.hpp>` from a
scratch file that lives *outside* this repo — while working a Taylor
problem in some other directory, say — you have three options, in the order
you should reach for them:

1. **Keep using `-I` explicitly, every time (do this first, for a long time).**
   A one-line function in `~/.zshrc`:
   ```bash
   nc++() { clang++ -std=c++17 -Wall -Wextra -I "$HOME/Github/newtonian-library/include" "$@"; }
   ```
   Then `nc++ scratch.cpp -o scratch` works from any directory. This is the
   safest option: the path is right there in your shell config and in every
   command you type, nothing is hidden, and if it breaks you can see
   exactly why.

2. **Symlink `include/newton` into `/usr/local/include` once the API has
   stabilized.** `/usr/local/include` is on Apple clang's built-in search
   path already — but on a fresh Mac the directory doesn't exist yet, so
   clang silently skips it (`ignoring nonexistent directory
   "/usr/local/include"` in the `-v` output above confirms this on this
   machine). Create it and link in:
   ```bash
   sudo mkdir -p /usr/local/include
   sudo ln -s "$PWD/include/newton" /usr/local/include/newton
   ```
   A **symlink** (not a copy) means edits to files under this repo's
   `include/newton/` are instantly visible system-wide with no reinstall
   step — exactly what you want while the headers are still changing. `make
   install` (`sudo make install`) does a `cp -R` instead of a symlink, which
   is a deliberate choice for a *frozen, versioned* install rather than a
   live one; prefer the symlink for your own day-to-day use and save `make
   install`/`uninstall` for "pin a stable copy independent of my working
   tree."

3. **`CPLUS_INCLUDE_PATH` (know it exists; use it only for quick, temporary
   experiments).**
   ```bash
   export CPLUS_INCLUDE_PATH="$HOME/Github/newtonian-library/include:$CPLUS_INCLUDE_PATH"
   ```
   clang and gcc both honor this colon-separated environment variable as an
   extra angle-bracket search path. It works, but it's *invisible* — a
   build that succeeds in your shell fails with no explanation for anyone
   (including future you, in a fresh terminal or a CI job) who doesn't have
   the same variable set. It also marks the directory as a *system* include
   directory, which quietly suppresses compiler warnings from inside it —
   fine for a third-party library, actively unhelpful for your own code that
   you want warnings from. Reach for option 1 or 2 instead for anything that
   needs to keep working tomorrow.

## Reading order for the rest of `notes/`

1. **`notes/00-start-here.md`** — this file. The map.
2. **`notes/01-cpp-for-this-library.md`** — every non-basic C++ feature the
   library's headers actually use (fold expressions, `inline constexpr`,
   template callables, and the rest), each with a real snippet and file
   path, so that reading the headers themselves afterward is smooth instead
   of stopping every few lines to look something up.
3. **`notes/02-physics-to-code.md`** — the bridge from a Taylor section to
   the function that implements it, chapter by chapter through ch. 1–5.
4. **`notes/03-numerical-methods.md`** — the four integrators
   (`euler_step`, `semi_implicit_euler_step`, `velocity_verlet_step`,
   `rk4_step`) in depth: what each one buys you and costs you, and how to
   pick.
5. **`notes/04-extending.md`** — a semester-length roadmap for growing the
   library past v0.1, including the v0.2 plan from `DESIGN.md`'s roadmap
   section (N-particle systems, adaptive step sizes, and more) that is
   deliberately *not* implemented here.
6. **`notes/05-library-reference.md`**: every function in the installed
   library with a tested usage example, how C++ libraries are traditionally
   built, a tour of the rest of the repo, and how the tests work. Keep it
   open while you write code that uses the library.
7. **`notes/06-fields-and-potentials.md`**: the physics, numerics, and C++
   behind `newton::fields` and `newton::potentials`, for when Taylor chapter 4
   is in front of you: how a numerical derivative and a line integral work,
   how `check_conservative` decides whether a force is conservative (and what
   it cannot prove), why a time-dependent potential makes the energy change,
   and how `energy_balance` keeps the books for drag. Read it alongside
   `examples/07_potential_from_force.cpp` and
   `examples/08_force_from_potential.cpp`.

1–2 are C++-focused and apply even if you're just skimming the source; 3–4
assume you're actively working a Taylor problem and want the paper-to-code
mapping; 5 is for when v0.1 starts feeling limiting; 6 is a lookup table to
keep open; 7 is the chapter 4 companion.
