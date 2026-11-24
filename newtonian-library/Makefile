# newton — Makefile
#
# A header-only library has no library to build: this file's only jobs are
# to compile each examples/*.cpp and tests/*.cpp to its own binary, run the
# tests, and clean up. See docs/research/cpp-library-design.md §6 for the
# tutorial this file is based on (make syntax, why `-I include`, etc).
#
# Usage:
#   make            # build every example into build/
#   make tests      # build every test into build/tests/ (but don't run them)
#   make check      # build AND run every test; nonzero exit if any fails
#   make clean      # delete build/
#   make debug      # rebuild everything with sanitizers + -O0 -g, then check
#   make install    # copy include/newton to $(PREFIX)/include (default /usr/local)
#   make uninstall  # remove it again
#   make -j8 check  # any target, built in parallel across 8 cores

# ---- configuration ----------------------------------------------------
# `?=` sets a variable only if it isn't already set in the environment, so
# `make CXX=g++-14` (or `CXX=g++-14 make`) overrides this without editing
# the file.
CXX      ?= clang++
CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -O2 -I include

# `?=` again: `make install PREFIX=$HOME/.local` works with no edit here.
PREFIX ?= /usr/local

# ---- file discovery -----------------------------------------------------
# $(wildcard ...) expands a glob into a list at the moment make reads this
# file — add examples/07_whatever.cpp tomorrow and it's picked up with no
# Makefile change. $(patsubst pattern,replacement,list) rewrites each
# element: examples/foo.cpp -> build/foo, tests/foo.cpp -> build/tests/foo.
EXAMPLE_SRCS := $(wildcard examples/*.cpp)
EXAMPLE_BINS := $(patsubst examples/%.cpp,build/%,$(EXAMPLE_SRCS))
TEST_SRCS    := $(wildcard tests/*.cpp)
TEST_BINS    := $(patsubst tests/%.cpp,build/tests/%,$(TEST_SRCS))

# Every binary depends on every header: because the library is header-only
# there is no separate compilation step where a header change would only
# affect the .o files that include it — touching ANY header must rebuild
# ANY binary, so we just list them all as a blanket prerequisite below.
# $(shell ...) runs a real shell command once, at Makefile-read time.
HEADERS := $(shell find include -name '*.hpp')

# ---- targets --------------------------------------------------------------
# .PHONY tells make these target names are not real files. Without this, a
# stray file literally named `clean` or `check` in the repo root would
# short-circuit its rule (make would see the file "up to date" and do
# nothing) — declaring them phony makes the recipe run every time regardless.
.PHONY: all examples tests check clean debug install uninstall

all: examples

examples: $(EXAMPLE_BINS)
tests:    $(TEST_BINS)

# Pattern rules: `%` matches the same stem in the target and the
# prerequisite. `$<` is the first prerequisite (the .cpp file), `$@` is the
# target (the binary path). `| build` / `| build/tests` are ORDER-ONLY
# prerequisites: the output directory must exist before the recipe runs, but
# its own timestamp changing (e.g. adding another file into it) must NOT
# force every binary in it to be relinked — that's exactly what a normal
# (non-order-only) prerequisite would do.
#
# make picks the more specific of two matching pattern rules by shortest
# stem, so `tests/vec3.cpp` matches `build/tests/%: tests/%.cpp` (stem
# "vec3") rather than `build/%: examples/%.cpp` (which would need a stem of
# "tests/vec3" and a nonexistent examples/tests/vec3.cpp) — this is verified
# behaviour, not a happy accident (docs/research/cpp-library-design.md §6).
build/%: examples/%.cpp $(HEADERS) | build
	$(CXX) $(CXXFLAGS) $< -o $@

build/tests/%: tests/%.cpp tests/test_util.hpp $(HEADERS) | build/tests
	$(CXX) $(CXXFLAGS) $< -o $@

# build/tests is nested inside build, so creating it also creates build/ —
# `mkdir -p` makes both existing-already and missing-parents non-errors.
build build/tests:
	mkdir -p $@

# check builds every test (via the `tests` target above) then runs each
# binary in turn. `$$t` is a literal `$t` handed to the shell (make eats a
# single `$`); `|| exit 1` stops the whole recipe — and fails the target,
# giving `make check` a nonzero exit code — the instant one test binary
# returns nonzero, per test_util.hpp's report() contract.
check: tests
	@for t in $(TEST_BINS); do ./$$t || exit 1; done
	@echo "All tests passed."

clean:
	rm -rf build

# Target-specific variable: this CXXFLAGS override is in effect only while
# building `debug` and everything `debug` depends on (transitively —
# `clean` and `check`, and everything `check` in turn needs: `tests`, then
# every build/tests/% pattern rule) — it does not leak into a plain `make
# check` run afterward. `-fsanitize=address,undefined` catches out-of-bounds
# access and undefined behaviour at ~2x runtime cost; `-g -O0` keeps stack
# traces and variable names readable if something does trip a sanitizer.
# `clean` first guarantees every binary is actually recompiled with the
# sanitizer flags rather than reusing an already-up-to-date -O2 build.
debug: CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -g -O0 -fsanitize=address,undefined -I include
debug: clean check

# Copies (not symlinks) include/newton into $(PREFIX)/include, so programs
# elsewhere on the machine can `#include <newton/newton.hpp>` with no `-I`
# flag. Run as `sudo make install` if $(PREFIX) needs root (the default
# /usr/local usually does on macOS). See
# docs/research/cpp-library-design.md §2 for the symlink alternative, better
# suited to active development of the library itself.
install:
	mkdir -p $(PREFIX)/include
	cp -R include/newton $(PREFIX)/include/

uninstall:
	rm -rf $(PREFIX)/include/newton
