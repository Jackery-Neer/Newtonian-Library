#!/usr/bin/env bash
# build.sh — compiles every example and every test, then runs the tests.
# A plain-shell alternative to `make`/`make check` that needs no tool beyond
# bash itself. If you have never read a .sh file before, read this one
# top to bottom — every non-obvious line is commented.

# The "shebang" line above tells the OS which interpreter runs this file.
# `/usr/bin/env bash` searches PATH for `bash` instead of hardcoding a path
# like `/bin/bash` — on macOS, `/bin/bash` is a very old bash (3.2, kept for
# licensing reasons), so `env bash` is what finds a newer one if you've
# installed one (e.g. via `brew install bash`).

set -euo pipefail
# set -e            exit this script immediately if any command fails (i.e.
#                   returns a nonzero exit status) instead of the default
#                   shell behaviour of plowing ahead to the next line.
# set -u            treat using an undefined variable as an error rather
#                   than silently substituting an empty string — catches a
#                   surprising number of typos (e.g. $CXXFALGS).
# set -o pipefail   in a pipeline `a | b`, the pipeline's exit status is
#                   normally just b's; pipefail makes it fail if ANY stage
#                   fails, not only the last one. (No pipelines below yet,
#                   but it's cheap insurance and standard practice.)
# These three together are sometimes called "bash strict mode."

cd "$(dirname "$0")"
# $0 is the path this script was invoked with (e.g. "./build.sh" or
# "/some/where/build.sh"). `dirname` strips the filename, leaving the
# directory the script lives in — the repo root, since build.sh sits there.
# `cd`-ing there means `./build.sh` behaves identically no matter what
# directory you were sitting in when you ran it. The double quotes matter:
# without them, a path containing a space would be split into multiple
# words and `cd` would see the wrong (or multiple) arguments.

CXX="${CXX:-clang++}"
# ${VAR:-default} means "use $VAR if it's set and non-empty, otherwise use
# `default`" — this does NOT triger `set -u`'s undefined-variable check,
# because the `:-` form is exactly the shell's built-in way of providing a
# fallback for a variable that might not exist. Run `CXX=g++-14 ./build.sh`
# to override it for one invocation without editing this file.

CXXFLAGS=(-std=c++17 -Wall -Wextra -Wpedantic -O2 -I include)
# A bash ARRAY (the parentheses make it one), not a plain string. Using an
# array — and later expanding it as "${CXXFLAGS[@]}" — keeps each flag as
# its own separate word even if a future flag ever contained a space; a
# plain string variable would need to be word-split by the shell, which is
# exactly the kind of fragile behaviour `set -u`-style "strict mode"
# scripts try to avoid elsewhere too.

mkdir -p build/tests
# -p: create build/ and build/tests/ together, and don't error if either
# (or both) already exist — the flag this repo's .gitignore relies on,
# since build/ is never checked into git.

echo "== compiling examples =="
for src in examples/*.cpp; do
  # Bash expands `examples/*.cpp` to the list of matching files BEFORE the
  # loop starts. If the directory genuinely has no .cpp files, the glob
  # doesn't match anything and bash (by default) leaves the literal text
  # "examples/*.cpp" as the one "file" — the next line guards against that
  # so an empty examples/ doesn't try to compile a file that doesn't exist.
  [ -e "$src" ] || continue

  # ${src%.cpp} strips a trailing ".cpp" suffix (the shortest match, from
  # the % operator) off $src, e.g. "examples/01_free_fall.cpp" becomes
  # "examples/01_free_fall". `basename` then strips the directory part too,
  # leaving just "01_free_fall" — so the binary lands directly in build/,
  # not build/examples/.
  out="build/$(basename "${src%.cpp}")"

  echo "  $src -> $out"
  # "${CXXFLAGS[@]}" expands the array to one shell word per element — the
  # quotes are what make this safe; "${CXXFLAGS[*]}" (with a star) would
  # instead glue everything into a single word, which is NOT what a
  # compiler command line wants. This is the array-equivalent of always
  # quoting a plain string variable.
  "$CXX" "${CXXFLAGS[@]}" "$src" -o "$out"
done

echo "== compiling and running tests =="
status=0
# Track failure in a variable rather than exiting immediately (which `set
# -e` would otherwise do on the first failing test) — the whole point of
# running a test SUITE is to see every failure in one pass, not just the
# first one. `"$out" || status=1` catches a nonzero exit from the test
# binary itself and records it, without tripping `set -e`, because the
# failure is the left side of an `||`, and `set -e` only fires on a bare
# failing command, not one that's part of a larger `&&`/`||` expression.
for src in tests/*.cpp; do
  [ -e "$src" ] || continue
  out="build/tests/$(basename "${src%.cpp}")"
  echo "  $src -> $out"
  "$CXX" "${CXXFLAGS[@]}" "$src" -o "$out"
  "$out" || status=1
done

exit "$status"
# A nonzero exit code here is what tells a human (or a CI system) at a
# glance that something failed, without having to scroll back through the
# output — `echo $?` right after running this script shows the same number.
