#pragma once
// tests/test_util.hpp — the library's whole test framework: three small
// functions and two macros. Not a framework in the Catch2/doctest sense —
// just enough to record a pass/fail, report where a failure happened, and
// produce a nonzero exit code so `make check` fails loudly. Every
// tests/test_*.cpp file includes this (quoted: it sits right next to them)
// alongside <newton/newton.hpp> (angle brackets: reached via -I include).
//
// C++ note: kept dependency-free on purpose — no #include of any newton/
// header — so this file works for testing anything, not just newton types.
// close_vec3() below is templated on any type with public .x/.y/.z members
// rather than hardcoded to newton::Vec3, for the same reason.

#include <cmath>
#include <cstdio>

namespace newton_test {

// C++17 inline variables: one shared counter per program, even though this
// header may be #included (separately) by several test .cpp translation
// units — each test file is its own executable, so each gets its own fresh
// counters (see core/constants.hpp's C++ note on why `inline` is what makes
// this safe rather than an ODR violation).
inline int checks_run = 0;
inline int failures = 0;

/// Floating-point comparison with a mixed absolute/relative tolerance.
/// Never compare computed doubles with `==` — 0.1 + 0.2 != 0.3 exactly,
/// because binary floating point cannot represent those decimals exactly.
///   - the absolute test (diff <= tol) lets values near zero compare
///     sensibly at all — a purely relative test is meaningless when the
///     expected value is 0 (everything within tol*0 = 0 would fail);
///   - the relative test (diff <= tol * scale) keeps large values from being
///     held to an impossible absolute precision — a double carries ~15-16
///     significant decimal digits, so a value near 1e6 simply cannot match
///     to an absolute 1e-12.
/// A value passes if EITHER test succeeds, so `tol` reads naturally as "this
/// many digits of agreement, in whichever sense makes sense for the size of
/// the numbers involved."
[[nodiscard]] inline bool close(double a, double b, double tol) {
    const double diff = std::fabs(a - b);
    if (diff <= tol) return true;  // absolute test
    const double scale = std::fmax(std::fabs(a), std::fabs(b));
    return diff <= tol * scale;  // relative test
}

/// Component-wise close() for any 3-component vector type (newton::Vec3, or
/// anything else with public .x/.y/.z) — every test file in this suite
/// compares positions/velocities/forces, so this is the one genuinely
/// necessary extension beyond the verified scalar harness.
template <class V>
[[nodiscard]] inline bool close_vec3(const V& a, const V& b, double tol) {
    return close(a.x, b.x, tol) && close(a.y, b.y, tol) && close(a.z, b.z, tol);
}

inline void record(bool ok, const char* expr, const char* file, int line) {
    ++checks_run;
    if (!ok) {
        ++failures;
        std::printf("  FAIL %s:%d  %s\n", file, line, expr);
    }
}

/// Call at the end of main() and return its value: 0 (success) iff every
/// check in this translation unit passed, nonzero otherwise — that's what
/// makes a failing test fail the `check` target and any CI that runs it.
[[nodiscard]] inline int report(const char* suite) {
    std::printf("%s: %d checks, %d failed\n", suite, checks_run, failures);
    return failures == 0 ? 0 : 1;
}

}  // namespace newton_test

// Macros, not functions, for two things only a macro can do:
//   1. #expr stringifies the source text of the condition, so a failure
//      prints exactly what was checked, e.g. "norm(a) == 3.0";
//   2. __FILE__ / __LINE__ expand at the *call site* — inside a function
//      they would always report test_util.hpp's own location instead of the
//      line in the test file that actually failed.
// The leading `::` qualifies from the global namespace so these expand
// correctly no matter what namespace surrounds the call site.
//
// Deliberately non-fatal: record() does not throw or abort, so a failed
// check does not stop the rest of the test file from running — one run of
// one binary tells you everything that broke, not just the first thing.
//
// C++ note: NT_CHECK is variadic (`...`/__VA_ARGS__) rather than a plain
// single-parameter macro. Reason: the preprocessor splits macro arguments on
// every top-level comma it sees, and it only understands *parentheses* as
// grouping — a brace-init list like `Vec3{1.0, 2.0, 3.0}`, which this test
// suite writes constantly, has top-level commas of its own. A one-parameter
// macro would misparse `NT_CHECK(v == Vec3{1.0, 2.0, 3.0})` as three
// arguments and fail to compile. `__VA_ARGS__` re-glues everything between
// the outer parens back together, commas included, so braces "just work".
#define NT_CHECK(...) ::newton_test::record((__VA_ARGS__), #__VA_ARGS__, __FILE__, __LINE__)

#define NT_CLOSE(a, b, tol)                                                       \
    ::newton_test::record(::newton_test::close((a), (b), (tol)), #a " ~= " #b, __FILE__, \
                          __LINE__)

#define NT_CLOSE_VEC3(a, b, tol)                                                       \
    ::newton_test::record(::newton_test::close_vec3((a), (b), (tol)), #a " ~= " #b, __FILE__, \
                          __LINE__)
