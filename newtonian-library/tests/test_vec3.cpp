// tests/test_vec3.cpp — arithmetic, dot/cross identities, norms, the
// normalized() zero-guard. Vec3 is the foundation every other header builds
// on, so this is the library's first smoke test: if this file fails,
// nothing downstream can be trusted either.
#include "test_util.hpp"        // quotes: sits next to this file
#include <newton/newton.hpp>    // angle brackets: reached via -I include

using newton::Vec3;
using newton::cross;
using newton::dot;
using newton::norm;
using newton::norm2;
using newton::normalized;

int main() {
    // --- construction & equality --------------------------------------
    // Default member initializers must zero every component (core/vec3.hpp's
    // whole "no NaN surprises" promise) — exact == is legitimate here
    // because nothing was computed, only assigned/defaulted.
    Vec3 zero{};
    NT_CHECK(zero == Vec3{0.0, 0.0, 0.0});
    NT_CHECK(zero.x == 0.0 && zero.y == 0.0 && zero.z == 0.0);

    Vec3 a{1.0, 2.0, 3.0};
    Vec3 b{4.0, 5.0, 6.0};
    NT_CHECK(a != b);
    NT_CHECK(a == Vec3{1.0, 2.0, 3.0});

    // --- arithmetic -----------------------------------------------------
    NT_CHECK(a + b == Vec3{5.0, 7.0, 9.0});
    NT_CHECK(b - a == Vec3{3.0, 3.0, 3.0});
    NT_CHECK(a * 2.0 == Vec3{2.0, 4.0, 6.0});
    NT_CHECK(2.0 * a == Vec3{2.0, 4.0, 6.0});       // both operand orders (vec3.hpp's whole point)
    NT_CHECK(a * 2.0 == 2.0 * a);
    NT_CHECK(Vec3{2.0, 4.0, 6.0} / 2.0 == a);
    NT_CHECK(-a == Vec3{-1.0, -2.0, -3.0});

    Vec3 c = a;
    c += b;
    NT_CHECK(c == Vec3{5.0, 7.0, 9.0});
    c -= b;
    NT_CHECK(c == a);
    c *= 3.0;
    NT_CHECK(c == Vec3{3.0, 6.0, 9.0});
    c /= 3.0;
    NT_CHECK(c == a);

    // --- dot product ------------------------------------------------------
    // Perpendicular unit vectors: a.b = 0.
    NT_CHECK(dot(Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}) == 0.0);
    // Parallel vectors: a.a = |a|^2.
    NT_CLOSE(dot(a, a), norm2(a), 1e-12);
    // Hand-computed: (1,2,3).(4,5,6) = 4+10+18 = 32.
    NT_CLOSE(dot(a, b), 32.0, 1e-12);
    NT_CHECK(dot(a, b) == dot(b, a));  // commutative

    // --- cross product ------------------------------------------------
    // Right-hand rule: x_hat cross y_hat = z_hat.
    NT_CHECK(cross(Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}) == Vec3{0.0, 0.0, 1.0});
    // Anti-commutative: a x b = -(b x a).
    NT_CHECK(cross(a, b) == -cross(b, a));
    // a x a = 0 for any a (parallel vectors span no area).
    NT_CHECK(cross(a, a) == Vec3{});
    // The defining geometric property: a x b is perpendicular to BOTH a and
    // b — this is the check DESIGN.md calls out explicitly, and it's the one
    // that actually verifies the formula (rather than just its symmetries).
    const Vec3 axb = cross(a, b);
    NT_CLOSE(dot(axb, a), 0.0, 1e-12);
    NT_CLOSE(dot(axb, b), 0.0, 1e-12);

    // --- norms ----------------------------------------------------------
    // The 1-2-2 "Pythagorean-style" triple: sqrt(1+4+4) = sqrt(9) = 3
    // exactly (9 is a perfect square, so IEEE-754 sqrt rounds to exactly
    // 3.0) — DESIGN.md's specific worked example, tight tolerance is fine.
    NT_CLOSE(norm(Vec3{1.0, 2.0, 2.0}), 3.0, 1e-12);
    NT_CLOSE(norm2(Vec3{1.0, 2.0, 2.0}), 9.0, 1e-12);
    NT_CHECK(norm(zero) == 0.0);  // exact == is fine: no arithmetic happened worth rounding
    // norm(v)^2 == norm2(v) for a non-trivial vector too.
    NT_CLOSE(norm(a) * norm(a), norm2(a), 1e-12);

    // --- normalized(): zero-guard and unit length ------------------------
    // Zero vector must return the zero vector, not NaN from 0/0.
    NT_CHECK(normalized(zero) == Vec3{});
    NT_CHECK(!std::isnan(normalized(zero).x));
    // Any nonzero vector normalizes to unit length...
    NT_CLOSE(norm(normalized(a)), 1.0, 1e-12);
    // ...and points the same direction (cross product with the original is zero).
    NT_CLOSE(norm2(cross(normalized(a), a)), 0.0, 1e-12);
    // An already-unit vector is (numerically) unchanged.
    const Vec3 unit_x{1.0, 0.0, 0.0};
    NT_CLOSE_VEC3(normalized(unit_x), unit_x, 1e-12);

    return newton_test::report("vec3");
}
