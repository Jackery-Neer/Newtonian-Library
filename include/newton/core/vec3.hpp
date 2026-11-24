#pragma once
// newton/core/vec3.hpp — the 3D vector type everything else is built on.
//
// Physics: positions, velocities, accelerations, and forces are all vectors
// in R^3 (Taylor §1.2). 2D problems (projectiles, oscillators) simply use
// z = 0 — one type, no special cases.
//
// C++ note: Vec3 is an *aggregate* — no constructors, all-public members with
// default member initializers. That keeps brace-init working (Vec3{1, 2, 3}),
// guarantees zero-initialization by default (no NaN surprises from
// uninitialized doubles), and lets structured bindings decompose it.
// Everything here is `constexpr` where the math allows (std::sqrt is not
// constexpr in C++17, so norm() and friends are plain inline).

#include <cmath>

namespace newton {

struct Vec3 {
    double x{0.0};
    double y{0.0};
    double z{0.0};

    // Compound assignment as members; the symmetric binary operators are free
    // functions below so that `2.0 * v` and `v * 2.0` both work.
    constexpr Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    constexpr Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    constexpr Vec3& operator*=(double s)      { x *= s;   y *= s;   z *= s;   return *this; }
    constexpr Vec3& operator/=(double s)      { x /= s;   y /= s;   z /= s;   return *this; }
};

// --- arithmetic ------------------------------------------------------------

[[nodiscard]] constexpr Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
[[nodiscard]] constexpr Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
[[nodiscard]] constexpr Vec3 operator*(Vec3 v, double s)      { return v *= s; }
[[nodiscard]] constexpr Vec3 operator*(double s, Vec3 v)      { return v *= s; }
[[nodiscard]] constexpr Vec3 operator/(Vec3 v, double s)      { return v /= s; }
[[nodiscard]] constexpr Vec3 operator-(const Vec3& v)         { return {-v.x, -v.y, -v.z}; }

// Exact equality. Legitimate for values that were *assigned*; for anything
// *computed*, compare with a tolerance instead (see tests/test_util.hpp).
[[nodiscard]] constexpr bool operator==(const Vec3& a, const Vec3& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
[[nodiscard]] constexpr bool operator!=(const Vec3& a, const Vec3& b) { return !(a == b); }

// --- products and norms ----------------------------------------------------

/// Dot product: a·b = |a||b|cosθ (Taylor §1.2).
[[nodiscard]] constexpr double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/// Cross product: a×b, right-hand rule. Used for angular momentum ℓ = r×p
/// (Taylor §3.4) and the magnetic force qv×B (Taylor §2.5).
[[nodiscard]] constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

/// Squared magnitude |v|². Prefer this over norm() when comparing lengths —
/// it avoids the sqrt entirely.
[[nodiscard]] constexpr double norm2(const Vec3& v) { return dot(v, v); }

/// Magnitude |v|.
[[nodiscard]] inline double norm(const Vec3& v) { return std::sqrt(norm2(v)); }

/// Unit vector v/|v|. Guarded: returns the zero vector when |v| is (near)
/// zero rather than dividing by zero and producing NaNs. Physics code that
/// needs v̂ (e.g. quadratic drag) should prefer the singularity-free algebraic
/// form where one exists: -c|v|·v instead of -c|v|²·v̂.
[[nodiscard]] inline Vec3 normalized(const Vec3& v) {
    const double n2 = norm2(v);
    if (n2 <= 0.0) return {};
    return v / std::sqrt(n2);
}

}  // namespace newton
