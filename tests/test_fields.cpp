// tests/test_fields.cpp — numerical gradient, Jacobian, curl, and divergence
// (newton::fields, Taylor §4.3-4.4): gradients against the analytic -F of the
// library's potentials, exact zeros for constant fields, the textbook curl and
// divergence examples, the vector-calculus identities, the spherical-polar
// projection, the classifier's relative_curl score, the known kink pitfall,
// and the atomic-scale failure that motivates `length_scale`.
//
// Every expected value comes from calculus (closed-form derivatives) or from
// hand calculation. Vectors are compared by NORM, not component by component:
// a per-component relative test falsely fails a tiny component of a huge
// vector (an Earth-orbit force has components 1e-3 of its length).
#include "test_util.hpp"        // quotes: sits next to this file
#include <newton/newton.hpp>    // angle brackets: reached via -I include
#include <newton/fields/derivatives.hpp>

#include <algorithm>
#include <cmath>

using newton::Vec3;
using newton::norm;
using newton::fields::Jacobian;
using newton::fields::curl;
using newton::fields::default_step;
using newton::fields::divergence;
using newton::fields::frobenius_norm;
using newton::fields::gradient;
using newton::fields::jacobian;

namespace {

// Relative error of vector a against reference b: |a - b| / max(|b|, tiny).
// `tiny` only stops a division by zero when b is the zero vector.
double rel_err(const Vec3& a, const Vec3& b) {
    return norm(a - b) / std::max(norm(b), 1e-300);
}

// Relative error of a Jacobian against a reference, in the Frobenius norm.
double rel_err(const Jacobian& a, const Jacobian& b) {
    const Jacobian d{a.dx - b.dx, a.dy - b.dy, a.dz - b.dz};
    return frobenius_norm(d) / std::max(frobenius_norm(b), 1e-300);
}

// ---- test potentials, each with its analytic gradient (so analytic F = -grad) ----
// Central potentials U(rho) have gradient U'(rho) * d / rho, d = r - center.

// Point-mass gravity: U = -GMm / rho.  grad U = +GMm d / rho^3.
struct GravityU {
    double GMm{1.0};
    double operator()(const Vec3& r) const { return -GMm / norm(r); }
    Vec3 grad(const Vec3& r) const { return (GMm / std::pow(norm(r), 3.0)) * r; }
};

// Morse: U = D (e^{-2x} - 2 e^{-x}),  x = (rho - r0)/S.  U'(rho) = (2D/S)(e^{-x} - e^{-2x}).
struct MorseU {
    double D{1.0}, r0{1.0}, S{0.5};
    double operator()(const Vec3& r) const {
        const double x = (norm(r) - r0) / S;
        return D * (std::exp(-2.0 * x) - 2.0 * std::exp(-x));
    }
    Vec3 grad(const Vec3& r) const {
        const double rho = norm(r);
        const double x = (rho - r0) / S;
        return ((2.0 * D / S) * (std::exp(-x) - std::exp(-2.0 * x)) / rho) * r;
    }
};

// Lennard-Jones: U = U0 [(r0/rho)^12 - 2 (r0/rho)^6].  U'(rho) = -(12 U0/rho)[(r0/rho)^12 - (r0/rho)^6].
struct LennardJonesU {
    double U0{1.0}, r0{1.0};
    double operator()(const Vec3& r) const {
        const double q = std::pow(r0 / norm(r), 6.0);
        return U0 * (q * q - 2.0 * q);
    }
    Vec3 grad(const Vec3& r) const {
        const double rho = norm(r);
        const double q = std::pow(r0 / rho, 6.0);
        return (-(12.0 * U0 / rho) * (q * q - q) / rho) * r;
    }
};

// Yukawa (screened Coulomb): U = -(k/rho) e^{-rho/a}.  U'(rho) = k e^{-rho/a} (1/rho^2 + 1/(a rho)).
struct YukawaU {
    double k{1.0}, a{1.0};
    double operator()(const Vec3& r) const { return -(k / norm(r)) * std::exp(-norm(r) / a); }
    Vec3 grad(const Vec3& r) const {
        const double rho = norm(r);
        const double dU = k * std::exp(-rho / a) * (1.0 / (rho * rho) + 1.0 / (a * rho));
        return (dU / rho) * r;
    }
};

// A non-central trig field: U = sin(1.3x) cos(0.7y) e^{0.4z}.
struct TrigU {
    double operator()(const Vec3& r) const {
        return std::sin(1.3 * r.x) * std::cos(0.7 * r.y) * std::exp(0.4 * r.z);
    }
    Vec3 grad(const Vec3& r) const {
        const double e = std::exp(0.4 * r.z);
        return {1.3 * std::cos(1.3 * r.x) * std::cos(0.7 * r.y) * e,
                -0.7 * std::sin(1.3 * r.x) * std::sin(0.7 * r.y) * e,
                0.4 * std::sin(1.3 * r.x) * std::cos(0.7 * r.y) * e};
    }
};

// Helpers that give a plain-function-pointer scalar field and a vector field,
// used below to show the templates accept any callable.
double plain_function_U(const Vec3& r) { return 2.0 * norm2(r); }  // U = 2 r^2, the k = 4 spring
Vec3 swirl(const Vec3& r) { return {-r.y, r.x, 0.0}; }

}  // namespace

int main() {
    using namespace newton::forces;
    using newton::norm2;

    // The relative tolerance for gradient/Jacobian vs analytic, every field in
    // the standard battery: measured worst case is ~1e-10, so 1e-9 has a
    // margin of 10x on the hardest field and ~300x on the easy ones.
    const double tol_fd = 1e-9;

    // ---------------------------------------------------------------------
    // default_step: the step rule h = 1e-4 * max(|r|, length_scale).
    // ---------------------------------------------------------------------
    NT_CLOSE(default_step(Vec3{0.0, 0.0, 0.0}), 1e-4, 1e-15);        // floored by 1 m at the origin
    NT_CLOSE(default_step(Vec3{3.0, 4.0, 0.0}), 5e-4, 1e-15);        // |r| = 5 m beats the floor
    NT_CLOSE(default_step(Vec3{7e6, 0.0, 0.0}), 700.0, 1e-12);       // orbital scale: h grows with |r|
    NT_CLOSE(default_step(Vec3{1e-10, 0.0, 0.0}, 1e-10), 1e-14, 1e-24);  // atomic scale, user floor

    // ---------------------------------------------------------------------
    // gradient vs analytic: the numerical -grad U must equal the library's F.
    // Here F is the library's own force functor, U is written independently.
    // ---------------------------------------------------------------------
    {
        // Spring k = 4: U = 2|r|^2, F = -4 r.  A polynomial: the 5-point stencil
        // is exact up to rounding, so a tighter tolerance (1e-11) applies.
        const Spring spring{4.0, {}, 0.0};
        const auto U = [](const Vec3& r) { return 2.0 * newton::norm2(r); };
        for (const Vec3 r : {Vec3{1.0, 2.0, 3.0}, Vec3{-0.3, 0.2, 0.5}, Vec3{2.0, -1.0, 0.0}}) {
            const Vec3 g = gradient(U, r);
            NT_CHECK(rel_err(-g, spring(0.0, r, Vec3{})) < 1e-11);
        }

        // Point gravity at lab scale: GM = 3, m = 2, so GMm = 6.  Includes a point
        // inside the 1 m floor (|r| = 0.5), where h = 1e-4 m, still much less than |r|.
        const PointGravity grav{3.0, 2.0, {}, 0.0};
        const GravityU Ug{6.0};
        for (const Vec3 r : {Vec3{1.0, 2.0, 2.0}, Vec3{0.3, -0.2, 0.3}, Vec3{-2.0, 0.5, 1.0}}) {
            const Vec3 g = gradient(Ug, r);
            NT_CHECK(rel_err(-g, grav(0.0, r, Vec3{})) < tol_fd);
            NT_CHECK(rel_err(g, Ug.grad(r)) < tol_fd);
        }

        // Earth-orbit scale: GM = 3.986e14 m^3/s^2, m = 1000 kg, |r| ~ 7e6 m.
        // h = 1e-4 * |r| = ~700 m; the force is ~8e3 N and has tiny components
        // (the y and z parts are small fractions of x), which is exactly why
        // the comparison is by norm.
        const PointGravity earth{3.986e14, 1000.0, {}, 0.0};
        const GravityU Ue{3.986e14 * 1000.0};
        for (const Vec3 r : {Vec3{7.0e6, 1.0e6, -2.0e6}, Vec3{6.6e6, 0.0, 1.0}, Vec3{-5.0e6, 4.5e6, 0.1e6}}) {
            const Vec3 g = gradient(Ue, r);
            NT_CHECK(rel_err(-g, earth(0.0, r, Vec3{})) < tol_fd);
        }

        // Yukawa, k = 1, a = 0.5, at |r| = 0.5..1.5.
        // Spec check from the physics notes: for a = 1 the radial force at r = 1 is
        // -2/e = -0.735759 N (negative = attractive).
        const YukawaU Uy{1.0, 0.5};
        for (const Vec3 r : {Vec3{0.5, 0.3, 0.1}, Vec3{1.0, 0.0, 0.0}, Vec3{-0.6, 0.9, 0.7}}) {
            NT_CHECK(rel_err(gradient(Uy, r), Uy.grad(r)) < tol_fd);
        }
        const YukawaU Uy1{1.0, 1.0};
        const Vec3 gy = gradient(Uy1, Vec3{1.0, 0.0, 0.0});
        NT_CLOSE(-gy.x, -2.0 / std::exp(1.0), 1e-9);  // F_r(1) = -0.735759 N
        NT_CLOSE(gy.y, 0.0, 1e-12);
        NT_CLOSE(gy.z, 0.0, 1e-12);

        // A trig field with no symmetry at all: the stencil must get every
        // component right, not just a radial one.
        const TrigU Ut;
        for (const Vec3 r : {Vec3{0.4, -0.2, 0.3}, Vec3{1.1, 0.8, -0.5}, Vec3{-0.7, 0.3, 0.9}}) {
            NT_CHECK(rel_err(gradient(Ut, r), Ut.grad(r)) < tol_fd);
        }
    }

    // ---------------------------------------------------------------------
    // Atomic scale: length_scale matters.
    // Morse: D = 1e-19 J, r0 = 1e-10 m, S = 3e-11 m.  Lennard-Jones:
    // U0 = 1.6e-21 J, r0 = 3.8e-10 m (argon-like).
    // ---------------------------------------------------------------------
    {
        const MorseU Um{1e-19, 1e-10, 3e-11};
        const Vec3 r{0.8e-10, 0.6e-10, 0.4e-10};  // rho = 1.077e-10 m: just outside r0
        const Vec3 g_exact = Um.grad(r);
        NT_CHECK(norm(g_exact) > 1e-10);          // sanity: a genuinely nonzero force, ~1e-9 N

        // With the atomic length scale the gradient is accurate.
        const Vec3 g_good = gradient(Um, r, default_step(r, 1e-10));
        NT_CHECK(rel_err(g_good, g_exact) < tol_fd);

        // With the default 1 m floor, h = 1e-4 m: all four stencil points are
        // ~1e-4 m away, where exp(-x) with x ~ 3e6 underflows to exactly 0, so U = 0 at every
        // point and the gradient is EXACTLY zero: a silent 100% error.
        // This is a regression test that documents the failure.
        const Vec3 g_bad = gradient(Um, r);
        NT_CHECK(g_bad == Vec3{});
        NT_CHECK(rel_err(g_bad, g_exact) == 1.0);

        // Lennard-Jones at atomic scale.  rho = 4e-10: slightly beyond r0.
        const LennardJonesU Ul{1.6e-21, 3.8e-10};
        const Vec3 rl{2.4e-10, -2.4e-10, 2.0e-10};  // rho = 4.0e-10 m
        NT_CHECK(rel_err(gradient(Ul, rl, default_step(rl, 1e-10)), Ul.grad(rl)) < tol_fd);
        // Lab-scale Lennard-Jones (everything of order 1): fine with the default step.
        const LennardJonesU U1{1.0, 1.0};
        NT_CHECK(rel_err(gradient(U1, Vec3{0.9, 0.5, 0.3}), U1.grad(Vec3{0.9, 0.5, 0.3})) < tol_fd);
        // Lab-scale Morse A = 1, r0 = 1, S = 0.5.
        const MorseU M1{1.0, 1.0, 0.5};
        NT_CHECK(rel_err(gradient(M1, Vec3{0.9, 0.5, 0.3}), M1.grad(Vec3{0.9, 0.5, 0.3})) < tol_fd);
    }

    // ---------------------------------------------------------------------
    // Overloads and callable kinds.
    // ---------------------------------------------------------------------
    {
        const Vec3 r{0.4, -0.2, 0.3};
        const TrigU Ut;
        // The no-h overload is exactly the explicit-h overload with default_step(r).
        NT_CHECK(gradient(Ut, r) == gradient(Ut, r, default_step(r)));
        // A function pointer, a const lambda, and a temporary lambda all work.
        NT_CLOSE_VEC3(gradient(plain_function_U, Vec3{1.0, 2.0, 3.0}), (Vec3{4.0, 8.0, 12.0}), 1e-10);
        const auto const_lambda = [](const Vec3& p) { return p.x * p.y * p.z; };
        NT_CLOSE_VEC3(gradient(const_lambda, Vec3{1.0, 2.0, 3.0}), (Vec3{6.0, 3.0, 2.0}), 1e-10);
        NT_CLOSE_VEC3(gradient([](const Vec3& p) { return p.x * p.x * p.x; }, Vec3{2.0, 0.0, 0.0}),
                      (Vec3{12.0, 0.0, 0.0}), 1e-9);
        // The Jacobian/curl/divergence no-h overloads agree with the explicit-h ones.
        const Vec3 q{0.7, -0.4, 0.2};
        NT_CHECK(curl(swirl, q) == curl(swirl, q, default_step(q)));
        NT_CHECK(divergence(swirl, q) == divergence(swirl, q, default_step(q)));
        NT_CHECK(rel_err(jacobian(swirl, q), jacobian(swirl, q, default_step(q))) == 0.0);
    }

    // ---------------------------------------------------------------------
    // Constant fields: every derivative is EXACTLY zero (not just small).
    // F(r + s) - F(r - s) is bit-for-bit 0, so ==, not a tolerance, is
    // legitimate here.
    // ---------------------------------------------------------------------
    {
        const ConstantForce cf{{1.0, 2.0, 3.0}};
        const UniformGravity ug{2.0, {0.0, -9.80665, 0.0}};
        const ConstantForce zero{};
        const auto field_cf = [&](const Vec3& r) { return cf(0.0, r, Vec3{}); };
        const auto field_ug = [&](const Vec3& r) { return ug(0.0, r, Vec3{}); };
        const auto field_zero = [&](const Vec3& r) { return zero(0.0, r, Vec3{}); };
        for (const Vec3 r : {Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 2.0, 3.0}, Vec3{-5e6, 1e-3, 40.0}}) {
            const Jacobian J1 = jacobian(field_cf, r);
            NT_CHECK(J1.dx == Vec3{} && J1.dy == Vec3{} && J1.dz == Vec3{});
            const Jacobian J2 = jacobian(field_ug, r);
            NT_CHECK(J2.dx == Vec3{} && J2.dy == Vec3{} && J2.dz == Vec3{});
            const Jacobian J3 = jacobian(field_zero, r);
            NT_CHECK(J3.dx == Vec3{} && J3.dy == Vec3{} && J3.dz == Vec3{});
            NT_CHECK(curl(field_cf, r) == Vec3{});
            NT_CHECK(curl(field_ug, r) == Vec3{});
            NT_CHECK(divergence(field_cf, r) == 0.0);
            NT_CHECK(divergence(field_ug, r) == 0.0);
        }
    }

    // ---------------------------------------------------------------------
    // Jacobian layout: columns.  F = (a x + b y, c x + d y + e z, f z + g x)
    // has dF/dx = (a, c, g), dF/dy = (b, d, 0), dF/dz = (0, e, f).
    // Linear, so the stencil is exact up to rounding.
    // ---------------------------------------------------------------------
    {
        const auto F = [](const Vec3& r) {
            return Vec3{1.0 * r.x + 2.0 * r.y, 3.0 * r.x + 4.0 * r.y + 5.0 * r.z, 6.0 * r.z + 7.0 * r.x};
        };
        const Jacobian J = jacobian(F, Vec3{0.3, -0.8, 1.2});
        NT_CLOSE_VEC3(J.dx, (Vec3{1.0, 3.0, 7.0}), 1e-10);
        NT_CLOSE_VEC3(J.dy, (Vec3{2.0, 4.0, 0.0}), 1e-10);
        NT_CLOSE_VEC3(J.dz, (Vec3{0.0, 5.0, 6.0}), 1e-10);
        NT_CLOSE(J.dx.y, 3.0, 1e-10);  // dFy/dx, as documented
        NT_CLOSE(frobenius_norm(J), std::sqrt(1.0 + 4.0 + 9.0 + 16.0 + 25.0 + 36.0 + 49.0), 1e-10);
        // The constexpr Jacobian-only overloads, on hand-written entries:
        // curl = (dy.z - dz.y, dz.x - dx.z, dx.y - dy.x).
        constexpr Jacobian K{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
        static_assert(curl(K) == Vec3{6.0 - 8.0, 7.0 - 3.0, 2.0 - 4.0}, "curl(Jacobian) is constexpr");
        static_assert(divergence(K) == 15.0, "divergence(Jacobian) is constexpr");
        NT_CHECK(curl(K) == (Vec3{-2.0, 4.0, -2.0}));
        NT_CHECK(divergence(K) == 15.0);
    }

    // ---------------------------------------------------------------------
    // Curl: the textbook examples (Taylor §4.4).  1e-10 absolute is
    // rounding only for these O(1) linear/quadratic fields.
    // ---------------------------------------------------------------------
    {
        const Vec3 r{0.7, -1.3, 0.4};

        // Swirl (-y, x, 0): curl = (0, 0, 2) everywhere.
        NT_CLOSE_VEC3(curl(swirl, r), (Vec3{0.0, 0.0, 2.0}), 1e-10);

        // Rigid rotation F = omega x r: curl = 2 omega.
        const Vec3 omega{0.3, -1.1, 0.8};
        const auto rigid = [&](const Vec3& p) { return cross(omega, p); };
        NT_CLOSE_VEC3(curl(rigid, r), 2.0 * omega, 1e-10);

        // (a y, b x, 0): curl = (0, 0, b - a); conservative iff a == b.
        {
            const double a = 3.0, b = 1.0;
            const auto F = [&](const Vec3& p) { return Vec3{a * p.y, b * p.x, 0.0}; };
            NT_CLOSE_VEC3(curl(F, r), (Vec3{0.0, 0.0, b - a}), 1e-10);
            const auto G = [&](const Vec3& p) { return Vec3{a * p.y, a * p.x, 0.0}; };
            NT_CLOSE_VEC3(curl(G, r), (Vec3{}), 1e-10);   // a == b: curl-free (G = grad(a x y))
        }

        // F = x r = (x^2, x y, x z): central but not spherically symmetric
        // (§4.8 counterexample). curl = (0, -z, y); at (1, 2, 3) that is (0, -3, 2).
        {
            const auto F = [](const Vec3& p) { return p.x * p; };
            NT_CLOSE_VEC3(curl(F, Vec3{1.0, 2.0, 3.0}), (Vec3{0.0, -3.0, 2.0}), 1e-9);
            NT_CLOSE_VEC3(curl(F, r), (Vec3{0.0, -r.z, r.y}), 1e-9);
            NT_CLOSE(divergence(F, Vec3{1.0, 2.0, 3.0}), 4.0, 1e-9);  // div (x r) = 4x
        }

        // (x y, y z, z x): curl = (-y, -z, -x), div = x + y + z.  At (1,2,3): (-2,-3,-1), 6.
        {
            const auto F = [](const Vec3& p) { return Vec3{p.x * p.y, p.y * p.z, p.z * p.x}; };
            NT_CLOSE_VEC3(curl(F, Vec3{1.0, 2.0, 3.0}), (Vec3{-2.0, -3.0, -1.0}), 1e-9);
            NT_CLOSE(divergence(F, Vec3{1.0, 2.0, 3.0}), 6.0, 1e-9);
        }

        // The classic textbook conservative field (2xy + z^3, x^2, 3xz^2): curl = 0.
        {
            const auto F = [](const Vec3& p) {
                return Vec3{2.0 * p.x * p.y + p.z * p.z * p.z, p.x * p.x, 3.0 * p.x * p.z * p.z};
            };
            NT_CLOSE_VEC3(curl(F, Vec3{1.0, 1.0, 1.0}), (Vec3{}), 1e-9);
            NT_CLOSE(divergence(F, Vec3{1.0, 1.0, 1.0}), 8.0, 1e-9);  // 2y + 6xz = 8
        }

        // Conservative library forces, curl relative to the size of the Jacobian.
        // The plain |curl| / ||J|| is ~1e-11 for these; tolerance 1e-9.
        {
            const PointGravity grav{3.0, 2.0, {}, 0.0};
            const Spring sp{4.0, {1.0, 0.0, 0.0}, 0.5};        // anchored, natural length 0.5
            const AnisotropicSpring asp{1.0, 2.0, 3.0};
            const auto fg = [&](const Vec3& p) { return grav(0.0, p, Vec3{}); };
            const auto fs = [&](const Vec3& p) { return sp(0.0, p, Vec3{}); };
            const auto fa = [&](const Vec3& p) { return asp(0.0, p, Vec3{}); };
            for (const Vec3 p : {Vec3{1.0, 2.0, 2.0}, Vec3{0.4, -0.3, 0.8}, Vec3{3.0, 1.0, 1.0}}) {
                NT_CHECK(norm(curl(fg, p)) / frobenius_norm(jacobian(fg, p)) < 1e-9);
                NT_CHECK(norm(curl(fs, p)) / frobenius_norm(jacobian(fs, p)) < 1e-9);
                NT_CHECK(norm(curl(fa, p)) / frobenius_norm(jacobian(fa, p)) < 1e-9);
            }
        }
    }

    // ---------------------------------------------------------------------
    // Divergence.
    // ---------------------------------------------------------------------
    {
        // Spring: div(-k r) = -3k.  k = 4: -12.
        const Spring sp{4.0, {}, 0.0};
        const auto fsp = [&](const Vec3& p) { return sp(0.0, p, Vec3{}); };
        NT_CLOSE(divergence(fsp, Vec3{1.0, -2.0, 0.5}), -12.0, 1e-10);
        NT_CLOSE(divergence(fsp, Vec3{0.0, 0.0, 0.0}), -12.0, 1e-10);

        // Anisotropic spring: -(kx + ky + kz) = -6.
        const AnisotropicSpring asp{1.0, 2.0, 3.0};
        const auto fas = [&](const Vec3& p) { return asp(0.0, p, Vec3{}); };
        NT_CLOSE(divergence(fas, Vec3{0.3, 0.2, 0.1}), -6.0, 1e-10);

        // Position field F = r: div r = 3.
        const auto id = [](const Vec3& p) { return p; };
        NT_CLOSE(divergence(id, Vec3{1.0, 2.0, 3.0}), 3.0, 1e-10);

        // Inverse-square field: div = 0 away from the center. Compare to the
        // size of the Jacobian, since the true answer is 0.
        const PointGravity grav{3.0, 2.0, {}, 0.0};
        const auto fg = [&](const Vec3& p) { return grav(0.0, p, Vec3{}); };
        for (const Vec3 p : {Vec3{1.0, 2.0, 2.0}, Vec3{0.3, -0.2, 0.3}, Vec3{-2.0, 0.5, 1.0}}) {
            NT_CHECK(std::fabs(divergence(fg, p)) / frobenius_norm(jacobian(fg, p)) < 1e-9);
        }

        // Softened gravity: div = -3 GMm eps^2 / (r^2 + eps^2)^(5/2).  At the
        // center this is -3 * 6 * 0.25 / 0.5^5 = -144 N/m  (eps = 0.5).
        // The center is the worst case: the field is smooth there, so the stencil works.
        const PointGravity soft{3.0, 2.0, {}, 0.5};
        const auto fsoft = [&](const Vec3& p) { return soft(0.0, p, Vec3{}); };
        NT_CLOSE(divergence(fsoft, Vec3{}), -144.0, 1e-8);

        // Swirl: div = 0 even though curl != 0 (divergence is not the
        // conservative-force test).
        NT_CLOSE(divergence(swirl, Vec3{0.3, 0.9, 0.2}), 0.0, 1e-10);
    }

    // ---------------------------------------------------------------------
    // Identities.  Both derivatives are numerical (a finite difference inside
    // a finite difference), so the inner error ~1e-12 is divided by h again
    // and the achievable agreement is ~1e-8: 1e-6 relative to the second
    // derivatives is the honest tolerance.  (Measured: 2.3e-8 and 6.2e-9.)
    // ---------------------------------------------------------------------
    {
        // curl(grad U) = 0 for any smooth U.
        const TrigU Ut;
        const YukawaU Uy{1.0, 0.5};
        const GravityU Ug{6.0};
        const auto g_trig = [&](const Vec3& p) { return gradient(Ut, p); };
        const auto g_yuk = [&](const Vec3& p) { return gradient(Uy, p); };
        const auto g_grav = [&](const Vec3& p) { return gradient(Ug, p); };
        for (const Vec3 p : {Vec3{0.4, -0.2, 0.3}, Vec3{1.1, 0.8, -0.5}}) {
            NT_CHECK(norm(curl(g_trig, p)) / frobenius_norm(jacobian(g_trig, p)) < 1e-6);
            NT_CHECK(norm(curl(g_yuk, p)) / frobenius_norm(jacobian(g_yuk, p)) < 1e-6);
            NT_CHECK(norm(curl(g_grav, p)) / frobenius_norm(jacobian(g_grav, p)) < 1e-6);
        }

        // div(curl A) = 0 for any smooth A.  A = (y z^2, x^3 + sin z, e^{x y}).
        const auto A = [](const Vec3& p) {
            return Vec3{p.y * p.z * p.z, p.x * p.x * p.x + std::sin(p.z), std::exp(p.x * p.y)};
        };
        // Analytic curl A = (x e^{xy} - cos z, 2 y z - y e^{xy}, 3 x^2 - z^2).
        const auto curlA_exact = [](const Vec3& p) {
            const double e = std::exp(p.x * p.y);
            return Vec3{p.x * e - std::cos(p.z), 2.0 * p.y * p.z - p.y * e, 3.0 * p.x * p.x - p.z * p.z};
        };
        const auto curlA_num = [&](const Vec3& p) { return curl(A, p); };
        for (const Vec3 p : {Vec3{0.5, 0.6, 0.7}, Vec3{-0.8, 0.4, 1.1}}) {
            NT_CHECK(rel_err(curl(A, p), curlA_exact(p)) < 1e-9);                      // curl A vs closed form
            NT_CHECK(std::fabs(divergence(A, p)) < 1e-10);                             // div A = 0 here (exact)
            NT_CHECK(std::fabs(divergence(curlA_num, p)) / frobenius_norm(jacobian(curlA_num, p)) < 1e-6);
        }
    }

    // ---------------------------------------------------------------------
    // Spherical projection (Taylor §4.8): the gradient of a function of r
    // alone is purely radial.  U = 1/rho: -grad U = r_hat / rho^2, so the
    // radial component is 1/rho^2 and the theta and phi components are 0.
    // The unit vectors are built by hand from the polar angles.
    // ---------------------------------------------------------------------
    {
        const auto U = [](const Vec3& p) { return 1.0 / norm(p); };
        for (const Vec3 p : {Vec3{1.0, 2.0, 2.0}, Vec3{-0.6, 0.3, 0.7}}) {
            const newton::coords::SphericalPolar s = newton::coords::to_spherical(p);
            const Vec3 rhat{std::sin(s.theta) * std::cos(s.phi), std::sin(s.theta) * std::sin(s.phi), std::cos(s.theta)};
            const Vec3 that{std::cos(s.theta) * std::cos(s.phi), std::cos(s.theta) * std::sin(s.phi), -std::sin(s.theta)};
            const Vec3 phat{-std::sin(s.phi), std::cos(s.phi), 0.0};
            const Vec3 F = -gradient(U, p);
            NT_CLOSE(dot(F, rhat), 1.0 / (s.r * s.r), 1e-9);           // radial: 1/rho^2
            NT_CLOSE(dot(F, that) * s.r * s.r, 0.0, 1e-9);             // tangential parts vanish
            NT_CLOSE(dot(F, phat) * s.r * s.r, 0.0, 1e-9);             // (scaled by rho^2: O(1) test)
        }
    }

    // ---------------------------------------------------------------------
    // detail::relative_curl: the classifier's unit-free curl score.
    // Sample at r, with ell = max(|r|, L) and h = 1e-4 ell, L = 1 m.
    // ---------------------------------------------------------------------
    {
        using newton::fields::detail::jacobian_with_error;
        using newton::fields::detail::relative_curl;
        const auto score = [](auto&& field, const Vec3& r) {
            const double ell = std::max(norm(r), 1.0);
            const auto s = jacobian_with_error(field, r, 1e-4 * ell);
            NT_CHECK(s.finite);
            return relative_curl(s, ell);
        };

        // Conservative fields score <= 1e-9 (measured: ~1e-14; the threshold the classifier uses is 1e-6).
        {
            const PointGravity grav{3.0, 2.0, {}, 0.0};
            const Spring sp{4.0, {1.0, 0.0, 0.0}, 0.5};
            const auto fg = [&](const Vec3& p) { return grav(0.0, p, Vec3{}); };
            const auto fs = [&](const Vec3& p) { return sp(0.0, p, Vec3{}); };
            const auto fc = [](const Vec3&) { return Vec3{1.0, 2.0, 3.0}; };
            const auto fu = [](const Vec3&) { return Vec3{}; };  // F = 0: the 0/0 edge case
            const YukawaU Uy{1.0, 0.5};
            const auto fy = [&](const Vec3& p) { return -Uy.grad(p); };
            for (const Vec3 p : {Vec3{1.0, 2.0, 2.0}, Vec3{0.4, -0.3, 0.8}, Vec3{-0.5, 0.9, 0.2}}) {
                NT_CHECK(score(fg, p) <= 1e-9);
                NT_CHECK(score(fs, p) <= 1e-9);
                NT_CHECK(score(fc, p) == 0.0);
                NT_CHECK(score(fu, p) == 0.0);
                NT_CHECK(score(fy, p) <= 1e-9);
            }
        }

        // Spring k = 4 plus a weak swirl delta (-y, x, 0), delta = 4e-4, so delta/k = 1e-4.
        // curl = (0, 0, 2 delta); ||J||_F ~ sqrt(3 k^2 + 2 delta^2), so the score is
        // ~ 2 delta / (sqrt(3) k) = 1.155 delta/k = 1.155e-4.  Allow x0.5 .. x2.
        {
            const double k = 4.0, delta = 4e-4;
            const auto F = [&](const Vec3& p) { return Vec3{-k * p.x - delta * p.y, -k * p.y + delta * p.x, -k * p.z}; };
            for (const Vec3 p : {Vec3{0.5, 0.2, -0.3}, Vec3{-0.8, 0.6, 0.1}}) {
                const double sc = score(F, p);
                NT_CHECK(sc > 0.5 * 1.155e-4 && sc < 2.0 * 1.155e-4);
                NT_CHECK(sc > 1e-6);  // far above the classifier's curl_tol, so it is flagged
            }
        }

        // Fields whose Jacobian is (almost) pure noise.  Plain |curl|/||J|| is meaningless
        // there (0/0, or noise/noise = O(1)); the F_max/ell term in relative_curl's
        // denominator is what keeps the score honest.
        {
            // (1) Hidden catastrophic cancellation, cross-coupled (research-differentiation §2, b6):
            // s = ((x+y+z+C) - C) - (x+y+z) is exactly 0 in real arithmetic but rounding noise in
            // doubles, and it feeds all three components, so the noise Jacobian has a nonzero curl.
            const double C = 1e4;
            const auto hidden = [&](const Vec3& p) {
                const double q = p.x + p.y + p.z;
                const double sn = ((q + C) - C) - q;
                return Vec3{1.0 + sn, 2.0 + sn, 3.0 + sn};
            };
            // p0 was found by scanning random points: there the noise happens to give a plain
            // |curl|/||J|| of 0.59, pure noise over noise, which would flag a conservative field.
            const Vec3 p0{0.16210661573163976, 0.99102391635793152, -0.1160165633639243};
            for (const Vec3 p : {p0, Vec3{0.31, -0.42, 0.17}, Vec3{0.1, 0.2, 0.3}}) {
                NT_CHECK(score(hidden, p) <= 1e-6);  // the classifier's curl_tol
            }
            // The plain score at p0 (0.59 on the machine this was written on) is deliberately
            // NOT asserted: it depends on that machine's exact rounding pattern, and another
            // compiler or CPU (e.g. one that fuses multiply-adds) scatters the noise differently,
            // so asserting it would test the hardware, not the library. What must hold on every
            // machine is the check above: relative_curl stays below curl_tol however the noise falls.
        }
        {
            // (2) A strong uniform field plus a weak swirl: F = (0, -19.6, 0) + delta (-y, x, 0),
            // delta = 1e-5 N/m.  The Jacobian is ONLY the swirl, so plain |curl|/||J|| = sqrt(2) = 1.41
            // however weak the swirl is.  Relative to the field's own scale |F|/ell = 19.6 the swirl
            // is ~1e-6: that is the score relative_curl must give (research-differentiation §2).
            const double delta = 1e-5;
            const auto F = [&](const Vec3& p) { return Vec3{-delta * p.y, -19.6 + delta * p.x, 0.0}; };
            for (const Vec3 p : {Vec3{0.5, 0.2, -0.3}, Vec3{-0.8, 0.6, 0.1}}) {
                const double sc = score(F, p);
                NT_CHECK(sc > 0.5e-6 && sc < 2e-6);
                NT_CHECK(norm(curl(F, p)) / frobenius_norm(jacobian(F, p)) > 1.0);  // the plain score
            }
        }

        // 1/r^2 field sampled 3 mm from its center, with L = 1 m so h/rho = 0.033: the plain
        // |curl|/||J|| is ~7e-6 there (a false alarm), but the stencil's own 4th-difference
        // error estimate is subtracted, and the score is exactly 0.
        {
            const PointGravity grav{3.0, 2.0, {}, 0.0};
            const auto fg = [&](const Vec3& p) { return grav(0.0, p, Vec3{}); };
            NT_CHECK(score(fg, Vec3{3e-3, 0.0, 0.0}) == 0.0);
            // Off the symmetry axis (|q| = 3.1 mm) the plain score is a genuine false alarm.
            const Vec3 q{1.2e-3, -1.6e-3, 2.4e-3};
            NT_CHECK(score(fg, q) == 0.0);
            // The same point with the plain criterion, to show why the error term exists.
            NT_CHECK(norm(curl(fg, q)) / frobenius_norm(jacobian(fg, q)) > 1e-6);
        }

        // Non-finite samples are flagged: PointGravity with a center exactly on a stencil point.
        {
            const auto singular = [](const Vec3& p) { return Vec3{1.0 / p.x, 0.0, 0.0}; };
            const Vec3 p{1e-4, 0.0, 0.0};  // h = 1e-4: the stencil hits x = 0 exactly
            const auto s = jacobian_with_error(singular, p, 1e-4);
            NT_CHECK(!s.finite);
        }
    }

    // ---------------------------------------------------------------------
    // Documented pitfall: a kink.  U = |x| has slope exactly 1 for x > 0, but
    // within 2h of x = 0 the 5-point stencil is wrong.  With the point exactly
    // one step from the kink (x = h) the stencil reads f(-h) = h, f(0) = 0,
    // f(2h) = 2h, f(3h) = 3h, so
    //   d1 = f(2h) - f(0) = 2h,  d2 = f(3h) - f(-h) = 2h,
    //   (8 d1 - d2) / (12 h) = 14/12 = 7/6:  a 17% overshoot of the true slope 1.
    // Three steps away it is exact again.
    // ---------------------------------------------------------------------
    {
        const auto absx = [](const Vec3& p) { return std::fabs(p.x); };
        const double h = 0.01;
        NT_CLOSE(gradient(absx, Vec3{h, 0.0, 0.0}, h).x, 7.0 / 6.0, 1e-9);   // inside 2h: wrong by 7/6
        NT_CLOSE(gradient(absx, Vec3{3.0 * h, 0.0, 0.0}, h).x, 1.0, 1e-9);  // outside 2h: exact
        NT_CLOSE(gradient(absx, Vec3{-3.0 * h, 0.0, 0.0}, h).x, -1.0, 1e-9);
    }

    // ---- line integrals (added in the next step) ----

    return newton_test::report("fields");
}
