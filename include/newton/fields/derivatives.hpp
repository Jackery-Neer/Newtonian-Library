#pragma once
// newton/fields/derivatives.hpp — numerical gradient, Jacobian, curl, and
// divergence of fields given as plain callables.
//
// Physics: Taylor §4.3-4.4 turns the two conditions for a conservative force
// into statements about derivatives of fields: F = -grad U (a force from a
// potential energy) and curl F = 0 (the local test for path independence).
// This file computes those derivatives numerically, for ANY scalar field
// U(r) or vector field F(r) you can write as a function — including fields
// with no tidy closed-form derivative, which is exactly when you need it.
//
// A *scalar field* is any callable   double U(const Vec3& r)    [J, say].
// A *vector field* is any callable   Vec3   F(const Vec3& r)    [N, say].
// Fields here depend on position only. A library Force F(t, r, v) becomes a
// field by freezing t and v inside a lambda:
//     [&](const Vec3& r) { return force(t, r, Vec3{}); }
//
// Numerical method: 5-point central differences (4th order), described in
// detail at detail::central_difference. The default step is
// h = 1e-4 * max(|r|, length_scale); see default_step() for why.
//
// C++ note: everything here is a function template over the callable's type,
// so a lambda is called directly (and can be inlined) with no std::function
// overhead. None of these functions is `noexcept`: they call user code, and
// user code may throw.

#include <algorithm>
#include <cmath>
#include <type_traits>

#include <newton/core/vec3.hpp>

namespace newton::fields {

// --- step size ---------------------------------------------------------------

/// Default finite-difference step h = 1e-4 * max(|r|, length_scale)   [m].
///
/// Physics: a derivative is a limit h -> 0, but a computer cannot take it
/// literally. The step must be small compared with the distance over which
/// the field changes appreciably (its "length scale"), or the difference
/// quotient averages over real structure and is wrong.
///
/// C++ note (the step-size trade-off): two errors compete. Truncation error
/// (the stencil's Taylor series is cut off) shrinks like h^4. Rounding error
/// (subtracting two nearly equal doubles loses digits) grows like eps/h, with
/// eps = 2.2e-16. Their sum is smallest around h ~ 1e-3..1e-4 times the
/// field's length scale; this library picks the smaller, safer 1e-4 because
/// real fields often vary faster than |r| (a Yukawa range shorter than |r|,
/// gravity at |r| < 1 m). Measured worst-case relative error for the
/// library's standard fields is about 3e-12.
///
/// Why `length_scale` exists: the step grows with |r| so it works from lab
/// scale to orbital scale, but near the origin |r| -> 0 would give h -> 0.
/// So |r| is floored by `length_scale`, 1 m by default (everyday mechanics).
/// At atomic scale that floor is a silent disaster. With a Morse potential of
/// width 3e-11 m, a 1 m floor gives h = 1e-4 m, far past where the potential
/// has decayed to exactly 0.0 (exp underflows), so U is 0 at all four stencil
/// points and the gradient comes out as EXACTLY zero: no error, no warning,
/// a molecule that feels no force. Pass the smallest distance over which your
/// field changes appreciably:
///     gradient(U, r, default_step(r, 1e-10));
///
/// @param r            evaluation point [m]
/// @param length_scale smallest distance over which the field changes
///                     appreciably [m]; default 1 m
[[nodiscard]] inline double default_step(const Vec3& r, double length_scale = 1.0) {
    return 1e-4 * std::max(norm(r), length_scale);
}

// --- Jacobian and the things you read off it ----------------------------------

/// The nine first partial derivatives of a vector field F, stored as COLUMNS:
///   dx = dF/dx = (dFx/dx, dFy/dx, dFz/dx),  likewise dy and dz.
/// So `J.dx.y` is dFy/dx, and `J.dy.x` is dFx/dy.
///
/// Physics: the matrix J_ij = dF_i/dx_j describes how F changes when you move
/// a little in each direction. Its antisymmetric part is the curl (rotation),
/// its trace is the divergence (expansion). For F = -grad U it is minus the
/// Hessian of U, which is symmetric (mixed partials commute): that symmetry
/// IS curl F = 0.
///
/// C++ note: columns, not rows, because the natural way to compute it is one
/// partial derivative at a time: "move along x, see how the whole vector F
/// responds" gives one Vec3, the column dF/dx. The textbook matrix is the
/// transpose of this layout (rows are the components of F); the members are
/// named for the variable you differentiate by, so you never need the
/// transpose in your head. Aggregate with defaults on every member, like Vec3.
struct Jacobian {
    Vec3 dx{};  // dF/dx [F-units per m]
    Vec3 dy{};  // dF/dy
    Vec3 dz{};  // dF/dz
};

/// Curl from a Jacobian: curl F = (dFz/dy - dFy/dz, dFx/dz - dFz/dx,
/// dFy/dx - dFx/dy)  [F-units per m].
///
/// Physics: the curl measures how much F circulates around a point: its
/// component along a unit vector n is the work per unit area, (closed-loop
/// integral of F.dr) / A, around a tiny loop with normal n (Stokes'
/// theorem). Taylor §4.4 shows that a force F(r) has path-independent work,
/// condition (ii) of §4.2, exactly when curl F = 0 throughout a simply
/// connected region, one where every loop can be shrunk to a point (removing
/// a point keeps 3D space simply connected; removing a line does not). For
/// a time-dependent F(r, t), take the curl with t frozen: curl F = 0 at every
/// instant is what a time-dependent potential U(r, t) needs (§4.5).
///
/// C++ note: this can be `constexpr` because it only does arithmetic on
/// Vec3's members, with no sqrt and no calls to user code. It is
/// distinguished from the curl(F, r) overloads below by arity.
[[nodiscard]] constexpr Vec3 curl(const Jacobian& J) {
    return {J.dy.z - J.dz.y, J.dz.x - J.dx.z, J.dx.y - J.dy.x};
}

/// Divergence from a Jacobian: div F = dFx/dx + dFy/dy + dFz/dz
/// [F-units per m].
///
/// Physics: the divergence measures how much F spreads out from a point: the
/// net outward flux (closed-surface integral of F.dA) per unit volume of a
/// tiny closed surface (the divergence theorem). It is NOT one of Taylor's
/// tests for a conservative force. F = -k r is conservative and has
/// div F = -3k, and the swirl (-y, x, 0) has div F = 0 yet is not
/// conservative. What divergence does is connect a field to its sources: for
/// gravity div g = -4 pi G rho, so PointGravity has div F = 0 everywhere
/// except at its center, while the softened version has
/// div F = -3 GMm eps^2 / (r^2 + eps^2)^(5/2) where its smoothed mass lives.
///
/// C++ note: `constexpr` for the same reason as curl(const Jacobian&).
[[nodiscard]] constexpr double divergence(const Jacobian& J) {
    return J.dx.x + J.dy.y + J.dz.z;
}

/// Frobenius norm of the Jacobian: sqrt of the sum of squares of all nine
/// entries, a single number for "how big are the derivatives here"
/// [F-units per m]. Used to judge whether a curl is large or small *relative
/// to* the field's own variation.
[[nodiscard]] inline double frobenius_norm(const Jacobian& J) {
    return std::sqrt(norm2(J.dx) + norm2(J.dy) + norm2(J.dz));
}

// --- the stencil (implementation details) ----------------------------------------

namespace detail {

/// r moved by s along one axis (0 = x, 1 = y, 2 = z); the other components are
/// untouched.
///
/// C++ note: takes r by value, so the caller's Vec3 is not modified: we
/// change our private copy and return it.
[[nodiscard]] constexpr Vec3 shifted(Vec3 r, int axis, double s) {
    if (axis == 0) {
        r.x += s;
    } else if (axis == 1) {
        r.y += s;
    } else {
        r.z += s;
    }
    return r;
}

/// 5-point (4th-order) central difference along one axis:
///     f' ~ [8 (f(x+h) - f(x-h)) - (f(x+2h) - f(x-2h))] / (12 h).
///
/// Physics/numerics: this is *Richardson extrapolation* of the ordinary
/// 2-point difference D(h) = [f(x+h) - f(x-h)] / (2h). D(h) = f' + c h^2 + ...;
/// the same formula with step 2h has four times the h^2 error, so
/// D(h) + [D(h) - D(2h)] / 3 cancels the h^2 term and leaves an error
/// of order h^4. Expanding that combination gives exactly the stencil above.
/// Cost: 4 evaluations of f per derivative.
///
/// C++ note: templated on the *value type* of f. `auto` return and the
/// operators on Vec3 mean one function serves double (a gradient component)
/// and Vec3 (a Jacobian column). `Field&` (not `Field&&`) because we only
/// call it; it binds to const callables too, by deduction (Field = const T).
template <class Field>
[[nodiscard]] auto central_difference(Field& f, const Vec3& r, int axis, double h) {
    const auto d1 = f(shifted(r, axis, h)) - f(shifted(r, axis, -h));
    const auto d2 = f(shifted(r, axis, 2.0 * h)) - f(shifted(r, axis, -2.0 * h));
    return (8.0 * d1 - d2) / (12.0 * h);
}

/// True if all three components are finite (not NaN, not +-infinity).
[[nodiscard]] inline bool is_finite(const Vec3& v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

/// What the conservative-force classifier needs at one sample point: the
/// Jacobian, an error estimate for it, and the field's size.
struct JacobianSample {
    Jacobian J{};          // 5-point Jacobian [F-units per m]
    double error{0.0};     // 4th-difference estimate of the numerical noise in |curl| [F-units per m]
    double F_max{0.0};     // largest |F| among the 13 samples [F-units]
    bool finite{true};     // false if any of the 13 F values was NaN or infinite
};

/// The 5-point Jacobian (12 evaluations of F) plus F(r) itself (1 more) and a
/// numerical-noise estimate from 4th differences.
///
/// Numerics: along each axis the 4th difference
///     F(x+2h) - 4F(x+h) + 6F(x) - 4F(x-h) + F(x-2h)
/// is exactly zero for any cubic polynomial, so whatever it measures is NOT
/// the smooth field: it is rounding noise, or a singularity the step cannot
/// resolve. A curl smaller than that (divided by h, to put it in the same
/// units as a derivative) could have been manufactured by the method itself,
/// so the classifier subtracts it. Being an even-order difference, it is
/// independent of the odd-order derivative we are estimating.
///
/// @param F field Vec3(const Vec3&)
/// @param r sample point [m]
/// @param h step [m]
template <class VectorField>
[[nodiscard]] JacobianSample jacobian_with_error(VectorField& F, const Vec3& r, double h) {
    const Vec3 F0 = F(r);
    JacobianSample s;
    s.F_max = norm(F0);
    s.finite = is_finite(F0);
    Vec3* const column[3] = {&s.J.dx, &s.J.dy, &s.J.dz};
    for (int i = 0; i < 3; ++i) {
        const Vec3 p1 = F(shifted(r, i, h));
        const Vec3 m1 = F(shifted(r, i, -h));
        const Vec3 p2 = F(shifted(r, i, 2.0 * h));
        const Vec3 m2 = F(shifted(r, i, -2.0 * h));
        s.finite = s.finite && is_finite(p1) && is_finite(m1) && is_finite(p2) && is_finite(m2);
        *column[i] = (8.0 * (p1 - m1) - (p2 - m2)) / (12.0 * h);
        s.error += norm(p2 - 4.0 * p1 + 6.0 * F0 - 4.0 * m1 + m2) / h;
        s.F_max = std::max({s.F_max, norm(p1), norm(m1), norm(p2), norm(m2)});
    }
    return s;
}

/// A dimensionless, unit-free "how much does F really circulate here?":
///     max(0, |curl J| - error) / (frobenius_norm(J) + F_max / ell),
/// and 0 when the numerator is 0.
///
/// Physics: |curl F| alone has the units of F per meter, so "is it small?"
/// has no answer until you say small compared with what. The denominator is
/// the field's own scale: its derivatives, plus |F|/ell, the variation you
/// would expect from a field of size |F| changing over the length ell. Every
/// term is linear in F, so the score does not depend on the units or the
/// overall size of F. Conservative fields score ~1e-12 or less; a swirl of
/// relative strength delta scores about delta. The classifier flags scores
/// above 1e-6.
///
/// @param s   sample from jacobian_with_error (should have s.finite == true)
/// @param ell the length over which F is expected to vary, ell = max(|r|, L) [m]
[[nodiscard]] inline double relative_curl(const JacobianSample& s, double ell) {
    const double excess = norm(curl(s.J)) - s.error;
    if (!(excess > 0.0)) return 0.0;  // also catches NaN: "no evidence of a curl"
    return excess / (frobenius_norm(s.J) + s.F_max / ell);
}

}  // namespace detail

// --- public numerical derivatives ----------------------------------------------------

/// Numerical gradient of a scalar field, grad U = (dU/dx, dU/dy, dU/dz), with
/// an explicit step.
///
/// Physics: the gradient ∇U = (∂U/∂x, ∂U/∂y, ∂U/∂z) points in the direction in
/// which the scalar field U increases fastest, and its length is that rate of
/// increase per meter (Taylor §4.3). In mechanics it turns a potential energy
/// into a force, F = -∇U: a conservative force points "downhill",
/// perpendicular to the surfaces of constant U. In spherical coordinates
/// ∇U = r̂ ∂U/∂r + θ̂ (1/r) ∂U/∂θ + φ̂ (1/(r sin θ)) ∂U/∂φ, so a U that depends
/// only on r gives a purely radial force (§4.8). To see it, project this
/// Cartesian result onto r̂, θ̂, φ̂, which you build at the point from
/// coords::to_spherical's θ and φ (the library has no unit-vector helper yet):
/// r̂ = (sinθ cosφ, sinθ sinφ, cosθ), θ̂ = (cosθ cosφ, cosθ sinφ, −sinθ),
/// φ̂ = (−sinφ, cosφ, 0). tests/test_fields.cpp does exactly this.
///
/// Pitfalls (finite differences cannot see them coming):
///  - Near a kink, e.g. U = |x| at x = 0, the 5-point stencil is wrong within
///    2h of the kink. The error is largest one step away (x = h), where the
///    stencil overshoots the true slope by a factor 7/6 (about 17%); closer
///    than h it undershoots (the result is 7x/(6h) for 0 <= x <= h, so it
///    returns 0 at the kink itself).
///  - Near a singularity, e.g. U = -GMm/r at r = 0, accuracy degrades
///    gradually: the relative error falls off like about 4 (h/rho)^4, where
///    rho is the distance to the singularity (measured for U = -1/rho:
///    1.7e-4 at rho = 10h, 2e-6 at 30h, 1.7e-8 at 100h). Even 100h away that
///    is far above the ~3e-12 quoted for smooth fields. Inside about 2h the
///    stencil straddles or touches the infinity and the answer is
///    meaningless. Near a singular point, pass a smaller h.
///  - A multi-valued U, e.g. U = -atan2(y, x), has a jump on a cut; the
///    gradient is right everywhere except within 2h of the cut.
///
/// C++ note: `ScalarField&&` is a *forwarding reference*: it accepts lambdas,
/// function objects, and function pointers, temporary or not, without
/// copying. We only call U, never move from it. The static_assert turns the
/// wall of template errors you would otherwise get for a wrong callable into
/// one readable sentence.
///
/// @param U scalar field double(const Vec3&) [J]
/// @param r point [m]
/// @param h finite-difference step [m]; see default_step()
/// @return grad U [J/m = N]
template <class ScalarField>
[[nodiscard]] Vec3 gradient(ScalarField&& U, const Vec3& r, double h) {
    static_assert(std::is_invocable_r_v<double, ScalarField&, const Vec3&>,
                  "gradient() needs a scalar field: a callable double U(const Vec3& r)");
    return {detail::central_difference(U, r, 0, h), detail::central_difference(U, r, 1, h),
            detail::central_difference(U, r, 2, h)};
}

/// gradient() with the default step h = default_step(r) (length scale 1 m).
///
/// C++ note: two overloads (with and without h) instead of one function with
/// "h <= 0 means choose automatically". A sentinel like that turns a typo or
/// a computed h that happens to be 0 into a silent change of behavior; an
/// overload makes "I want the default" something you can see at the call site.
/// For atomic-scale fields pass your own step: gradient(U, r, default_step(r, 1e-10)).
///
/// @param U scalar field double(const Vec3&) [J]
/// @param r point [m]
/// @return grad U [N]
template <class ScalarField>
[[nodiscard]] Vec3 gradient(ScalarField&& U, const Vec3& r) {
    return gradient(U, r, default_step(r));
}

/// Numerical Jacobian of a vector field (all nine first partial derivatives,
/// 12 evaluations of F), with an explicit step. Column layout: see Jacobian.
///
/// Physics: how F changes when you move a little along x, y, and z. For a
/// force from a potential, F = -grad U, it is minus the Hessian of U.
///
/// Pitfalls: as for gradient(): the result is unreliable within 2h of a kink
/// or a singularity of F.
///
/// C++ note: pass a Force by freezing t and v in a lambda:
///     jacobian([&](const Vec3& r) { return force(t, r, Vec3{}); }, r)
///
/// @param F vector field Vec3(const Vec3&) [N]
/// @param r point [m]
/// @param h finite-difference step [m]
/// @return J, with J.dx.y = dFy/dx etc. [N/m]
template <class VectorField>
[[nodiscard]] Jacobian jacobian(VectorField&& F, const Vec3& r, double h) {
    static_assert(std::is_invocable_r_v<Vec3, VectorField&, const Vec3&>,
                  "this needs a vector field: a callable Vec3 F(const Vec3& r). To use a Force "
                  "F(t, r, v), freeze t and v: [&](const Vec3& r) { return force(t, r, Vec3{}); }");
    return {detail::central_difference(F, r, 0, h), detail::central_difference(F, r, 1, h),
            detail::central_difference(F, r, 2, h)};
}

/// jacobian() with the default step h = default_step(r).
/// @param F vector field Vec3(const Vec3&) [N]
/// @param r point [m]
/// @return J [N/m]
template <class VectorField>
[[nodiscard]] Jacobian jacobian(VectorField&& F, const Vec3& r) {
    return jacobian(F, r, default_step(r));
}

/// Numerical curl of a vector field at a point, with an explicit step.
/// See curl(const Jacobian&) for the physics.
///
/// Pitfall: a numerical curl of a conservative field is not exactly zero, only
/// ~1e-11 of the field's variation (rounding), so test "is it small
/// compared with the Jacobian" rather than "== 0". (Fields that are exactly
/// constant do give exactly 0.)
///
/// @param F vector field Vec3(const Vec3&) [N]
/// @param r point [m]
/// @param h finite-difference step [m]
/// @return curl F [N/m]
template <class VectorField>
[[nodiscard]] Vec3 curl(VectorField&& F, const Vec3& r, double h) {
    return curl(jacobian(F, r, h));
}

/// curl() with the default step h = default_step(r).
/// @param F vector field Vec3(const Vec3&) [N]
/// @param r point [m]
/// @return curl F [N/m]
template <class VectorField>
[[nodiscard]] Vec3 curl(VectorField&& F, const Vec3& r) {
    return curl(jacobian(F, r, default_step(r)));
}

/// Numerical divergence of a vector field at a point, with an explicit step.
/// See divergence(const Jacobian&) for the physics.
///
/// @param F vector field Vec3(const Vec3&) [N]
/// @param r point [m]
/// @param h finite-difference step [m]
/// @return div F [N/m]
template <class VectorField>
[[nodiscard]] double divergence(VectorField&& F, const Vec3& r, double h) {
    return divergence(jacobian(F, r, h));
}

/// divergence() with the default step h = default_step(r).
/// @param F vector field Vec3(const Vec3&) [N]
/// @param r point [m]
/// @return div F [N/m]
template <class VectorField>
[[nodiscard]] double divergence(VectorField&& F, const Vec3& r) {
    return divergence(jacobian(F, r, default_step(r)));
}

}  // namespace newton::fields
