#pragma once
// newton/fields/line_integral.hpp — adaptive numerical line integrals ∫ F·dr
// along segments, polylines and circles, plus the circulation ∮ F·dr.
//
// Physics: the work a force does along a path C is W = ∫_C F·dr (Taylor
// §4.1). Taylor's second condition for a conservative force (§4.2) is that W
// between two points does not depend on the path, i.e. that the circulation
// ∮ F·dr around every closed loop is zero. This file turns both statements
// into numbers you can compute for ANY force field you can write as a
// function, however messy, with no closed-form integral needed.
//
// A *vector field* is any callable   Vec3 F(const Vec3& r)   [N, say].
// Fields here depend on position only. A library Force F(t, r, v) becomes a
// field by freezing t and v inside a lambda:
//     [&](const Vec3& r) { return force(t, r, Vec3{}); }
//
// Numerical method: global adaptive Gauss-Kronrod quadrature (G7-K15), the
// QUADPACK "QAG" strategy. See detail::integrate_adaptive for the algorithm
// and for why it is "global" rather than recursive.
//
// C++ note: everything public is a function template over the callable's
// type, so a lambda is called directly (and can be inlined) with no
// std::function overhead. None of these functions is `noexcept`: they call
// user code, and user code may throw. This header does not include (or need)
// derivatives.hpp: the two fields/ headers are independent.

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include <newton/core/constants.hpp>
#include <newton/core/vec3.hpp>

namespace newton::fields {

/// The result of one adaptive line integral.
///
/// C++ note: a numerical integral can fail in ways a plain `double` return
/// cannot express (a singularity on the path, a budget that ran out). Rather
/// than throw, or silently return a wrong number, the routine returns this
/// small aggregate: the answer, how much to trust it, and a `converged`
/// flag. Check `converged` (or use the plain-double overloads, which turn a
/// failure into NaN so it cannot be mistaken for a result). Every member has
/// a default initializer so `LineIntegral{}` is a well-defined "failed, zero"
/// value.
struct LineIntegral {
    double value{0.0};        // ∫ F·dr [J if F is a force]; NaN if F was non-finite on the path
    double error{0.0};        // estimated absolute error, |K15 - G7| summed over pieces [same unit]
    double scale{0.0};        // S = ∫ |F| |dr|: the work if F always pointed along the path [same unit]
    int    evaluations{0};    // number of field evaluations used
    bool   converged{false};  // false: non-finite F met, piece limit hit, pieces collapsed, or bad arguments
};

namespace detail {

// ---- Gauss-Kronrod constants (QUADPACK qk15) -----------------------------
//
// C++ note: namespace-scope `inline constexpr` arrays are the C++17 way to
// have one shared table in a header: `constexpr` makes them compile-time
// constants and `inline` guarantees a single definition program-wide, even
// though many .cpp files include this header (see core/constants.hpp).
//
// Physics/numerics: Gauss-Legendre with n points integrates polynomials up to
// degree 2n-1 exactly. A Gauss-Kronrod pair adds n+1 extra points to an
// n-point Gauss rule so that ALL 2n+1 points integrate degree 3n+1 exactly (3n+2 when n is odd)
// while still containing the Gauss points. With n = 7: 15 field evaluations
// give two answers, K15 (exact for polynomials up to degree 23, i.e. 3n+2 for odd n) and G7 (degree 13), and
// their difference is a free error estimate.
//
// These are the nodes on [-1, 1], positive half only (the rule is symmetric
// about 0). kGkNodes[7] = 0 is the centre. The odd-indexed nodes [1], [3],
// [5] are also the 7-point Gauss nodes (with the centre, node [7]).
inline constexpr double kGkNodes[8] = {
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.0};
// Kronrod weights, one per node above.
inline constexpr double kKronrodWeights[8] = {
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
// Gauss weights for the shared nodes kGkNodes[1], [3], [5], [7].
inline constexpr double kGaussWeights[4] = {
    0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

// ---- limits ---------------------------------------------------------------
// At most this many pieces: 15 + 30 * (200 - 1) = 5985 field evaluations.
inline constexpr int kMaxPieces = 200;
// A piece narrower than this fraction of the whole parameter range "is a
// point": splitting it further cannot help (see integrate_adaptive).
inline constexpr double kMinPieceFraction = 1e-12;
// The integral of |F||dr| must itself be resolved to this fraction of S before
// we declare convergence. Why: the work integrand F·dr is often ODD about a
// singular point (gravity's F·d changes sign through the center), and an odd
// divergence cancels exactly in both K15 and G7, so the work error estimate
// alone can read 0 for a divergent integral. |F||dr| is EVEN (never negative),
// cannot cancel, and so exposes the divergence.
// Trade-off, stated honestly: this is not free. |F||dr| can have a kink even
// on a perfectly smooth path (a spring segment through its anchor has |F| = k|x|,
// a V), and K15 - G7 then exceeds 1e-6 S on the first piece, so such a path
// needs extra splits. Measured: a segment that crosses a spring's anchor takes
// 45-105 evaluations instead of 15, and a spring with a natural length (|F| has
// a kink wherever the spring is relaxed) up to ~600. Results stay accurate;
// only the cost rises. We accept that, because it is exactly what makes a
// symmetric path through a singularity report failure instead of a principal value.
inline constexpr double kScaleErrorTol = 1e-6;

/// What a path integrand returns at parameter s: F·(dr/ds) and |F||dr/ds|.
struct PathSample {
    double work_rate{0.0};   // F · dr/ds  [N·m per unit s]
    double scale_rate{0.0};  // |F| |dr/ds| [N·m per unit s]
};

/// One Gauss-Kronrod piece of the parameter interval [a, b].
struct Piece {
    double a{0.0}, b{0.0};   // parameter interval
    double value{0.0};       // K15 estimate of ∫ F·dr over the piece
    double error{0.0};       // |K15 - G7|
    double scale{0.0};       // K15 estimate of ∫ |F||dr| over the piece
    double scale_error{0.0}; // |K15 - G7| of the |F||dr| integrand over the piece
    bool   finite{true};     // false if any of the 15 samples was NaN or infinite
};

/// Apply the G7-K15 pair to f on [a, b]: 15 evaluations of the integrand.
template <class Integrand>
[[nodiscard]] Piece gk15_piece(Integrand& f, double a, double b) {
    const double mid = 0.5 * (a + b);
    const double half = 0.5 * (b - a);
    const PathSample c = f(mid);
    double kronrod = kKronrodWeights[7] * c.work_rate;
    double gauss = kGaussWeights[3] * c.work_rate;
    double scale = kKronrodWeights[7] * c.scale_rate;
    double scale_gauss = kGaussWeights[3] * c.scale_rate;
    bool finite = std::isfinite(c.work_rate) && std::isfinite(c.scale_rate);
    for (int j = 0; j < 7; ++j) {
        const PathSample lo = f(mid - half * kGkNodes[j]);
        const PathSample hi = f(mid + half * kGkNodes[j]);
        // Summing first is a cheap way to test four values at once: any NaN
        // or infinity poisons the sum.
        finite = finite && std::isfinite(lo.work_rate + hi.work_rate + lo.scale_rate + hi.scale_rate);
        kronrod += kKronrodWeights[j] * (lo.work_rate + hi.work_rate);
        scale += kKronrodWeights[j] * (lo.scale_rate + hi.scale_rate);
        if (j % 2 == 1) {
            gauss += kGaussWeights[j / 2] * (lo.work_rate + hi.work_rate);
            scale_gauss += kGaussWeights[j / 2] * (lo.scale_rate + hi.scale_rate);
        }
    }
    return {a, b, half * kronrod, std::fabs(half * (kronrod - gauss)), half * scale,
            std::fabs(half * (scale - scale_gauss)), finite};
}

/// Global adaptive Gauss-Kronrod quadrature of a path integrand over [a, b].
///
/// Algorithm: start with one piece. Add up the pieces' errors; if the total is
/// small enough, stop. Otherwise split the ONE piece with the largest error
/// into two halves and repeat. "Small enough" means
///     error <= max(rel_tol * |value|, scale_tol * scale)
/// The first term is the usual relative goal. The second is a roundoff floor:
/// when the true answer is 0 (a conservative loop) or tiny compared with the
/// integrand's size, no routine can reach a purely relative goal, because the
/// positive and negative contributions cancel to within ε·S (S = ∫|F||dr|),
/// so we accept an error of a few times 1e-14·S instead.
///
/// C++ note: "global" adaptivity beats the classic recursive "split every
/// piece, halve the tolerance" scheme. With recursion, once the per-piece
/// tolerance falls below roundoff, every piece near a singularity fails and
/// BOTH its children recurse: we measured 24-44 million evaluations on a
/// path ending near a point mass before adding a budget. Here work is spent
/// only where the error actually is, and the piece budget (kMaxPieces) caps
/// the total cost.
///
/// Failure modes, all reported as converged = false:
///  - the integrand was NaN or infinite anywhere sampled (value = NaN);
///  - the piece limit was reached;
///  - the worst piece is narrower than kMinPieceFraction * (b - a): further
///    splitting only chases a singularity (the "pieces collapsed" case).
///
/// Convergence also requires the |F||dr| integral to be resolved (its own
/// |K15 - G7| summed over pieces <= kScaleErrorTol * S). This closes a trap: a
/// divergent integrand that is exactly antisymmetric about a quadrature node,
/// such as PointGravity along a segment through its center at the exact
/// midpoint, cancels in K15 and G7 alike, so the work error estimate is 0 and
/// the integral would "converge" to its principal value. |F||dr| is even about
/// the singularity, cannot cancel, and keeps the pieces splitting until the
/// minimum width ends the hunt with converged = false. While the work error is
/// already acceptable, the piece with the largest |F||dr| error is split next.
/// Cost: this check can add evaluations even on smooth paths when |F||dr| has
/// a kink (a spring segment through its anchor: 45-105 evaluations instead of
/// 15, up to ~600 for a spring with a natural length); see
/// kScaleErrorTol.
template <class Integrand>
[[nodiscard]] LineIntegral integrate_adaptive(Integrand&& f, double a, double b, double rel_tol,
                                              double scale_tol) {
    const double min_width = kMinPieceFraction * std::fabs(b - a);
    std::vector<Piece> pieces{gk15_piece(f, a, b)};
    LineIntegral out;
    out.evaluations = 15;
    for (;;) {
        double value = 0.0, error = 0.0, scale = 0.0, scale_error = 0.0;
        bool finite = true;
        std::size_t worst = 0, worst_scale = 0;
        for (std::size_t i = 0; i < pieces.size(); ++i) {
            value += pieces[i].value;
            error += pieces[i].error;
            scale += pieces[i].scale;
            scale_error += pieces[i].scale_error;
            finite = finite && pieces[i].finite;
            if (pieces[i].error > pieces[worst].error) worst = i;
            if (pieces[i].scale_error > pieces[worst_scale].scale_error) worst_scale = i;
        }
        out.value = value;
        out.error = error;
        out.scale = scale;
        if (!finite) {  // F was NaN/inf somewhere on the path: no trustworthy number exists
            out.value = std::numeric_limits<double>::quiet_NaN();
            out.error = std::numeric_limits<double>::infinity();
            out.converged = false;
            return out;
        }
        const bool work_ok = error <= std::fmax(rel_tol * std::fabs(value), scale_tol * scale);
        const bool scale_ok = scale_error <= kScaleErrorTol * scale;
        if (work_ok && scale_ok) {
            out.converged = true;
            return out;
        }
        // Refine where the unmet criterion is worst.
        const Piece w = pieces[work_ok ? worst_scale : worst];
        const std::size_t split = work_ok ? worst_scale : worst;
        if (static_cast<int>(pieces.size()) >= kMaxPieces || (w.b - w.a) < min_width) {
            out.converged = false;  // budget exhausted or pieces collapsed onto a point
            return out;
        }
        const double m = 0.5 * (w.a + w.b);
        pieces[split] = gk15_piece(f, w.a, m);
        pieces.push_back(gk15_piece(f, m, w.b));
        out.evaluations += 30;
    }
}

/// Orthonormal (u, w) perpendicular to the unit normal n with u × w = n, so
/// that θ ↦ cos θ u + sin θ w runs counterclockwise seen from the tip of n
/// (the right-hand rule).
struct CircleFrame {
    Vec3 u{}, w{};
};

/// Build the frame for any nonzero normal (it is normalized here, so a
/// non-unit n is fine). u comes from crossing n with the
/// coordinate axis it is LEAST aligned with, so the cross product is never
/// small (a fixed axis would break down when n is parallel to it).
[[nodiscard]] inline CircleFrame circle_frame(const Vec3& normal) {
    const Vec3 n = normalized(normal);
    const double ax = std::fabs(n.x), ay = std::fabs(n.y), az = std::fabs(n.z);
    const Vec3 e = (ax <= ay && ax <= az) ? Vec3{1.0, 0.0, 0.0}
                                          : (ay <= az ? Vec3{0.0, 1.0, 0.0} : Vec3{0.0, 0.0, 1.0});
    const Vec3 u = normalized(cross(n, e));
    return {u, cross(n, u)};
}

}  // namespace detail

/// ∫ F·dr along the straight segment a → b, with full diagnostics.
///
/// Physics: the work done by F along a path C from r₁ to r₂ is W = ∫_C F·dr
/// (Taylor §4.1). When F is the net force, the work-KE theorem says W = ΔT.
/// In general W depends on the path. Condition (ii) of §4.2 says a force is
/// conservative exactly when it does not, so comparing two polylines with the
/// same endpoints is the direct numerical test of path independence. This
/// integrates a field F(r) along straight segments. For a force that depends
/// on v (drag), the work depends on how fast the path is traversed, not just
/// on its shape: use analysis::work_done on a trajectory instead.
///
/// The segment is parametrized r(s) = a + s (b - a), s ∈ [0, 1], so
/// ∫ F·dr = ∫₀¹ F(r(s)) · (b - a) ds, which is a one-dimensional integral.
///
/// C++ note: `VectorField&&` is a forwarding reference: it accepts lambdas,
/// function objects and function pointers, temporary or not, without copying
/// (we only call F, never move from it). The static_assert turns the wall of
/// template errors you would get from a wrong callable into one sentence.
/// The result is a struct, not a bare double, so a failed integral (for
/// instance a path through PointGravity's center, where F is infinite)
/// cannot be mistaken for an answer.
///
/// @param F         vector field Vec3 F(const Vec3& r) [N]
/// @param a         start point [m]
/// @param b         end point [m]
/// @param rel_tol   relative accuracy goal for the value (dimensionless)
/// @param scale_tol roundoff floor relative to S = ∫|F||dr| (dimensionless)
/// @return value [J], error estimate [J], scale S [J], evaluation count, converged flag
template <class VectorField>
[[nodiscard]] LineIntegral line_integral_detailed(VectorField&& F, const Vec3& a, const Vec3& b,
                                                  double rel_tol = 1e-12, double scale_tol = 1e-14) {
    static_assert(std::is_invocable_r_v<Vec3, VectorField&, const Vec3&>,
                  "line_integral() needs a vector field: a callable Vec3 F(const Vec3& r). To use a "
                  "Force F(t, r, v), freeze t and v: [&](const Vec3& r) { return force(t, r, Vec3{}); }");
    const Vec3 d = b - a;
    const double len = norm(d);
    auto integrand = [&](double s) {
        const Vec3 f = F(a + s * d);
        return detail::PathSample{dot(f, d), norm(f) * len};
    };
    return detail::integrate_adaptive(integrand, 0.0, 1.0, rel_tol, scale_tol);
}

/// ∫ F·dr along the straight segment a → b, as a plain number [J].
/// Returns NaN if the integral did not converge (see line_integral_detailed),
/// so a failure can never pass for a result.
/// @param F vector field Vec3 F(const Vec3& r) [N]
/// @param a start point [m]
/// @param b end point [m]
template <class VectorField>
[[nodiscard]] double line_integral(VectorField&& F, const Vec3& a, const Vec3& b) {
    const LineIntegral r = line_integral_detailed(F, a, b);
    return r.converged ? r.value : std::numeric_limits<double>::quiet_NaN();
}

/// ∫ F·dr along the polyline path[0] → path[1] → ... (the sum of its
/// straight segments), as a plain number [J].
///
/// Physics: any path can be approximated by a polyline, and work is additive
/// over the pieces of a path. Two polylines with the same endpoints give the
/// same number for a conservative F (Taylor §4.2 condition ii); for the swirl
/// (-y, x, 0) from (1,0,0) to (0,1,0), the direct segment gives 1 J, while
/// the corner route through (1,1,0) gives 2 J.
///
/// Returns NaN if any segment fails to converge, and 0 for fewer than two
/// points (no segments, no work).
/// @param F    vector field Vec3 F(const Vec3& r) [N]
/// @param path corner points in order [m]
template <class VectorField>
[[nodiscard]] double line_integral(VectorField&& F, const std::vector<Vec3>& path) {
    double total = 0.0;
    for (std::size_t i = 0; i + 1 < path.size(); ++i) {
        const LineIntegral seg = line_integral_detailed(F, path[i], path[i + 1]);
        if (!seg.converged) return std::numeric_limits<double>::quiet_NaN();
        total += seg.value;
    }
    return total;
}

/// Circulation ∮ F·dr once around a circle, with full diagnostics.
///
/// Physics: the circulation is the net work F does on a particle carried once
/// around a closed loop. A conservative force does zero work around every
/// closed loop, which is just condition (ii) of §4.2 restated, and by Stokes'
/// theorem ∮F·dr = ∫(∇×F)·n̂ dA over any surface spanning the loop (§4.4):
/// that is how the library's conservative-force classifier uses small loops
/// to look for curl. A nonzero circulation around a loop on which the curl is
/// zero means the loop encircles a "hole": a line where F is singular, like
/// the magnetic field around a current-carrying wire, whose circulation is
/// μ₀I (Ampère's law) although its curl is zero everywhere off the wire. For
/// the hole field (-y, x, 0)/(x² + y²) the circulation around the z axis is
/// exactly 2π; for the swirl (-y, x, 0) a circle of radius R gives 2πR².
///
/// The circle is traversed counterclockwise when seen from the tip of
/// `normal` (the right-hand rule: curl the fingers of your right hand along
/// the loop and your thumb points along `normal`). Flipping the normal flips
/// the sign of the result.
///
/// C++ note: the circle is r(θ) = center + R (cos θ u + sin θ w) with
/// (u, w, n̂) a right-handed orthonormal frame built by detail::circle_frame,
/// then integrated over θ ∈ [0, 2π] by the same adaptive routine as a segment.
/// A zero (or non-finite) normal or a radius that is not positive and finite
/// does not describe a circle: the result is `converged = false` with zero
/// evaluations, rather than a quiet wrong answer.
///
/// @param F         vector field Vec3 F(const Vec3& r) [N]
/// @param center    center of the circle [m]
/// @param normal    normal to the plane of the circle (any nonzero length) [dimensionless direction]
/// @param radius    radius, must be > 0 [m]
/// @param rel_tol   relative accuracy goal (dimensionless)
/// @param scale_tol roundoff floor relative to ∮|F||dr| (dimensionless)
template <class VectorField>
[[nodiscard]] LineIntegral circulation_detailed(VectorField&& F, const Vec3& center, const Vec3& normal,
                                                double radius, double rel_tol = 1e-12,
                                                double scale_tol = 1e-14) {
    static_assert(std::is_invocable_r_v<Vec3, VectorField&, const Vec3&>,
                  "circulation() needs a vector field: a callable Vec3 F(const Vec3& r). To use a "
                  "Force F(t, r, v), freeze t and v: [&](const Vec3& r) { return force(t, r, Vec3{}); }");
    const double nlen = norm(normal);
    if (!(nlen > 0.0) || !std::isfinite(nlen) || !(radius > 0.0) || !std::isfinite(radius)) {
        return {};  // all-default LineIntegral: converged = false
    }
    const detail::CircleFrame fr = detail::circle_frame(normal / nlen);
    auto integrand = [&](double th) {
        const double c = std::cos(th), s = std::sin(th);
        const Vec3 f = F(center + radius * (c * fr.u + s * fr.w));
        const Vec3 dr = radius * (-s * fr.u + c * fr.w);  // dr/dθ
        return detail::PathSample{dot(f, dr), norm(f) * radius};
    };
    return detail::integrate_adaptive(integrand, 0.0, 2.0 * constants::pi, rel_tol, scale_tol);
}

/// Circulation ∮ F·dr around a circle as a plain number [J]; NaN if the
/// integral did not converge or the circle is invalid (see circulation_detailed).
/// @param F      vector field Vec3 F(const Vec3& r) [N]
/// @param center center of the circle [m]
/// @param normal normal to the circle's plane, right-hand rule gives the direction of travel
/// @param radius radius, must be > 0 [m]
template <class VectorField>
[[nodiscard]] double circulation(VectorField&& F, const Vec3& center, const Vec3& normal, double radius) {
    const LineIntegral r = circulation_detailed(F, center, normal, radius);
    return r.converged ? r.value : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace newton::fields
