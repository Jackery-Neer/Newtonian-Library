#pragma once
// newton/potentials/potential_from_force.hpp — from a conservative force to
// its potential energy: U(r) = U_ref - ∫ F·dr (Taylor §4.2, §4.5).
//
// Physics: Taylor defines the potential energy of a conservative force as
// minus the work it does carrying the particle from a fixed reference point
// r_ref to r:
//     U(r) = U_ref - ∫_{r_ref}^{r} F(r')·dr'.
// Path independence (condition (ii) of §4.2) is what makes this a function of
// r alone. This file turns that definition into code: potential_from_force()
// builds U for a force you already know is conservative, and find_potential()
// does the whole chapter-4 job in one call: classify the force, and if a
// potential exists, hand it back, together with where it is zero.
//
// Only differences of U are physical: U_ref is the arbitrary additive
// constant ("the zero of the potential"). Conventions differ by force:
// uniform gravity mgy is zero at y = 0, a spring ½k|r|² at r = 0, but
// gravity and Coulomb are zero at INFINITY, which no finite integral can
// start from (see the U(∞) = 0 recipe at potential_from_force).
//
// C++ note: the potential is returned as a *function object* (a small struct
// with operator()), the same style as the forces: a value you can pass to
// analysis::energy_drift, force_from_potential, or call yourself.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <newton/core/vec3.hpp>
#include <newton/fields/line_integral.hpp>
#include <newton/potentials/conservative.hpp>

namespace newton::potentials {

namespace detail {

/// The result of one recovered line integral: the work, and the scale
/// S = ∫|F||dr| of the path that produced it (NaN / infinity if none did).
struct PathWork {
    double value{std::numeric_limits<double>::quiet_NaN()};
    double scale{std::numeric_limits<double>::infinity()};
};

/// Clamps a point into a box, component by component.
[[nodiscard]] inline Vec3 clamp_to_box(const Vec3& p, const Box& b) {
    return {std::clamp(p.x, b.lo.x, b.hi.x), std::clamp(p.y, b.lo.y, b.hi.y), std::clamp(p.z, b.lo.z, b.hi.z)};
}

/// ∫_a^b F·dr for a field Vec3(const Vec3&), with recovery from the two ways a
/// straight-line integral goes wrong near a singularity. value is NaN if no path works.
///
/// Numerics (research-integration.md Q4): the straight segment is tried first
/// at rel_tol 1e-12. If it fails to converge, OR converges but is
/// ill-conditioned (S = ∫|F||dr| more than 1e3 times |W|: the answer is a
/// small difference of huge numbers, and the quadrature's own error estimate
/// cannot see the rounding of the sample positions), four DETOURS are tried:
/// two-segment paths a -> p -> b, with p the midpoint pushed sideways by
/// 0.25·|b - a| along +u, +w, -u, -w (u, w ⊥ the segment), and, if a region is
/// given, clamped into it (the box is convex, so both legs then stay inside).
/// The converged candidate with the smallest S wins. If nothing converged the
/// value is NaN: never a silently wrong number.
///
/// Physics: a detour gives the same work only if F is path-independent on
/// the region swept by the detour. For a conservative force with a point
/// singularity (gravity, Coulomb) that holds, because removing a point leaves
/// 3D space simply connected. Around a singular LINE it does not: detours on
/// opposite sides of the hole field differ by exactly 2π. So detours are
/// legitimate only after check_conservative() has cleared the region, and only
/// inside it: that is why find_potential exists, and why it hands the box to
/// the potential.
template <class VectorField>
[[nodiscard]] PathWork work_between_detailed(const VectorField& F, const Vec3& a, const Vec3& b,
                                             const std::optional<Box>& region = std::nullopt) {
    const fields::LineIntegral straight = fields::line_integral_detailed(F, a, b, 1e-12, 1e-14);
    const bool ill_conditioned = straight.converged && straight.scale > 1e3 * std::fabs(straight.value);
    if (straight.converged && !ill_conditioned) return {straight.value, straight.scale};

    PathWork best;
    if (straight.converged) best = {straight.value, straight.scale};
    const fields::detail::CircleFrame frame = fields::detail::circle_frame(normalized(b - a));
    const Vec3 sides[4] = {frame.u, frame.w, -frame.u, -frame.w};
    for (const Vec3& side : sides) {
        Vec3 p = 0.5 * (a + b) + (0.25 * norm(b - a)) * side;
        if (region) p = clamp_to_box(p, *region);
        const fields::LineIntegral first = fields::line_integral_detailed(F, a, p, 1e-12, 1e-14);
        const fields::LineIntegral second = fields::line_integral_detailed(F, p, b, 1e-12, 1e-14);
        if (first.converged && second.converged && first.scale + second.scale < best.scale) {
            best = {first.value + second.value, first.scale + second.scale};
            if (!ill_conditioned) break;  // a plain failure: the first good detour is enough
        }
    }
    return best;
}

/// ∫_a^b F·dr with recovery; NaN if no path converged. See work_between_detailed.
template <class VectorField>
[[nodiscard]] double work_between(const VectorField& F, const Vec3& a, const Vec3& b,
                                  const std::optional<Box>& region = std::nullopt) {
    return work_between_detailed(F, a, b, region).value;
}

}  // namespace detail

/// A potential energy defined by a line integral of a force from a reference
/// point: U(r, t) = U_ref - ∫_{r_ref}^{r} F(r', t, v = 0)·dr'.
///
/// Physics: Taylor §4.2 (static) and §4.5 (time-dependent: t is held FIXED
/// during the integral, so each instant has its own U(·, t) with F = -∇U
/// there). The reference value is imposed at every instant, U(r_ref, t) = U_ref:
/// that is one particular gauge, because U → U + c(t) leaves the force
/// unchanged for ANY function c(t). Textbook potentials often use a different
/// gauge (for a spring whose support moves, ½k|r - s(t)|² differs from this
/// by ½k|s(t)|²). Consequence: E(t) = T + U(r, t) depends on the gauge, but
/// the identity E(t) - E(0) = ∫∂U/∂t dt holds in every gauge; compare
/// differences in one gauge, or check that identity.
///
/// Both call forms work: `U(t, r)` for a time-dependent force, and `U(r)`
/// (t = 0) for a static one. analysis::total_energy/energy_drift pick the
/// (t, r) form automatically.
///
/// C++ note: a struct with operator() is a *function object*. It stores a
/// copy of the force, so it is self-contained and can outlive the original.
/// An aggregate with a default for every member (`Force force{}` is only
/// instantiated if you leave the member out, so lambda forces work too).
/// `std::optional<Box> region` says "a Box, or nothing" in the type.
///
/// Cost: one U evaluation is an adaptive quadrature, 15 force evaluations for
/// a smooth force, ~1e3-4e3 near a singularity, and up to 9 × 5985 = 53,865 in
/// the worst case (the straight path plus 4 detours of 2 legs, each capped at
/// 200 pieces); a hole-field graze measured 13,425. analysis::energy_balance
/// treats this potential as time-dependent (it always accepts (t, r)), which
/// costs 4 extra integrals per sample even for a static force: wrap it,
/// `[&](const Vec3& r) { return U(r); }`, to get the static path.
/// `operator()` is const: evaluating U changes nothing.
template <class Force>
struct LineIntegralPotential {
    Force  force{};       // the force; any accepted shape (see conservative.hpp)
    Vec3   r_ref{};       // reference point: U(r_ref) = U_ref [m]
    double U_ref{0.0};    // [J]; for a time-dependent force, U(r_ref, t) = U_ref at every t (a gauge choice)
    /// The region in which the potential was certified (find_potential sets it
    /// to CheckOptions::region). If set, the detours of the integration are
    /// clamped into it, so they never leave the tested box. Empty (the
    /// default of potential_from_force) means "no region": detours are unclamped.
    std::optional<Box> region{};

    /// U(r, t) = U_ref - ∫_{r_ref}^{r} F(r', t, v=0)·dr'  [J]. NaN if no integration
    /// path converged (a singularity on every path): never a silently wrong number.
    /// @param t time at which the force is frozen [s]
    /// @param r where to evaluate the potential [m]
    [[nodiscard]] double operator()(double t, const Vec3& r) const {
        static_assert(detail::shape_v<Force> != detail::Shape::None, "see detail::kShapeMessage");
        const auto f = detail::view_force(force);
        auto field = [&f, t](const Vec3& x) { return f(t, x, Vec3{}); };
        return U_ref - detail::work_between(field, r_ref, r, region);
    }

    /// U(r), evaluated at t = 0: the natural call for a static force.
    /// @param r where to evaluate the potential [m]
    [[nodiscard]] double operator()(const Vec3& r) const { return (*this)(0.0, r); }
};

/// Builds the potential energy of a force you know is conservative.
///
/// Physics: U(r) = U_ref - ∫_{r_ref}^{r} F·dr' (Taylor §4.2), computed along the
/// straight segment from r_ref to r, detouring around a singular point when
/// the segment hits one: allowed because removing a point leaves a 3D region
/// simply connected. This function does NOT check that F is conservative; use
/// find_potential for that. For a time-dependent force, t is held fixed
/// during the integral (§4.5), which gives U(r_ref, t) = U_ref at every t.
///
/// The default r_ref = origin reproduces the textbook forms mgy, ½k|r|² and
/// -F·r exactly. But the origin is the singular point of every central force
/// centered there: U(r) from r_ref = origin is NaN for PointGravity, Coulomb
/// and friends.
///
/// The U(∞) = 0 recipe: gravity, Coulomb, Yukawa, Morse and Lennard-Jones
/// potentials are conventionally zero at infinity, which a finite integral cannot
/// start from. Pass any regular r_ref together with U_ref = U_textbook(r_ref),
/// for example for PointGravity centered at c:
///     potential_from_force(gravity, r_ref, -GM*m / norm(r_ref - c))
/// and U then matches -GMm/r inside the box you tested, to ~1e-12 (it is a line
/// integral, so it is certified only there). For a softened PointGravity, whose
/// potential is the Plummer form -GMm/sqrt(|r - c|² + ε²), use
///     U_ref = -GM*m / sqrt(norm2(r_ref - c) + eps*eps).
///
/// Reliability: each U evaluation is an adaptive quadrature (15 force
/// evaluations for a smooth force, ~1.3e4 in a bad graze, at most 53,865). The result is certified only where F is path-independent and
/// the integration stays inside the region you tested: outside it, a field like
/// the hole (-y, x, 0)/(x² + y²) gives different U along different paths, by
/// 2π. Use find_potential, and stay inside its box.
///
/// @param F     force, any accepted shape (copied)
/// @param r_ref reference point [m], where U = U_ref
/// @param U_ref the value of U at r_ref [J]: your choice of the arbitrary additive constant
template <class Force>
[[nodiscard]] LineIntegralPotential<Force> potential_from_force(Force F, const Vec3& r_ref = {},
                                                                double U_ref = 0.0) {
    static_assert(detail::shape_v<Force> != detail::Shape::None, "see detail::kShapeMessage");
    return LineIntegralPotential<Force>{std::move(F), r_ref, U_ref};
}

/// What find_potential() returns: the verdict, the potential, and where its
/// zero is.
///
/// C++ note: a class template because the potential stores the force's type.
/// `U` is meaningful only if has_potential(): for a force with no potential
/// energy it is still filled in (so the struct is always usable) but its
/// values mean nothing. `note` is a string literal (or empty): no allocation.
template <class Force>
struct PotentialResult {
    ConservativeReport report{};
    LineIntegralPotential<Force> U{};   // meaningful only when has_potential()
    Vec3 r_ref_used{};                  // where U = U_ref actually is [m]
    bool reference_ok{false};           // r_ref_used is inside the box and passed the corner checks
    bool reference_moved{false};        // true if the default reference was singular and an alternative was picked
    const char* note{""};               // why reference_ok is false, if the verdict allowed a potential

    /// A potential exists AND its reference point works: U can be trusted.
    [[nodiscard]] bool has_potential() const { return report.has_potential() && reference_ok; }
};

namespace detail {

/// The 8 corners of a box, in the order (lo/hi) with x fastest.
[[nodiscard]] inline std::vector<Vec3> box_corners(const Box& b) {
    std::vector<Vec3> c;
    c.reserve(8);
    for (int j = 0; j < 8; ++j) {
        c.push_back({(j & 1) ? b.hi.x : b.lo.x, (j & 2) ? b.hi.y : b.lo.y, (j & 4) ? b.hi.z : b.lo.z});
    }
    return c;
}

[[nodiscard]] inline bool inside_box(const Vec3& p, const Box& b) {
    return b.lo.x <= p.x && p.x <= b.hi.x && b.lo.y <= p.y && p.y <= b.hi.y && b.lo.z <= p.z && p.z <= b.hi.z;
}

/// The off-axis waypoint of the path-independence check for the segment
/// ref -> corner: the segment's midpoint pushed sideways by 0.3·|corner - ref|
/// along a direction ⊥ to the segment and generic (the cross product with a
/// fixed, irrational-looking vector g), clamped into the box.
///
/// Why generic: a numerical-gradient potential is wrong inside a core of a few
/// step lengths around its singularity, by an amount that depends on the
/// direction. Paths that all approach along the cube's symmetry axes (the
/// diagonals, the coordinate axes) see the same wrong core, and agree with each
/// other; a generic detour does not.
[[nodiscard]] inline Vec3 check_waypoint(const Vec3& ref, const Vec3& corner, const Box& b) {
    const Vec3 d = corner - ref;
    const Vec3 g1 = normalized(Vec3{0.37, 0.61, 0.70});
    const Vec3 g2 = normalized(Vec3{0.81, -0.29, 0.52});
    Vec3 n = cross(d, g1);
    if (norm(n) < 1e-3 * norm(d)) n = cross(d, g2);  // d nearly parallel to g1: use the second vector
    return clamp_to_box(ref + 0.5 * d + (0.3 * norm(d)) * normalized(n), b);
}

/// True if, at every frozen time, for every corner of the box: the work from
/// `ref` along the recovered straight path and along an independent detour
/// both converge and AGREE to 1e-9·max(S).
///
/// Physics/numerics: this is a spot check of path independence (condition (ii)
/// of Taylor §4.2) at the reference. Convergence alone proves little: a finite
/// F(ref) proves nothing (PointGravity returns {} exactly at its centre, so
/// F looks fine while the integral that starts there diverges), and a
/// numerical gradient (force_from_potential) of a singular U is FINITE
/// near the singularity, so integrals from it converge, to a wrong, direction-
/// dependent core value. Two different paths to the same corner expose both.
template <class Force>
[[nodiscard]] bool reference_is_valid(const Force& F, const Vec3& ref, const std::vector<double>& times,
                                      const std::vector<Vec3>& corners, const Box& box) {
    const auto f = view_force(F);
    const std::optional<Box> region = box;
    for (const double t : times) {
        auto field = [&f, t](const Vec3& x) { return f(t, x, Vec3{}); };
        for (const Vec3& corner : corners) {
            if (norm(corner - ref) == 0.0) continue;  // the work from a point to itself is 0
            const PathWork straight = work_between_detailed(field, ref, corner, region);
            if (!std::isfinite(straight.value)) return false;
            const Vec3 w = check_waypoint(ref, corner, box);
            const fields::LineIntegral first = fields::line_integral_detailed(field, ref, w, 1e-12, 1e-14);
            const fields::LineIntegral second = fields::line_integral_detailed(field, w, corner, 1e-12, 1e-14);
            if (!first.converged || !second.converged) return false;
            const double detour = first.value + second.value;
            const double scale = std::max(straight.scale, first.scale + second.scale);
            if (!(std::fabs(straight.value - detour) <= 1e-9 * scale)) return false;
        }
    }
    return true;
}

}  // namespace detail

/// The headline function: is this force conservative, and if so what is its
/// potential energy?
///
/// Physics: Taylor's chapter-4 question in one call: is this force
/// conservative (§4.2, §4.4), and if so, what is its potential energy? It runs
/// check_conservative and, when a potential exists (U(r) for a conservative
/// force, or U(r, t) for a force that is curl-free at every sampled instant,
/// §4.5), returns it as a line integral from r_ref. Otherwise the report says
/// why not. The potential is certified only inside the checked region: outside
/// it, a field like the hole (-y, x, 0)/(x² + y²) can give different U for
/// different paths; the returned potential therefore keeps the box and never
/// lets its detours leave it. Gauge: for a time-dependent force,
/// U(r_ref, t) = U_ref at every t, an arbitrary choice of the function c(t)
/// that U → U + c(t) leaves free (see LineIntegralPotential).
///
/// The reference point (research-integration.md Q6):
///  - If you pass no r_ref: the origin when it lies in the box (inclusive), so
///    that mgy, ½k|r|² and -F·r come out in textbook form; otherwise the box
///    centre.
///  - It is then VALIDATED: from it to each of the 8 box corners, at t_min and
///    (for a time-dependent force) at every time slice, the work along the
///    straight path and along an independent off-axis detour must both converge
///    and agree to 1e-9. This catches a singular reference (PointGravity at its
///    centre) and a reference inside the wrong core of a numerically
///    differentiated potential (find_potential(force_from_potential(-GMm/r))).
///  - If validation fails and you did NOT choose r_ref, the point is moved: the
///    candidates r_ref + ½(corner_j - r_ref), j = 1..8, then the same points
///    shifted by a generic offset (so that no path to a corner runs through the
///    old reference) are tried in turn and the first valid one is used;
///    `reference_moved` is then true and `r_ref_used` tells you where
///    U = U_ref actually is. Shift U yourself with U_ref if you need a textbook zero.
///  - A r_ref you chose is NEVER moved. If it is outside the box, or fails
///    validation, `reference_ok` is false and `note` / describe(result) say why.
///    Choose a regular point inside the box.
///  - If the verdict has no potential, validation is skipped (`reference_ok`
///    stays false).
///
/// C++ note: `std::optional<Vec3>` says "a Vec3, or nothing" in the type. It
/// replaces a sentinel such as "{0,0,0} means unset", which would make the
/// origin impossible to request. Test it with has_value(), read it with
/// value() / value_or(default). std::nullopt is the empty one.
///
/// @param F      force, any accepted shape (copied into the result)
/// @param opt    region [m], time window [s], tolerances (see CheckOptions)
/// @param r_ref  reference point [m], inside the box; empty = choose automatically
/// @param U_ref  value of U at the reference point [J]
template <class Force>
[[nodiscard]] PotentialResult<Force> find_potential(Force F, const CheckOptions& opt = {},
                                                    std::optional<Vec3> r_ref = std::nullopt,
                                                    double U_ref = 0.0) {
    static_assert(detail::shape_v<Force> != detail::Shape::None, "see detail::kShapeMessage");
    const ConservativeReport report = check_conservative(F, opt);

    const bool caller_chose = r_ref.has_value();
    const Box box = opt.region;
    const bool origin_in_box = detail::inside_box(Vec3{}, box);
    const Vec3 ref = r_ref.value_or(origin_in_box ? Vec3{} : 0.5 * (box.lo + box.hi));

    bool ok = false;
    bool moved = false;
    const char* note = "";
    Vec3 used = ref;
    if (report.has_potential()) {
        if (!detail::inside_box(ref, box)) {
            note = "the reference point is outside the tested box, where U is not certified";
        } else {
            // t_min always, then every golden-ratio slice for a time-dependent force.
            std::vector<double> times{opt.t_min};
            if (report.depends_on_time) {
                const std::vector<double> slices = detail::slice_times(opt, true);
                times.insert(times.end(), slices.begin(), slices.end());
            }
            const std::vector<Vec3> corners = detail::box_corners(box);
            ok = detail::reference_is_valid(F, ref, times, corners, box);
            if (!ok && !caller_chose) {
                std::vector<Vec3> candidates;
                for (const Vec3& corner : corners) candidates.push_back(ref + 0.5 * (corner - ref));
                // The same points, nudged by a generic vector: paths from a half-corner to the opposite
                // corner all run through `ref`, which is the very point that just failed.
                const Vec3 nudge = (0.2 * detail::box_length_scale(box)) * normalized(Vec3{0.37, 0.61, 0.70});
                for (const Vec3& corner : corners) {
                    candidates.push_back(detail::clamp_to_box(ref + 0.5 * (corner - ref) + nudge, box));
                }
                for (const Vec3& candidate : candidates) {
                    if (detail::reference_is_valid(F, candidate, times, corners, box)) {
                        used = candidate;
                        ok = true;
                        moved = true;
                        break;
                    }
                }
            }
            if (!ok) {
                note = caller_chose ? "an integral from the chosen reference point to a corner of the box "
                                      "failed or depends on the path (a singularity on the way)"
                                    : "no usable reference point was found in the box";
            }
        }
    }
    LineIntegralPotential<Force> U{std::move(F), used, U_ref, box};
    return PotentialResult<Force>{report, std::move(U), used, ok, moved, note};
}

/// describe() for a PotentialResult: the verdict, then what happened to the
/// reference point.
///
/// @param r a result from find_potential
template <class Force>
[[nodiscard]] std::string describe(const PotentialResult<Force>& r) {
    std::string s = describe(r.report);
    if (!r.report.has_potential()) return s;
    char buf[480];
    if (r.reference_ok && r.reference_moved) {
        std::snprintf(buf, sizeof buf,
                      " The default reference point was singular (or inside the unreliable core of a "
                      "numerical derivative), so U = U_ref was moved to (%.4g, %.4g, %.4g) m.",
                      r.r_ref_used.x, r.r_ref_used.y, r.r_ref_used.z);
        s += buf;
    } else if (!r.reference_ok) {
        std::snprintf(buf, sizeof buf,
                      " But the reference point (%.4g, %.4g, %.4g) m is not usable: %s, so U is not "
                      "trustworthy. Choose a regular reference point inside the box.",
                      r.r_ref_used.x, r.r_ref_used.y, r.r_ref_used.z,
                      (r.note != nullptr && r.note[0] != '\0') ? r.note : "validation failed");
        s += buf;
    }
    return s;
}

}  // namespace newton::potentials
