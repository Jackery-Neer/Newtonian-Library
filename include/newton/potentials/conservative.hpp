#pragma once
// newton/potentials/conservative.hpp — "Is this force conservative?"
// (Taylor §4.2-4.5), answered by experiment.
//
// Physics: Taylor §4.2 gives two conditions for a force to be conservative:
//   (i)  F depends only on the particle's position r (not on v, not on t);
//   (ii) the work W = ∫F·dr between two points is the same along every path
//        (equivalently: ∮F·dr = 0 around every closed loop).
// Only then does a potential energy U(r) exist, and E = T + U is conserved.
// §4.4 turns condition (ii) into a local test, curl F = 0 (valid on a simply
// connected region), and §4.5 relaxes (i) a little: a force that depends on t
// but is curl-free at every instant has a potential U(r, t), though the
// energy then changes at the rate dE/dt = ∂U/∂t.
//
// check_conservative() applies those conditions, in that order, to ANY force
// you can write as a function, by sampling it: it probes for v and t
// dependence, measures the curl at random points, and integrates F around
// closed loops. The verdict is *evidence, not proof*: a nonzero curl that is
// confirmed at two step sizes and by a small loop (Stokes), or a loop that does
// nonzero work, is strong evidence a force is NOT conservative; sampling can
// only ever gather evidence that it IS. The report says how much evidence
// (sample counts, worst numbers).
//
// Accepted force shapes (detected at compile time, checked in this order):
//     Vec3 F(double t, const Vec3& r, const Vec3& v)   the library's Force
//     Vec3 F(double t, const Vec3& r)                  a time-dependent field
//     Vec3 F(const Vec3& r)                            a static field
//
// C++ note: the shape detection uses `if constexpr` and
// std::is_invocable_r_v, see detail::shape_v. A callable that matches none of
// the three is rejected with a readable static_assert, not a page of template
// errors. Nothing here is `noexcept`: it all calls user code.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <newton/core/vec3.hpp>
#include <newton/fields/derivatives.hpp>
#include <newton/fields/line_integral.hpp>

namespace newton::potentials {

// --- options and result types ----------------------------------------------------

/// An axis-aligned box [lo, hi]: the region of space in which a force is
/// tested [m].
///
/// Physics: "conservative" is a statement about a REGION. The hole field
/// (-y, x, 0)/(x² + y²) is conservative in a box that avoids the z axis and
/// not conservative in a box that contains it (the Stokes argument behind
/// Taylor §4.4 needs a simply connected region: one where every loop can be
/// shrunk to a point).
///
/// C++ note: an aggregate with a default for every member, so `Box{}` is the
/// cube [-1, 1]³ and `Box{{0, 0, 0}, {2, 2, 2}}` is a custom one.
struct Box {
    Vec3 lo{-1.0, -1.0, -1.0};  // corner with the smallest x, y, z [m]
    Vec3 hi{1.0, 1.0, 1.0};     // corner with the largest x, y, z [m]
};

/// Knobs for check_conservative(), find_potential() and probe_dependence().
/// The defaults were chosen from measurements on a 30-force battery (the
/// loop and curl margins are quoted per member below); most users only set
/// `region` and, for a time-dependent force, the time window.
///
/// C++ note: an aggregate with a default member initializer on every member,
/// so you can write `CheckOptions opt; opt.region = Box{...};` and everything
/// you do not touch keeps a sensible default.
struct CheckOptions {
    /// Where F is tested [m]. It must contain the motion you care about, and
    /// must have positive size in all three directions (a flat, inverted or
    /// non-finite box gives Inconclusive). Its half-largest side L is also the
    /// *length scale* of the test: finite-difference steps are
    /// 1e-4·max(|r - box centre|, L), measured from the BOX CENTRE, so a
    /// molecule-sized box (L ~ 1e-10 m) gets molecule-sized steps wherever it
    /// sits. (Limit: far from the origin the coordinates themselves are only
    /// resolved to ~1e-16·|r|, so a box much smaller than 1e-10·|r| cannot be
    /// resolved at all.)
    Box    region{};
    double t_min{0.0};   // [s] start of the window for the t-probe and the time slices
    double t_max{10.0};  // [s] end of that window
    /// [m/s] The v-probe uses velocities whose components lie in [-v_max, v_max];
    /// every probe sample also compares against v = 0 (where speed thresholds
    /// and v/|v| factors behave specially). Measured reach: a speed threshold
    /// is found reliably up to about 1.4·v_max (the corners of the cube), and
    /// not beyond. Likewise a switch-on time very close to t_max is found only
    /// some of the time.
    double v_max{10.0};
    /// The curl, loop and time-slice counts must be positive. A callable that
    /// can depend on time or velocity also needs positive dependence_samples.
    /// check_conservative() returns Inconclusive rather than quietly
    /// certifying a force after a caller has disabled required evidence.
    int    dependence_samples{64};  // probe samples (v and t together)
    int    curl_samples{256};       // curl sample points, per time slice
    int    time_slices{5};          // instants tested if F depends on t (a static force uses one slice at t_min)
    /// How many of the loop-grid sizes {1, 2, 3, 5} to use (1 to 4; default 4).
    /// In each of 3 parallel planes per axis (the box centre and centre ± half
    /// a side/2) a g×g grid of rectangular cell loops tiles the box's cross
    /// section, for each size g used: 3 axes × 3 planes × (1 + 4 + 9 + 25)
    /// = 351 loops per time slice at the default, 4 segment integrals each.
    /// The sizes are coprime so that a singular line cannot sit on the cell
    /// edges of every grid at once.
    int    loop_grid_levels{4};
    /// A sample's relative curl (see fields::detail::relative_curl) above this
    /// is a candidate for NonzeroCurl. Measured on tests/test_potentials.cpp's
    /// battery: the worst conservative field scores 2.7e-13 (the error term
    /// scores pure noise and near-singular samples as 0), and the weakest
    /// non-conservative field it must catch, (y, (1 + 1e-4) x, 0), scores 4.9e-5.
    double curl_tol{1e-6};
    /// A loop's |∮F·dr| / ∮|F||dr| above this is PathDependent. Measured:
    /// worst conservative loop 3.3e-15, weakest non-conservative one 5.0e-5.
    double loop_tol{1e-8};
    /// Seed of the deterministic random generator (std::mt19937_64, with a hand-
    /// written conversion to [0, 1), because std::uniform_real_distribution
    /// is allowed to differ between standard libraries and we want the same
    /// verdict on every machine).
    std::uint64_t seed{20261124};
};

/// What changing t or v did to F. See probe_dependence().
struct Dependence {
    bool on_time{false};      // F(t1, r, v) != F(t2, r, v) for some r, v
    bool on_velocity{false};  // F(t, r, v1) != F(t, r, v2) for some t, r
};

/// The answer to "is this force conservative?", ordered as Taylor's tests are
/// applied (§4.2 condition (i) first, then (ii)).
///
/// C++ note: `enum class` (a *scoped* enumeration): the names live inside
/// `Verdict::`, so `Verdict::Conservative` cannot clash with another enum's
/// `Conservative`, and a Verdict never converts silently to an int. Use
/// to_string() when you want to print one.
enum class Verdict {
    /// F = F(r), curl F ~ 0 and closed loops do no work in this region: U(r) exists and E = T + U is conserved (§4.2-4.4).
    Conservative,
    /// Depends on t, not v, and conservative at each of the sampled frozen instants: U(r, t) exists with F = -grad U, but dE/dt = ∂U/∂t (§4.5). Only `time_slices` instants are tested: a curl that exists only between them is missed (increase time_slices).
    TimeDependentPotential,
    /// Depends on v: condition (i) fails, no potential energy of any kind (§4.2). See ConservativeReport::velocity_power.
    VelocityDependent,
    /// curl F != 0 at a sampled point (and instant), confirmed at two step sizes and by a small loop (Stokes): the work is path-dependent (§4.4).
    NonzeroCurl,
    /// curl F ~ 0 at every sample, yet a closed loop does nonzero work: a singular LINE through the region (a "hole"), or a curl too weak or too concentrated for point samples to see.
    PathDependent,
    /// Too many non-finite force values, the loop integrals did not converge, or the region is degenerate: no verdict can be given (see ConservativeReport::note).
    Inconclusive,
};

/// For a VelocityDependent force: what does its velocity-dependent PART do to
/// the energy? Taylor §4.1: the power of a force is P = F·v.
///
/// Physics: the label describes F_v = F(t, r, v) - F(t, r, 0), the part that
/// exists only because the particle moves (see velocity_dependent_part), and
/// is judged from the sign of cos = v·F_v / (|v||F_v|) at the probe samples.
/// It does NOT say what the whole force does. In particular, for drag in a
/// wind, F = -b(v - w), the constant push b·w belongs to the "at rest" part
/// (it is booked as a formal potential -b w·r) and the label is Dissipative
/// even though the true drag power -b(v - w)·v can be positive; pass such a
/// force whole as F_nc to analysis::energy_balance (research-physics.md §8(b)).
enum class VelocityPower {
    NotApplicable,  // the verdict is not VelocityDependent
    NoWork,         // F_v ⊥ v at every usable probe (magnetic force q v×B): on its own it does no work, so it never changes the kinetic energy
    Dissipative,    // F_v·v < 0 at every usable probe (drag, kinetic friction): it removes energy
    Feeding,        // F_v·v > 0 at every usable probe (negative damping +b v): it supplies energy
    Mixed,          // both signs genuinely occur (or F_v is partly perpendicular, partly not)
    Unknown,        // no usable cosine: F_v was not finite at rest (friction -mu v/|v| is NaN at v = 0) or too small to measure
};

/// Everything check_conservative() learned. The verdict first; the rest is the
/// evidence behind it, so a student can judge how much to trust it.
struct ConservativeReport {
    Verdict       verdict{Verdict::Inconclusive};
    VelocityPower velocity_power{VelocityPower::NotApplicable};
    bool   depends_on_time{false};
    bool   depends_on_velocity{false};
    /// Worst fields::detail::relative_curl over the samples [dimensionless],
    /// or 0 if no curl was computed (a velocity-dependent force stops at the
    /// probe). Candidates that failed the h/2 or Stokes confirmation are
    /// excluded (they are counted in unconfirmed_curl_candidates). For
    /// NonzeroCurl this is the confirmed sample that showed it; the search
    /// stops there.
    double max_relative_curl{0.0};
    Vec3   worst_curl_point{};   // [m] where max_relative_curl occurred
    double worst_curl_time{0.0}; // [s] the frozen time of that sample
    /// Worst |∮F·dr| / ∮|F||dr| over the loops that converged [dimensionless].
    double max_loop_ratio{0.0};
    Vec3   worst_loop_center{};   // [m] centre of the loop with the worst ratio: where to look for a hole
    double worst_loop_time{0.0};  // [s] the frozen time of that loop
    int    samples_tested{0};    // curl samples evaluated (all slices)
    int    skipped_samples{0};   // curl or probe samples skipped: a non-finite force value (a singularity)
    int    unconfirmed_curl_candidates{0};  // samples above curl_tol that failed the h/2 check (artifacts, ignored)
    int    loops_tested{0};
    int    loops_skipped{0};     // loop integrals that did not converge: inconclusive, never a failure
    const char* note{""};        // why the verdict is Inconclusive (a string literal); empty otherwise

    /// True for Conservative and TimeDependentPotential: some U exists.
    [[nodiscard]] bool has_potential() const {
        return verdict == Verdict::Conservative || verdict == Verdict::TimeDependentPotential;
    }
    /// True for Conservative only: U(r) exists and E = T + U is conserved.
    [[nodiscard]] bool is_conservative() const { return verdict == Verdict::Conservative; }
};

/// Short name of a verdict, e.g. "conservative".
[[nodiscard]] inline const char* to_string(Verdict v) {
    switch (v) {
        case Verdict::Conservative:           return "conservative";
        case Verdict::TimeDependentPotential: return "time-dependent potential";
        case Verdict::VelocityDependent:      return "velocity-dependent";
        case Verdict::NonzeroCurl:            return "nonzero curl";
        case Verdict::PathDependent:          return "path-dependent";
        case Verdict::Inconclusive:           return "inconclusive";
    }
    return "unknown";  // unreachable for a valid enum; keeps -Wreturn-type quiet
}

/// Short name of a velocity-power label, e.g. "dissipative".
[[nodiscard]] inline const char* to_string(VelocityPower p) {
    switch (p) {
        case VelocityPower::NotApplicable: return "not applicable";
        case VelocityPower::NoWork:        return "no work";
        case VelocityPower::Dissipative:   return "dissipative";
        case VelocityPower::Feeding:       return "feeding";
        case VelocityPower::Mixed:         return "mixed";
        case VelocityPower::Unknown:       return "unknown";
    }
    return "unknown";
}

/// One or two student-facing sentences explaining the verdict, with the
/// evidence (worst point and time where relevant).
///
/// Physics: the wording follows Taylor's two conditions (§4.2), the curl test
/// (§4.4) and time-dependent potentials (§4.5).
[[nodiscard]] inline std::string describe(const ConservativeReport& r) {
    char buf[480];
    std::string s;
    switch (r.verdict) {
        case Verdict::Conservative:
            std::snprintf(buf, sizeof buf,
                          "conservative: F depends only on r, its curl is zero and closed loops do no "
                          "work in this region, so U(r) exists and E = T + U is conserved (Taylor "
                          "sec. 4.2-4.4). Evidence, not proof: worst relative curl %.2g over %d samples, "
                          "worst loop ratio %.2g over %d loops.",
                          r.max_relative_curl, r.samples_tested, r.max_loop_ratio,
                          r.loops_tested - r.loops_skipped);
            s = buf;
            break;
        case Verdict::TimeDependentPotential:
            std::snprintf(buf, sizeof buf,
                          "not conservative (F depends on t), but curl-free at every sampled instant: "
                          "U(r,t) exists with F = -grad U, and E = T + U changes at the rate "
                          "dE/dt = partial U/partial t, the change of U at fixed r (Taylor sec. 4.5). "
                          "Worst relative curl %.2g.",
                          r.max_relative_curl);
            s = buf;
            break;
        case Verdict::VelocityDependent:
            s = "not conservative: F depends on v, so no potential energy exists (Taylor sec. 4.2, "
                "condition (i)). About its velocity-dependent PART F_v = F(v) - F(0):";
            switch (r.velocity_power) {
                case VelocityPower::NoWork:
                    s += " it is perpendicular to v (like q v x B), so it does no work, and on its own "
                         "it never changes the kinetic energy.";
                    break;
                case VelocityPower::Dissipative:
                    s += " it opposes v (like drag or friction), removing energy at the rate -F_v.v.";
                    break;
                case VelocityPower::Feeding:
                    s += " it points along v (like negative damping), so it supplies energy at the "
                         "rate F_v.v.";
                    break;
                case VelocityPower::Mixed:
                    s += " it both opposes and aids the motion in different places, so it can remove "
                         "or supply energy.";
                    break;
                case VelocityPower::Unknown:
                case VelocityPower::NotApplicable:
                    s += " its effect on the energy could not be determined (no usable sample: F was "
                         "not finite at rest, or F_v was too small to measure).";
                    break;
            }
            s += " (This describes the velocity-dependent part only: whatever F does at rest, e.g. "
                 "the constant push b w of drag in a wind, is not included.)";
            break;
        case Verdict::NonzeroCurl:
            std::snprintf(buf, sizeof buf,
                          "not conservative: curl F is nonzero at (%.4g, %.4g, %.4g) m [t = %.4g s] "
                          "(relative curl %.2g), so the work depends on the path (Taylor sec. 4.4).",
                          r.worst_curl_point.x, r.worst_curl_point.y, r.worst_curl_point.z,
                          r.depends_on_time ? r.worst_curl_time : 0.0, r.max_relative_curl);
            s = buf;
            if (!r.depends_on_time) {  // a static field has no time to report
                const std::size_t k = s.find(" [t = ");
                const std::size_t e = s.find("s]", k);
                if (k != std::string::npos && e != std::string::npos) s.erase(k, e + 2 - k);
            }
            break;
        case Verdict::PathDependent:
            std::snprintf(buf, sizeof buf,
                          "not conservative in this region: curl F ~ 0 at every sample, yet a closed "
                          "loop does nonzero work (ratio %.2g, loop centred at (%.4g, %.4g, %.4g) m). "
                          "Either F is singular on a line through the region (like the magnetic field "
                          "around a wire), or it has a curl too weak or too concentrated for point "
                          "samples to see.",
                          r.max_loop_ratio, r.worst_loop_center.x, r.worst_loop_center.y,
                          r.worst_loop_center.z);
            s = buf;
            break;
        case Verdict::Inconclusive:
            std::snprintf(buf, sizeof buf,
                          "no verdict: %s (%d samples skipped, %d of %d loops skipped).",
                          (r.note != nullptr && r.note[0] != '\0')
                              ? r.note
                              : "F was not finite at too many sample points, or an integral did not converge",
                          r.skipped_samples, r.loops_skipped, r.loops_tested);
            s = buf;
            break;
    }
    return s;
}

// --- implementation details --------------------------------------------------------

namespace detail {

/// Which signature a force callable has. See shape_v.
enum class Shape { None, StaticField, TimeField, Force };

/// Compile-time detection of the callable's shape, in the spec's order.
///
/// C++ note: std::is_invocable_r_v<R, F, Args...> is true when `F` can be
/// called with arguments of types Args... and the result converts to R. It is
/// a *variable template* (the _v suffix), evaluated at compile time. We test
/// `const F&` because the classifier only ever calls the force through a
/// const reference (a `mutable` lambda would be rejected, on purpose: a force
/// that changes between calls would break the whole "same input, same output"
/// premise). A force that is a plain `Vec3 f(const Vec3&)` *function* works
/// too: pass its address or wrap it in a lambda.
template <class F>
inline constexpr Shape shape_v =
    std::is_invocable_r_v<Vec3, const F&, double, const Vec3&, const Vec3&> ? Shape::Force
    : std::is_invocable_r_v<Vec3, const F&, double, const Vec3&>            ? Shape::TimeField
    : std::is_invocable_r_v<Vec3, const F&, const Vec3&>                    ? Shape::StaticField
                                                                            : Shape::None;

inline constexpr const char* kShapeMessage =
    "this needs a force callable as Vec3 F(double t, const Vec3& r, const Vec3& v) (the library's "
    "Force), or Vec3 F(double t, const Vec3& r) (a time-dependent field), or Vec3 F(const Vec3& r) "
    "(a static field), and it must be callable on a const object (no `mutable` lambdas)";

/// Wraps a force of ANY accepted shape as the 3-argument Force concept
/// (t, r, v), ignoring the arguments the original does not take.
///
/// C++ note: `F` may be a value type (the wrapper owns a copy: see as_force)
/// or a const reference type (the wrapper only looks at the caller's object:
/// see view_force). A reference member is fine in an aggregate like this; the
/// price is that the wrapper must not outlive what it refers to. (Unlike the
/// library's other aggregates, `f` has no default member initializer: a
/// reference member cannot have a sensible default, and a lambda force has no
/// default constructor.)
template <class F>
struct AsForce {
    F f;
    [[nodiscard]] Vec3 operator()([[maybe_unused]] double t, const Vec3& r,
                                  [[maybe_unused]] const Vec3& v) const {
        constexpr Shape s = shape_v<std::remove_cv_t<std::remove_reference_t<F>>>;
        static_assert(s != Shape::None, "see detail::kShapeMessage");
        if constexpr (s == Shape::Force) {
            return f(t, r, v);
        } else if constexpr (s == Shape::TimeField) {
            return f(t, r);
        } else {
            return f(r);
        }
    }
};

/// An owning wrapper: a copy of `f` adapted to (t, r, v). Used for the
/// callables the split helpers return, which must keep their force alive.
template <class F>
[[nodiscard]] AsForce<F> as_force(F f) {
    return AsForce<F>{std::move(f)};
}

/// A non-owning wrapper around the caller's force; free to construct.
template <class F>
[[nodiscard]] AsForce<const F&> view_force(const F& f) {
    return AsForce<const F&>{f};
}

/// The deterministic random source. std::mt19937_64's output sequence is
/// fully specified by the C++ standard, so it is bit-for-bit the same on
/// every compiler; the conversion to doubles is ours for the same reason.
class Rng {
public:
    explicit Rng(std::uint64_t seed) : engine_(seed) {}

    /// Uniform in [0, 1): the top 53 bits of a 64-bit draw times 2^-53. A
    /// double has 53 significant bits, so every value is hit exactly once;
    /// `0x1.0p-53` is a hexadecimal floating-point literal (C++17) for 2^-53.
    [[nodiscard]] double uniform() { return static_cast<double>(engine_() >> 11) * 0x1.0p-53; }
    [[nodiscard]] double uniform(double a, double b) { return a + (b - a) * uniform(); }
    [[nodiscard]] Vec3 in_box(const Box& b) {
        const double x = uniform(b.lo.x, b.hi.x);
        const double y = uniform(b.lo.y, b.hi.y);
        const double z = uniform(b.lo.z, b.hi.z);
        return {x, y, z};
    }
    [[nodiscard]] Vec3 in_cube(double half_side) {
        const double x = uniform(-half_side, half_side);
        const double y = uniform(-half_side, half_side);
        const double z = uniform(-half_side, half_side);
        return {x, y, z};
    }

private:
    std::mt19937_64 engine_;
};

/// The golden-ratio conjugate, (sqrt(5) - 1)/2. See slice_times().
inline constexpr double kGoldenConjugate = 0.6180339887498949;

/// The instants at which F is frozen and tested: one at t_min for a force that
/// does not depend on time, otherwise K golden-ratio offsets
///     t_k = t_min + (k + 0.618...)/K · (t_max - t_min),  k = 0..K-1.
///
/// Physics: the induced swirl sin(t)(-y, x, 0) has curl 2 sin(t) ẑ, which is
/// exactly zero at t = 0, π, 2π, ...  A test that only looks at t = 0 (or at
/// "nice" round instants) would call it conservative. An irrational offset
/// inside each of K equal sub-windows never lands on such an instant and, being
/// deterministic, gives the same answer every run.
[[nodiscard]] inline std::vector<double> slice_times(const CheckOptions& opt, bool time_dependent) {
    if (!time_dependent) return {opt.t_min};
    const int K = std::max(opt.time_slices, 1);
    std::vector<double> t;
    t.reserve(static_cast<std::size_t>(K));
    for (int k = 0; k < K; ++k) {
        t.push_back(opt.t_min + (static_cast<double>(k) + kGoldenConjugate) / static_cast<double>(K) *
                                    (opt.t_max - opt.t_min));
    }
    return t;
}

/// The length scale L of a box: half its largest side [m] (1 m if the box is degenerate).
[[nodiscard]] inline double box_length_scale(const Box& b) {
    const double side = std::max({b.hi.x - b.lo.x, b.hi.y - b.lo.y, b.hi.z - b.lo.z});
    return side > 0.0 ? 0.5 * side : 1.0;
}

/// Component `axis` (0 = x, 1 = y, 2 = z) of a vector, as a reference.
[[nodiscard]] inline double& component(Vec3& v, int axis) {
    return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
}
[[nodiscard]] inline double component(const Vec3& v, int axis) {
    return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
}

/// The loop-grid sizes: g×g cells per plane. Pairwise coprime, so a singular
/// line cannot lie on the cell edges of every grid at once.
inline constexpr int kLoopGridSizes[4] = {1, 2, 3, 5};

/// ∮F·dr around one rectangular cell of a plane, counterclockwise about +axis a:
/// the plane is x_a = mid_a + offset, the cell is [i0, i1] × [j0, j1] in the
/// in-plane coordinates (i, j) with (a, i, j) a cyclic permutation of (x, y, z),
/// so that (i-direction) × (j-direction) = +a.
///
/// The four straight sides are integrated separately with
/// fields::line_integral_detailed(rel_tol 1e-12, scale_tol): the values and
/// scales (S = ∫|F||dr|) add, and the loop is `converged` only if every side
/// is (a side that touches a singularity, or meets a non-finite F, makes the
/// loop inconclusive: skipped, never evidence).
template <class VectorField>
[[nodiscard]] fields::LineIntegral cell_circulation(VectorField& field, const Vec3& mid, int a, int i, int j,
                                                    double offset, double i0, double i1, double j0, double j1,
                                                    double scale_tol) {
    auto corner = [&](double ci, double cj) {
        Vec3 p = mid;
        component(p, a) += offset;
        component(p, i) = ci;
        component(p, j) = cj;
        return p;
    };
    const Vec3 c[4] = {corner(i0, j0), corner(i1, j0), corner(i1, j1), corner(i0, j1)};
    fields::LineIntegral total;
    total.converged = true;
    for (int k = 0; k < 4; ++k) {
        const fields::LineIntegral side = fields::line_integral_detailed(field, c[k], c[(k + 1) % 4], 1e-12, scale_tol);
        total.value += side.value;
        total.error += side.error;
        total.scale += side.scale;
        total.evaluations += side.evaluations;
        total.converged = total.converged && side.converged;
    }
    return total;
}

/// What the dependence probe found.
struct ProbeResult {
    Dependence    dep{};
    VelocityPower power{VelocityPower::NotApplicable};  // meaningful if dep.on_velocity
    int           samples{0};  // probe samples attempted (0 for a position-only signature)
    int           skipped{0};  // samples whose base force value was not finite
};

/// Threshold on |cos| below which F_v is "perpendicular to v" (and on cos
/// beyond ±threshold for "opposing" / "aiding" v).
///
/// Numerics: F_v is a SUBTRACTION, F(v) - F(0), so its rounding error is
/// ~ eps·max(|F(v)|, |F(0)|) however small F_v itself is. The cosine
/// divides by |F_v|, so its noise is about that error over |F_v|. A sample is
/// therefore used only if the noise estimate 2·eps·max(|F(v)|, |F(0)|)/|F_v|
/// is below kCosineTol, i.e. |F_v| >= 2·eps·max(|F(v)|, |F(0)|)/kCosineTol
/// (about 4.4e-7 of |F|);
/// a force whose velocity part is smaller than that (a tiny magnetic part on a
/// huge E field) gives Unknown rather than a wrong label.
inline constexpr double kCosineTol = 1e-9;
inline constexpr double kEps = 2.220446049250313e-16;  // double precision machine epsilon

/// The probe, on a force already adapted to (t, r, v).
///
/// Each sample draws r in the box, two times t1, t2 in [t_min, t_max], and a
/// velocity va in the cube [-v_max, v_max]³. It evaluates base = F(t1, r, va),
/// then:
///  - time:     compares with F(t2, r, va);
///  - velocity: compares with F(t1, r, 0) (the "at rest" force) and, on even
///              samples, with F(t1, r, vb) for a second random velocity vb.
/// Any exact difference between two FINITE values proves dependence. A
/// sample whose base value is not finite is skipped (a singularity or a NaN).
///
/// The same data give the velocity power: with F_v = F(t1, r, va) - F(t1, r, 0),
/// the cosine of the angle between va and F_v. Magnetic force: F_v ⊥ va, so
/// cos = 0. Drag or friction: F_v opposes va, cos = -1. Negative damping
/// (+b v): cos = +1. A sample with no usable cosine (F(0) not finite, as for
/// friction written -mu v/|v|; or F_v below the rounding noise, see
/// kCosineTol) is left out; if no sample is usable the label is Unknown.
template <class AdaptedForce>
[[nodiscard]] ProbeResult run_probe(const AdaptedForce& F, Shape shape, const CheckOptions& opt, Rng& rng) {
    ProbeResult out;
    if (shape == Shape::StaticField) return out;  // position-only by signature: nothing to probe
    const bool can_depend_on_v = (shape == Shape::Force);
    int cosines = 0;
    int perpendicular = 0;
    int opposing = 0;
    int aiding = 0;
    for (int i = 0; i < opt.dependence_samples; ++i) {
        // Draw every random number first, in a fixed order, so that the stream
        // (hence everything after the probe) does not depend on what F does.
        const Vec3 r = rng.in_box(opt.region);
        const double t1 = rng.uniform(opt.t_min, opt.t_max);
        const double t2 = rng.uniform(opt.t_min, opt.t_max);
        const Vec3 va = rng.in_cube(opt.v_max);
        const Vec3 vb = rng.in_cube(opt.v_max);
        ++out.samples;
        const Vec3 base = F(t1, r, va);
        if (!fields::detail::is_finite(base)) {
            ++out.skipped;
            continue;
        }
        const Vec3 later = F(t2, r, va);
        if (fields::detail::is_finite(later) && later != base) out.dep.on_time = true;
        if (!can_depend_on_v) continue;

        const Vec3 at_rest = F(t1, r, Vec3{});
        if (fields::detail::is_finite(at_rest)) {
            if (at_rest != base) out.dep.on_velocity = true;
            const Vec3 fv = base - at_rest;
            const double fv_size = norm(fv);
            const double noise = 2.0 * kEps * std::max(norm(base), norm(at_rest));
            if (norm(va) > 0.0 && fv_size > 0.0 && noise < kCosineTol * fv_size) {
                const double c = dot(va, fv) / (norm(va) * fv_size);
                ++cosines;
                if (std::fabs(c) <= kCosineTol) {
                    ++perpendicular;
                } else if (c < 0.0) {
                    ++opposing;
                } else {
                    ++aiding;
                }
            }
        }
        if (i % 2 == 0) {
            const Vec3 other = F(t1, r, vb);
            if (fields::detail::is_finite(other) && other != base) out.dep.on_velocity = true;
        }
    }
    if (out.dep.on_velocity) {
        if (cosines == 0) {
            out.power = VelocityPower::Unknown;  // dependence proven, but no usable angle
        } else if (perpendicular == cosines) {
            out.power = VelocityPower::NoWork;
        } else if (opposing == cosines) {
            out.power = VelocityPower::Dissipative;
        } else if (aiding == cosines) {
            out.power = VelocityPower::Feeding;
        } else {
            out.power = VelocityPower::Mixed;  // both signs genuinely occur
        }
    }
    return out;
}

}  // namespace detail

// --- the probe --------------------------------------------------------------------------

/// Does F change when only t, or only v, changes?
///
/// Physics: condition (i) of Taylor §4.2 requires F to depend on position
/// alone: not on velocity (drag §2.2-2.4, the magnetic force §2.5) and not
/// explicitly on time (a driving force §5.5, a time-dependent potential
/// §4.5). The probe holds r fixed and changes only t, then only v; any change
/// in F proves that dependence. No finite set of probes can prove
/// independence: a force that switches on at t = 100 s looks static on
/// [0, 10] s, and a drag that only acts above v_max is invisible.
///
/// "Change" means an exact difference between two finite values, and -0.0 equals
/// +0.0, so the question asked is what the force DOES, not what its signature
/// allows: LorentzForce with B = 0 computes q(E + v × 0), which is bitwise
/// independent of v, and so is reported as position-only.
///
/// C++ note: a callable that takes no v (or no t, or neither) is detected at
/// compile time and is not probed for what it cannot depend on.
///
/// @param F   force, any accepted shape (see the top of this file)
/// @param opt region [m], time window [s], v_max [m/s], seed, sample count
template <class Force>
[[nodiscard]] Dependence probe_dependence(const Force& F, const CheckOptions& opt = {}) {
    static_assert(detail::shape_v<Force> != detail::Shape::None, "see detail::kShapeMessage");
    detail::Rng rng(opt.seed);
    return detail::run_probe(detail::view_force(F), detail::shape_v<Force>, opt, rng).dep;
}

// --- the classifier ----------------------------------------------------------------------

/// Classifies a force: conservative, time-dependent potential, velocity-
/// dependent, nonzero curl, path-dependent, or inconclusive. Deterministic:
/// the same force and options always give the same report.
///
/// Physics: applies Taylor's two conditions (§4.2) in order. Condition (i),
/// F = F(r), comes first: velocity dependence rules out a potential energy of
/// any kind, while time dependence alone still allows a time-dependent U(r, t)
/// (§4.5). Then condition (ii), path independence, is tested locally by
/// curl F = 0 at sample points (§4.4) and globally by the work around closed
/// loops, which catches "holes" where zero curl is not enough (Stokes' theorem
/// needs a simply connected region). The verdict applies to the given box and
/// time window and is evidence, not proof.
///
/// The steps, in this order:
///  0. A degenerate region (flat, inverted or non-finite box), or invalid
///     options that disable a required source of evidence, is Inconclusive.
///  1. Probe v (and t). F depends on v: VelocityDependent, plus velocity_power
///     from the same probe data. Stop.
///  2. Curl at `curl_samples` points, deterministic and the same in every time
///     slice, with t frozen at each slice. For each point the Jacobian comes
///     from fields::detail::jacobian_with_error with h = 1e-4·max(|r - c|, L)
///     (c = box centre, L = half the largest side); a sample with
///     relative_curl > curl_tol is a *candidate*. A candidate counts only if
///     (a) the curl recomputed with h/2 agrees to 10%: a true curl does not
///     depend on h, while truncation artifacts do; and (b) a small Stokes loop
///     agrees: the circulation around a circle of radius 0.05·L centred on the
///     sample, normal along the curl, divided by πR², must match |curl| to
///     25%. (b) rejects coherent rounding noise, which passes (a) but does not
///     integrate to a circulation. The first confirmed candidate gives
///     NonzeroCurl. Stop.
///  3. Closed loops. In each of 3 parallel planes per axis (the box centre
///     and centre ± half a side/2), grids of g×g rectangular cells with
///     g = 1, 2, 3, 5 tile the box's cross-section; each cell boundary is a
///     counterclockwise loop of 4 straight segments, integrated with
///     fields::line_integral_detailed. By Stokes' theorem ∮F·dr = ∫(curl F)·n dA, so
///     a loop sees a curl too weak or too concentrated for point samples, and
///     (the reason they are here) the delta-function curl of a hole. A loop that
///     converged and has |∮F·dr| / ∮|F||dr| > loop_tol gives PathDependent. A loop
///     any of whose segments did not converge (it touched a singularity) is
///     skipped and counted, never treated as evidence either way. All loops
///     are run, so max_loop_ratio is the true worst.
///     Cost: 351 loops (1404 segment integrals, ~15-1000 force evaluations
///     each) per time slice at the default loop_grid_levels = 4.
///  4. Otherwise: TimeDependentPotential if F depends on t, else Conservative.
///  Inconclusive instead, when no positive evidence was found and any of these
///  holds: the probe skipped every sample; more than half the curl samples
///  of a slice were skipped; every loop of a slice was skipped.
///
/// Why rectangles in several planes: a singular line (a "hole") is seen by a
/// loop only if it pierces the loop's disc. A straight line pierces the
/// cross-section of the box in some axis direction, but not necessarily near
/// the centre of the box. A fixed set of circles in the three centre planes
/// missed 12.7% of random oblique lines, and antiparallel wire pairs whose
/// circulations cancel in the one big loop. Cell loops cancel a wire pair only
/// if both wires sit in the same cell, which the finer grids and the shifted
/// planes avoid.
///
/// Known limits (documented, measured):
///  - Sampling cannot prove conservativeness. A line that crosses the box
///    only near a corner, or runs exactly along the cell edges of every grid
///    in every plane, can be missed.
///  - The loop ratio divides by ∮|F||dr|, so a huge uniform field added to a
///    circulating one dilutes it. Measured for the hole field plus a uniform
///    (F0, 0, 0): ratio 3.9e-6 at F0 = 1e6 N, 3.9e-8 at 1e8 N (still caught),
///    and 3.9e-9 at 1e9 N, which falls below loop_tol and is reported
///    Conservative.
///  - A force whose own evaluation is noisy at the 1e-9 relative level or
///    worse may be misjudged: the noise can survive the curl confirmation (it
///    is then usually rejected by the Stokes loop, but not always), and noisy
///    loops do not converge and are skipped. Differentiation cannot rescue it.
///  - Time dependence outside [t_min, t_max] and velocity dependence outside
///    about 1.4·v_max are invisible (see probe_dependence), and only
///    `time_slices` instants are tested: a curl that exists only between
///    them is missed (increase time_slices).
///  - Fields with a kink or singular point inside the box are differentiated
///    wrongly within 2h of it; such samples are caught by the error estimate,
///    the h/2 check and the Stokes check, and the loops are the global backstop.
///  - A box much smaller than 1e-10·|r| (far from the origin) cannot be
///    resolved in double precision.
///
/// C++ note: takes the force by const reference and never copies it; a
/// non-const-callable force (`mutable` lambda) is rejected at compile time.
///
/// @param F   force, any accepted shape
/// @param opt region [m], time window [s], tolerances, sample counts, seed
template <class Force>
[[nodiscard]] ConservativeReport check_conservative(const Force& F, const CheckOptions& opt = {}) {
    constexpr detail::Shape shape = detail::shape_v<Force>;
    static_assert(shape != detail::Shape::None, "see detail::kShapeMessage");
    const auto force = detail::view_force(F);
    detail::Rng rng(opt.seed);
    ConservativeReport rep;

    // ---- step 0: the region and configuration must support a real check --------------------
    const Vec3 lo = opt.region.lo;
    const Vec3 hi = opt.region.hi;
    if (!(fields::detail::is_finite(lo) && fields::detail::is_finite(hi) && hi.x > lo.x && hi.y > lo.y &&
          hi.z > lo.z)) {
        rep.verdict = Verdict::Inconclusive;
        rep.note = "the region is degenerate (flat, inverted or not finite), so no volume was tested";
        return rep;
    }
    const bool bad_window = !(std::isfinite(opt.t_min) && std::isfinite(opt.t_max) && opt.t_max >= opt.t_min);
    const bool bad_counts = opt.curl_samples <= 0 || opt.time_slices <= 0 ||
                            opt.loop_grid_levels < 1 || opt.loop_grid_levels > 4;
    const bool needs_dependence_probe = shape != detail::Shape::StaticField;
    const bool bad_probe = needs_dependence_probe && opt.dependence_samples <= 0;
    const bool bad_velocity_range = shape == detail::Shape::Force &&
                                    !(std::isfinite(opt.v_max) && opt.v_max > 0.0);
    const bool bad_tolerances = !(std::isfinite(opt.curl_tol) && opt.curl_tol > 0.0 &&
                                  std::isfinite(opt.loop_tol) && opt.loop_tol > 0.0);
    if (bad_window || bad_counts || bad_probe || bad_velocity_range || bad_tolerances) {
        rep.verdict = Verdict::Inconclusive;
        rep.note = "check options disable required evidence or contain an invalid range or tolerance";
        return rep;
    }

    // ---- step 1: probe v (and t) -------------------------------------------------------
    const detail::ProbeResult probe = detail::run_probe(force, shape, opt, rng);
    rep.depends_on_time = probe.dep.on_time;
    rep.depends_on_velocity = probe.dep.on_velocity;
    rep.skipped_samples = probe.skipped;
    if (probe.dep.on_velocity) {
        rep.verdict = Verdict::VelocityDependent;
        rep.velocity_power = probe.power;
        return rep;
    }
    if (probe.samples > 0 && probe.skipped == probe.samples) {
        rep.verdict = Verdict::Inconclusive;  // we never saw a finite force value
        rep.note = "F was not finite at any probe sample";
        return rep;
    }

    // The same sample points are used in every time slice.
    const double L = detail::box_length_scale(opt.region);
    const Vec3 mid = 0.5 * (lo + hi);
    const Vec3 half = 0.5 * (hi - lo);
    std::vector<Vec3> points;
    points.reserve(static_cast<std::size_t>(std::max(opt.curl_samples, 0)));
    for (int i = 0; i < opt.curl_samples; ++i) points.push_back(rng.in_box(opt.region));
    const std::vector<double> times = detail::slice_times(opt, rep.depends_on_time);

    // ---- step 2: curl at the samples -----------------------------------------------------
    bool too_many_curl_skips = false;
    for (const double t : times) {
        // Freeze t, and v = 0, to get a plain field; the callables below take it by reference.
        auto field = [&force, t](const Vec3& r) { return force(t, r, Vec3{}); };
        int skipped_in_slice = 0;
        for (const Vec3& r : points) {
            const double ell = std::max(norm(r - mid), L);  // distance from the BOX centre, not the origin
            const double h = 1e-4 * ell;
            const fields::detail::JacobianSample s = fields::detail::jacobian_with_error(field, r, h);
            if (!s.finite) {
                ++rep.skipped_samples;
                ++skipped_in_slice;
                continue;
            }
            ++rep.samples_tested;
            const double rel = fields::detail::relative_curl(s, ell);
            if (rel > opt.curl_tol) {
                // Candidate. Confirmation (a): recompute at h/2. Compare curls, not scores.
                const fields::detail::JacobianSample s2 = fields::detail::jacobian_with_error(field, r, 0.5 * h);
                const Vec3 c1 = fields::curl(s.J);
                const Vec3 c2 = fields::curl(s2.J);
                bool confirmed = s2.finite && norm(c2 - c1) <= 0.1 * norm(c1);
                if (confirmed) {
                    // Confirmation (b): Stokes. The circulation around a small circle centred here,
                    // normal along the curl, must be |curl|·πR². rel_tol 1e-3: a real circulation is
                    // large against its quadrature error; noise is not (and fails to converge).
                    const double R = 0.05 * L;
                    const fields::LineIntegral loop = fields::circulation_detailed(field, r, c1 / norm(c1), R, 1e-3, 1e-14);
                    const double pi_R2 = 3.14159265358979323846 * R * R;
                    confirmed = loop.converged && std::fabs(loop.value / pi_R2 - norm(c1)) <= 0.25 * norm(c1);
                }
                if (confirmed) {
                    rep.verdict = Verdict::NonzeroCurl;
                    rep.max_relative_curl = rel;  // exceeds every earlier (non-candidate) score
                    rep.worst_curl_point = r;
                    rep.worst_curl_time = t;
                    return rep;
                }
                ++rep.unconfirmed_curl_candidates;
            } else if (rel > rep.max_relative_curl) {
                rep.max_relative_curl = rel;
                rep.worst_curl_point = r;
                rep.worst_curl_time = t;
            }
        }
        if (2 * skipped_in_slice > opt.curl_samples) too_many_curl_skips = true;
    }

    // ---- step 3: closed loops --------------------------------------------------------------
    const int sizes_used = std::clamp(opt.loop_grid_levels, 0, 4);
    bool a_slice_lost_all_loops = false;
    bool path_dependent = false;
    for (const double t : times) {
        auto field = [&force, t](const Vec3& r) { return force(t, r, Vec3{}); };
        int tested_in_slice = 0;
        int skipped_in_slice = 0;
        for (int a = 0; a < 3; ++a) {  // planes with normal axis a; in-plane axes i, j, with (a, i, j) cyclic
            const int i = (a + 1) % 3;
            const int j = (a + 2) % 3;
            const double ext_i = detail::component(half, i);
            const double ext_j = detail::component(half, j);
            const double ext_a = detail::component(half, a);
            for (int plane = 0; plane < 3; ++plane) {
                const double offset = plane == 0 ? 0.0 : (plane == 1 ? -0.5 * ext_a : 0.5 * ext_a);
                for (int size = 0; size < sizes_used; ++size) {
                    const int g = detail::kLoopGridSizes[size];
                    for (int gi = 0; gi < g; ++gi) {
                        for (int gj = 0; gj < g; ++gj) {
                            const double i0 = detail::component(mid, i) - ext_i + 2.0 * ext_i * gi / g;
                            const double i1 = detail::component(mid, i) - ext_i + 2.0 * ext_i * (gi + 1) / g;
                            const double j0 = detail::component(mid, j) - ext_j + 2.0 * ext_j * gj / g;
                            const double j1 = detail::component(mid, j) - ext_j + 2.0 * ext_j * (gj + 1) / g;
                            const fields::LineIntegral I = detail::cell_circulation(field, mid, a, i, j, offset, i0, i1, j0, j1,
                                                                                    1e-2 * opt.loop_tol);
                            ++rep.loops_tested;
                            ++tested_in_slice;
                            if (!I.converged) {
                                ++rep.loops_skipped;
                                ++skipped_in_slice;
                                continue;
                            }
                            const double ratio = I.scale > 0.0 ? std::fabs(I.value) / I.scale : 0.0;
                            if (ratio > rep.max_loop_ratio) {
                                rep.max_loop_ratio = ratio;
                                Vec3 c = mid;
                                detail::component(c, a) += offset;
                                detail::component(c, i) = 0.5 * (i0 + i1);
                                detail::component(c, j) = 0.5 * (j0 + j1);
                                rep.worst_loop_center = c;
                                rep.worst_loop_time = t;
                            }
                            if (ratio > opt.loop_tol) path_dependent = true;
                        }
                    }
                }
            }
        }
        if (tested_in_slice > 0 && skipped_in_slice == tested_in_slice) a_slice_lost_all_loops = true;
    }

    // ---- steps 4 and 5: the verdict --------------------------------------------------------
    // Positive evidence (a confirmed curl, a converged loop with work) beats Inconclusive: a
    // force can be shown non-conservative even if parts of the box were not evaluable.
    if (path_dependent) {
        rep.verdict = Verdict::PathDependent;
    } else if (too_many_curl_skips || a_slice_lost_all_loops) {
        rep.verdict = Verdict::Inconclusive;
        rep.note = too_many_curl_skips ? "F was not finite at more than half the curl samples of a time slice"
                                       : "every loop integral of a time slice failed to converge (a singularity on each)";
    } else {
        rep.verdict = rep.depends_on_time ? Verdict::TimeDependentPotential : Verdict::Conservative;
    }
    return rep;
}

// --- splitting at rest ----------------------------------------------------------------------------

/// The part of F that acts on a particle at rest: (t, r, v) -> F(t, r, 0).
///
/// Physics: any force splits as F(t, r, v) = F(t, r, 0) + [F(t, r, v) - F(t, r, 0)]:
/// what it does to a particle at rest, plus what it does only because the
/// particle moves. If the rest part is curl-free it has a potential U (§4.4),
/// or U(r, t) (§4.5), and then d(T + U)/dt = F_v·v (+ ∂U/∂t), where F_v is
/// the velocity-dependent part. This is exact bookkeeping, NOT Taylor's
/// physical classification: he sorts forces by what they are (gravity and
/// springs conservative; friction and drag not).
///
/// Caveat, drag in a wind: F = -b(v - w) splits into the constant push
/// F(t, r, 0) = b w, booked as a "potential" -b w·r, and -b v, booked as
/// dissipation. The ledger balances, but U = -b w·r is not stored energy (the
/// wind supplies it), and -b|v|² is not the heat produced (that is
/// -b|v - w|²). For such a force pass the whole drag as F_nc to
/// analysis::energy_balance instead.
///
/// C++ note: a *factory returning a function object*: it returns a lambda
/// that captured (copied) the adapted force. The return type is `auto`
/// because a lambda's type has no name; callers use `auto` too. The lambda
/// has the 3-argument Force signature whatever shape F had, so it plugs into
/// simulate(), check_conservative(), find_potential() and so on.
///
/// @param F force, any accepted shape (copied)
template <class Force>
[[nodiscard]] auto velocity_independent_part(Force F) {
    static_assert(detail::shape_v<Force> != detail::Shape::None, "see detail::kShapeMessage");
    return [f = detail::as_force(std::move(F))](double t, const Vec3& r, const Vec3& /*v*/) -> Vec3 {
        return f(t, r, Vec3{});
    };
}

/// The part of F that exists only because the particle moves:
/// (t, r, v) -> F(t, r, v) - F(t, r, 0). See velocity_independent_part for the
/// physics and the wind-drag caveat.
///
/// Physics: for the magnetic force this is all of it, q v × B, and P = F_v·v = 0:
/// it does no work, so T + U_rest is conserved. For drag, F_v·v < 0 and energy
/// is dissipated. check_conservative() reports these as velocity_power.
///
/// Numerics: computed as a subtraction, so rounding in the result scales with
/// |F(t, r, 0)|, not with |F_v|: for a force that is mostly position-dependent
/// the "perpendicular" test of a magnetic part is good to ~1e-16·|F|/|F_v|.
///
/// @param F force, any accepted shape (copied)
template <class Force>
[[nodiscard]] auto velocity_dependent_part(Force F) {
    static_assert(detail::shape_v<Force> != detail::Shape::None, "see detail::kShapeMessage");
    return [f = detail::as_force(std::move(F))](double t, const Vec3& r, const Vec3& v) -> Vec3 {
        return f(t, r, v) - f(t, r, Vec3{});
    };
}

}  // namespace newton::potentials
