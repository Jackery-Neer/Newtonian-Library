#pragma once
// newton/analysis/energy.hpp — energy diagnostics for a simulated trajectory.
//
// Physics: mechanical energy E = T + U is conserved for any force derivable
// from a time-independent potential, F = -grad(U) (Taylor §4.1-4.4). No
// integrator step function in this library ever evaluates U — a numerical
// integrator only ever sees forces — so watching E(t) computed *after the
// fact* from a recorded Trajectory is the standard, integrator-agnostic way
// to grade how well a run conserved energy without touching the stepper's
// internals. This is the single most useful integrator diagnostic (research
// numerical-integration.md §4.3): a symplectic method (semi-implicit Euler,
// velocity Verlet) shows bounded oscillation around 0 drift with no trend;
// RK4 shows a small monotone drift; explicit Euler on an oscillator shows
// unbounded exponential growth. See energy_drift() below.
//
// Non-conservative forces (drag, drives) change E_mech, and that is not a
// failure: Delta(T + U) = W_nc, the work the non-conservative forces do
// (Taylor §4.1-4.2). power(), work_done() and energy_balance() below keep
// that ledger from the recorded samples, so a damped or driven run can be
// audited just as a conservative one is by energy_drift().
//
// C++ note: PotentialFn is a template callable (`double U(const Vec3& r)`
// or, for a time-dependent potential, `double U(double t, const Vec3& r)`),
// the same "template, not std::function" contract the Force concept uses
// (DESIGN.md). These functions typically run once per sample in a
// post-processing pass rather than inside a hot integration loop, but
// matching the rest of the library's calling convention costs nothing and
// keeps the two contracts easy to remember together.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include <newton/core/particle.hpp>
#include <newton/core/trajectory.hpp>
#include <newton/core/vec3.hpp>

namespace newton::analysis {

/// Kinetic energy T = 1/2 m |v|^2 (Taylor §4.1).
/// @param mass particle mass [kg]
/// @param v    velocity [m/s]
[[nodiscard]] inline double kinetic_energy(double mass, const Vec3& v) {
    return 0.5 * mass * norm2(v);
}

namespace detail {

/// True when U can be called as U(t, r): a time-dependent potential U(r, t)
/// (Taylor §4.5).
///
/// C++ note: std::is_invocable_v asks the compiler "would this call compile?"
/// without making it. Paired with `if constexpr` (below) the unused branch is
/// not even instantiated, so a plain U(r) lambda never sees a 2-argument call.
template <class PotentialFn>
inline constexpr bool is_time_dependent_potential_v =
    std::is_invocable_v<PotentialFn&, double, const Vec3&>;

/// U evaluated at one instant, preferring U(t, r) when U accepts it.
template <class PotentialFn>
[[nodiscard]] double potential_at(PotentialFn& U, double t, const Vec3& r) {
    if constexpr (is_time_dependent_potential_v<PotentialFn>) {
        return U(t, r);
    } else {
        (void)t;  // a static potential ignores the time
        return U(r);
    }
}

/// A time step for a centred finite-difference stencil. Start with a small
/// fraction of the locally recorded time gap, then round *up* to a multiple
/// of the local floating-point spacing at t. This keeps t +/- d and
/// t +/- 2d distinct representable times without letting an absolute
/// timestamp turn a millisecond-scale stencil into a huge one.
[[nodiscard]] inline double time_derivative_step(double t, double local_gap) {
    const double desired = 1e-2 * local_gap;
    const double forward_ulp = std::nextafter(t, std::numeric_limits<double>::infinity()) - t;
    const double backward_ulp = t - std::nextafter(t, -std::numeric_limits<double>::infinity());
    const double ulp = std::max(forward_ulp, backward_ulp);

    // A finite trajectory timestamp has a positive neighbouring spacing.
    // Keep the gap-based choice as a harmless fallback for malformed input
    // rather than manufacturing a NaN derivative here.
    if (!(ulp > 0.0) || !std::isfinite(ulp)) return desired;

    // Snapping to the coarser neighbouring ULP also keeps the negative
    // offsets exact at a power-of-two binade boundary, where the spacing
    // immediately below t is half the spacing immediately above it.
    const double multiples = std::ceil(std::max(desired, ulp) / ulp);
    if (!std::isfinite(multiples)) return desired;
    return multiples * ulp;
}

/// Cumulative integral of samples f_i taken at times t_i (strictly
/// increasing, but NOT necessarily evenly spaced): out[0] = 0 and
/// out[i] = integral of f from t_0 to t_i. Used for work (f = power) and for
/// the dU/dt term of the energy ledger.
///
/// Method (research: composite Simpson on pairs of intervals, non-uniform):
///   - Even i: out[i] = out[i-2] + Simpson over samples (i-2, i-1, i):
///       H/6 [ (2 - h2/h1) f0 + H^2/(h1 h2) f1 + (2 - h1/h2) f2 ],
///       h1 = t1 - t0, h2 = t2 - t1, H = h1 + h2   (equal gaps: h/3 [1 4 1]).
///   - Odd i >= 3: out[i] = out[i-1] + the last interval integrated under the
///     parabola through samples (i-2, i-1, i): trapezoid minus a curvature
///     term, h2 (f1 + f2)/2 - C h2^3 / 6, C = the parabola's second divided
///     difference.
///   - i == 1: the same, using the parabola through samples 0, 1, 2.
///   - Two samples only: the trapezoid. One: out = {0}.
/// Every piece takes arbitrary gaps, so a short final interval (from
/// simulate_until) or a record_every that does not divide n is fine.
/// Accuracy O(h^4) like Simpson's rule. A zero or negative gap (a duplicated
/// time stamp) cannot be fitted by a parabola, so that piece falls back to
/// the trapezoid.
[[nodiscard]] inline std::vector<double> cumulative_integral(const std::vector<double>& t,
                                                             const std::vector<double>& f) {
    const std::size_t n = t.size();
    std::vector<double> out(n, 0.0);
    if (n < 2) return out;

    auto trapezoid = [&](std::size_t a, std::size_t b) {
        return 0.5 * (t[b] - t[a]) * (f[a] + f[b]);
    };
    if (n == 2) {
        out[1] = trapezoid(0, 1);
        return out;
    }

    // Integral over the single interval [t_b, t_c] of the parabola through
    // samples a, b, c (in time order a < b < c) taken as the *last* interval
    // (b -> c), or as the *first* interval (a -> b) when `first` is set.
    // Both are "trapezoid + C * integral of (x - t1)(x - t2)" = trapezoid -
    // C h^3 / 6, where (t1, t2) is the interval being integrated.
    auto one_interval = [&](std::size_t a, std::size_t b, std::size_t c, bool first) {
        const double h1 = t[b] - t[a];
        const double h2 = t[c] - t[b];
        const std::size_t lo = first ? a : b;
        const std::size_t hi = first ? b : c;
        if (!(h1 > 0.0) || !(h2 > 0.0)) return trapezoid(lo, hi);
        const double C = (h1 * (f[c] - f[b]) - h2 * (f[b] - f[a])) / (h1 * h2 * (h1 + h2));
        const double h = t[hi] - t[lo];
        return trapezoid(lo, hi) - C * h * h * h / 6.0;
    };

    for (std::size_t i = 1; i < n; ++i) {
        if (i == 1) {
            out[1] = one_interval(0, 1, 2, /*first=*/true);
        } else if (i % 2 == 0) {
            const double h1 = t[i - 1] - t[i - 2];
            const double h2 = t[i] - t[i - 1];
            if (!(h1 > 0.0) || !(h2 > 0.0)) {
                out[i] = out[i - 1] + trapezoid(i - 1, i);
            } else if (h2 < 1e-3 * h1) {
                // Lopsided pair (simulate_until lands an event just past a
                // step boundary): the Simpson weights H^2/(h1 h2) and
                // (2 - h1/h2) are ~h1/h2 and cancel catastrophically, losing
                // ~eps*h1/h2 of accuracy. Integrate the tiny last interval
                // from the parabola through the last three samples instead,
                // which is well conditioned for any h2.
                out[i] = out[i - 1] + one_interval(i - 2, i - 1, i, /*first=*/false);
            } else {
                const double H = h1 + h2;
                out[i] = out[i - 2] + H / 6.0 *
                    ((2.0 - h2 / h1) * f[i - 2] + H * H / (h1 * h2) * f[i - 1] +
                     (2.0 - h1 / h2) * f[i]);
            }
        } else {
            out[i] = out[i - 1] + one_interval(i - 2, i - 1, i, /*first=*/false);
        }
    }
    return out;
}

}  // namespace detail

/// Total mechanical energy E = T + U at one instant.
///
/// Physics: with a time-dependent potential U(r, t) (Taylor §4.5) E is no
/// longer constant even when nothing dissipates: dE/dt = dU/dt|_r, the
/// power supplied by whoever moves the potential. Also, E(t) then depends on
/// the additive "gauge" function c(t) in U; only differences made with the
/// identity E(t) - E(0) = integral of dU/dt are gauge-free.
///
/// C++ note: U is called as U(s.t, s.r) when it accepts (double, Vec3),
/// otherwise as U(s.r) (`if constexpr` + std::is_invocable_v), so every
/// existing static-potential call keeps compiling unchanged.
/// @param mass particle mass [kg]
/// @param s    state (r, v, t) to evaluate
/// @param U    potential energy: double U(const Vec3& r) [J] or
///             double U(double t, const Vec3& r) [J], t in [s], r in [m]
template <class PotentialFn>
[[nodiscard]] double total_energy(double mass, const State& s, PotentialFn&& U) {
    return kinetic_energy(mass, s.v) + detail::potential_at(U, s.t, s.r);
}

/// Relative energy drift over a trajectory: (E(t) - E(0)) / |E(0)|, one entry
/// per sample in traj.states, same length and order — so drift[i] lines up
/// with traj.states[i] for plotting drift(t).
///
/// Physics: for a static, conservative U this is pure integrator error. If U
/// is time-dependent, U(t, r), the drift ALSO contains real physics, the
/// integral of dU/dt (Taylor §4.5), and non-conservative forces such as drag
/// change E for real (§4.2). To separate physics from numerics use
/// energy_balance(), which subtracts those terms.
/// @param traj trajectory to diagnose; traj.mass is used as the particle mass
/// @param U    potential energy: double U(const Vec3& r) [J] or
///             double U(double t, const Vec3& r) [J]
template <class PotentialFn>
[[nodiscard]] std::vector<double> energy_drift(const Trajectory& traj, PotentialFn&& U) {
    std::vector<double> drift;
    drift.reserve(traj.size());
    if (traj.empty()) return drift;

    const double E0 = total_energy(traj.mass, traj.states.front(), U);
    // Guard the (legitimate) E0 == 0 case — e.g. a particle released from
    // rest where U was zeroed at the release point — so a genuinely-zero
    // reference energy can't turn every later sample's drift into +-inf/NaN
    // by dividing by zero. Mirrors the zero-guard pattern in
    // Vec3::normalized() (core/vec3.hpp).
    const double scale = (std::fabs(E0) > 0.0) ? std::fabs(E0) : 1.0;

    for (const State& s : traj.states) {
        drift.push_back((total_energy(traj.mass, s, U) - E0) / scale);
    }
    return drift;
}

/// Instantaneous power P = F . v delivered by a force to the particle [W]
/// (Taylor §4.1).
///
/// Physics: dT/dt = F_net . v, and the same dot product for one force alone is
/// the rate at which that force does work. A magnetic force q v x B is
/// perpendicular to v, so P = 0 exactly; drag F = -b v gives P = -b v^2 <= 0.
/// @param F force: Vec3 F(double t, const Vec3& r, const Vec3& v) [N]
/// @param s state (t [s], r [m], v [m/s]) at which to evaluate it
template <class Force>
[[nodiscard]] double power(Force&& F, const State& s) {
    return dot(F(s.t, s.r, s.v), s.v);
}

/// Cumulative work W[i] = integral from t_0 to t_i of F . v dt, one entry per
/// sample, aligned with traj.states (W[0] = 0) [J].
///
/// Physics: the work a force does along the actual motion is
/// W(t) = integral of F . dr = integral of F(t', r, v) . v dt' (Taylor §4.1),
/// the time integral of the power P = F . v. Unlike fields::line_integral,
/// this follows the recorded trajectory in time, so it works for every force,
/// including velocity-dependent ones (drag, magnetic) and time-dependent ones
/// (drives). For the net force it reproduces the work-KE theorem, W = delta T.
/// The magnetic force always gives W = 0 (it is perpendicular to v) and drag
/// always gives W <= 0.
///
/// C++ note: power is evaluated exactly at each sample, then integrated with
/// a non-uniform composite Simpson rule (see detail::cumulative_integral), so
/// it handles unevenly spaced samples: simulate_until's short last interval,
/// or a record_every that does not divide the step count. Accuracy is
/// O((record_every * dt)^4) plus the trajectory's own integrator error: for
/// bookkeeping record every step (record_every = 1) and resolve about 20 or
/// more samples per oscillation period. With 0 samples the result is empty,
/// with 1 it is {0}, with 2 it is the trapezoid.
/// @param traj trajectory (positions [m], velocities [m/s], times [s])
/// @param F    force: Vec3 F(double t, const Vec3& r, const Vec3& v) [N]
template <class Force>
[[nodiscard]] std::vector<double> work_done(const Trajectory& traj, Force&& F) {
    std::vector<double> t;
    std::vector<double> P;
    t.reserve(traj.size());
    P.reserve(traj.size());
    for (const State& s : traj.states) {
        t.push_back(s.t);
        P.push_back(power(F, s));
    }
    return detail::cumulative_integral(t, P);
}

/// Residual of the work-energy ledger at every sample, same length and order
/// as traj.states:
///   [ (T+U)_i - (T+U)_0 - W_nc[i] - integral of dU/dt ] / scale
/// where W_nc comes from work_done(traj, F_nc).
///
/// Physics: the work-energy theorem with non-conservative forces: if the force
/// you simulated is F = -grad U + F_nc, then delta(T + U) = W_nc, where W_nc
/// is the work done by F_nc (Taylor §4.1-4.2; add the integral of dU/dt if U
/// depends on t, §4.5). This returns the residual of that ledger at every
/// sample. It should stay at integrator plus quadrature error for any run
/// (damped, dragged or driven), because it balances energy rather than
/// assuming it is conserved. A residual that grows means U and F_nc do not
/// add up to the force you actually simulated (a forgotten or double-counted
/// force), or that dt is too large. Requirement: -grad U + F_nc must equal the
/// force that was simulated.
///
/// C++ note: for a time-dependent U(t, r) the dU/dt term is estimated at each
/// sample with a 5-point central difference in t. Its step is 1% of the
/// smaller adjacent sample gap, rounded up only as needed to make the
/// neighbouring timestamps representable in `double`; it therefore follows
/// the time scale the samples resolve, not the run length or the absolute
/// start time. It is integrated with the same cumulative Simpson rule as the
/// work. For a static U that term is skipped entirely. The scale is
/// max over samples of (T_i, |U_i - U_0|, |W_nc[i]|, |integral dU/dt|_i), or
/// 1 J if all vanish: the largest energy that crosses the ledger. It is
/// gauge-invariant (only differences of U enter), and it contains neither
/// |E_0| nor |E_i - E_0|, so a correct run with E_0 = 0 reads ~0 while a wrong
/// U or a forgotten force reads ~1 at any unit scale.
///
/// Sign convention: W_nc < 0 means F_nc removes energy (drag dissipates,
/// -b v^2 for linear drag, -c|v|^3 for quadratic drag). The residual is
/// (energy change) - (booked work), so a NEGATIVE residual means energy was
/// lost that F_nc does not book (e.g. a forgotten drag), a positive one that
/// energy appeared from nowhere (e.g. a forgotten driving force).
///
/// Accuracy: the residual contains the quadrature error, O((record_every *
/// dt)^4), PLUS the integrator's own energy error: RK4 O(dt^4), velocity
/// Verlet O(dt^2), (semi-implicit) Euler O(dt). Record every step.
/// @param traj trajectory; traj.mass is the particle mass [kg]
/// @param U    potential energy: double U(const Vec3& r) or
///             double U(double t, const Vec3& r) [J]
/// @param F_nc the non-conservative part: Vec3 F(double t, const Vec3& r,
///             const Vec3& v) [N]; pass a zero force when there is none
template <class PotentialFn, class Force>
[[nodiscard]] std::vector<double> energy_balance(const Trajectory& traj, PotentialFn&& U,
                                                 Force&& F_nc) {
    std::vector<double> residual;
    if (traj.empty()) return residual;
    const std::size_t n = traj.size();
    residual.reserve(n);

    const std::vector<double> W = work_done(traj, F_nc);

    // Cumulative integral of dU/dt along the path (zeros for a static U).
    std::vector<double> dUdt_int(n, 0.0);
    if constexpr (detail::is_time_dependent_potential_v<PotentialFn>) {
        std::vector<double> t(n);
        std::vector<double> dUdt(n);
        for (std::size_t i = 0; i < n; ++i) t[i] = traj.states[i].t;
        const double span = t.back() - t.front();
        const double t_span = (span > 0.0) ? span : 1.0;
        for (std::size_t i = 0; i < n; ++i) {
            const Vec3& r = traj.states[i].r;
            // Step follows the time scale the samples resolve (the doc asks
            // for >= 20 samples per period of U). `time_derivative_step`
            // only enlarges it enough to survive floating-point timestamp
            // resolution; tying it to |t| would make the stencil's
            // (omega delta)^4 error grow with the absolute start time.
            double gap = t_span;
            if (i > 0 && t[i] - t[i - 1] > 0.0) gap = std::min(gap, t[i] - t[i - 1]);
            if (i + 1 < n && t[i + 1] - t[i] > 0.0) gap = std::min(gap, t[i + 1] - t[i]);
            const double d = detail::time_derivative_step(t[i], gap);
            dUdt[i] = (-U(t[i] + 2.0 * d, r) + 8.0 * U(t[i] + d, r) - 8.0 * U(t[i] - d, r) +
                       U(t[i] - 2.0 * d, r)) / (12.0 * d);
        }
        dUdt_int = detail::cumulative_integral(t, dUdt);
    }

    const double E0 = total_energy(traj.mass, traj.states.front(), U);
    const double U0 = detail::potential_at(U, traj.states.front().t, traj.states.front().r);

    // Scale = the largest energy that actually crosses the ledger:
    //   max_i of ( T_i, |U_i - U_0|, |W_nc[i]|, |integral dU/dt|_i ).
    // It is built only from kinetic energy and from DIFFERENCES of U, so it is
    // gauge-invariant (adding a constant to U changes neither the residual nor
    // the scale). It deliberately does NOT contain |E_0| or |E_i - E_0|: when
    // nothing else crosses the ledger, E_i - E_0 IS the residual, and putting
    // it in the scale would cancel a real error to "1" (a correct free fall
    // would read as total failure) while a large constant in U would swamp a
    // forgotten drag.
    double scale = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const State& s = traj.states[i];
        scale = std::max({scale, kinetic_energy(traj.mass, s.v),
                          std::fabs(detail::potential_at(U, s.t, s.r) - U0), std::fabs(W[i]),
                          std::fabs(dUdt_int[i])});
    }
    if (!(scale > 0.0)) scale = 1.0;  // everything zero: 1 J, as in energy_drift

    for (std::size_t i = 0; i < n; ++i) {
        const double dE = total_energy(traj.mass, traj.states[i], U) - E0;
        residual.push_back((dE - W[i] - dUdt_int[i]) / scale);
    }
    return residual;
}

/// U(r) = m g r.y: potential energy for newton::forces::UniformGravity's
/// default convention (g points in -y — this library's 2D "y is up" rule,
/// see core/constants.hpp). Pass g as a positive magnitude (e.g.
/// constants::g_earth), not a signed vector component: F = -dU/dr then comes
/// out to (0, -m g, 0), matching that force.
/// @param mass particle mass [kg]
/// @param g    gravitational acceleration magnitude [m/s^2]
[[nodiscard]] inline auto uniform_gravity_U(double mass, double g) {
    return [mass, g](const Vec3& r) { return mass * g * r.y; };
}

/// U(r) = 1/2 k |r|^2: potential energy for an isotropic spring anchored at
/// the origin (newton::forces::Spring with anchor = {} and
/// natural_length = 0 — the reduced form F = -k r, Taylor §5.1).
/// @param k spring constant [N/m]
[[nodiscard]] inline auto spring_U(double k) {
    return [k](const Vec3& r) { return 0.5 * k * norm2(r); };
}

}  // namespace newton::analysis
