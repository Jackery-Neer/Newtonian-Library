#pragma once
// newton/analytic/drag.hpp — closed-form drag solutions (Taylor §2.2-2.4).
//
// Physics: linear drag (F = -b v) is a *linear* ODE, so the full vector
// motion — including the coupling to gravity — has an exact solution; it is
// the single richest 2D validation oracle in the book (research map §8, row
// 6). Quadratic drag (F = -c|v|v) is only integrable in the two decoupled
// 1D cases below (pure horizontal coast, pure vertical fall from rest) —
// the general 2D quadratic-drag trajectory has NO closed form (that's the
// whole reason this library's numeric integrators exist; see
// forces/drag.hpp and docs/research/taylor-ch1-7-map.md §2).
//
// Numerics: the textbook formulas below are written with tau = m/b, which is
// HUGE when the drag is small, multiplied by (1 - e^(-t/tau)), which is TINY.
// Huge x tiny is the classic recipe for catastrophic cancellation. Every
// function here is therefore rewritten in terms of small dimensionless
// ratios (u = b t/m, eps = b x/(m vx0), z = g t/v_ter) and the Taylor-series
// functions of them (phi1, phi2, psi, ...) in `detail::` — the very same
// small-drag expansions Taylor §2.4 uses to do the physics. As a result
// every function is exact (not merely close) at b = 0 / c = 0 and nothing
// divides by b.
//
// C++ note: every function here calls std::exp/log/tanh/cosh, none of which
// are constexpr in C++17, so nothing in this file is constexpr (matches
// core/vec3.hpp's norm()/normalized() precedent). LinearDragParams::tau()
// and v_ter() are pure arithmetic and ARE constexpr.

#include <cmath>
#include <limits>

#include <newton/core/constants.hpp>
#include <newton/core/vec3.hpp>

namespace newton::analytic {

namespace detail {

// Physics + C++ note (the cancellation, once, for the whole file): take
// m = 1 kg, a ball launched at 30 m/s at 45 degrees, t = 3 s, so the vacuum
// height is y = 19.5097 m. The textbook y(t) = (vy0 + v_ter) tau (1 - e^(-t/tau))
// - v_ter t at b = 1e-10 kg/s has tau = 1e10 s and 1 - e^(-3e-10) = 3e-10:
// the computer forms 1 - (1 - 3e-10) and the subtraction keeps only the
// ~6 digits double stores in 1 - 3e-10, then multiplying by tau = 1e10
// blows that error up: the answer comes out as 24 KILOMETRES. At b = 1e-14
// it is -2.35e12 m, and at b = 0 it is 0 x infinity = NaN.
//
// The cure is two-fold. (1) std::expm1(x) and std::log1p(x) are library
// functions that compute e^x - 1 and ln(1 + x) WITHOUT the cancellation: they
// stay accurate when x is tiny, where forming exp(x) - 1 or log(1 + x)
// yourself would throw away every digit. (2) Where even those leave a
// subtraction of nearly equal numbers (the *second* difference in phi2 and
// psi below) we use Taylor series, which has no subtraction to lose digits in.
// The series are not a numerical trick bolted onto the physics: they ARE the
// small-drag expansion of Taylor §2.4 (expand ln(1 - eps) and ln(1 + eps)
// for small eps = "drag is a small correction"). Same math, two jobs.

/// phi1(u) = (1 - e^(-u))/u, the "how much of a full time step survives drag"
/// factor; phi1(0) = 1 (no drag: the full t), phi1(u) -> 1/u for large u
/// (the particle stops after a distance v0 tau).
///
/// Physics: with u = t/tau, the drag-modified displacement is v0 t phi1(u)
/// (Taylor §2.2: x(t) = vx0 tau (1 - e^(-t/tau)) = vx0 t phi1(t/tau)).
/// C++ note: -expm1(-u) is accurate for tiny u (no 1 - e^-u cancellation), so
/// only u == 0 (0/0) needs a special case.
/// @param u dimensionless time b t / m [-], >= 0
[[nodiscard]] inline double phi1(double u) {
    if (u == 0.0) return 1.0;
    return -std::expm1(-u) / u;
}

/// phi2(u) = (u - 1 + e^(-u))/u^2 = sum_{n>=0} (-u)^n/(n+2)!
///         = 1/2 - u/6 + u^2/24 - u^3/120 + ...; phi2(0) = 1/2.
///
/// Physics: the gravity term of linear-drag motion is a t^2 phi2(u) — for
/// u -> 0 this is the familiar (1/2) a t^2 of vacuum free fall, and the -u/6
/// term is the first drag correction (Taylor's expansion of e^(-t/tau)).
/// C++ note: for u < 1 we sum the series (terms shrink by <= u/(n+2) each, so
/// 20 terms leave a truncation error ~1/22! ~ 1e-21, far below roundoff, and
/// being alternating-but-rapidly-decreasing it loses no digits). The direct
/// form (u - 1 + e^-u)/u^2 subtracts nearly equal numbers for small u and
/// has relative error ~2 eps_machine/u: ~2e-16 at u = 1, which is why that
/// is the crossover. Above 1 we use the equivalent (1 - phi1(u))/u, which has
/// no cancellation there (phi1 <= 0.64) and cannot overflow u*u.
/// @param u dimensionless time b t / m [-], >= 0
[[nodiscard]] inline double phi2(double u) {
    if (u < 1.0) {
        double term = 0.5;  // n = 0: 1/2!
        double sum = term;
        for (int n = 1; n < 20; ++n) {
            term *= -u / static_cast<double>(n + 2);
            sum += term;
        }
        return sum;
    }
    return (1.0 - phi1(u)) / u;
}

/// psi(eps) = -(eps + ln(1 - eps))/eps^2 = sum_{n>=2} eps^(n-2)/n
///          = 1/2 + eps/3 + eps^2/4 + eps^3/5 + ...; psi(0) = 1/2.
///
/// Physics: THE small-drag function of linear drag (Taylor §2.4). The exact
/// trajectory y(x) has the term v_ter tau ln(1 - eps), eps = x/(vx0 tau);
/// expanding ln(1 - eps) = -eps - eps^2/2 - eps^3/3 - ... turns y(x) into
/// the vacuum parabola plus corrections, y = (vy0/vx0) x - (g x^2/vx0^2) psi(eps).
/// The leading 1/2 is the parabola, eps/3 the first drag correction. Only
/// defined for eps < 1 (x short of the vertical asymptote x = vx0 tau).
/// C++ note: for |eps| < 1/2 we evaluate the series by Horner's rule with 56
/// terms (truncation ~ 2^-56/58 ~ 2e-19). The direct form has error
/// ~2 eps_machine/eps ~ 4e-16 at eps = 1/2 and rapidly worse below; the
/// series has none. Returns NaN for eps >= 1.
/// @param eps dimensionless x b/(m vx0) [-], must be < 1
[[nodiscard]] inline double psi(double eps) {
    if (!(eps < 1.0)) return std::numeric_limits<double>::quiet_NaN();
    if (std::fabs(eps) < 0.5) {
        double sum = 0.0;
        for (int n = 55; n >= 0; --n) sum = sum * eps + 1.0 / static_cast<double>(n + 2);
        return sum;
    }
    return -(eps + std::log1p(-eps)) / (eps * eps);
}

/// ln cosh(z) without overflow or cancellation (even in z).
///
/// C++ note: std::cosh(z) overflows to infinity for |z| > ~710, so a plain
/// log(cosh(z)) is infinite there although ln cosh(z) is only ~|z|. Two
/// identities fix it: cosh z = 1 + 2 sinh^2(z/2) feeds log1p, which stays
/// accurate for small z (a plain log(cosh z) would lose every digit when
/// cosh z ~ 1), and for |z| > 10, ln cosh z = |z| - ln 2 + log1p(e^(-2|z|)),
/// which involves only e^(-2|z|) <= 2e-9 and cannot overflow.
/// @param z argument [-]
[[nodiscard]] inline double ln_cosh(double z) {
    const double a = std::fabs(z);
    if (a > 10.0) {
        constexpr double ln2 = 0.693147180559945309417232121458;
        return a - ln2 + std::log1p(std::exp(-2.0 * a));
    }
    const double sh = std::sinh(0.5 * a);
    return std::log1p(2.0 * sh * sh);
}

/// lambda(z) = 2 ln cosh(z)/z^2 = 1 - z^2/6 + 2 z^4/45 - ...; lambda(0) = 1.
///
/// Physics: the quadratic-drag fall distance is (1/2) g t^2 lambda(g t/v_ter)
/// (Taylor §2.4). Small z is "terminal speed far away": lambda -> 1 recovers
/// the vacuum (1/2) g t^2, and -z^2/6 is the first drag correction. Large z
/// is "terminal speed reached": lambda ~ 2(z - ln 2)/z^2, so the distance
/// becomes v_ter t - (v_ter^2/g) ln 2.
/// C++ note: below |z| = 1e-3 the 3-term series is used (next term
/// 17 z^6/1260 ~ 1e-20), which also covers z = 0 and avoids 0/0 when z^2
/// underflows.
/// @param z dimensionless g t/v_ter [-] (0 when v_ter is infinite)
[[nodiscard]] inline double lambda(double z) {
    const double a = std::fabs(z);
    if (a < 1e-3) {
        const double z2 = z * z;
        return 1.0 - z2 / 6.0 + 2.0 * z2 * z2 / 45.0;
    }
    return 2.0 * (ln_cosh(a) / a) / a;  // two divisions: a*a could overflow
}

/// tanh(z)/z = 1 - z^2/3 + 2 z^4/15 - 17 z^6/315 + ...; value 1 at z = 0.
///
/// Physics: the quadratic-drag fall speed is g t tanh(z)/z with z = g t/v_ter;
/// the 1 is free fall v = g t, and -z^2/3 is the first drag correction.
/// C++ note: below |z| = 1e-3 the 3-term series is used (next term
/// ~5e-20), which also makes z = 0 exact instead of 0/0.
/// @param z dimensionless g t/v_ter [-]
[[nodiscard]] inline double tanh_over_z(double z) {
    if (std::fabs(z) < 1e-3) {
        const double z2 = z * z;
        return 1.0 - z2 / 3.0 + 2.0 * z2 * z2 / 15.0;
    }
    return std::tanh(z) / z;
}

/// ell(eps) = ln(1 + eps)/eps = 1 - eps/2 + eps^2/3 - eps^3/4 + ...;
/// ell(0) = 1. Defined for eps > -1.
///
/// Physics: quadratic-drag horizontal distance x = v0 t ell(c v0 t/m)
/// (Taylor §2.4): 1 is coasting with no drag, -eps/2 the first slowdown
/// correction — Taylor's expansion of ln(1 + eps).
/// C++ note: for |eps| < 0.01 we sum 10 terms by Horner (truncation
/// ~ 1e-20/11); log1p is accurate beyond that.
/// @param eps dimensionless c v0 t/m [-]
[[nodiscard]] inline double ln1p_over(double eps) {
    if (std::fabs(eps) < 0.01) {
        double sum = 0.0;
        for (int n = 9; n >= 0; --n) {
            const double coef = ((n % 2 == 0) ? 1.0 : -1.0) / static_cast<double>(n + 1);
            sum = sum * eps + coef;
        }
        return sum;
    }
    return std::log1p(eps) / eps;
}

}  // namespace detail

/// Parameters for the linear-drag closed forms below: m ẍ = -m g ŷ - b v,
/// i.e. gravity along -y at magnitude g (matching forces::UniformGravity's
/// default 2D convention, "y is up") plus linear drag -b v (Taylor §2.2).
///
/// Physics: b = 0 (the default) is the vacuum limit, and every function in
/// this file handles it exactly.
struct LinearDragParams {
    double m{1.0};                 ///< mass [kg]
    double b{0.0};                 ///< linear drag coefficient [kg/s], >= 0; 0 = no drag
    double g{constants::g_earth};  ///< gravitational acceleration [m/s^2]

    /// Time constant τ = m/b [s] — the e-folding time for velocity
    /// relaxation toward the drift solution (Taylor §2.2). Returns +∞ at
    /// b = 0 (IEEE division of a positive number by zero: "no drag, never
    /// relaxes"). It is for display and for asymptotes like x_max = vx0 τ;
    /// none of the functions in this file multiply it by something that is
    /// zero, they use the ratios b t/m instead.
    [[nodiscard]] constexpr double tau() const { return m / b; }

    /// Terminal (vertical, free-fall) speed v_ter = m g / b [m/s]
    /// (Taylor §2.2) — the speed at which drag exactly balances gravity.
    /// Returns +∞ at b = 0 (vacuum has no terminal speed); see tau() for
    /// why the closed forms never multiply it by zero.
    [[nodiscard]] constexpr double v_ter() const { return m * g / b; }
};

/// Velocity at time t under gravity + linear drag, from initial velocity v0
/// (Taylor §2.2-2.3). Vector form of both 1D results at once: writing
/// a = (0, -g, 0) (matching LinearDragParams's -y gravity convention),
/// the ODE dv/dt = a - v/τ integrates to
///   v(t) = v0 e^(-u) + a t φ1(u),   u = b t/m = t/τ,  φ1(u) = (1 - e^(-u))/u,
/// which is the textbook v0 e^(-t/τ) + a τ (1 - e^(-t/τ)) rewritten so that
/// τ·(1 - e^(-u)) = t·φ1(u) never multiplies a huge number by a tiny one.
/// The x-component reduces to Taylor's vx(t) = vx0 e^(-t/τ); the y-component
/// to vy(t) = -v_ter + (vy0+v_ter) e^(-t/τ), using v_ter = g τ.
///
/// Physics: φ1(0) = 1, so at b = 0 this is exactly the vacuum v0 + a t.
/// The small-drag expansion e^(-u) = 1 - u + ..., φ1 = 1 - u/2 + ... is the
/// Taylor-series idea of §2.4 and is what makes the numerics safe.
/// C++ note: see detail:: above for the full story of why τ(1 - e^(-t/τ))
/// is a numerical trap (24 km instead of 19.5 m at b = 1e-10) and how
/// std::expm1 (e^x - 1 without cancellation) avoids it.
/// @param p  drag/gravity parameters
/// @param v0 initial velocity [m/s]
/// @param t  elapsed time [s]
[[nodiscard]] inline Vec3 linear_drag_velocity(const LinearDragParams& p, Vec3 v0, double t) {
    const double u = p.b * t / p.m;
    const Vec3 a_grav{0.0, -p.g, 0.0};
    return v0 * std::exp(-u) + a_grav * (t * detail::phi1(u));
}

/// Position at time t under gravity + linear drag, from initial position r0
/// and velocity v0 (Taylor §2.2-2.3). Obtained by integrating
/// linear_drag_velocity() in closed form:
///   r(t) = r0 + v0 t φ1(u) + a t² φ2(u),   u = b t/m,
///   φ1(u) = (1 - e^(-u))/u,   φ2(u) = (u - 1 + e^(-u))/u²,
/// the textbook r0 + (v0 - a τ) τ (1 - e^(-t/τ)) + a τ t regrouped.
/// Components reduce to Taylor's x(t) = vx0 τ (1-e^(-t/τ)) and
/// y(t) = -v_ter t + (vy0+v_ter) τ (1-e^(-t/τ)).
///
/// Physics: φ1(0) = 1 and φ2(0) = 1/2, so at b = 0 this is exactly the
/// vacuum r0 + v0 t - (1/2) g t² ŷ. For small u, φ2 = 1/2 - u/6 + ..., so the
/// first drag correction to the height is +g t² u/6 (gravity acts for less
/// time-integrated effect as the ball slows).
/// C++ note: φ2 is computed from its Taylor series for u < 1 because the
/// direct form u - 1 + e^(-u) subtracts nearly equal numbers; see detail::.
/// @param p  drag/gravity parameters
/// @param r0 initial position [m]
/// @param v0 initial velocity [m/s]
/// @param t  elapsed time [s]
[[nodiscard]] inline Vec3 linear_drag_position(const LinearDragParams& p, Vec3 r0, Vec3 v0, double t) {
    const double u = p.b * t / p.m;
    const Vec3 a_grav{0.0, -p.g, 0.0};
    return r0 + v0 * (t * detail::phi1(u)) + a_grav * (t * t * detail::phi2(u));
}

/// The exact 2D trajectory y(x) for gravity + linear drag (Taylor §2.3),
/// for a launch from the origin:
///   y(x) = ((vy0+v_ter)/vx0) x + v_ter τ ln(1 - x/(vx0 τ))
///        = (vy0/vx0) x - (g x²/vx0²) ψ(ε),   ε = x/(vx0 τ) = b x/(m vx0),
///   ψ(ε) = -(ε + ln(1-ε))/ε² = 1/2 + ε/3 + ε²/4 + ...
/// eliminating time entirely — a rare case where a 2D drag problem has a
/// closed-form shape, and the best 2D integrator oracle in the book
/// (research map §8, row 6).
///
/// Physics: this is exactly Taylor's small-drag method (§2.4). Expanding
/// ln(1 - ε) = -ε - ε²/2 - ε³/3 - ... cancels the two leading terms of the
/// textbook formula, leaving the vacuum parabola (ψ = 1/2) plus
/// corrections: ψ ≈ 1/2 + ε/3 gives
///   y ≈ y_vacuum - g x³/(3 vx0³ τ),
/// i.e. drag always pulls the trajectory below the parabola, by an amount
/// growing like x³. At b = 0 the result IS the vacuum parabola, exactly.
/// Only valid for x < vx0·τ (ε < 1): as t -> infinity, x(t) -> vx0·τ (see
/// linear_drag_position) — the horizontal position has a genuine finite
/// asymptote, so y(x) has a real vertical asymptote there, not a numerical
/// artifact. For ε >= 1 (a point the projectile never reaches) this
/// function returns NaN.
/// C++ note: the textbook form multiplies the huge v_ter τ by the tiny
/// ln(1 - ε) and adds a huge (vy0+v_ter) x/vx0 — two enormous terms
/// cancelling to a modest answer. detail::psi() uses the ln(1 - ε) Taylor
/// series for small ε (so ε + ln(1-ε) is never formed) and log1p otherwise.
/// @param p   drag/gravity parameters
/// @param vx0 initial horizontal speed [m/s], assumed > 0
/// @param vy0 initial vertical speed [m/s]
/// @param x   horizontal position [m], must satisfy x < vx0 * p.tau()
[[nodiscard]] inline double linear_drag_trajectory_y(const LinearDragParams& p, double vx0, double vy0,
                                                      double x) {
    const double eps = p.b * x / (p.m * vx0);
    return (vy0 / vx0) * x - (p.g * x * x / (vx0 * vx0)) * detail::psi(eps);
}

/// Horizontal range on level ground (y returns to 0) under gravity + linear
/// drag, for a launch from the origin (Taylor §2.4). Setting y(R) = 0 in
/// linear_drag_trajectory_y() and dividing by R > 0 gives
///   R = R_vac / (2 ψ(b R/(m vx0))),   R_vac = 2 vx0 vy0/g,
/// an implicit equation (R appears on both sides) with no closed-form
/// solution; Taylor solves it by expanding in small drag — see
/// linear_drag_range_small_drag().
///
/// Physics: since ψ >= 1/2, R <= R_vac: drag always shortens the range.
/// C++ note: the root of f(R) = vy0/vx0 - (g R/vx0²) ψ(bR/(m vx0)) is found
/// by BISECTION on (0, min(R_vac, vx0 τ)): f(0⁺) > 0, and f <= 0 at the upper
/// end (at R_vac because ψ >= 1/2; at vx0 τ because ψ -> ∞ at the vertical
/// asymptote), and f is strictly decreasing (R ψ(bR/..) increases), so the
/// root is unique. Bisection is slow-but-sure; it halves the bracket until
/// the two ends are adjacent doubles, i.e. full double precision.
/// Returns 0 if vy0 <= 0 or vx0 <= 0 (the projectile never rises), NaN if
/// b < 0, and exactly R_vac at b = 0.
/// @param p   drag/gravity parameters
/// @param vx0 initial horizontal speed [m/s]
/// @param vy0 initial vertical speed [m/s]
[[nodiscard]] inline double linear_drag_range(const LinearDragParams& p, double vx0, double vy0) {
    if (p.b < 0.0) return std::numeric_limits<double>::quiet_NaN();
    if (!(vy0 > 0.0) || !(vx0 > 0.0)) return 0.0;
    const double r_vac = 2.0 * vx0 * vy0 / p.g;
    if (p.b == 0.0) return r_vac;

    const double k = p.b / (p.m * vx0);  // eps = k R
    const double x_max = 1.0 / k;        // = vx0 tau, where eps = 1
    auto f = [&](double R) {
        const double eps = k * R;
        if (!(eps < 1.0)) return -std::numeric_limits<double>::infinity();
        return vy0 / vx0 - (p.g * R / (vx0 * vx0)) * detail::psi(eps);
    };
    double lo = 0.0;
    double hi = std::fmin(r_vac, x_max);
    for (int iter = 0; iter < 1200; ++iter) {  // 1200 >> needed; loop exits on adjacent doubles
        const double mid = 0.5 * (lo + hi);
        if (mid <= lo || mid >= hi) break;
        if (f(mid) > 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

/// Taylor's first-order small-drag range (§2.4):
///   R ≈ R_vac (1 - (4/3) η),   R_vac = 2 vx0 vy0/g,   η = vy0/v_ter = vy0 b/(m g).
///
/// Physics: substitute ψ ≈ 1/2 + ε/3 in R = R_vac/(2ψ) with ε = 2ηR/R_vac
/// and expand to first order in η. The NEXT term, from ψ ≈ 1/2 + ε/3 + ε²/4,
/// is
///   R ≈ R_vac (1 - (4/3) η + (14/9) η² - ...),
/// so the error of this first-order formula scales like (14/9) η² R_vac.
/// It is good when η << 1, i.e. when the launch's vertical speed is well
/// below the terminal speed (drag changes the range by a few percent).
/// C++ note: η is computed as vy0 b/(m g), never vy0/v_ter, so b = 0 gives
/// η = 0 and the result is exactly R_vac instead of 0·∞ trouble.
/// @param p   drag/gravity parameters
/// @param vx0 initial horizontal speed [m/s]
/// @param vy0 initial vertical speed [m/s]
[[nodiscard]] inline double linear_drag_range_small_drag(const LinearDragParams& p, double vx0, double vy0) {
    const double r_vac = 2.0 * vx0 * vy0 / p.g;
    const double eta = vy0 * p.b / (p.m * p.g);
    return r_vac * (1.0 - (4.0 / 3.0) * eta);
}

/// Quadratic-drag terminal speed for vertical fall: v_ter = sqrt(m g / c)
/// (Taylor §2.4) — the speed at which quadratic drag c v² balances gravity.
/// Returns +∞ at c = 0 (no drag, no terminal speed).
/// @param m mass [kg]
/// @param c quadratic drag coefficient [kg/m]
/// @param g gravitational acceleration [m/s^2]
[[nodiscard]] inline double quadratic_drag_terminal_speed(double m, double c, double g) {
    return std::sqrt(m * g / c);
}

/// Fall speed under quadratic drag from rest:
///   v(t) = v_ter tanh(g t/v_ter) = g t · tanh(z)/z,   z = g t/v_ter
/// (Taylor §2.4) — the exact solution to m v̇ = m g - c v²; monotonically
/// approaches v_ter as t -> infinity (tanh saturates to 1).
///
/// Physics: the tanh(z)/z form shows the physics at a glance: tanh(z)/z =
/// 1 - z²/3 + 2z⁴/15 - ... so v ≈ g t (free fall) with first drag correction
/// -g t z²/3. At v_ter = ∞ (c = 0), z = 0 and v = g t exactly.
/// C++ note: v_ter·tanh(g t/v_ter) is ∞·0 = NaN when v_ter is infinite; the
/// rewrite never forms that product. detail::tanh_over_z uses the series
/// for small z. Signature unchanged.
/// @param v_ter terminal speed [m/s], from quadratic_drag_terminal_speed()
/// @param g     gravitational acceleration [m/s^2]
/// @param t     elapsed time since release [s]
[[nodiscard]] inline double quadratic_drag_fall_speed(double v_ter, double g, double t) {
    const double z = g * t / v_ter;
    return g * t * detail::tanh_over_z(z);
}

/// Fall distance under quadratic drag from rest:
///   y(t) = (v_ter²/g) ln cosh(g t/v_ter) = (1/2) g t² λ(z),
///   z = g t/v_ter,   λ(z) = 2 ln cosh(z)/z² = 1 - z²/6 + ...
/// (Taylor §2.4) — the time integral of quadratic_drag_fall_speed().
///
/// Physics: at c = 0 (z = 0) this is exactly the vacuum (1/2) g t². For
/// z >> 1 the ball is at terminal speed and y ≈ v_ter t - (v_ter²/g) ln 2:
/// the ln 2 is the distance "lost" while still speeding up.
/// C++ note: cosh(z) overflows once z > ~710 (so does the naive formula
/// for any long fall), and for tiny z ln cosh z cancels; detail::ln_cosh
/// and detail::lambda handle both ends. Signature unchanged.
/// @param v_ter terminal speed [m/s], from quadratic_drag_terminal_speed()
/// @param g     gravitational acceleration [m/s^2]
/// @param t     elapsed time since release [s]
[[nodiscard]] inline double quadratic_drag_fall_distance(double v_ter, double g, double t) {
    const double z = g * t / v_ter;
    return 0.5 * g * t * (t * detail::lambda(z));
}

/// Horizontal speed under quadratic drag alone (no gravity component along
/// the motion): v(t) = v0 / (1 + c v0 t / m) (Taylor §2.4) — power-law
/// decay, in contrast to linear drag's exponential decay. Never reaches
/// zero in finite time, but the associated distance (see
/// quadratic_drag_horizontal_distance()) grows without bound
/// (logarithmically) rather than approaching a finite asymptote the way
/// linear drag's x_∞ = vx0 τ does. Exact at c = 0 (gives v0).
/// @param m  mass [kg]
/// @param c  quadratic drag coefficient [kg/m]
/// @param v0 initial speed [m/s]
/// @param t  elapsed time [s]
[[nodiscard]] inline double quadratic_drag_horizontal_speed(double m, double c, double v0, double t) {
    return v0 / (1.0 + c * v0 * t / m);
}

/// Horizontal distance under quadratic drag alone (Taylor §2.4):
///   x(t) = (m/c) ln(1 + ε) = v0 t · ln(1+ε)/ε,   ε = c v0 t/m.
/// The time integral of quadratic_drag_horizontal_speed().
///
/// Physics: this is the ln(1 + ε) side of Taylor's small-drag method:
/// ln(1+ε)/ε = 1 - ε/2 + ε²/3 - ε³/4 + ... so x ≈ v0 t (1 - ε/2), the
/// coasting distance minus a first drag correction. At c = 0 (ε = 0) this
/// is exactly v0 t. Unlike linear drag, x grows forever (like ln t).
/// C++ note: (m/c) ln(1+ε) is ∞·0 at c = 0 and cancels for tiny c; the
/// factored form uses the series for |ε| < 0.01 and log1p otherwise
/// (log1p = ln(1 + x) computed without the cancellation).
/// @param m  mass [kg]
/// @param c  quadratic drag coefficient [kg/m]
/// @param v0 initial speed [m/s]
/// @param t  elapsed time [s]
[[nodiscard]] inline double quadratic_drag_horizontal_distance(double m, double c, double v0, double t) {
    const double eps = c * v0 * t / m;
    return v0 * t * detail::ln1p_over(eps);
}

}  // namespace newton::analytic
