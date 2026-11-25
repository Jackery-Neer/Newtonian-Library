// tests/test_analytic.cpp — internal consistency of newton::analytic's
// closed forms, checked against each other and against known textbook
// identities (not against the numerical integrators — that's
// test_numeric_vs_analytic.cpp).
#include <cmath>
#include <limits>

#include "test_util.hpp"
#include <newton/newton.hpp>

using namespace newton;

int main() {
    // ======================================================================
    // SHM: x(0) = x0, v(0) = v0 exactly (the whole point of the closed-form
    // x0*cos + (v0/w0)*sin construction over an amplitude/phase round trip
    // — see analytic/oscillator.hpp's C++ note on SHM::x()).
    // ======================================================================
    {
        for (double omega0 : {0.5, 1.0, 3.7}) {
            for (double x0 : {-2.0, 0.0, 1.5}) {
                for (double v0 : {-1.0, 0.0, 4.0}) {
                    analytic::SHM shm{omega0, x0, v0};
                    NT_CLOSE(shm.x(0.0), x0, 1e-12);
                    NT_CLOSE(shm.v(0.0), v0, 1e-12);
                    // Energy computed via omega0/x0/v0 must equal the direct
                    // 1/2 m (v0^2 + omega0^2 x0^2) formula at t=0, and stay
                    // constant at a later time too (SHM is conservative).
                    const double m = 2.0;
                    const double E0 = 0.5 * m * (v0 * v0 + omega0 * omega0 * x0 * x0);
                    NT_CLOSE(shm.energy(m), E0, 1e-9);
                    const double t = 1.3;
                    const double Et = 0.5 * m * (shm.v(t) * shm.v(t) + omega0 * omega0 * shm.x(t) * shm.x(t));
                    NT_CLOSE(Et, E0, 1e-9);
                }
            }
        }
        // period() is independent of amplitude (defining feature of SHM).
        analytic::SHM small{2.0, 0.001, 0.0};
        analytic::SHM large{2.0, 1000.0, 0.0};
        NT_CLOSE(small.period(), large.period(), 1e-12);
        NT_CLOSE(small.period(), 2.0 * constants::pi / 2.0, 1e-12);
    }

    // ======================================================================
    // Damped oscillator: x(0)=x0, v(0)=v0 in ALL four regimes, and the three
    // formula branches are continuous in t across the beta==omega0 boundary
    // that regime() routes around (research: the underdamped/overdamped
    // formulas individually blow up there; the physics does not).
    // ======================================================================
    {
        const double omega0 = 2.0;
        const double x0 = 1.5;
        const double v0 = -0.7;

        analytic::DampedOscillator under{omega0, 0.5, x0, v0};    // beta < omega0
        analytic::DampedOscillator crit{omega0, omega0, x0, v0};  // beta == omega0
        analytic::DampedOscillator over{omega0, 3.0, x0, v0};     // beta > omega0

        NT_CHECK(under.regime() == analytic::Regime::Underdamped);
        NT_CHECK(crit.regime() == analytic::Regime::CriticallyDamped);
        NT_CHECK(over.regime() == analytic::Regime::Overdamped);
        // beta == 0 is the distinct Undamped regime, not folded into Underdamped.
        analytic::DampedOscillator undamped{omega0, 0.0, x0, v0};
        NT_CHECK(undamped.regime() == analytic::Regime::Undamped);

        for (const analytic::DampedOscillator* d : {&under, &crit, &over, &undamped}) {
            NT_CLOSE(d->x(0.0), x0, 1e-12);
            NT_CLOSE(d->v(0.0), v0, 1e-12);
        }

        // Continuity across the critical boundary: as beta -> omega0 from
        // both sides, x(t) and v(t) converge to the critically-damped values
        // at a fixed later time t. Test at beta = omega0 +- a small delta
        // well outside regime()'s own kRelEps routing window (so we're
        // genuinely exercising the Underdamped/Overdamped formula branches,
        // not the CriticallyDamped one), yet small enough that the physical
        // continuity is visible.
        const double t = 0.8;
        const double x_crit = crit.x(t);
        analytic::DampedOscillator near_under{omega0, omega0 - 1e-4, x0, v0};
        analytic::DampedOscillator near_over{omega0, omega0 + 1e-4, x0, v0};
        NT_CHECK(near_under.regime() == analytic::Regime::Underdamped);
        NT_CHECK(near_over.regime() == analytic::Regime::Overdamped);
        NT_CLOSE(near_under.x(t), x_crit, 1e-3);
        NT_CLOSE(near_over.x(t), x_crit, 1e-3);
    }

    // ======================================================================
    // Driven oscillator: amplitude() peaks at resonant_omega() (Taylor
    // §5.6) — sample on both sides and check the resonance value is the
    // local maximum, and that it exceeds amplitude() far from resonance.
    // ======================================================================
    {
        analytic::DrivenOscillator drv{2.0, 0.3, 1.0, 2.0};  // omega set per-sample below
        const double omega_res = drv.resonant_omega();
        NT_CHECK(!std::isnan(omega_res));  // beta=0.3 < omega0/sqrt(2)=1.414: a real peak exists

        drv.omega = omega_res;
        const double A_res = drv.amplitude();
        drv.omega = omega_res - 0.05;
        const double A_below = drv.amplitude();
        drv.omega = omega_res + 0.05;
        const double A_above = drv.amplitude();

        NT_CHECK(A_res > A_below);
        NT_CHECK(A_res > A_above);

        // Far from resonance the amplitude is much smaller.
        drv.omega = omega_res * 5.0;
        const double A_far = drv.amplitude();
        NT_CHECK(A_res > 10.0 * A_far);

        // Phase passes through pi/2 exactly at omega = omega0 (Taylor §5.6),
        // independent of beta.
        analytic::DrivenOscillator at_omega0{2.0, 0.3, 1.0, 2.0};
        NT_CLOSE(at_omega0.phase(), constants::pi / 2.0, 1e-9);

        // Quality factor Q = omega0/(2 beta).
        NT_CLOSE(drv.quality_factor(), 2.0 / (2.0 * 0.3), 1e-12);
        // FWHM = 2 beta.
        NT_CLOSE(drv.fwhm(), 0.6, 1e-12);
    }

    // ======================================================================
    // Vacuum projectile: range at 45 degrees is the textbook v0^2/g maximum
    // (Taylor §1.6), and position()/time_of_flight()/apex_height() agree
    // with each other internally.
    // ======================================================================
    {
        const double v0 = 25.0;
        const double g = constants::g_earth;
        NT_CLOSE(analytic::range(v0, constants::pi / 4.0, g), v0 * v0 / g, 1e-12);

        // 45 degrees is the maximum: any other angle gives a smaller range.
        NT_CHECK(analytic::range(v0, constants::pi / 4.0, g) > analytic::range(v0, constants::pi / 6.0, g));
        NT_CHECK(analytic::range(v0, constants::pi / 4.0, g) > analytic::range(v0, constants::pi / 3.0, g));

        // position() at the time_of_flight() lands back at y = 0, at x = range().
        const double theta = 0.7;  // rad, an arbitrary non-special angle
        const double T = analytic::time_of_flight(v0, theta, g);
        const Vec3 landing = analytic::position(T, v0, theta, g);
        NT_CLOSE(landing.x, analytic::range(v0, theta, g), 1e-9);
        NT_CLOSE(landing.y, 0.0, 1e-9);

        // apex_height() matches position() evaluated at the apex time T/2
        // (symmetric parabola in vacuum).
        const Vec3 apex_pos = analytic::position(T / 2.0, v0, theta, g);
        NT_CLOSE(apex_pos.y, analytic::apex_height(v0, theta, g), 1e-9);
    }

    // ======================================================================
    // Linear drag, small-b regime (Taylor §2.4's small-drag expansion, and
    // the bug it fixes): the old closed forms multiplied tau = m/b (huge) by
    // 1 - e^(-t/tau) (tiny) and returned 24 km at b = 1e-10, NaN at b = 0.
    // Every check here uses ONLY quantities that cannot be wrong by
    // cancellation: the vacuum formula plus the analytic small-drag series.
    // ======================================================================
    {
        const double g = constants::g_earth;
        const double theta = 45.0 * constants::deg;
        const double vx0 = 30.0 * std::cos(theta), vy0 = 30.0 * std::sin(theta);
        const double t = 3.0;
        const double y_vac = vy0 * t - 0.5 * g * t * t;  // 19.5097 m
        const double x_vac = vx0 * t;
        const Vec3 v0{vx0, vy0, 0.0};

        // b = 0 (and the default-constructed params): exactly the vacuum formulas.
        for (const analytic::LinearDragParams p : {analytic::LinearDragParams{}, analytic::LinearDragParams{1.0, 0.0, g},
                                                   analytic::LinearDragParams{3.0, 0.0, g}}) {
            const Vec3 r = analytic::linear_drag_position(p, Vec3{}, v0, t);
            NT_CHECK(std::isfinite(r.x) && std::isfinite(r.y));
            NT_CLOSE(r.x, x_vac, 1e-15);
            NT_CLOSE(r.y, y_vac, 1e-15);
            const Vec3 v = analytic::linear_drag_velocity(p, v0, t);
            NT_CLOSE(v.x, vx0, 1e-15);
            NT_CLOSE(v.y, vy0 - g * t, 1e-15);
            // y(x): the vacuum parabola.
            NT_CLOSE(analytic::linear_drag_trajectory_y(p, vx0, vy0, 45.0),
                     vy0 / vx0 * 45.0 - g * 45.0 * 45.0 / (2.0 * vx0 * vx0), 1e-14);
            // Range: exactly R_vac.
            NT_CLOSE(analytic::linear_drag_range(p, vx0, vy0), 2.0 * vx0 * vy0 / g, 1e-15);
            NT_CLOSE(analytic::linear_drag_range_small_drag(p, vx0, vy0), 2.0 * vx0 * vy0 / g, 1e-15);
        }
        // Non-trivial r0 passes through untouched.
        {
            const analytic::LinearDragParams p{1.0, 0.0, g};
            const Vec3 r = analytic::linear_drag_position(p, Vec3{1.0, 2.0, 3.0}, v0, t);
            NT_CLOSE(r.x, 1.0 + x_vac, 1e-15);
            NT_CLOSE(r.y, 2.0 + y_vac, 1e-15);
            NT_CLOSE(r.z, 3.0, 1e-15);
        }

        // Tiny b: vacuum + first and second order in u = b t/m (from
        // phi1 = 1 - u/2 + u^2/6, phi2 = 1/2 - u/6 + u^2/24 — the e^(-u)
        // Taylor series), to 1e-12 RELATIVE. The b = 1e-10 case is the old
        // 24-km disaster; b = 1e-14 was -2.35e12 m.
        for (double b : {1e-6, 1e-10, 1e-14, 1e-300}) {
            const analytic::LinearDragParams p{1.0, b, g};
            const double u = b * t / p.m;
            const Vec3 r = analytic::linear_drag_position(p, Vec3{}, v0, t);
            const double y_exp = y_vac + u * (g * t * t / 6.0 - vy0 * t / 2.0) +
                                 u * u * (vy0 * t / 6.0 - g * t * t / 24.0);
            const double x_exp = x_vac * (1.0 - u / 2.0 + u * u / 6.0);
            NT_CLOSE(r.y, y_exp, 1e-12);
            NT_CLOSE(r.x, x_exp, 1e-12);
            // The textbook value: 19.5097 m (NOT 24405.8 m).
            NT_CHECK(std::fabs(r.y - 19.5097) < 1e-4);

            // Velocity: v0 e^-u + a t phi1(u) -> vacuum to O(u).
            const Vec3 v = analytic::linear_drag_velocity(p, v0, t);
            NT_CLOSE(v.x, vx0 * (1.0 - u + u * u / 2.0), 1e-12);
            NT_CLOSE(v.y, vy0 * (1.0 - u + u * u / 2.0) - g * t * (1.0 - u / 2.0 + u * u / 6.0), 1e-12);

            // y(x) at x = 45 m: eps = b x/(m vx0); y = y_vac(x) - g x^2/vx0^2 (eps/3 + eps^2/4 + eps^3/5).
            const double x = 45.0;
            const double eps = b * x / (p.m * vx0);
            const double yx_vac = vy0 / vx0 * x - g * x * x / (2.0 * vx0 * vx0);
            const double yx_exp = yx_vac - g * x * x / (vx0 * vx0) * (eps / 3.0 + eps * eps / 4.0 + eps * eps * eps / 5.0);
            NT_CLOSE(analytic::linear_drag_trajectory_y(p, vx0, vy0, x), yx_exp, 1e-12);
            // ... and the leading drag correction -g x^3/(3 vx0^3 tau) quoted in the doc comment.
            NT_CLOSE(analytic::linear_drag_trajectory_y(p, vx0, vy0, x) - yx_vac,
                     -g * x * x * x / (3.0 * vx0 * vx0 * vx0 * p.tau()), 1e-5 + 3.0 * eps);
        }

        // y(x) is consistent with x(t), y(t) at moderate drag (the two
        // closed forms are independent derivations; this is a real check).
        {
            const analytic::LinearDragParams p{0.5, 0.2, g};
            for (double tt : {0.5, 1.0, 2.0, 4.0}) {
                const Vec3 r = analytic::linear_drag_position(p, Vec3{}, v0, tt);
                NT_CLOSE(analytic::linear_drag_trajectory_y(p, vx0, vy0, r.x), r.y, 1e-11);
            }
        }

        // x beyond the vertical asymptote vx0*tau: NaN, as documented.
        {
            const analytic::LinearDragParams p{1.0, 0.5, g};
            const double x_max = vx0 * p.tau();
            NT_CHECK(std::isnan(analytic::linear_drag_trajectory_y(p, vx0, vy0, x_max)));
            NT_CHECK(std::isnan(analytic::linear_drag_trajectory_y(p, vx0, vy0, 1.5 * x_max)));
            NT_CHECK(std::isfinite(analytic::linear_drag_trajectory_y(p, vx0, vy0, 0.99 * x_max)));
        }

        // tau()/v_ter() are +infinity at b = 0 and nothing multiplies them by 0.
        NT_CHECK(std::isinf(analytic::LinearDragParams{}.tau()));
        NT_CHECK(std::isinf(analytic::LinearDragParams{}.v_ter()));
    }

    // ======================================================================
    // The series/direct crossovers of phi2 (u = 1), psi (eps = 1/2), ln_cosh
    // (z = 10), lambda and tanh/z (1e-3), ln(1+eps)/eps (0.01): continuous to
    // ~1e-15 relative, and each agrees with an independent high-term sum.
    // ======================================================================
    {
        auto phi2_series = [](double u) {  // 60 terms, valid for moderate u
            double term = 0.5, sum = 0.5;
            for (int n = 1; n < 60; ++n) {
                term *= -u / static_cast<double>(n + 2);
                sum += term;
            }
            return sum;
        };
        auto psi_series = [](double e) {  // 3000 terms, valid for e < ~0.95
            double sum = 0.0;
            for (int n = 2999; n >= 0; --n) sum = sum * e + 1.0 / static_cast<double>(n + 2);
            return sum;
        };

        NT_CLOSE(analytic::detail::phi2(0.0), 0.5, 1e-16);
        NT_CLOSE(analytic::detail::phi1(0.0), 1.0, 1e-16);
        for (double u : {1e-12, 1e-6, 0.01, 0.1, 0.5, 0.999999, 1.0, 1.000001, 2.0, 5.0}) {
            NT_CLOSE(analytic::detail::phi2(u), phi2_series(u), 1e-14);
            // phi1: series 1 - u/2 + u^2/6 - u^3/24 for tiny u, textbook form where it is accurate.
            if (u < 1e-3) {
                NT_CLOSE(analytic::detail::phi1(u), 1.0 - u / 2.0 + u * u / 6.0 - u * u * u / 24.0, 1e-15);
            } else {
                NT_CLOSE(analytic::detail::phi1(u), (1.0 - std::exp(-u)) / u, 1e-13);
            }
        }
        // Across the threshold (series below u = 1, direct above).
        NT_CLOSE(analytic::detail::phi2(std::nextafter(1.0, 0.0)), analytic::detail::phi2(1.0), 2e-15);
        NT_CLOSE(analytic::detail::phi2(std::nextafter(1.0, 2.0)), analytic::detail::phi2(1.0), 2e-15);
        // Large u: phi2 ~ 1/u and no overflow even where u*u would overflow.
        NT_CLOSE(analytic::detail::phi2(1e200), 1e-200, 1e-14);

        NT_CLOSE(analytic::detail::psi(0.0), 0.5, 1e-16);
        for (double e : {1e-12, 1e-6, 0.01, 0.2, 0.4999, 0.5, 0.5001, 0.7, 0.9, -0.3, -0.7}) {
            NT_CLOSE(analytic::detail::psi(e), psi_series(e), 1e-14);
        }
        NT_CLOSE(analytic::detail::psi(std::nextafter(0.5, 0.0)), analytic::detail::psi(0.5), 2e-15);
        NT_CLOSE(analytic::detail::psi(std::nextafter(0.5, 1.0)), analytic::detail::psi(0.5), 2e-15);
        NT_CHECK(std::isnan(analytic::detail::psi(1.0)));
        NT_CHECK(std::isnan(analytic::detail::psi(2.0)));
        // psi >= 1/2 for eps >= 0 (it is why drag always shortens the range).
        for (double e : {0.0, 1e-9, 0.1, 0.5, 0.9, 0.999}) NT_CHECK(analytic::detail::psi(e) >= 0.5);

        // ln cosh across z = 10, and against log(cosh) where that is safe.
        for (double z : {0.0, 1e-8, 0.5, 3.0, 9.999, 10.0, 10.001, 15.0}) {
            if (z > 0.0) NT_CLOSE(analytic::detail::ln_cosh(z), std::log(std::cosh(z)), 1e-14);
            NT_CLOSE(analytic::detail::ln_cosh(-z), analytic::detail::ln_cosh(z), 1e-16);
        }
        NT_CLOSE(analytic::detail::ln_cosh(std::nextafter(10.0, 0.0)), analytic::detail::ln_cosh(10.0), 2e-15);
        NT_CLOSE(analytic::detail::ln_cosh(std::nextafter(10.0, 11.0)), analytic::detail::ln_cosh(10.0), 2e-15);
        NT_CLOSE(analytic::detail::ln_cosh(1e5), 1e5 - std::log(2.0), 1e-15);  // no overflow

        // lambda and tanh/z across 1e-3, and at 0.
        NT_CLOSE(analytic::detail::lambda(0.0), 1.0, 1e-16);
        NT_CLOSE(analytic::detail::tanh_over_z(0.0), 1.0, 1e-16);
        NT_CLOSE(analytic::detail::lambda(std::nextafter(1e-3, 0.0)), analytic::detail::lambda(1e-3), 2e-15);
        NT_CLOSE(analytic::detail::tanh_over_z(std::nextafter(1e-3, 0.0)), analytic::detail::tanh_over_z(1e-3), 2e-15);
        for (double z : {1e-5, 1e-3, 0.01, 0.5, 2.0, -2.0}) {
            NT_CLOSE(analytic::detail::tanh_over_z(z), std::tanh(z) / z, 1e-14);
            // log(cosh z) itself cancels for small z, so only compare with it where it is accurate.
            if (std::fabs(z) >= 0.5) {
                NT_CLOSE(analytic::detail::lambda(z), 2.0 * std::log(std::cosh(z)) / (z * z), 1e-13);
            } else {
                const double z2 = z * z;
                NT_CLOSE(analytic::detail::lambda(z),
                         1.0 - z2 / 6.0 + 2.0 * z2 * z2 / 45.0 - 17.0 * z2 * z2 * z2 / 1260.0 + 62.0 * z2 * z2 * z2 * z2 / 14175.0,
                         1e-15);
            }
        }
        // ln(1+eps)/eps across 0.01.
        NT_CLOSE(analytic::detail::ln1p_over(0.0), 1.0, 1e-16);
        NT_CLOSE(analytic::detail::ln1p_over(std::nextafter(0.01, 0.0)), analytic::detail::ln1p_over(0.01), 2e-15);
        for (double e : {1e-9, 0.005, 0.01, 0.5, 3.0, -0.5, -0.005}) {
            NT_CLOSE(analytic::detail::ln1p_over(e), std::log1p(e) / e, 1e-14);
        }
    }

    // ======================================================================
    // linear_drag_range (Taylor §2.4): the root of y(R) = 0, found by
    // bisection; self-consistent, shortened by drag, and matched by the
    // small-drag expansion R_vac (1 - 4 eta/3 + 14 eta^2/9 - ...).
    // ======================================================================
    {
        const double g = constants::g_earth;
        const double vx0 = 15.0, vy0 = 20.0, m = 0.7;
        const double r_vac = 2.0 * vx0 * vy0 / g;

        double r_prev = r_vac;
        for (double b : {1e-12, 1e-6, 1e-3, 0.05, 0.3, 1.0}) {
            const analytic::LinearDragParams p{m, b, g};
            const double R = analytic::linear_drag_range(p, vx0, vy0);
            NT_CHECK(R > 0.0 && R <= r_prev);  // drag shortens the range, monotonically in b
            r_prev = R;
            NT_CHECK(R < vx0 * p.tau());  // stays inside the vertical asymptote
            // y(R) = 0, and R = R_vac/(2 psi): the defining relations.
            const double eps = b * R / (m * vx0);
            NT_CLOSE(R, r_vac / (2.0 * analytic::detail::psi(eps)), 1e-13);
            NT_CHECK(std::fabs(analytic::linear_drag_trajectory_y(p, vx0, vy0, R)) < 1e-9 * (1.0 + R));
        }
        // Overwhelming drag (b = 5, 50: tau = 0.14 s, 0.014 s): the ball stops
        // within a hair of the asymptote vx0 tau, which the root approaches
        // (to double precision at b = 50) but never passes.
        for (double b : {5.0, 50.0}) {
            const analytic::LinearDragParams p{m, b, g};
            const double R = analytic::linear_drag_range(p, vx0, vy0);
            NT_CHECK(R <= vx0 * p.tau() && R > 0.99 * vx0 * p.tau() && R < r_prev);
            r_prev = R;
        }
        // b -> 0 limit is continuous, vy0 <= 0 gives 0.
        NT_CLOSE(analytic::linear_drag_range(analytic::LinearDragParams{m, 1e-300, g}, vx0, vy0), r_vac, 1e-15);
        NT_CHECK(analytic::linear_drag_range(analytic::LinearDragParams{m, 0.1, g}, vx0, 0.0) == 0.0);
        NT_CHECK(analytic::linear_drag_range(analytic::LinearDragParams{m, 0.1, g}, vx0, -3.0) == 0.0);
        NT_CHECK(analytic::linear_drag_range(analytic::LinearDragParams{m, 0.0, g}, vx0, -3.0) == 0.0);

        // Small-drag expansion: eta = vy0 b/(m g). The first-order error
        // scales like eta^2 (halving eta quarters it), and the measured
        // coefficient (R - R1)/(R_vac eta^2) tends to Taylor's 14/9.
        auto err_first_order = [&](double eta) {
            const double b = eta * m * g / vy0;
            const analytic::LinearDragParams p{m, b, g};
            const double R = analytic::linear_drag_range(p, vx0, vy0);
            return R - analytic::linear_drag_range_small_drag(p, vx0, vy0);
        };
        const double e1 = err_first_order(2e-3), e2 = err_first_order(1e-3);
        NT_CLOSE(e1 / e2, 4.0, 5e-3);
        NT_CLOSE(e2 / (r_vac * 1e-3 * 1e-3), 14.0 / 9.0, 1e-2);
        NT_CLOSE(err_first_order(1e-4) / (r_vac * 1e-8), 14.0 / 9.0, 2e-3);
        // The second-order formula is far better than first-order at eta = 0.02.
        {
            const double eta = 0.02, b = eta * m * g / vy0;
            const analytic::LinearDragParams p{m, b, g};
            const double R = analytic::linear_drag_range(p, vx0, vy0);
            const double R2 = r_vac * (1.0 - 4.0 / 3.0 * eta + 14.0 / 9.0 * eta * eta);
            NT_CHECK(std::fabs(R - R2) < 0.05 * std::fabs(R - analytic::linear_drag_range_small_drag(p, vx0, vy0)));
        }
    }

    // ======================================================================
    // Quadratic drag closed forms (Taylor §2.4), rewritten with
    // z = g t/v_ter: exact at c = 0, finite at huge z, textbook values
    // elsewhere.
    // ======================================================================
    {
        const double g = constants::g_earth;
        const double inf = std::numeric_limits<double>::infinity();

        // c = 0: v_ter = infinity, and fall speed/distance are the vacuum g t, g t^2/2.
        NT_CHECK(std::isinf(analytic::quadratic_drag_terminal_speed(1.0, 0.0, g)));
        for (double t : {0.0, 0.5, 3.0, 100.0}) {
            NT_CLOSE(analytic::quadratic_drag_fall_speed(inf, g, t), g * t, 1e-15);
            NT_CLOSE(analytic::quadratic_drag_fall_distance(inf, g, t), 0.5 * g * t * t, 1e-15);
            NT_CLOSE(analytic::quadratic_drag_horizontal_distance(2.0, 0.0, 12.0, t), 12.0 * t, 1e-15);
        }

        // Moderate z vs the textbook formulas (accurate there).
        for (double v_ter : {1.0, 5.0, 40.0}) {
            for (double t : {0.1, 0.7, 2.0, 10.0}) {
                const double z = g * t / v_ter;
                if (z > 30.0) continue;
                NT_CLOSE(analytic::quadratic_drag_fall_speed(v_ter, g, t), v_ter * std::tanh(z), 1e-13);
                NT_CLOSE(analytic::quadratic_drag_fall_distance(v_ter, g, t),
                         v_ter * v_ter / g * std::log(std::cosh(z)), 1e-13);
            }
        }

        // Tiny z: first-order correction -z^2/3 in the speed, -z^2/6 in the distance factor.
        {
            const double v_ter = 1e4, t = 2.0, z = g * t / v_ter;
            NT_CLOSE(analytic::quadratic_drag_fall_speed(v_ter, g, t), g * t * (1.0 - z * z / 3.0 + 2.0 * z * z * z * z / 15.0), 1e-15);
            NT_CLOSE(analytic::quadratic_drag_fall_distance(v_ter, g, t), 0.5 * g * t * t * (1.0 - z * z / 6.0 + 2.0 * z * z * z * z / 45.0), 1e-15);
        }

        // No overflow at g t/v_ter = 1000 (the old cosh() overflowed past ~710):
        // d = v_ter t - (v_ter^2/g) ln 2 and v = v_ter.
        for (double zbig : {1000.0, 1e5}) {
            const double v_ter = 2.0, t = zbig * v_ter / g;
            const double d = analytic::quadratic_drag_fall_distance(v_ter, g, t);
            NT_CHECK(std::isfinite(d));
            NT_CLOSE(d, v_ter * t - v_ter * v_ter / g * std::log(2.0), 1e-14);
            NT_CLOSE(analytic::quadratic_drag_fall_speed(v_ter, g, t), v_ter, 1e-14);
        }

        // Horizontal distance (m/c) ln(1 + eps): vs the textbook form where it is accurate,
        // and its series limit; and consistency with the speed (x' = v) by finite difference.
        {
            const double m = 2.0, c = 0.3, v0 = 10.0;
            for (double t : {0.05, 0.5, 5.0, 50.0}) {
                NT_CLOSE(analytic::quadratic_drag_horizontal_distance(m, c, v0, t),
                         m / c * std::log(1.0 + c * v0 * t / m), 1e-13);
            }
            const double t = 3.0, h = 1e-5;
            const double dxdt = (analytic::quadratic_drag_horizontal_distance(m, c, v0, t + h) -
                                 analytic::quadratic_drag_horizontal_distance(m, c, v0, t - h)) / (2.0 * h);
            NT_CLOSE(dxdt, analytic::quadratic_drag_horizontal_speed(m, c, v0, t), 1e-8);
            // Tiny c: v0 t (1 - eps/2 + eps^2/3).
            const double ct = 1e-12, eps = ct * v0 * 4.0 / m;
            NT_CLOSE(analytic::quadratic_drag_horizontal_distance(m, ct, v0, 4.0), v0 * 4.0 * (1.0 - eps / 2.0 + eps * eps / 3.0), 1e-15);
        }
    }

    return newton_test::report("analytic");
}
