// tests/test_analytic.cpp — internal consistency of newton::analytic's
// closed forms, checked against each other and against known textbook
// identities (not against the numerical integrators — that's
// test_numeric_vs_analytic.cpp).
#include <cmath>

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

    return newton_test::report("analytic");
}
