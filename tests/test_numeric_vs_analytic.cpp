// tests/test_numeric_vs_analytic.cpp — the library's soul: does the general
// numerical machinery (forces + integrators + simulate()) actually reproduce
// every closed-form solution newton::analytic knows about? Where
// test_analytic.cpp checks the formulas are internally consistent and
// test_integrators.cpp checks the integrators' abstract numerical
// properties, this file is the one that proves the two halves of the
// library agree with each other on real physical setups.
#include <cmath>

#include "test_util.hpp"
#include <newton/newton.hpp>

using namespace newton;

int main() {
    // ======================================================================
    // Linear drag: the exact 2D trajectory y(x) (Taylor §2.3), the richest
    // 2D validation oracle in the book — see analytic/drag.hpp.
    // ======================================================================
    {
        const double m = 1.0, b = 0.5, g = constants::g_earth;
        analytic::LinearDragParams params{m, b, g};
        const double vx0 = 10.0, vy0 = 5.0;
        // x -> vx0*tau = 20 m as t -> infinity; running to t=1.5s keeps x well
        // inside that asymptote (x(1.5) ~ 10.6 m) so y(x) stays finite.

        Particle p{m, Vec3{}, Vec3{vx0, vy0, 0.0}};
        auto force = forces::sum(forces::UniformGravity{m}, forces::LinearDrag{b});
        SimOptions opt;
        opt.t0 = 0.0;
        opt.t_end = 1.5;
        opt.dt = 1e-4;
        opt.method = Method::RK4;
        const Trajectory traj = simulate(p, force, opt);

        double max_err = 0.0;
        for (const State& s : traj.states) {
            if (s.r.x <= 0.0) continue;  // y(x) is undefined/trivial at x=0 (log(1)=0, skip the singular start)
            const double y_exact = analytic::linear_drag_trajectory_y(params, vx0, vy0, s.r.x);
            max_err = std::fmax(max_err, std::fabs(s.r.y - y_exact));
        }
        NT_CHECK(max_err < 1e-6);  // DESIGN.md: rtol 1e-6 with RK4 dt=1e-4 (actual error ~1e-13)
    }

    // ======================================================================
    // Quadratic drag: vertical fall from rest, the one nonlinear-drag case
    // with a closed form (Taylor §2.4) — tanh speed profile and its
    // terminal-velocity limit.
    // ======================================================================
    {
        const double m = 1.0, c = 0.1, g = constants::g_earth;
        const double v_ter = analytic::quadratic_drag_terminal_speed(m, c, g);

        Particle p{m, Vec3{}, Vec3{}};  // dropped from rest
        auto force = forces::sum(forces::UniformGravity{m}, forces::QuadraticDrag{c});
        SimOptions opt;
        opt.t0 = 0.0;
        opt.t_end = 3.0;  // ~3 characteristic times v_ter/g; enough to approach v_ter, not enough to saturate to roundoff
        opt.dt = 1e-3;
        opt.method = Method::RK4;
        const Trajectory traj = simulate(p, force, opt);

        double max_err_v = 0.0, max_err_y = 0.0;
        for (const State& s : traj.states) {
            const double speed_exact = analytic::quadratic_drag_fall_speed(v_ter, g, s.t);
            const double dist_exact = analytic::quadratic_drag_fall_distance(v_ter, g, s.t);
            max_err_v = std::fmax(max_err_v, std::fabs(-s.v.y - speed_exact));  // falling: v.y < 0
            max_err_y = std::fmax(max_err_y, std::fabs(-s.r.y - dist_exact));   // falling: r.y < 0
        }
        NT_CHECK(max_err_v < 1e-6);
        NT_CHECK(max_err_y < 1e-6);

        // Terminal-velocity approach: after ~3 characteristic times, speed is
        // within 1% of v_ter (tanh(3) = 0.995, so this is the physics, not
        // just integrator slop).
        const double final_speed = -traj.back().v.y;
        NT_CHECK(std::fabs(final_speed - v_ter) / v_ter < 0.01);
    }

    // ======================================================================
    // Underdamped oscillator vs closed form (Taylor §5.4) — exercises a
    // velocity-dependent force (LinearDrag) through RK4.
    // ======================================================================
    {
        const double m = 1.0, k = 4.0, b = 0.4;
        const double omega0 = std::sqrt(k / m);
        const double beta = b / (2.0 * m);
        const double omega1 = std::sqrt(omega0 * omega0 - beta * beta);
        const double T1 = 2.0 * constants::pi / omega1;

        Particle p{m, Vec3{1.0, 0.0, 0.0}, Vec3{}};
        auto force = forces::sum(forces::Spring{k, Vec3{}, 0.0}, forces::LinearDrag{b});
        SimOptions opt;
        opt.t0 = 0.0;
        opt.t_end = 3.0 * T1;
        opt.dt = T1 / 1000.0;
        opt.method = Method::RK4;
        const Trajectory traj = simulate(p, force, opt);

        analytic::DampedOscillator ref{omega0, beta, 1.0, 0.0};
        double max_err = 0.0;
        for (const State& s : traj.states) {
            max_err = std::fmax(max_err, std::fabs(s.r.x - ref.x(s.t)));
        }
        NT_CHECK(max_err < 1e-6);
    }

    // ======================================================================
    // Driven damped oscillator: steady-state amplitude within 0.5% (Taylor
    // §5.5-5.6) — this is the test that catches a stage-time bug in RK4
    // (research §4.1(d)): an autonomous-force test can't see "t_n used where
    // t_n + dt/2 belongs", but a time-dependent drive can.
    // ======================================================================
    {
        const double m = 1.0, k = 4.0, b = 0.4, F0 = 1.0, omega = 1.0;
        const double omega0 = std::sqrt(k / m);
        const double beta = b / (2.0 * m);
        const double f0 = F0 / m;

        Particle p{m, Vec3{}, Vec3{}};
        auto force = forces::sum(forces::Spring{k, Vec3{}, 0.0}, forces::LinearDrag{b},
                                  forces::SinusoidalDrive{F0, omega, 0.0, Vec3{1.0, 0.0, 0.0}});
        const double T_drive = 2.0 * constants::pi / omega;
        SimOptions opt;
        opt.t0 = 0.0;
        opt.t_end = 60.0;  // >> 1/beta = 5 s: the transient has long decayed
        opt.dt = T_drive / 200.0;
        opt.method = Method::RK4;
        const Trajectory traj = simulate(p, force, opt);

        analytic::DrivenOscillator ref{omega0, beta, f0, omega};
        const double amp_measured = analysis::measure_amplitude(traj, 40.0);  // 8 damping times in
        const double amp_exact = ref.amplitude();
        NT_CHECK(std::fabs(amp_measured - amp_exact) / amp_exact < 0.005);  // DESIGN.md: within 0.5%
    }

    // ======================================================================
    // Rocket: systems::simulate_rocket's variable-mass velocity vs
    // analytic::rocket_velocity_gravity (Tsiolkovsky + constant gravity,
    // Taylor §3.2) throughout the burn, including the exact burnout instant.
    // ======================================================================
    {
        systems::RocketParams rp{1000.0, 200.0, 2500.0, 50.0, constants::g_earth};
        const double t_burnout = (rp.m0 - rp.m_dry) / rp.burn_rate;
        const auto samples = systems::simulate_rocket(rp, 0.01, t_burnout);

        double max_err = 0.0;
        bool found_burnout_sample = false;
        for (const auto& s : samples) {
            const double v_exact = analytic::rocket_velocity_gravity(rp.v_ex, rp.m0, rp.burn_rate, rp.g, s.t);
            max_err = std::fmax(max_err, std::fabs(s.v - v_exact));
            if (std::fabs(s.t - t_burnout) < 1e-9) found_burnout_sample = true;
        }
        NT_CHECK(max_err < 1e-6);
        // DESIGN.md requires the burnout instant to be an exact sample (RK4
        // stages land the last burn-phase step precisely on it), not merely
        // straddled by the nearest grid point.
        NT_CHECK(found_burnout_sample);
    }

    // ======================================================================
    // Magnetic force: |v| conservation to < 1e-12 and the cyclotron period
    // (Taylor §2.5) — the magnetic force does no work, so |v| is exactly
    // constant along the TRUE trajectory; here we check how tightly RK4
    // reproduces that invariant at a small dt.
    // ======================================================================
    {
        const double m = 1.0, q = 1.0;
        const Vec3 B{0.0, 0.0, 2.0};
        const double T_c = analytic::cyclotron_period(m, q, B.z);

        Particle p{m, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 3.0, 0.0}};  // v0 entirely perpendicular to B
        forces::MagneticForce force{q, B};
        SimOptions opt;
        opt.t0 = 0.0;
        opt.t_end = T_c;  // exactly one full revolution
        opt.dt = T_c / 1000.0;
        opt.method = Method::RK4;
        const Trajectory traj = simulate(p, force, opt);

        const double v0_norm = norm(traj.front().v);
        double max_rel_dev = 0.0;
        for (const State& s : traj.states) {
            max_rel_dev = std::fmax(max_rel_dev, std::fabs(norm(s.v) - v0_norm) / v0_norm);
        }
        NT_CHECK(max_rel_dev < 1e-12);  // DESIGN.md's explicit number

        // Cyclotron period check: after exactly one period, velocity (not
        // just its magnitude) has returned to its starting value.
        NT_CLOSE_VEC3(traj.back().v, traj.front().v, 1e-9);
        NT_CLOSE_VEC3(traj.back().r, traj.front().r, 1e-9);

        // cyclotron_radius() matches the geometric radius: the farthest the
        // particle gets from its starting center is 2r (diameter), reached
        // at the half-period point.
        const double r_expected = analytic::cyclotron_radius(m, v0_norm, q, B.z);
        NT_CLOSE(r_expected, m * v0_norm / std::fabs(q * B.z), 1e-12);
    }

    // ======================================================================
    // simulate_until: event-refined projectile range matches the vacuum
    // analytic range to 1e-9 relative — the whole point of bisecting on dt'
    // instead of accepting an O(dt) crossing (research §7.2/simulate.hpp).
    // ======================================================================
    {
        const double v0 = 20.0, theta = constants::pi / 6.0, g = constants::g_earth;
        Particle p{1.0, Vec3{}, Vec3{v0 * std::cos(theta), v0 * std::sin(theta), 0.0}};
        forces::UniformGravity force{1.0};  // vacuum: no drag
        SimOptions opt;
        opt.t0 = 0.0;
        opt.t_end = 10.0;
        opt.dt = 0.01;  // deliberately coarse: refinement, not fine dt, should carry the accuracy
        opt.method = Method::RK4;
        const Trajectory traj = simulate_until(p, force, opt, events::ground_impact());

        const double range_exact = analytic::range(v0, theta, g);
        const double T_exact = analytic::time_of_flight(v0, theta, g);
        NT_CLOSE(traj.back().r.x, range_exact, 1e-9);
        NT_CLOSE(traj.back().t, T_exact, 1e-9);
        // Sanity: at dt=0.01 the naive "first sign change" answer would be
        // off from the exact range by O(dt * vx) ~ 0.01*17 ~ 0.17 m — three
        // orders of magnitude looser than the 1e-9 this test demands, so
        // this really is checking that refinement happened, not just that
        // RK4 is accurate.
        NT_CHECK(std::fabs(traj.back().r.x - range_exact) < 1e-6);
    }

    // ======================================================================
    // Linear drag, small-b regime and range (Taylor §2.4): RK4 against the
    // rewritten closed forms. At b = 1e-10 the old tau*(1 - e^(-t/tau))
    // form returned 24 km for a 19.5 m flight; RK4 never divides by b, so it
    // is the independent witness that the new answer is right.
    // ======================================================================
    {
        const double g = constants::g_earth;
        const double m = 1.0;
        const double theta = 45.0 * constants::deg;
        const double vx0 = 30.0 * std::cos(theta), vy0 = 30.0 * std::sin(theta);
        for (double b : {0.0, 1e-10, 1e-6, 0.05}) {
            const analytic::LinearDragParams params{m, b, g};
            Particle p{m, Vec3{}, Vec3{vx0, vy0, 0.0}};
            auto force = forces::sum(forces::UniformGravity{m}, forces::LinearDrag{b});
            SimOptions opt;
            opt.t0 = 0.0;
            opt.t_end = 3.0;
            opt.dt = 1e-3;
            opt.method = Method::RK4;
            const Trajectory traj = simulate(p, force, opt);
            const Vec3 r = analytic::linear_drag_position(params, Vec3{}, Vec3{vx0, vy0, 0.0}, 3.0);
            const Vec3 v = analytic::linear_drag_velocity(params, Vec3{vx0, vy0, 0.0}, 3.0);
            NT_CLOSE(traj.back().r.x, r.x, 1e-10);
            NT_CLOSE(traj.back().r.y, r.y, 1e-10);
            NT_CLOSE(traj.back().v.x, v.x, 1e-10);
            NT_CLOSE(traj.back().v.y, v.y, 1e-10);
        }
    }

    // linear_drag_range vs RK4 + event-bisected ground impact (simulate_until
    // with events::ground_impact), several drag strengths.
    {
        const double g = constants::g_earth;
        const double m = 0.145;
        const double v0 = 30.0, theta = 45.0 * constants::deg;
        const double vx0 = v0 * std::cos(theta), vy0 = v0 * std::sin(theta);
        for (double b : {1e-9, 1e-3, 0.0145, 0.1, 0.5}) {
            Particle p{m, Vec3{}, Vec3{vx0, vy0, 0.0}};
            auto force = forces::sum(forces::UniformGravity{m}, forces::LinearDrag{b});
            SimOptions opt;
            opt.t0 = 0.0;
            opt.t_end = 20.0;
            opt.dt = 1e-3;
            opt.method = Method::RK4;
            const Trajectory traj = simulate_until(p, force, opt, events::ground_impact(0.0));
            const double R = analytic::linear_drag_range(analytic::LinearDragParams{m, b, g}, vx0, vy0);
            NT_CHECK(std::fabs(traj.back().r.x - R) <= 1e-9 * R);
            NT_CHECK(std::fabs(traj.back().r.y) < 1e-9);
        }
    }

    // Quadratic drag: horizontal coast x(t) = (m/c) ln(1 + c v0 t/m) and its
    // speed v0/(1 + c v0 t/m), RK4 with no gravity; and the c -> 0 limit.
    {
        const double m = 2.0, v0 = 10.0;
        for (double c : {0.0, 1e-9, 0.05, 0.3}) {
            Particle p{m, Vec3{}, Vec3{v0, 0.0, 0.0}};
            auto force = forces::QuadraticDrag{c};
            SimOptions opt;
            opt.t0 = 0.0;
            opt.t_end = 5.0;
            opt.dt = 1e-3;
            opt.method = Method::RK4;
            const Trajectory traj = simulate(p, force, opt);
            NT_CLOSE(traj.back().r.x, analytic::quadratic_drag_horizontal_distance(m, c, v0, 5.0), 1e-10);
            NT_CLOSE(traj.back().v.x, analytic::quadratic_drag_horizontal_speed(m, c, v0, 5.0), 1e-10);
        }
    }

    return newton_test::report("numeric_vs_analytic");
}
