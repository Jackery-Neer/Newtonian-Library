// examples/02_projectile_drag.cpp — Taylor's baseball (Taylor §2.1-2.4).
//
// Physical situation: the flagship demonstration of Chapter 2. A baseball
// (diameter D = 7 cm, mass m = 0.145 kg) launched at 30 m/s at 45 degrees,
// compared under three force models sharing the same initial conditions:
//   vacuum      — gravity only, exact parabola (analytic/projectile.hpp)
//   linear      — gravity + forces::LinearDrag (Stokes-like, F = -b v)
//   quadratic   — gravity + forces::QuadraticDrag (F = -c|v|v, the
//                 physically correct regime for a ball-sized object moving
//                 through air — see forces/drag.hpp)
//
// Analytic expectation: only the vacuum case has a closed form here
// (analytic::range/time_of_flight); quadratic drag in 2D has NO closed form
// at all (Taylor's own justification for numerical integration, §2.4) — that
// is the whole point of this example. What CAN be checked quantitatively is
// that quadratic drag collapses the range dramatically relative to vacuum
// (Taylor quotes vacuum ~92 m falling to roughly half that for this ball),
// and that the optimal launch angle under quadratic drag sits measurably
// below the vacuum-optimal 45 degrees.
//
// What to look at: build/projectile_drag.csv holds all three trajectories
// (plot x vs y for each `series` value: 0 = vacuum, 1 = linear, 2 =
// quadratic) — the quadratic curve should visibly fall short and steepen
// late in flight, since drag grows with the square of speed. The angle
// sweep at the bottom should report an optimum below 45 degrees.

#include <cmath>
#include <cstdio>
#include <vector>

#include <newton/newton.hpp>

// Leaf translation unit (own main(), nothing includes this file) — see
// 01_free_fall.cpp for why `using namespace newton;` is fine here and never
// in a header.
using namespace newton;

namespace {
void row(const char* label, double value) { std::printf("%-28s % .6g\n", label, value); }

// Launches the baseball under `force` and returns the trajectory up to
// ground impact (y = 0), refined to the integrator's own accuracy by
// simulate_until (simulate.hpp) — shared by the CSV runs and the angle
// sweep below.
template <class Force>
Trajectory launch(double mass, double v0, double theta, Force&& force, double dt) {
    const Particle p{mass, Vec3{}, Vec3{v0 * std::cos(theta), v0 * std::sin(theta), 0.0}};
    SimOptions opt;
    opt.t0 = 0.0;
    opt.t_end = 20.0;  // generous upper bound; ground_impact() stops us well before this
    opt.dt = dt;
    opt.method = Method::RK4;
    return simulate_until(p, force, opt, events::ground_impact(0.0));
}
}  // namespace

int main() {
    const double D = 0.07;                       // baseball diameter [m]
    const double mass = 0.145;                    // baseball mass [kg]
    const double v0 = 30.0;                        // launch speed [m/s]
    const double theta = 45.0 * constants::deg;    // launch angle [rad]
    const double g = constants::g_earth;

    const double b = forces::linear_drag_sphere_air(D);     // [kg/s]
    const double c = forces::quadratic_drag_sphere_air(D);  // [kg/m]

    const forces::UniformGravity gravity{mass};
    const auto& vacuum_force = gravity;
    const auto linear_force = forces::sum(gravity, forces::LinearDrag{b});
    const auto quadratic_force = forces::sum(gravity, forces::QuadraticDrag{c});

    const Trajectory vac_traj = launch(mass, v0, theta, vacuum_force, 1e-4);
    const Trajectory lin_traj = launch(mass, v0, theta, linear_force, 1e-4);
    const Trajectory quad_traj = launch(mass, v0, theta, quadratic_force, 1e-4);

    std::printf("=== range at 45 deg, three drag models ===\n");
    row("vacuum range [m]", vac_traj.back().r.x);
    row("analytic vacuum range [m]", analytic::range(v0, theta, g));
    row("linear-drag range [m]", lin_traj.back().r.x);
    row("quadratic-drag range [m]", quad_traj.back().r.x);
    row("quadratic/vacuum ratio", quad_traj.back().r.x / vac_traj.back().r.x);
    std::printf("\n");

    // One CSV, three series: column `series` tags each row 0/1/2 so a
    // plotting tool can pick out vacuum/linear/quadratic (scripts/plot_csv.py;
    // write_csv's generic column overload just writes whatever numbers it is
    // handed, see io/csv.hpp).
    std::vector<std::vector<double>> rows;
    auto append_series = [&rows](double series_id, const Trajectory& traj) {
        for (const State& s : traj.states) rows.push_back({series_id, s.t, s.r.x, s.r.y});
    };
    append_series(0.0, vac_traj);
    append_series(1.0, lin_traj);
    append_series(2.0, quad_traj);
    const bool wrote_csv = io::write_csv("build/projectile_drag.csv", {"series", "t", "x", "y"}, rows);
    std::printf("%-28s %s\n", "wrote projectile_drag.csv", wrote_csv ? "ok" : "FAILED");
    std::printf("\n");

    // --- optimal launch angle under quadratic drag ----------------------------
    // No closed form exists (see file header), so find it the way Taylor's
    // computer problems ask you to: a brute-force sweep. 0.25 deg resolution
    // is coarser than a golden-section search but reads far more plainly.
    std::printf("=== optimal launch angle, quadratic drag ===\n");
    double best_theta_deg = 0.0;
    double best_range = -1.0;
    for (double theta_deg = 20.0; theta_deg <= 60.0; theta_deg += 0.25) {
        const double th = theta_deg * constants::deg;
        const Trajectory t = launch(mass, v0, th, forces::sum(gravity, forces::QuadraticDrag{c}), 1e-3);
        const double r = t.back().r.x;
        if (r > best_range) {
            best_range = r;
            best_theta_deg = theta_deg;
        }
    }
    row("optimal angle [deg]", best_theta_deg);
    row("optimal range [m]", best_range);
    row("range at 45 deg [m]", quad_traj.back().r.x);
    std::printf(
        "(optimal angle should be BELOW 45 deg: drag grows with the square\n"
        "of speed, and speed is largest early in a high, floaty trajectory --\n"
        "flattening the launch trades a bit of hang time for a lot less drag.)\n");

    return 0;
}
