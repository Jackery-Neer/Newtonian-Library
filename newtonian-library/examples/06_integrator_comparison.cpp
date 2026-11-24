// examples/06_integrator_comparison.cpp — the money demo: four integrators,
// one oscillator, 200 periods (docs/research/numerical-integration.md §2, §4.3).
//
// Physical situation: undamped SHM (forces::Spring, F = -k x) integrated by
// all four methods in newton::integrators, at the SAME step size
// dt = T/100 (100 steps per period — a perfectly reasonable-looking choice),
// run for 200 periods. Energy is conserved by the true dynamics, so
// analysis::energy_drift(traj, spring_U(k)) — the relative deviation of
// numeric E(t) from E(0) — isolates each integrator's own error with no
// analytic trajectory solution required.
//
// Analytic expectation (see each integrator's own file header for the
// derivation):
//   Euler              — energy multiplies by (1+(omega dt)^2) EVERY step:
//                         unconditional exponential growth, unusable for
//                         oscillators at any dt.
//   SemiImplicitEuler  — symplectic: energy is trapped in a band of relative
//                         half-width ~ omega*dt/2, oscillating forever,
//                         never trending.
//   VelocityVerlet     — also symplectic, one order more accurate: band of
//                         relative half-width ~ (omega*dt)^2/4.
//   RK4                — NOT symplectic: a slow, monotone energy decay,
//                         dE/E ~ -(omega*dt)^6/72 per step — tiny at this
//                         dt, but a real secular trend, unlike the other two.
//
// What to look at: build/energy_drift.csv has one relative-drift column per
// method (t, E_euler, E_semi, E_verlet, E_rk4) — Euler's column should
// dwarf the other three by many orders of magnitude; semi/verlet should
// look like flat, bounded noise bands; RK4 should look like a nearly flat
// line with a barely-visible downward slope.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

#include <newton/newton.hpp>

// Leaf translation unit — see 01_free_fall.cpp for why this is fine here.
using namespace newton;

namespace {
void row(const char* label, double value) { std::printf("%-28s % .6g\n", label, value); }
}  // namespace

int main() {
    const double mass = 1.0;
    const double omega0 = 1.0;  // natural frequency [rad/s]
    const double k = mass * omega0 * omega0;
    const double period = 2.0 * constants::pi / omega0;

    const Particle particle{mass, Vec3{1.0, 0.0, 0.0}, Vec3{}};
    const forces::Spring spring{k};

    SimOptions opt;
    opt.t0 = 0.0;
    opt.t_end = 200.0 * period;
    opt.dt = period / 100.0;  // the classic "looks fine" step size
    opt.record_every = 10;    // thin the CSV; first & last samples always kept

    struct MethodCase {
        const char* name;
        Method method;
    };
    const MethodCase cases[4] = {
        {"Euler", Method::Euler},
        {"SemiImplicitEuler", Method::SemiImplicitEuler},
        {"VelocityVerlet", Method::VelocityVerlet},
        {"RK4", Method::RK4},
    };

    std::vector<Trajectory> trajs;
    std::vector<std::vector<double>> drifts;
    std::printf("=== final relative energy error, dt = T/100, 200 periods ===\n");
    for (const MethodCase& c : cases) {
        opt.method = c.method;
        Trajectory traj = simulate(particle, spring, opt);
        std::vector<double> drift = analysis::energy_drift(traj, analysis::spring_U(k));

        row(c.name, drift.back());
        trajs.push_back(std::move(traj));
        drifts.push_back(std::move(drift));
    }
    std::printf("\n");

    std::printf("=== bounded vs. trending: max |drift| over the run ===\n");
    for (std::size_t m = 0; m < 4; ++m) {
        const double max_abs = std::fabs(*std::max_element(
            drifts[m].begin(), drifts[m].end(), [](double a, double b) { return std::fabs(a) < std::fabs(b); }));
        row(cases[m].name, max_abs);
    }
    std::printf(
        "(Euler's max IS its final value -- monotone growth. Semi-implicit\n"
        "and Verlet's max should be well above their final value if the run\n"
        "ends near a zero crossing of the oscillation -- that's the bounded\n"
        "wobble, not drift. RK4's max should be close to its final value too\n"
        "-- monotone decay, just far smaller in magnitude than Euler's.)\n\n");

    // All four trajectories share the same t0/dt/t_end/record_every, so they
    // share the same recorded time grid and the columns line up by index.
    std::vector<std::vector<double>> rows;
    for (std::size_t i = 0; i < trajs[0].size(); ++i) {
        rows.push_back({trajs[0].states[i].t, drifts[0][i], drifts[1][i], drifts[2][i], drifts[3][i]});
    }
    const bool wrote_csv =
        io::write_csv("build/energy_drift.csv", {"t", "E_euler", "E_semi", "E_verlet", "E_rk4"}, rows);
    std::printf("%-28s %s\n", "wrote energy_drift.csv", wrote_csv ? "ok" : "FAILED");

    return 0;
}
