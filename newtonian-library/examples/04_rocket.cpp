// examples/04_rocket.cpp — vertical rocket launch under gravity (Taylor §3.2).
//
// Physical situation: a single-stage rocket burns propellant at a constant
// rate, launched straight up against uniform gravity. This is the one
// example in the library that does NOT use Particle/simulate() — a burning
// rocket has time-varying mass, and "F = ma" only handles that correctly
// once the thrust term v_ex * burn_rate is included explicitly (see
// systems/rocket.hpp's file header for why d(mv)/dt is the wrong equation
// to integrate).
//
// Analytic expectation (analytic/rocket.hpp): while burning, the closed-form
// velocity v(t) = v_ex ln(m0/(m0 - burn_rate t)) - g t (exact integral of
// the ODE with constant burn rate and constant g) should match
// systems::simulate_rocket's numeric v(t) essentially to roundoff — this is
// RK4 on a smooth scalar ODE at dt << burn time. That includes the very last
// burn sample, t = t_burnout, and the reason it does is a lesson in itself:
// thrust switches off discontinuously at burnout, and an earlier version of
// systems/rocket.hpp decided burning-vs-coasting by comparing t >= t_burnout
// inside the derivative function. The burn phase's final RK4 step evaluates
// its k4 stage exactly AT t_burnout, so that comparison silently handed k4
// the coast-phase (thrust = 0) acceleration — one polluted stage cost the
// whole step ~1e-5 relative error. The fix (see rocket_accel()): the caller
// tells the derivative which phase it is integrating via an explicit
// `burning` flag, so the burn phase's own stages always use the one-sided
// "still firing" limit and the kink lives BETWEEN phases, never inside a
// step. Integrating across a discontinuity in a(t) degrades RK4's order;
// splitting the run at the discontinuity keeps every step smooth.
//
// What to look at: the numeric/analytic velocity table below during the
// burn, then the reported burnout speed and altitude — the instant thrust
// cuts off and the rocket begins an unpowered coast. No closed form is
// offered for altitude or for the coast phase (see
// rocket_velocity_gravity's doc comment), so those are numeric-only.

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
    systems::RocketParams params;
    params.m0 = 1.0e5;          // wet mass [kg]
    params.m_dry = 4.0e4;       // dry mass [kg]
    params.v_ex = 3000.0;       // exhaust speed [m/s]
    params.burn_rate = 1000.0;  // [kg/s]
    params.g = constants::g_earth;

    const double t_burnout = (params.m0 - params.m_dry) / params.burn_rate;
    row("burn duration [s]", t_burnout);

    const double dt = 0.01;
    const double t_end = t_burnout + 30.0;  // burn, then 30s of unpowered coast
    const std::vector<systems::RocketSample> samples = systems::simulate_rocket(params, dt, t_end);

    // The burn phase is re-indexed on its own t0 + i*h grid (systems/rocket.hpp),
    // and dt divides t_burnout exactly here, so sample index i corresponds to
    // t = i * dt. Times are recomputed fresh each iteration (t = 10*k, never
    // accumulated by +=) for the same reason simulate() does (research §7.2).
    std::printf("\n=== velocity during burn: numeric vs analytic ===\n");
    std::printf("%-10s %-16s %-16s %-16s\n", "t [s]", "numeric v", "analytic v", "rel err");
    for (int k = 0; 10.0 * k <= t_burnout + 1e-9; ++k) {
        const double t = 10.0 * k;
        const std::size_t i = static_cast<std::size_t>(std::llround(t / dt));
        const double v_num = samples[i].v;
        const double v_ana = analytic::rocket_velocity_gravity(params.v_ex, params.m0, params.burn_rate,
                                                                 params.g, samples[i].t);
        // At t = 0 both velocities are exactly 0, so a relative error would be
        // 0/0; fall back to a scale of 1 there (absolute error, also 0).
        const double scale = std::fabs(v_ana) > 0.0 ? std::fabs(v_ana) : 1.0;
        const double rel_err = std::fabs(v_num - v_ana) / scale;
        std::printf("%-10.4g %-16.8g %-16.8g %-16.3g\n", samples[i].t, v_num, v_ana, rel_err);
    }
    std::printf(
        "(rel err ~1e-15 at every row INCLUDING t = t_burnout -- the burn\n"
        "phase's final step is kept smooth by the explicit burning flag in\n"
        "systems/rocket.hpp; see the file header for the story.)\n");

    // Locate the burnout sample: systems::simulate_rocket guarantees it is
    // included exactly (the burn phase's last step lands precisely on
    // t_burnout — see the file header on why that matters).
    const auto burnout_it = std::find_if(samples.begin(), samples.end(),
                                          [t_burnout](const systems::RocketSample& s) {
                                              return std::fabs(s.t - t_burnout) < 1e-9;
                                          });

    std::printf("\n=== burnout ===\n");
    if (burnout_it != samples.end()) {
        row("burnout time [s]", burnout_it->t);
        row("numeric burnout v [m/s]", burnout_it->v);
        row("numeric burnout altitude [m]", burnout_it->y);
        row("analytic burnout v [m/s]", analytic::rocket_velocity_gravity(params.v_ex, params.m0,
                                                                            params.burn_rate, params.g, t_burnout));
    }
    row("Tsiolkovsky delta-v (no gravity) [m/s]",
        analytic::tsiolkovsky_delta_v(params.v_ex, params.m0, params.m_dry));
    row("gravity loss at burnout [m/s]", params.g * t_burnout);

    const systems::RocketSample& final_sample = samples.back();
    std::printf("\n=== end of run (%.0f s after burnout) ===\n", t_end - t_burnout);
    row("final time [s]", final_sample.t);
    row("final velocity [m/s]", final_sample.v);
    row("final altitude [m]", final_sample.y);

    return 0;
}
