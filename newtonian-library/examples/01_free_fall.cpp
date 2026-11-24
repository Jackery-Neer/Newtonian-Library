// examples/01_free_fall.cpp — free fall and the vacuum projectile (Taylor §1.6).
//
// Physical situation: two of the simplest force problems in the book, both
// under gravity alone (forces::UniformGravity) with no drag.
//   (1) A particle dropped from rest at height h0 — 1D free fall.
//   (2) A particle launched at 45 degrees above the horizontal — the classic
//       vacuum projectile, whose range/apex/time-of-flight are the first
//       closed-form results in the course (Taylor §1.6, analytic/projectile.hpp).
//
// Analytic expectation: RK4 is exact (to roundoff, ~1e-12 relative) for
// constant acceleration at ANY dt — there is no truncation error to speak
// of, because the true solution is a degree-2 polynomial in t and RK4
// integrates such polynomials exactly (see integrators/rk4.hpp and
// docs/research/numerical-integration.md §4.1(a)). So the numeric-vs-analytic
// agreement below is really testing the *event-crossing refinement*
// (simulate_until's bisection, simulate.hpp) as much as the integrator:
// landing on the impact time/position to ~1e-9 relative error, not merely
// O(dt), is the point of that machinery.
//
// What to look at: every numeric/analytic pair below should agree to
// roughly 1e-9 relative or better. If a pair only agrees to O(dt)
// (~1e-4 at the dt used here), the event refinement is broken, not the
// integrator.

#include <cmath>
#include <cstdio>

#include <newton/newton.hpp>

// Leaf translation unit: this file has its own main() and nothing includes
// it, so an unqualified `newton::` is safe here and reads far better in a
// lesson meant to be skimmed — never do this at namespace scope in a header.
using namespace newton;

namespace {
void row(const char* label, double value) { std::printf("%-28s % .6g\n", label, value); }
}  // namespace

int main() {
    // --- (1) free fall from rest --------------------------------------------
    std::printf("=== free fall from rest ===\n");
    const double h0 = 100.0;              // drop height [m]
    const double g = constants::g_earth;  // [m/s^2]

    // Mass is irrelevant to the trajectory here: UniformGravity::operator()
    // returns mass * g, and the integrator immediately divides by that same
    // mass to get acceleration, so a = g regardless of mass (Galileo's
    // equivalence of free-fall acceleration). Picked 2.0 kg just so the
    // cancellation is visible if you print momentum(drop).
    const Particle drop{/*mass=*/2.0, /*r=*/Vec3{0.0, h0, 0.0}, /*v=*/Vec3{}};
    const forces::UniformGravity gravity{drop.mass};

    SimOptions drop_opt;
    drop_opt.t0 = 0.0;
    drop_opt.t_end = 10.0;  // generous upper bound; the event stops us long before this
    drop_opt.dt = 1e-3;
    drop_opt.method = Method::RK4;

    const Trajectory drop_traj = simulate_until(drop, gravity, drop_opt, events::ground_impact(0.0));
    const State& impact = drop_traj.back();

    const double t_impact_analytic = std::sqrt(2.0 * h0 / g);  // 0 = h0 - 1/2 g t^2
    const double v_impact_analytic = -g * t_impact_analytic;   // v = -g t, downward

    row("numeric impact time [s]", impact.t);
    row("analytic impact time [s]", t_impact_analytic);
    row("numeric impact v.y [m/s]", impact.v.y);
    row("analytic impact v.y [m/s]", v_impact_analytic);
    std::printf("\n");

    // --- (2) vacuum projectile at 45 degrees ---------------------------------
    std::printf("=== vacuum projectile launched at 45 deg ===\n");
    const double v0 = 20.0;                        // launch speed [m/s]
    const double theta = 45.0 * constants::deg;     // constants::deg converts deg -> rad

    const Particle shot{/*mass=*/1.0, /*r=*/Vec3{},
                         /*v=*/Vec3{v0 * std::cos(theta), v0 * std::sin(theta), 0.0}};
    const forces::UniformGravity shot_gravity{shot.mass};

    SimOptions shot_opt;
    shot_opt.t0 = 0.0;
    shot_opt.t_end = 10.0;
    shot_opt.dt = 1e-4;
    shot_opt.method = Method::RK4;

    // Demonstrates simulate_until(): the run stops the instant r.y crosses
    // back through 0, refined by bisection to the integrator's own accuracy
    // (simulate.hpp) rather than landing O(dt) short/long the way a naive
    // "stop at the first negative y" would.
    const Trajectory range_traj = simulate_until(shot, shot_gravity, shot_opt, events::ground_impact(0.0));
    const State& landing = range_traj.back();

    // A second run stopped at the apex (v.y == 0) gives the peak height
    // directly, without scanning the recorded samples for a maximum.
    const Trajectory apex_traj = simulate_until(shot, shot_gravity, shot_opt, events::apex());
    const State& apex_state = apex_traj.back();

    const double range_analytic = analytic::range(v0, theta, g);
    const double apex_analytic = analytic::apex_height(v0, theta, g);
    const double tof_analytic = analytic::time_of_flight(v0, theta, g);

    row("numeric range [m]", landing.r.x);
    row("analytic range [m]", range_analytic);
    row("numeric apex height [m]", apex_state.r.y);
    row("analytic apex height [m]", apex_analytic);
    row("numeric time of flight [s]", landing.t);
    row("analytic time of flight [s]", tof_analytic);

    // Energy/speed sanity check that doesn't touch the analytic module at
    // all: on level ground with no drag, impact speed must equal launch
    // speed exactly (same height, conservative force).
    row("launch speed [m/s]", v0);
    row("impact speed [m/s]", norm(landing.v));

    return 0;
}
