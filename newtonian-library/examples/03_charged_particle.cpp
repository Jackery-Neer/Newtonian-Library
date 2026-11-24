// examples/03_charged_particle.cpp — helix in a uniform magnetic field
// (Taylor §2.5, §2.7).
//
// Physical situation: a charge q moving in a uniform magnetic field B feels
// F = q v x B (forces::MagneticForce) — always perpendicular to v, so it
// does no work: |v| is exactly conserved along the TRUE trajectory, and the
// motion is a helix (uniform circular motion transverse to B, uniform drift
// along B). A bonus section adds a uniform E field (forces::LorentzForce)
// and looks at the resulting E x B drift.
//
// Analytic expectation (analytic/magnetic.hpp): cyclotron_omega/_radius/
// _period give the exact transverse geometry. |v| should be conserved to
// close to machine precision — any deviation is pure integrator error,
// since the true dynamics conserve it exactly (see
// forces/electromagnetic.hpp).
//
// What to look at: build/helix.csv is one pure-B helix (5 gyrations while
// drifting along z) — plot x vs y for the circle, or t vs z for the linear
// drift. The |v| drift printed below should be many orders of magnitude
// smaller than any other number in this file.

#include <algorithm>
#include <cmath>
#include <cstdio>

#include <newton/newton.hpp>

// Leaf translation unit — see 01_free_fall.cpp for why this is fine here.
using namespace newton;

namespace {
void row(const char* label, double value) { std::printf("%-28s % .6g\n", label, value); }
}  // namespace

int main() {
    const double q = 1.0;         // charge [C]
    const double m = 1.0;         // mass [kg]
    const Vec3 B{0.0, 0.0, 1.0};  // uniform field along z [T]

    const double v_perp = 2.0;  // transverse speed [m/s]
    const double v_par = 0.5;   // speed along B [m/s]
    const Particle particle{m, Vec3{}, Vec3{v_perp, 0.0, v_par}};
    const forces::MagneticForce lorentz_b{q, B};

    const double omega = analytic::cyclotron_omega(q, B.z, m);
    const double period = analytic::cyclotron_period(m, q, B.z);
    const double radius_analytic = analytic::cyclotron_radius(m, v_perp, q, B.z);

    std::printf("=== helix: cyclotron motion ===\n");
    row("cyclotron omega [rad/s]", omega);
    row("cyclotron period [s]", period);
    row("analytic radius [m]", radius_analytic);

    const int n_periods = 5;
    SimOptions opt;
    opt.t0 = 0.0;
    opt.t_end = n_periods * period;
    opt.dt = period / 2000.0;  // 2000 steps/gyration: plenty for RK4 to be near-exact
    opt.method = Method::RK4;
    opt.record_every = 2;

    const Trajectory traj = simulate(particle, lorentz_b, opt);

    // |v| conservation: the sharpest possible check on this force, since the
    // TRUE dynamics conserve it exactly (F.v = q(v x B).v = 0 identically) —
    // any nonzero drift here is pure integrator error, not physics.
    const double v0_mag = norm(particle.v);
    double max_speed_drift = 0.0;
    for (const State& s : traj.states) {
        max_speed_drift = std::max(max_speed_drift, std::fabs(norm(s.v) - v0_mag) / v0_mag);
    }
    row("max relative |v| drift", max_speed_drift);

    // Numeric radius: average distance from the transverse guiding center,
    // estimated as the mean (x, y) over an integer number of full gyrations
    // (the circular part averages to exactly the center; see file header).
    double cx = 0.0, cy = 0.0;
    for (const State& s : traj.states) {
        cx += s.r.x;
        cy += s.r.y;
    }
    cx /= static_cast<double>(traj.size());
    cy /= static_cast<double>(traj.size());
    double radius_numeric = 0.0;
    for (const State& s : traj.states) radius_numeric += std::hypot(s.r.x - cx, s.r.y - cy);
    radius_numeric /= static_cast<double>(traj.size());
    row("numeric radius [m]", radius_numeric);

    // Numeric period: mean spacing between successive local maxima of r.x
    // (analysis::find_peaks needs no re-centering — a 3-point local-maximum
    // test is unaffected by a constant offset).
    const analysis::PeakList peaks = analysis::find_peaks(traj);
    double period_numeric = 0.0;
    if (peaks.t.size() >= 2) {
        period_numeric = (peaks.t.back() - peaks.t.front()) / static_cast<double>(peaks.t.size() - 1);
    }
    row("numeric period [s]", period_numeric);
    row("numeric v.z (axial) [m/s]", traj.back().v.z);
    row("analytic v.z (axial) [m/s]", v_par);
    std::printf("\n");

    const bool wrote_csv = io::write_csv("build/helix.csv", traj);
    std::printf("%-28s %s\n", "wrote helix.csv", wrote_csv ? "ok" : "FAILED");
    std::printf("\n");

    // --- bonus: crossed E x B drift (Taylor §2.7) ------------------------------
    std::printf("=== bonus: E x B drift ===\n");
    const Vec3 E{0.0, 1.0, 0.0};  // uniform electric field, perpendicular to B [V/m]
    const forces::LorentzForce lorentz_eb{q, E, B};

    const Particle drifting{m, Vec3{}, Vec3{v_perp, 0.0, 0.0}};
    const SimOptions drift_opt = opt;  // same fine timestep, same duration
    const Trajectory drift_traj = simulate(drifting, lorentz_eb, drift_opt);

    // The guiding-center drift is independent of q and m; averaging v(t)
    // over an integer number of full gyrations cancels the circular part and
    // leaves just the drift (same averaging trick as the radius above).
    Vec3 mean_v{};
    for (const State& s : drift_traj.states) mean_v += s.v;
    mean_v = mean_v / static_cast<double>(drift_traj.size());

    const Vec3 drift_analytic = cross(E, B) / norm2(B);
    row("numeric drift v.x [m/s]", mean_v.x);
    row("numeric drift v.y [m/s]", mean_v.y);
    row("analytic drift v.x [m/s]", drift_analytic.x);
    row("analytic drift v.y [m/s]", drift_analytic.y);

    return 0;
}
