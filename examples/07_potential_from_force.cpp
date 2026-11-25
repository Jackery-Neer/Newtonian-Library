// examples/07_potential_from_force.cpp — "Is this force conservative, and if so,
// what is its potential energy?" (Taylor §4.1-4.5, §4.8; the energy ledger of §5.4).
//
// Physical situation: you are handed a force and asked the central question of
// chapter 4. If the work it does depends only on the endpoints, a potential
// energy U exists, F = -grad U, and T + U is conserved. If not, there is no U.
// This lesson walks the question in seven steps, each one using a piece of the
// library you have not met in the earlier examples:
//   (1) the same endpoints along two paths     fields::line_integral   (§4.1-4.2)
//   (2) the local test, curl F = 0             fields::curl            (§4.4)
//   (3) a classification table                 potentials::check_conservative
//   (4) building U from F                      potentials::find_potential
//   (5) the hole field: a curl-free force that is NOT conservative
//   (6) a time-dependent U: E is not constant, but dE/dt = ∂U/∂t       (§4.5)
//   (7) a non-conservative force: gravity + drag, with a closed ledger (§4.1-4.2, §5.4)
//
// Analytic expectation (each line printed below is "numeric vs expected"):
//   (1) Spring{4}: W = -6 J on both paths. Swirl (-y, x, 0): 2 J and 4 J. The
//       difference is the circulation around the triangle, curl_z * area = 2 * 1 m^2.
//   (2) curl of (a y, b x, 0) is (0, 0, b - a); Spring and PointGravity are curl-free.
//   (4) U = 1/2 k r^2, U = m g y, and U = -GMm/r once you pin U(infinity) = 0.
//   (5) the hole field has circulation 2 pi around its axis, so its "U" = -atan2(y, x)
//       exists only in a box that does not enclose the axis.
//   (6) E(t) = 2 - 2 cos t - t sin t for a unit rotating field released from rest.
//   (7) W_drag(1) = -94.684 J, W_drag(2) = -126.762 J, W_drag(impact) = -145.879 J.
//
// What to look at: every "rel err" below should be tiny (1e-9 or better; most
// are near roundoff, and the RK4 runs agree to ~1e-12). One pair is
// deliberately LARGE: in (6) energy_drift is order 10 while energy_balance is
// ~1e-11, which is the point of that section: the energy that "drifts" is real
// physics (∂U/∂t), and the balance accounts for it.
//
// Error column: "rel err" is |numeric - expected| / |expected|; where the
// expected value is 0 it is the absolute error instead.
//
// The new headers below are also reachable through <newton/newton.hpp>; they're
// listed explicitly so you can see where each name comes from (#pragma once
// makes the double include harmless).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include <newton/newton.hpp>
#include <newton/fields/derivatives.hpp>
#include <newton/fields/line_integral.hpp>
#include <newton/potentials/conservative.hpp>
#include <newton/potentials/force_from_potential.hpp>
#include <newton/potentials/potential_from_force.hpp>

// Leaf translation unit — see 01_free_fall.cpp for why this is fine here.
using namespace newton;

namespace {

double rel_err(double numeric, double expected) {
    const double d = std::fabs(numeric - expected);
    return std::fabs(expected) > 1e-9 ? d / std::fabs(expected) : d;
}

// One self-explanatory line: what, numeric, expected, error.
void compare(const char* label, double numeric, double expected) {
    std::printf("  %-34s numeric % .9g   expected % .9g   rel err %.1e\n", label, numeric, expected,
                rel_err(numeric, expected));
}

// Vectors are compared by norm, never component by component.
void compare_vec(const char* label, const Vec3& n, const Vec3& e) {
    const double scale = norm(e) > 1e-9 ? norm(e) : 1.0;
    std::printf("  %-34s numeric (% .6f, % .6f, % .6f)   expected (% .1f, % .1f, % .1f)   err %.1e\n",
                label, n.x, n.y, n.z, e.x, e.y, e.z, norm(n - e) / scale);
}

void heading(const char* text) { std::printf("\n=== %s ===\n", text); }

double max_abs(const std::vector<double>& v) {
    // NaN-aware: std::max(m, fabs(NaN)) silently drops the NaN, so an all-NaN residual
    // would print as 0 ("balanced"). Any NaN makes the result NaN instead.
    double m = 0.0;
    for (double x : v) {
        if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
        m = std::max(m, std::fabs(x));
    }
    return m;
}

// A force that is zero: "no non-conservative force" for energy_balance.
const auto no_force = [](double, const Vec3&, const Vec3&) { return Vec3{}; };

}  // namespace

int main() {
    const double pi = constants::pi;

    // ------------------------------------------------------------------------------------
    // (1) Same endpoints, two paths (Taylor §4.1-4.2).
    //
    // Condition (ii) of §4.2 says W(a -> b) must not depend on the path. We test it
    // directly: fields::line_integral(F, polyline) adds up F . dr along straight legs.
    // C++ note: a plain lambda taking a position is all the integrator needs. It is
    // called directly, with no std::function, so it can be inlined.
    // ------------------------------------------------------------------------------------
    heading("(1) work between the same endpoints along two paths");
    const Vec3 a{1.0, 0.0, 0.0}, b{0.0, 2.0, 0.0};
    const std::vector<Vec3> straight = {a, b};
    const std::vector<Vec3> via_corner = {a, Vec3{1.0, 2.0, 0.0}, b};  // up, then left

    const forces::Spring spring{4.0};
    const auto spring_field = [spring](const Vec3& r) { return spring(0.0, r, Vec3{}); };
    const auto swirl = [](const Vec3& r) { return Vec3{-r.y, r.x, 0.0}; };

    std::printf("  a = (1, 0, 0) m, b = (0, 2, 0) m; path 1 straight, path 2 via (1, 2, 0)\n");
    // Spring: W = -(U(b) - U(a)) = -(1/2*4*4 - 1/2*4*1) = -6 J on every path.
    compare("spring W, straight [J]", fields::line_integral(spring_field, straight), -6.0);
    compare("spring W, via corner [J]", fields::line_integral(spring_field, via_corner), -6.0);
    // Swirl, straight leg: F.dr = 2 ds. Via the corner: 2 J up the first leg + 2 J along the second.
    compare("swirl W, straight [J]", fields::line_integral(swirl, straight), 2.0);
    compare("swirl W, via corner [J]", fields::line_integral(swirl, via_corner), 4.0);
    std::printf("  -> equal work on two paths is consistent with a conservative force; unequal work\n"
                "     on two paths PROVES the swirl is not. Their difference, 4 - 2 = 2 J, is the work\n"
                "     around the closed triangle: curl_z * area = 2 * 1 m^2 (Stokes' theorem).\n");

    // ------------------------------------------------------------------------------------
    // (2) The local test: curl (Taylor §4.4).
    // ------------------------------------------------------------------------------------
    heading("(2) the curl test at a point");
    const Vec3 p{0.7, -0.4, 0.3};
    compare_vec("curl of swirl (-y, x, 0)", fields::curl(swirl, p), Vec3{0.0, 0.0, 2.0});
    compare_vec("curl of spring", fields::curl(spring_field, p), Vec3{});
    std::printf("  -> the swirl has curl 2 z-hat everywhere; the spring has none. Zero curl is the\n"
                "     local test; section (5) shows why it is not the whole story.\n");

    // ------------------------------------------------------------------------------------
    // (3) A classification table (Taylor §4.2, §4.5).
    //
    // check_conservative samples the force on a box (default [-1, 1]^3 m, t in [0, 10] s,
    // |v| components up to 10 m/s): does F depend on v or t? what is the curl? do closed
    // loops do work? It reports evidence, not a proof.
    // ------------------------------------------------------------------------------------
    heading("(3) classification table (check_conservative, default box [-1,1]^3 m)");
    using potentials::to_string;
    const auto ab = [](double A, double B) { return [A, B](const Vec3& r) { return Vec3{A * r.y, B * r.x, 0.0}; }; };
    const auto hole = [](const Vec3& r) {  // (-y, x, 0)/(x^2 + y^2)
        const double s = r.x * r.x + r.y * r.y;
        return Vec3{-r.y / s, r.x / s, 0.0};
    };

    struct Row {
        const char* name;
        potentials::ConservativeReport rep;
        potentials::Verdict expected;
    };
    using V = potentials::Verdict;
    const std::vector<Row> table = {
        {"UniformGravity", potentials::check_conservative(forces::UniformGravity{}), V::Conservative},
        {"Spring{4}", potentials::check_conservative(spring), V::Conservative},
        {"PointGravity (center in box)", potentials::check_conservative(forces::PointGravity{3.0, 2.0}), V::Conservative},
        {"LinearDrag{0.5}", potentials::check_conservative(forces::LinearDrag{0.5}), V::VelocityDependent},
        {"MagneticForce", potentials::check_conservative(forces::MagneticForce{}), V::VelocityDependent},
        {"SinusoidalDrive", potentials::check_conservative(forces::SinusoidalDrive{2.0, 3.0}), V::TimeDependentPotential},
        {"(a y, b x, 0), a = b = 3", potentials::check_conservative(ab(3.0, 3.0)), V::Conservative},
        {"(a y, b x, 0), a = 3, b = 1", potentials::check_conservative(ab(3.0, 1.0)), V::NonzeroCurl},
        {"hole field, box holds z axis", potentials::check_conservative(hole), V::PathDependent},
    };
    std::printf("  %-30s %-26s %-15s %-9s %s\n", "force", "verdict", "v-power", "max curl", "as expected?");
    for (const Row& row : table) {
        // A velocity-dependent force stops at the probe: no curl was computed, so print "-", not 0.
        char curl[16] = "-";
        if (row.rep.verdict != V::VelocityDependent) std::snprintf(curl, sizeof curl, "%.1e", row.rep.max_relative_curl);
        std::printf("  %-30s %-26s %-15s %-9s %s\n", row.name, to_string(row.rep.verdict),
                    to_string(row.rep.velocity_power), curl, row.expected == row.rep.verdict ? "yes" : "NO (!)");
    }
    std::printf("  describe() turns a report into words. For a few of them:\n");
    for (std::size_t i : {std::size_t{3}, std::size_t{4}, std::size_t{7}, std::size_t{8}}) {
        std::printf("    %s: %s\n", table[i].name, potentials::describe(table[i].rep).c_str());
    }
    std::printf("  -> drag and magnetic force are both VelocityDependent (no U exists), but their\n"
                "     v-power labels differ: drag 'dissipative' removes energy, magnetic 'no work'\n"
                "     only turns the particle. The hole field has zero curl yet fails the loop test.\n");

    // ------------------------------------------------------------------------------------
    // (4) find_potential: from F to U (Taylor §4.2, §4.8).
    //
    // U(r) = U_ref - integral from r_ref to r of F . dr, along a path the library picks.
    // The default reference is the origin (so you get mgy and 1/2 k r^2 in textbook form);
    // U_ref is the arbitrary additive constant: only differences of U are physical.
    // ------------------------------------------------------------------------------------
    heading("(4) find_potential: U from F");
    const Vec3 pts[] = {{1.0, 2.0, 2.0}, {0.3, -0.5, 0.8}, {-0.9, 0.4, 0.1}};

    const auto spring_U = potentials::find_potential(spring);
    std::printf("  Spring{4}: U should be 1/2 k |r|^2 (reference: origin, U = 0)\n");
    for (const Vec3& q : pts) {
        compare("U_spring(r) [J]", spring_U.U(q), 0.5 * 4.0 * norm2(q));
    }

    const double m_g = 2.0;
    const auto gravity_U = potentials::find_potential(forces::UniformGravity{m_g});
    std::printf("  UniformGravity{m = 2 kg}: U should be m g y\n");
    for (const Vec3& q : pts) {
        compare("U_gravity(r) [J]", gravity_U.U(q), m_g * constants::g_earth * q.y);
    }

    // Point gravity: the textbook U = -GMm/r is zero at INFINITY, but the line integral
    // can only start from a finite, regular point. The recipe: pick any such r_ref and
    // tell the library the textbook value there, U_ref = -GMm/|r_ref|. Then U(infinity) = 0
    // comes along for free, because the textbook formula and this U differ by a constant
    // that this choice makes zero.
    const double GMm = 6.0;  // GM = 3 m^3/s^2, m = 2 kg
    const forces::PointGravity gravity{3.0, 2.0};
    potentials::CheckOptions big;
    big.region = potentials::Box{Vec3{-3.0, -3.0, -3.0}, Vec3{3.0, 3.0, 3.0}};
    const Vec3 r_ref{2.0, 0.0, 0.0};
    const auto pg = potentials::find_potential(gravity, big, r_ref, -GMm / norm(r_ref));
    std::printf("  PointGravity, r_ref = (2,0,0), U_ref = -GMm/|r_ref| = -3 J (the U(inf) = 0 recipe):\n");
    const Vec3 gpts[] = {{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {1.0, 2.0, 2.0}, {-1.0, 0.0, 0.0}, {0.3, 0.2, -0.1}};
    for (const Vec3& q : gpts) {
        compare("U(r) vs -GMm/|r| [J]", pg.U(q), -GMm / norm(q));
    }
    std::printf("  (U(-1,0,0) is reached by detouring around the singular center: any path inside the\n"
                "   tested box gives the same answer, because the force is conservative there.)\n");

    // The default reference is the origin, which is the SINGULAR point of PointGravity.
    // Rather than fail, find_potential moves the reference to a regular point and says so.
    const auto moved = potentials::find_potential(gravity);
    std::printf("  default reference for PointGravity: reference_moved = %s, r_ref_used = (%.2f, %.2f, %.2f), "
                "reference_ok = %s\n", moved.reference_moved ? "true" : "false", moved.r_ref_used.x,
                moved.r_ref_used.y, moved.r_ref_used.z, moved.reference_ok ? "true" : "false");
    std::printf("  -> the origin is where this force blows up (1/r^2), so U = 0 cannot sit there; the\n"
                "     library picked another point and put U = 0 THERE. Only differences are physical:\n");
    compare("U(1,0,0) - U(0,2,0) [J]", moved.U(Vec3{1.0, 0.0, 0.0}) - moved.U(Vec3{0.0, 2.0, 0.0}),
            -GMm / 1.0 + GMm / 2.0);

    // ------------------------------------------------------------------------------------
    // (5) The hole field: curl F = 0 and still not conservative (Taylor §4.4, simply connected).
    //
    // F = (-y, x, 0)/(x^2 + y^2) has zero curl everywhere except ON the z axis, where it
    // blows up. Around the axis the circulation is 2 pi, not 0. Zero curl implies path
    // independence only on a region with no holes through which a loop can be threaded;
    // Stokes' theorem needs the surface spanning the loop to lie inside the region where
    // curl is zero, and a loop around the axis has no such surface. Same as B around a wire.
    // ------------------------------------------------------------------------------------
    heading("(5) the hole field in two boxes");
    const auto in_axis_box = potentials::find_potential(hole);  // default box: contains the z axis
    std::printf("  box [-1,1]^3 (holds the axis): verdict = %s, has_potential = %s\n",
                to_string(in_axis_box.report.verdict), in_axis_box.has_potential() ? "true" : "false");
    compare("circulation around the axis", fields::circulation(hole, Vec3{}, Vec3{0.0, 0.0, 1.0}, 0.5), 2.0 * pi);

    potentials::CheckOptions side;  // a box to the right of the axis: x in [1, 3]
    side.region = potentials::Box{Vec3{1.0, -1.0, -1.0}, Vec3{3.0, 1.0, 1.0}};
    const auto off_axis = potentials::find_potential(hole, side, Vec3{2.0, 0.0, 0.0}, 0.0);
    std::printf("  box [1,3]x[-1,1]^2 (axis outside): verdict = %s, has_potential = %s\n",
                to_string(off_axis.report.verdict), off_axis.has_potential() ? "true" : "false");
    for (const Vec3& q : {Vec3{2.0, 1.0, 0.0}, Vec3{1.5, -0.7, 0.3}, Vec3{3.0, 1.0, -1.0}}) {
        compare("U(r) vs -atan2(y, x) [J]", off_axis.U(q), -std::atan2(q.y, q.x));
    }
    std::printf("  -> U = -atan2(y, x) is the polar angle (a 'potential' with a branch cut). It is a\n"
                "     fine single-valued U in the right-hand box and impossible around the axis.\n");

    // ------------------------------------------------------------------------------------
    // (6) A time-dependent potential (Taylor §4.5).
    //
    // F = F0 (cos t, sin t, 0), a uniform field rotating at 1 rad/s, F0 = 1 N, m = 1 kg,
    // released from rest at the origin. At every instant it is curl-free, so
    // find_potential returns a potential U(t, r) whose zero is pinned at r_ref = 0 for
    // all t: U = -(x cos t + y sin t). But E = T + U is not constant:
    //     dE/dt = ∂U/∂t at fixed r  (the power supplied by whoever turns the field).
    // Closed form from x = 1 - cos t, y = t - sin t:  E(t) = 2 - 2 cos t - t sin t.
    // ------------------------------------------------------------------------------------
    heading("(6) time-dependent U: the rotating uniform field, from rest");
    const auto rotating = [](double t, const Vec3&) { return Vec3{std::cos(t), std::sin(t), 0.0}; };
    potentials::CheckOptions wide;
    wide.region = potentials::Box{Vec3{-10.0, -10.0, -10.0}, Vec3{10.0, 10.0, 10.0}};
    const auto rot = potentials::find_potential(rotating, wide);
    std::printf("  verdict: %s (depends on t: %s)\n", to_string(rot.report.verdict),
                rot.report.depends_on_time ? "yes" : "no");

    SimOptions rot_opt;
    rot_opt.t_end = 3.5 * pi;
    rot_opt.dt = pi / 2000.0;  // pi/2, pi, 3pi/2, 2pi, 7pi/2 all land exactly on samples
    rot_opt.method = Method::RK4;
    const auto rot_force = [rotating](double t, const Vec3& r, const Vec3&) { return rotating(t, r); };
    const Trajectory rot_traj = simulate(Particle{1.0, Vec3{}, Vec3{}}, rot_force, rot_opt);

    const struct { const char* name; std::size_t index; double fraction; } samples[] = {
        {"E(pi/2)", 1000, 0.5}, {"E(pi)", 2000, 1.0}, {"E(3pi/2)", 3000, 1.5},
        {"E(2pi)", 4000, 2.0}, {"E(7pi/2)", 7000, 3.5}};
    for (const auto& smp : samples) {
        const State& s = rot_traj.states[smp.index];
        const double t = smp.fraction * pi;
        compare(smp.name, analysis::total_energy(1.0, s, rot.U), 2.0 - 2.0 * std::cos(t) - t * std::sin(t));
    }
    const std::vector<double> rot_balance = analysis::energy_balance(rot_traj, rot.U, no_force);
    const std::vector<double> rot_drift = analysis::energy_drift(rot_traj, rot.U);
    std::printf("  max |energy_balance residual| = %.1e   (accounts for ∂U/∂t)\n", max_abs(rot_balance));
    std::printf("  max |energy_drift|            = %.1e   (counts ∂U/∂t as 'error')\n", max_abs(rot_drift));
    std::printf("  -> E(t) is NOT constant, and that is not an integrator failure: it follows the closed form\n"
                "     to ~1e-14. The balance subtracts the integral of ∂U/∂t and is left with ~0. (E(2pi) = 0\n"
                "     because the particle is momentarily at rest again, at r = (0, 2 pi), where F is\n"
                "     perpendicular to r so U = 0; the expected value shows the roundoff of sin(2 pi).)\n");

    // ------------------------------------------------------------------------------------
    // (7) A non-conservative force: gravity + linear drag (Taylor §2.2-2.4, §4.1-4.2).
    //
    // m = 1 kg, b = 0.5 kg/s, launched at 20 m/s and 45 degrees. check_conservative says
    // VelocityDependent / dissipative: no U for the whole force. But F = F(t, r, 0) +
    // [F - F(t, r, 0)] splits it into a part that has a U (gravity) and a part that does
    // not (drag). The split is bookkeeping, not a claim about what the force "really"
    // is: it separates F at v = 0 from the rest. (Drag relative to a moving wind would
    // put a nonzero force into the "at rest" part: the split is by v, not by cause.)
    // The ledger: T + U changes by exactly the work of the drag, W_drag = integral F_drag . v dt.
    // ------------------------------------------------------------------------------------
    heading("(7) gravity + linear drag: classification, split, and the energy ledger");
    const forces::UniformGravity weight{1.0};
    const forces::LinearDrag drag{0.5};
    const auto both = forces::sum(weight, drag);
    const potentials::ConservativeReport whole = potentials::check_conservative(both);
    std::printf("  whole force: %s, v-power %s\n", to_string(whole.verdict), to_string(whole.velocity_power));
    std::printf("    %s\n", potentials::describe(whole).c_str());

    const auto at_rest = potentials::velocity_independent_part(both);
    const auto moving = potentials::velocity_dependent_part(both);
    const potentials::ConservativeReport rest_report = potentials::check_conservative(at_rest);
    const potentials::ConservativeReport moving_report = potentials::check_conservative(moving);
    std::printf("  velocity_independent_part: %s;  velocity_dependent_part: %s, v-power %s\n",
                to_string(rest_report.verdict), to_string(moving_report.verdict),
                to_string(moving_report.velocity_power));
    const auto weight_U = potentials::find_potential(at_rest);
    compare("U of the at-rest part, y = 3 [J]", weight_U.U(Vec3{0.0, 3.0, 0.0}), 1.0 * constants::g_earth * 3.0);

    const double speed = 20.0, c45 = std::cos(pi / 4.0);
    SimOptions flight_opt;
    flight_opt.t_end = 10.0;  // an upper bound: the event stops the run at impact
    flight_opt.dt = 1e-3;
    const Trajectory flight = simulate_until(Particle{1.0, Vec3{}, Vec3{speed * c45, speed * c45, 0.0}}, both,
                                             flight_opt, events::ground_impact(0.0));
    const std::vector<double> W_drag = analysis::work_done(flight, moving);
    const std::vector<double> ledger = analysis::energy_balance(flight, weight_U.U, moving);

    std::printf("  flight: %zu samples, impact at t = %.9f s\n", flight.size(), flight.states.back().t);
    compare("t_impact [s]", flight.states.back().t, 2.4110895210);
    compare("W_drag(1 s) [J]", W_drag[1000], -94.6841394160);
    compare("W_drag(2 s) [J]", W_drag[2000], -126.7624560775);
    compare("W_drag(impact) [J]", W_drag.back(), -145.8787809024);
    const double E0 = analysis::total_energy(1.0, flight.states.front(), weight_U.U);
    const double E_end = analysis::total_energy(1.0, flight.states.back(), weight_U.U);
    compare("E(impact) - E(0) vs W_drag [J]", E_end - E0, W_drag.back());
    std::printf("  max |energy_balance(traj, U_gravity, drag)| = %.1e\n", max_abs(ledger));
    std::printf("  -> W_drag is the work along the REAL path (a line integral of F . dr could not\n"
                "     compute it without knowing the path). It is negative, and equals the energy CHANGE E(t) - E(0) (the energy lost is -W_drag): E is 200 J\n"
                "     at launch, %.3f J at impact, and the ledger T + U - W_drag closes to roundoff (above).\n", E_end);

    return 0;
}
