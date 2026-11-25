// tests/test_potentials.cpp — newton::potentials: is a force conservative
// (Taylor §4.2-4.5), what is its potential energy, and what force does a
// potential energy give?
//
// Every expected value comes from physics: closed forms (mgy, ½k r², -F·r,
// -GMm/r, Morse, Yukawa, atan2 for the hole field), the exact energy of a
// particle in a rotating uniform field, and hand-computed curls and loops.
// Vectors are compared by norm, never component by component. Run with
// `--table` to print the verdict table for the whole battery:
//     ./build/tests/test_potentials --table
#include "test_util.hpp"        // quotes: sits next to this file
#include <newton/newton.hpp>    // angle brackets: reached via -I include
#include <newton/potentials/conservative.hpp>
#include <newton/potentials/force_from_potential.hpp>
#include <newton/potentials/potential_from_force.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <limits>
#include <random>
#include <string>
#include <vector>

using newton::Particle;
using newton::SimOptions;
using newton::Trajectory;
using newton::Vec3;
using newton::constants::pi;
using namespace newton::potentials;  // fine in a test .cpp; never in a header

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Vectors are compared by NORM: a per-component test falsely fails the small
// components of a large vector.
double rel_err(const Vec3& a, const Vec3& b) {
    return norm(a - b) / std::max(norm(b), 1e-300);
}

// Largest |x|, NaN-aware: std::max(m, fabs(NaN)) silently drops the NaN, so an
// all-NaN residual would read as 0 ("balanced"). Here any NaN makes the result NaN,
// and `NaN < tol` is false, so the check fails loudly.
double max_abs(const std::vector<double>& v) {
    double m = 0.0;
    for (const double x : v) {
        if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
        m = std::max(m, std::fabs(x));
    }
    return m;
}

// A straight singular line through p along unit vector d: the field of a thin wire,
// d × q/|q|² with q the part of (r - p) perpendicular to d. Curl-free off the line,
// circulation 2π around it.
auto wire_field(const Vec3& p, const Vec3& d) {
    return [p, d](const Vec3& r) {
        const Vec3 u = r - p;
        const Vec3 q = u - dot(u, d) * d;
        return cross(d, q) / norm2(q);
    };
}

// --- the homework / battery fields (research-physics.md §4) ----------------------

// Position-only fields, written with the SHORT signature Vec3(const Vec3&).
const auto hole = [](const Vec3& r) {  // (-y, x, 0)/(x²+y²): curl-free off the z axis, ∮ = 2π around it
    const double s = r.x * r.x + r.y * r.y;
    return Vec3{-r.y / s, r.x / s, 0.0};
};
auto ab_field(double a, double b) {  // (a y, b x, 0): curl (0, 0, b - a)
    return [a, b](const Vec3& r) { return Vec3{a * r.y, b * r.x, 0.0}; };
}
const auto classic = [](const Vec3& r) {  // (2xy + z³, x², 3xz²): the textbook "is it conservative?" field
    return Vec3{2.0 * r.x * r.y + r.z * r.z * r.z, r.x * r.x, 3.0 * r.x * r.z * r.z};
};
const auto xy_yz_zx = [](const Vec3& r) { return Vec3{r.x * r.y, r.y * r.z, r.z * r.x}; };
const auto x_times_r = [](const Vec3& r) { return r.x * r; };  // central but not spherically symmetric (§4.8)
const auto thick_wire = [](const Vec3& r) {  // the hole field with its singularity smoothed over w = 1e-4 m
    const double s = r.x * r.x + r.y * r.y + 1e-8;
    return Vec3{-r.y / s, r.x / s, 0.0};
};
const auto coulomb2 = [](const Vec3& r) {  // kqQ = 2 (repulsive): U = 2/r
    const double n = norm(r);
    return (2.0 / (n * n * n)) * r;
};
const auto yukawa = [](const Vec3& r) {  // U = -(k/r) e^{-r/a}, k = a = 1: F = -e^{-r}(1/r² + 1/r) r̂
    const double n = norm(r);
    return (-std::exp(-n) * (1.0 / (n * n) + 1.0 / n) / n) * r;
};
double morse_U(const Vec3& r) {  // A = 1, r0 = 1, S = 0.5: U = A[e^{-2(r-r0)/S} - 2e^{-(r-r0)/S}]
    const double e = std::exp(-(norm(r) - 1.0) / 0.5);
    return e * e - 2.0 * e;
}
const auto morse = [](const Vec3& r) {  // F = (2A/S)(e^{-2x} - e^{-x}) r̂, x = (r - r0)/S
    const double e = std::exp(-(norm(r) - 1.0) / 0.5);
    return (4.0 * (e * e - e)) * normalized(r);
};
const auto lennard_jones = [](const Vec3& r) {  // U0 = r0 = 1: F = (12/r)[(1/r)^12 - (1/r)^6] r̂
    const double n = norm(r);
    const double s6 = std::pow(1.0 / n, 6.0);
    return (12.0 / n * (s6 * s6 - s6)) * (r / n);
};

// Time-dependent fields, written with the (t, r) signature.
const auto pulsing_spring = [](double t, const Vec3& r) { return -4.0 * (1.0 + 0.5 * std::sin(t)) * r; };
const auto moving_support = [](double t, const Vec3& r) {  // -k(r - A cos(wt) x̂), k = 4, A = 0.5, w = 3
    return -4.0 * (r - Vec3{0.5 * std::cos(3.0 * t), 0.0, 0.0});
};
const auto rotating_field = [](double t, const Vec3&) { return Vec3{std::cos(t), std::sin(t), 0.0}; };  // F0 = 1 N, w = 1 rad/s
const auto decaying_coulomb = [](double t, const Vec3& r) {  // kqQ0 e^{-t/tau} r̂/r², kqQ0 = 1, tau = 2
    const double n = norm(r);
    return (std::exp(-t / 2.0) / (n * n * n)) * r;
};
const auto induced_swirl = [](double t, const Vec3& r) { return std::sin(t) * Vec3{-r.y, r.x, 0.0}; };

// A battery row: name, expected verdict, expected velocity power, expected depends_on_time, and
// how to run the classifier (a lambda so each row can pick its own box and its own force shape).
struct Row {
    const char* name;
    Verdict verdict;
    VelocityPower power;
    bool on_time;
    std::function<ConservativeReport()> run;
};

template <class F>
std::function<ConservativeReport()> classify(F f, CheckOptions opt = {}) {
    return [f, opt]() { return check_conservative(f, opt); };
}

CheckOptions box(const Vec3& lo, const Vec3& hi) {
    CheckOptions o;
    o.region = Box{lo, hi};
    return o;
}

}  // namespace

int main(int argc, char** argv) {
    const bool print_table = (argc > 1 && std::strcmp(argv[1], "--table") == 0);
    namespace forces = newton::forces;
    namespace analysis = newton::analysis;
    const Vec3 origin{};

    // --- the deterministic generator ---------------------------------------------
    {
        detail::Rng a(12345);
        std::mt19937_64 raw(12345);
        bool same = true, in_range = true;
        for (int i = 0; i < 100; ++i) {
            const double u = a.uniform();
            const double expect = static_cast<double>(raw() >> 11) * 0x1.0p-53;
            same = same && (u == expect);
            in_range = in_range && u >= 0.0 && u < 1.0;
        }
        NT_CHECK(same);
        NT_CHECK(in_range);
        // The standard fixes mt19937_64: its 10000th output (default seed) is a published value.
        std::mt19937_64 def;
        def.discard(9999);
        NT_CHECK(def() == 9981545732273789042ULL);
    }

    // --- golden-ratio time slices: never "nice" instants ----------------------------
    {
        CheckOptions o;
        o.t_min = 0.0;
        o.t_max = 2.0 * pi;
        const std::vector<double> ts = detail::slice_times(o, true);
        NT_CHECK(ts.size() == 5);
        NT_CLOSE(ts[0], (0.0 + 0.6180339887498949) / 5.0 * 2.0 * pi, 1e-14);
        NT_CLOSE(ts[4], (4.0 + 0.6180339887498949) / 5.0 * 2.0 * pi, 1e-14);
        bool away = true;
        for (const double t : ts) away = away && std::fabs(t - o.t_min) > 0.1;
        NT_CHECK(away);
        const std::vector<double> one = detail::slice_times(o, false);
        NT_CHECK(one.size() == 1 && one[0] == o.t_min);
    }

    // --- probe_dependence for every library force (research-physics.md §4) -------------
    {
        auto probe = [](const auto& F) { return probe_dependence(F); };
        auto expect = [&](const Dependence& d, bool t, bool v) { return d.on_time == t && d.on_velocity == v; };
        NT_CHECK(expect(probe(forces::UniformGravity{2.0}), false, false));
        NT_CHECK(expect(probe(forces::PointGravity{3.0, 2.0}), false, false));
        NT_CHECK(expect(probe(forces::Spring{4.0}), false, false));
        NT_CHECK(expect(probe(forces::AnisotropicSpring{1, 2, 3}), false, false));
        NT_CHECK(expect(probe(forces::ConstantForce{{1, 2, 3}}), false, false));
        NT_CHECK(expect(probe(forces::LorentzForce{2.0, {0, 0, 5}, {}}), false, false));  // B = 0: position-only in behaviour
        NT_CHECK(expect(probe(forces::MagneticForce{1.0, {}}), false, false));            // B = 0: the zero force
        NT_CHECK(expect(probe(forces::LinearDrag{0.5}), false, true));
        NT_CHECK(expect(probe(forces::QuadraticDrag{0.5}), false, true));
        NT_CHECK(expect(probe(forces::MagneticForce{}), false, true));
        NT_CHECK(expect(probe(forces::LorentzForce{2.0, {0, 0, 5}, {0, 0, 1}}), false, true));
        NT_CHECK(expect(probe(forces::SinusoidalDrive{2.0, 3.0, 0.5}), true, false));
        NT_CHECK(expect(probe(forces::sum(forces::Spring{4.0}, forces::LinearDrag{0.4}, forces::SinusoidalDrive{2.0, 3.0})),
                        true, true));
        NT_CHECK(expect(probe(hole), false, false));            // short signature: position-only by type
        NT_CHECK(expect(probe(pulsing_spring), true, false));   // (t, r) signature: t only
        // A step that switches on at t = 5 s is inside the window; at t = 100 s it is invisible (documented limit).
        auto step5 = [](double t, const Vec3&) { return Vec3{t > 5.0 ? 1.0 : 0.0, 0.0, 0.0}; };
        auto step100 = [](double t, const Vec3&) { return Vec3{t > 100.0 ? 1.0 : 0.0, 0.0, 0.0}; };
        NT_CHECK(probe(step5).on_time);
        NT_CHECK(!probe(step100).on_time);
        // v/|v| is NaN at v = 0 (half the samples compare against v = 0): dependence is still found.
        auto unit_v = [](double, const Vec3&, const Vec3& v) { return -v / norm(v); };
        NT_CHECK(probe(unit_v).on_velocity);
        // A drag that only acts faster than 12 m/s: the cube corners (|v| up to 17) reach it.
        auto fast_drag = [](double, const Vec3&, const Vec3& v) { return norm(v) > 12.0 ? -v : Vec3{}; };
        NT_CHECK(probe(fast_drag).on_velocity);
    }

    // --- the verdict battery ---------------------------------------------------------------
    const CheckOptions hole_b = box({1, -1, -1}, {3, 1, 1});
    const CheckOptions lj_box = box({0.5, -1, -1}, {2.5, 1, 1});
    const VelocityPower NA = VelocityPower::NotApplicable;
    const Verdict C = Verdict::Conservative, T = Verdict::TimeDependentPotential, V = Verdict::VelocityDependent,
                  N = Verdict::NonzeroCurl, P = Verdict::PathDependent;
    const std::vector<Row> rows = {
        {"UniformGravity{2}", C, NA, false, classify(forces::UniformGravity{2.0})},
        {"PointGravity, centre in box", C, NA, false, classify(forces::PointGravity{3.0, 2.0})},
        {"PointGravity, softened 0.5", C, NA, false, classify(forces::PointGravity{3.0, 2.0, {}, 0.5})},
        {"PointGravity, centre outside", C, NA, false, classify(forces::PointGravity{1.0, 1.0, {3, 0, 0}})},
        {"LinearDrag{0.5}", V, VelocityPower::Dissipative, false, classify(forces::LinearDrag{0.5})},
        {"QuadraticDrag{0.5}", V, VelocityPower::Dissipative, false, classify(forces::QuadraticDrag{0.5})},
        {"Spring{4}", C, NA, false, classify(forces::Spring{4.0})},
        {"Spring{4, (0.2,0.1,0), L=0.5}", C, NA, false, classify(forces::Spring{4.0, {0.2, 0.1, 0.0}, 0.5})},
        {"AnisotropicSpring{1,2,3}", C, NA, false, classify(forces::AnisotropicSpring{1, 2, 3})},
        {"SinusoidalDrive{2,3,0.5}", T, NA, true, classify(forces::SinusoidalDrive{2.0, 3.0, 0.5})},
        {"ConstantForce{1,2,3}", C, NA, false, classify(forces::ConstantForce{{1, 2, 3}})},
        {"MagneticForce{1,(0,0,1)}", V, VelocityPower::NoWork, false, classify(forces::MagneticForce{})},
        {"LorentzForce, B = 0", C, NA, false, classify(forces::LorentzForce{2.0, {0, 0, 5}, {}})},
        {"LorentzForce, B != 0", V, VelocityPower::NoWork, false, classify(forces::LorentzForce{2.0, {0, 0, 5}, {0, 0, 1}})},
        {"sum(UniformGravity, Spring)", C, NA, false, classify(forces::sum(forces::UniformGravity{2.0}, forces::Spring{4.0}))},
        {"sum(UniformGravity, LinearDrag)", V, VelocityPower::Dissipative, false,
         classify(forces::sum(forces::UniformGravity{2.0}, forces::LinearDrag{0.5}))},
        {"sum(Spring, LinearDrag, Drive)", V, VelocityPower::Dissipative, true,
         classify(forces::sum(forces::Spring{4.0}, forces::LinearDrag{0.4}, forces::SinusoidalDrive{2.0, 3.0}))},
        {"sum(Spring, MagneticForce)", V, VelocityPower::NoWork, false, classify(forces::sum(forces::Spring{4.0}, forces::MagneticForce{}))},
        {"sum(UG, Spring{1e-9})", C, NA, false, classify(forces::sum(forces::UniformGravity{2.0}, forces::Spring{1e-9}))},
        {"(3y, 3x, 0)", C, NA, false, classify(ab_field(3.0, 3.0))},
        {"(3y, 1x, 0)", N, NA, false, classify(ab_field(3.0, 1.0))},
        {"(y, (1+1e-4)x, 0)", N, NA, false, classify(ab_field(1.0, 1.0 + 1e-4))},
        {"(2xy+z^3, x^2, 3xz^2)", C, NA, false, classify(classic)},
        {"(xy, yz, zx)", N, NA, false, classify(xy_yz_zx)},
        {"Coulomb kqQ=2", C, NA, false, classify(coulomb2)},
        {"Yukawa k=a=1", C, NA, false, classify(yukawa)},
        {"Morse A=1 r0=1 S=.5", C, NA, false, classify(morse)},
        {"Lennard-Jones, box x in [.5,2.5]", C, NA, false, classify(lennard_jones, lj_box)},
        {"x r (central, not symmetric)", N, NA, false, classify(x_times_r)},
        {"hole, box contains the axis", P, NA, false, classify(hole)},
        {"hole, box x in [1,3]", C, NA, false, classify(hole, hole_b)},
        {"zero force", C, NA, false, classify([](double, const Vec3&, const Vec3&) { return Vec3{}; })},
        {"pulsing spring", T, NA, true, classify(pulsing_spring)},
        {"moving support", T, NA, true, classify(moving_support)},
        {"rotating uniform field", T, NA, true, classify(rotating_field)},
        {"decaying Coulomb", T, NA, true, classify(decaying_coulomb)},
        {"induced swirl sin(t)(-y,x,0)", N, NA, true, classify(induced_swirl)},
        {"hole + (1e6,0,0)", P, NA, false, classify([](const Vec3& r) { return hole(r) + Vec3{1e6, 0.0, 0.0}; })},
    };
    if (print_table) {
        std::printf("%-34s %-26s %-12s %-9s %-10s %-9s %-9s %-7s %-6s %7s\n", "force", "verdict", "power", "max curl",
                    "loop ratio", "samples", "skipped", "unconf", "lskip", "ms");
    }
    for (const Row& row : rows) {
        const auto t0 = std::chrono::steady_clock::now();
        const ConservativeReport rep = row.run();
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        if (print_table) {
            std::printf("%-34s %-26s %-12s %-9.1e %-10.1e %-9d %-9d %-7d %-6d %7.2f\n", row.name, to_string(rep.verdict),
                        to_string(rep.velocity_power), rep.max_relative_curl, rep.max_loop_ratio, rep.samples_tested,
                        rep.skipped_samples, rep.unconfirmed_curl_candidates, rep.loops_skipped, ms);
        }
        if (rep.verdict != row.verdict || rep.velocity_power != row.power || rep.depends_on_time != row.on_time) {
            std::printf("  BATTERY MISMATCH: %s -> %s / %s / t=%d\n", row.name, to_string(rep.verdict),
                        to_string(rep.velocity_power), rep.depends_on_time ? 1 : 0);
        }
        newton_test::record(rep.verdict == row.verdict, row.name, __FILE__, __LINE__);
        newton_test::record(rep.velocity_power == row.power, row.name, __FILE__, __LINE__);
        newton_test::record(rep.depends_on_time == row.on_time, row.name, __FILE__, __LINE__);
        NT_CHECK(rep.has_potential() == (row.verdict == C || row.verdict == T));
        NT_CHECK(rep.is_conservative() == (row.verdict == C));
    }

    // --- thick-wire vortex (w = 1e-4 m): physically nonconservative, caught by the curl OR the loops ----
    // research-physics.md §4 row 20 expects PathDependent (the curl at (0.5,0,0) is only 4e-8). But 256
    // random samples include a few with rho < 0.1 m, where the curl 2w²/(rho²+w²)² is large and REAL
    // (confirmed at h/2), so with the default seed the classifier reports NonzeroCurl (relative curl
    // 2.3e-6). Both verdicts mean "no potential"; see the report to the parent.
    {
        const ConservativeReport r = check_conservative(thick_wire);
        NT_CHECK(r.verdict == N || r.verdict == P);
        NT_CHECK(!r.has_potential());
        NT_CLOSE(newton::fields::circulation(thick_wire, origin, Vec3{0, 0, 1}, 0.5), 2.0 * pi * (1.0 - 4e-8), 1e-8);
        // A box far from the axis sees a curl too weak to matter: conservative. (Close to the axis,
        // e.g. x in [0.3, 3], the loops DO see the real circulation, ~2e-8 of |F||dr| > loop_tol.)
        NT_CHECK(check_conservative(thick_wire, box({1.5, -1, -1}, {3.5, 1, 1})).verdict == C);
    }

    // --- evidence in the report -----------------------------------------------------------------
    {
        const ConservativeReport r = check_conservative(hole);
        // The big square loop around the axis: ∮F·dr = 2π, ∮|F||dr| = 8 asinh(1) a / a  (|F| = 1/ρ on a square of half-side a)
        NT_CLOSE(r.max_loop_ratio, 2.0 * pi / (8.0 * std::asinh(1.0)), 1e-9);
        NT_CHECK(r.max_relative_curl < 1e-6);   // the curl cannot see the hole
        const ConservativeReport s = check_conservative(induced_swirl);
        NT_CHECK(s.worst_curl_time >= 0.0 && s.worst_curl_time <= 10.0);
        NT_CHECK(std::fabs(s.worst_curl_point.x) <= 1.0 && std::fabs(s.worst_curl_point.y) <= 1.0);
        NT_CHECK(s.max_relative_curl > 1e-6);
        // Even in a window whose endpoints are all zeros of sin(t): golden-ratio slices still see the curl.
        CheckOptions o;
        o.t_min = 0.0;
        o.t_max = 2.0 * pi;
        NT_CHECK(check_conservative(induced_swirl, o).verdict == N);
        const ConservativeReport g = check_conservative(forces::PointGravity{3.0, 2.0});
        NT_CHECK(g.loops_tested == 351);  // 3 axes × 3 planes × (1 + 4 + 9 + 25) cells
        NT_CHECK(g.samples_tested + g.skipped_samples == 256);
        NT_CHECK(g.max_loop_ratio < 1e-8);
    }

    // --- inconclusive: a force that is not finite --------------------------------------------------
    {
        auto nan3 = [](double, const Vec3&, const Vec3&) { return Vec3{kNaN, kNaN, kNaN}; };
        const ConservativeReport a = check_conservative(nan3);
        NT_CHECK(a.verdict == Verdict::Inconclusive);
        NT_CHECK(a.skipped_samples == 64);
        const ConservativeReport b = check_conservative([](const Vec3&) { return Vec3{kNaN, 0.0, 0.0}; });
        NT_CHECK(b.verdict == Verdict::Inconclusive);
        auto mostly_nan = [](const Vec3& r) { return r.x > -0.5 ? Vec3{kNaN, kNaN, kNaN} : Vec3{-r.x, 0.0, 0.0}; };
        NT_CHECK(check_conservative(mostly_nan).verdict == Verdict::Inconclusive);
        // A small bad patch (about 10% of the box) is skipped and counted, and the verdict stands.
        auto patchy = [](const Vec3& r) { return r.x > 0.8 ? Vec3{kNaN, 0.0, 0.0} : Vec3{-4.0 * r.x, -4.0 * r.y, -4.0 * r.z}; };
        const ConservativeReport c = check_conservative(patchy);
        NT_CHECK(c.verdict == Verdict::Conservative);
        NT_CHECK(c.skipped_samples > 0);
        NT_CHECK(c.loops_skipped > 0);
    }

    // --- describe / to_string -----------------------------------------------------------------
    {
        NT_CHECK(std::string(to_string(Verdict::PathDependent)) == "path-dependent");
        NT_CHECK(std::string(to_string(VelocityPower::NoWork)) == "no work");
        for (const Row& row : rows) {
            const std::string d = describe(row.run());
            NT_CHECK(d.size() > 40);
        }
        NT_CHECK(describe(check_conservative(forces::LinearDrag{0.5})).find("opposes v") != std::string::npos);
        NT_CHECK(describe(check_conservative(forces::MagneticForce{})).find("no work") != std::string::npos);
        NT_CHECK(describe(check_conservative(induced_swirl)).find("path") != std::string::npos);
        // §4.5: the rate is the PARTIAL derivative at fixed r, not dU/dt along the motion.
        const std::string td = describe(check_conservative(pulsing_spring));
        NT_CHECK(td.find("partial U/partial t") != std::string::npos);
        NT_CHECK(td.find("dU/dt") == std::string::npos);
        // PathDependent allows both explanations
        const std::string pd = describe(check_conservative(hole));
        NT_CHECK(pd.find("singular on a line") != std::string::npos && pd.find("too weak or too concentrated") != std::string::npos);
        // a static field's curl sentence prints no time
        NT_CHECK(describe(check_conservative(xy_yz_zx)).find("[t =") == std::string::npos);
        NT_CHECK(describe(check_conservative(induced_swirl)).find("[t =") != std::string::npos);
    }

    // --- potentials vs the textbook ----------------------------------------------------------------
    {
        const double mg = 2.0 * newton::constants::g_earth;
        const auto g = find_potential(forces::UniformGravity{2.0});
        NT_CHECK(g.has_potential() && !g.reference_moved && g.r_ref_used == origin);
        const Vec3 a{1, 2, 3}, b{-1, -1, 0};
        NT_CLOSE(g.U(a), mg * a.y, 1e-12);                 // U = mgy, zero at y = 0
        NT_CLOSE(g.U(a) - g.U(b), 58.8399, 1e-9);          // the research value 3·19.6133
        const auto sp = find_potential(forces::Spring{4.0});
        NT_CLOSE(sp.U(Vec3{1, 2, 2}), 2.0 * 9.0, 1e-12);   // ½k r² = 2 r²
        NT_CLOSE(sp.U(Vec3{1, 2, 2}) - sp.U(Vec3{0, 0, 1}), 16.0, 1e-12);
        const auto cf = find_potential(forces::ConstantForce{{1, 2, 3}});
        NT_CLOSE(cf.U(Vec3{1, 1, 1}), -6.0, 1e-12);        // U = -F·r
        NT_CLOSE(cf.U(Vec3{1, 1, 1}) - cf.U(origin), -6.0, 1e-12);
        const auto lz = find_potential(forces::LorentzForce{2.0, {0, 0, 5}, {}});
        NT_CLOSE(lz.U(Vec3{0, 0, 1}), -10.0, 1e-12);       // U = -qE·r
        // anchored spring: U = 2(d - 0.5)² up to a constant
        const auto as = find_potential(forces::Spring{4.0, {1, 0, 0}, 0.5}, box({-3, -3, -3}, {3, 3, 3}));
        NT_CHECK(as.has_potential());
        NT_CLOSE(as.U(Vec3{3, 0, 0}) - as.U(Vec3{1, 1, 0}), 4.0, 1e-9);
        NT_CLOSE(as.U(Vec3{3, 0, 0}) - as.U(Vec3{1, 0, 0}), 4.0, 1e-9);
        // the classic textbook field: U = -(x²y + xz³)
        const auto cl = find_potential(classic);
        NT_CLOSE(cl.U(Vec3{1, 1, 1}) - cl.U(origin), -2.0, 1e-10);
    }
    {
        // PointGravity with an explicit reference and the U(∞) = 0 recipe U_ref = -GMm/|r_ref|.
        const CheckOptions big = box({-3, -3, -3}, {3, 3, 3});
        const auto pg = find_potential(forces::PointGravity{3.0, 2.0}, big, Vec3{2, 0, 0}, -6.0 / 2.0);
        NT_CHECK(pg.has_potential() && !pg.reference_moved);
        const Vec3 pts[] = {{1, 0, 0}, {0, 2, 0}, {1, 2, 2}, {0, 0, 1}, {-2, 0, 0}, {-1, 0, 0}, {0.3, 0.2, -0.1}};
        for (const Vec3& p : pts) {
            NT_CLOSE(pg.U(p), -6.0 / norm(p), 1e-9);  // includes paths straight through the centre (detours)
        }
        NT_CLOSE(pg.U(Vec3{1, 0, 0}) - pg.U(Vec3{0, 2, 0}), -3.0, 1e-9);
        NT_CLOSE(pg.U(Vec3{1, 2, 2}) - pg.U(Vec3{0, 0, 1}), 4.0, 1e-9);
        NT_CLOSE(pg.U(Vec3{1, 0, 0}) - pg.U(Vec3{-1, 0, 0}), 0.0, 1e-9);
        NT_CHECK(std::isnan(pg.U(origin)));  // U = -infinity at the centre: NaN, never a wrong number

        // Default reference = origin = the singular point: find_potential moves it.
        const auto mv = find_potential(forces::PointGravity{3.0, 2.0});
        NT_CHECK(mv.report.has_potential());
        NT_CHECK(mv.reference_moved && mv.reference_ok && mv.has_potential());
        NT_CHECK(mv.r_ref_used != origin);
        NT_CLOSE(mv.U(mv.r_ref_used), 0.0, 1e-12);
        NT_CLOSE(mv.U(Vec3{1, 0, 0}) - mv.U(Vec3{0, 2, 0}), -3.0, 1e-9);
        NT_CLOSE(mv.U(Vec3{1, 2, 2}) - mv.U(Vec3{0, 0, 1}), 4.0, 1e-9);
        NT_CHECK(describe(mv).find("moved") != std::string::npos);

        // A caller's singular reference is never moved.
        const auto bad = find_potential(forces::PointGravity{3.0, 2.0}, CheckOptions{}, origin);
        NT_CHECK(bad.report.has_potential());
        NT_CHECK(!bad.reference_ok && !bad.reference_moved && !bad.has_potential());
        NT_CHECK(bad.r_ref_used == origin);
        NT_CHECK(describe(bad).find("not usable") != std::string::npos);

        // Softened gravity is regular at the centre, so the default reference works: Plummer U.
        const auto sg = find_potential(forces::PointGravity{3.0, 2.0, {}, 0.5});
        NT_CHECK(sg.has_potential() && !sg.reference_moved);
        NT_CLOSE(sg.U(Vec3{1, 0, 0}) - sg.U(origin), 6.0 / 0.5 - 6.0 / std::sqrt(1.25), 1e-10);  // 6.633437 J (negated: U rises outward)
    }
    {
        // Morse: U(∞) = 0 via r_ref = (r0,0,0), U_ref = -A.
        const CheckOptions big = box({-3, -3, -3}, {3, 3, 3});
        const auto m = find_potential(morse, big, Vec3{1, 0, 0}, -1.0);
        NT_CHECK(m.has_potential());
        const Vec3 pts[] = {{0.8, 0, 0}, {1, 0, 0}, {1.5, 0, 0}, {3, 0, 0}, {0, 2, 0}, {1, 1, 1}, {-1.5, 0, 0}};
        for (const Vec3& p : pts) NT_CLOSE(m.U(p), morse_U(p), 1e-9);
        NT_CLOSE(m.U(Vec3{1.5, 0, 0}) - m.U(Vec3{1, 0, 0}), 0.3995764, 1e-7);

        // Yukawa: U = -e^{-r}/r.
        const auto y = find_potential(yukawa, big, Vec3{1, 0, 0}, -std::exp(-1.0));
        NT_CHECK(y.has_potential());
        NT_CLOSE(y.U(Vec3{1, 0, 0}) - y.U(Vec3{0, 2, 0}), -0.3002118, 1e-6);
        for (const Vec3& p : {Vec3{0, 2, 0}, Vec3{0.5, 0.5, 0.5}, Vec3{-1, 1, 2}}) {
            NT_CLOSE(y.U(p), -std::exp(-norm(p)) / norm(p), 1e-9);
        }
        // Coulomb with the same recipe: U = 2/r.
        const auto co = find_potential(coulomb2, big, Vec3{1, 0, 0}, 2.0);
        NT_CLOSE(co.U(Vec3{1, 0, 0}) - co.U(Vec3{0, 0, 2}), 1.0, 1e-9);
    }
    {
        // The hole field: PathDependent in [-1,1]^3, conservative with U = -atan2(y, x) in [1,3]x[-1,1]².
        const auto h1 = find_potential(hole);
        NT_CHECK(h1.report.verdict == P && !h1.has_potential() && !h1.reference_ok);
        const auto h2 = find_potential(hole, hole_b, Vec3{2, 0, 0}, 0.0);
        NT_CHECK(h2.has_potential());
        NT_CLOSE(h2.U(Vec3{2, 1, 0}) - h2.U(Vec3{2, -1, 0}), -2.0 * std::atan(0.5), 1e-10);  // -0.9272952
        for (const Vec3& p : {Vec3{2, 1, 0}, Vec3{1.5, -0.7, 0.3}, Vec3{3, 1, -1}}) {
            NT_CLOSE(h2.U(p), -std::atan2(p.y, p.x), 1e-10);
        }
        NT_CLOSE(newton::fields::circulation(hole, Vec3{2, 0, 0}, Vec3{0, 0, 1}, 0.5), 0.0, 1e-9);  // a loop beside the axis
        NT_CLOSE(newton::fields::circulation(hole, origin, Vec3{0, 0, 1}, 0.5), 2.0 * pi, 1e-9);        // a loop around it
    }

    // --- time-dependent potentials --------------------------------------------------------------------
    {
        const CheckOptions big = box({-10, -10, -10}, {10, 10, 10});
        // pulsing spring: U = 2(1 + ½ sin t)|r|²; at t = π/2, ΔU((1,0,0), 0) = 3 J
        const auto pu = find_potential(pulsing_spring);
        NT_CHECK(pu.report.verdict == T && pu.has_potential());
        NT_CLOSE(pu.U(pi / 2.0, Vec3{1, 0, 0}) - pu.U(pi / 2.0, origin), 3.0, 1e-10);
        // moving support: gauge U(0, t) = 0 gives U = 2|r|² - 2 cos(3t) x
        const auto ms = find_potential(moving_support);
        NT_CHECK(ms.report.verdict == T && ms.has_potential());
        for (const double t : {0.0, 0.7, 2.3}) {
            for (const Vec3& p : {Vec3{0.5, 0.2, -0.3}, Vec3{-0.9, 0.9, 0.1}}) {
                NT_CLOSE(ms.U(t, p), 2.0 * norm2(p) - 2.0 * std::cos(3.0 * t) * p.x, 1e-10);
            }
        }
        // SinusoidalDrive (a Force-shaped t-dependent force): U = -2cos(3t + 0.5) x; at t = 1: 1.872913 J
        const auto dr = find_potential(forces::SinusoidalDrive{2.0, 3.0, 0.5});
        NT_CLOSE(dr.U(1.0, Vec3{1, 0, 0}) - dr.U(1.0, origin), -2.0 * std::cos(3.5), 1e-10);
        NT_CLOSE(dr.U(1.0, Vec3{1, 0, 0}) - dr.U(1.0, origin), 1.872913, 1e-6);

        // Rotating uniform field, released from rest at the origin: exact E(t) and the gauge identity.
        const auto rot = find_potential(rotating_field, big);
        NT_CHECK(rot.report.verdict == T && rot.has_potential());
        SimOptions so;
        so.t_end = 11.0;
        so.dt = 2e-3;
        const Trajectory traj = newton::simulate(Particle{1.0, {}, {}}, [](double t, const Vec3& r, const Vec3&) { return rotating_field(t, r); }, so);
        double worst = 0.0, scale = 0.0;
        for (const auto& s : traj.states) {
            const double E = analysis::total_energy(1.0, s, rot.U);
            const double exact = 2.0 - 2.0 * std::cos(s.t) - s.t * std::sin(s.t);  // F0²/(mω²) = 1
            worst = std::max(worst, std::fabs(E - exact));
            scale = std::max(scale, std::fabs(exact));
        }
        NT_CHECK(worst < 1e-8 * scale);
        // E(π) = 4 J: pick the sample nearest π
        std::size_t k = 0;
        for (std::size_t i = 0; i < traj.size(); ++i) {
            if (std::fabs(traj.states[i].t - pi) < std::fabs(traj.states[k].t - pi)) k = i;
        }
        const double tk = traj.states[k].t;
        NT_CLOSE(analysis::total_energy(1.0, traj.states[k], rot.U), 2.0 - 2.0 * std::cos(tk) - tk * std::sin(tk), 1e-8);
        NT_CLOSE(tk, pi, 2e-3);
        // E(t) - E(0) = ∫∂U/∂t dt: the energy_balance residual (with no non-conservative force) vanishes.
        const std::vector<double> res = analysis::energy_balance(traj, rot.U, forces::ConstantForce{});
        NT_CHECK(max_abs(res) < 1e-8);
        // ... while the plain drift is NOT small: the energy really changes (E(π) = 4 J from E(0) = 0).
        const std::vector<double> drift = analysis::energy_drift(traj, rot.U);
        NT_CHECK(max_abs(drift) > 1.0);
    }

    // --- force_from_potential ---------------------------------------------------------------------------------
    {
        // U = 2|r|²  =>  F = -4 r
        auto spring_U = [](const Vec3& r) { return 2.0 * norm2(r); };
        const auto Fs = force_from_potential(spring_U);
        for (const Vec3& p : {Vec3{1, 2, 3}, Vec3{-0.3, 0.1, 0.2}, Vec3{1e3, 0, -2e3}}) {
            NT_CHECK(rel_err(Fs(0.0, p, Vec3{}), -4.0 * p) < 1e-9);
        }
        // Point gravity at lab and orbital scale: U = -GMm/r
        auto grav_U = [](const Vec3& r) { return -6.0 / norm(r); };
        const auto Fg = force_from_potential(grav_U);
        const forces::PointGravity pg{3.0, 2.0};
        for (const Vec3& p : {Vec3{1, 0, 0}, Vec3{0.3, -2, 1}, Vec3{40, 3, -5}}) {
            NT_CHECK(rel_err(Fg(0.0, p, Vec3{}), pg(0.0, p, Vec3{})) < 1e-9);
        }
        // Morse at atomic scale needs length_scale ~ 1e-10; the default 1 m gives EXACTLY zero.
        const double A = 1e-19, r0 = 1.2e-10, S = 3e-11;
        auto atomic_U = [=](const Vec3& r) {
            const double e = std::exp(-(norm(r) - r0) / S);
            return A * (e * e - 2.0 * e);
        };
        const Vec3 p{1.3e-10, 0.4e-10, 0.2e-10};
        const double e = std::exp(-(norm(p) - r0) / S);
        const Vec3 exact = ((2.0 * A / S) * (e * e - e)) * normalized(p);
        const auto Fa = force_from_potential(atomic_U, 1e-10);
        NT_CHECK(rel_err(Fa(0.0, p, Vec3{}), exact) < 1e-7);
        NT_CHECK(norm(exact) > 1e-12);  // a real atomic-scale force (~1e-9 N)
        const auto Fbad = force_from_potential(atomic_U);  // length_scale = 1 m
        NT_CHECK(norm(Fbad(0.0, p, Vec3{})) == 0.0);
        // U(t, r) is preferred when callable: pulsing spring U = 2(1 + ½ sin t)|r|² => F = -4(1 + ½ sin t) r
        auto pulse_U = [](double t, const Vec3& r) { return 2.0 * (1.0 + 0.5 * std::sin(t)) * norm2(r); };
        const auto Fp = force_from_potential(pulse_U);
        for (const double t : {0.0, 1.3, 4.0}) {
            const Vec3 p2{0.5, -0.4, 0.9};
            NT_CHECK(rel_err(Fp(t, p2, Vec3{}), -4.0 * (1.0 + 0.5 * std::sin(t)) * p2) < 1e-9);
        }
        // The derived forces classify as the physics says.
        NT_CHECK(check_conservative(Fs).verdict == C);
        NT_CHECK(check_conservative(Fg).verdict == C);
        NT_CHECK(check_conservative(Fp).verdict == T);
        NT_CHECK(check_conservative(Fp).depends_on_time);
        NT_CHECK(check_conservative(force_from_potential([](const Vec3& r) { return std::sin(r.x) * r.y + r.z * r.z * r.x; })).verdict == C);
        // Kink pitfall: U = |x| overshoots by 7/6 at x = h (documented).
        auto vee = [](const Vec3& r) { return std::fabs(r.x); };
        const double h = 1e-4;
        NT_CLOSE(-force_from_potential(vee)(0.0, Vec3{h, 0, 0}, Vec3{}).x, 7.0 / 6.0, 1e-9);
    }

    // --- round trips -----------------------------------------------------------------------------------------------
    {
        // F -> U -> F
        const auto F = forces::sum(forces::UniformGravity{2.0}, forces::Spring{4.0});
        const auto res = find_potential(F);
        const auto back = force_from_potential(res.U);
        for (const Vec3& p : {Vec3{0.5, -0.3, 0.8}, Vec3{-0.9, 0.2, 0.1}, Vec3{0.1, 0.1, -0.7}}) {
            NT_CHECK(rel_err(back(0.0, p, Vec3{}), F(0.0, p, Vec3{})) < 1e-7);
        }
        // U -> F -> U: compare differences (the constant is free)
        auto U0 = [](const Vec3& r) { return r.x * r.x * r.y + std::sin(r.z) * r.x + 2.0 * norm2(r); };
        const auto F0 = force_from_potential(U0);
        const auto r2 = find_potential(F0, CheckOptions{}, origin, U0(origin));
        NT_CHECK(r2.has_potential());
        for (const Vec3& p : {Vec3{0.5, -0.3, 0.8}, Vec3{-0.9, 0.2, 0.1}, Vec3{0.7, 0.7, -0.7}}) {
            NT_CLOSE(r2.U(p) - r2.U(Vec3{0.1, 0.2, 0.3}), U0(p) - U0(Vec3{0.1, 0.2, 0.3}), 1e-8);
        }
    }

    // --- split helpers --------------------------------------------------------------------------------------------------
    {
        const auto gd = forces::sum(forces::UniformGravity{2.0}, forces::LinearDrag{0.5});
        const auto rest = velocity_independent_part(gd);
        const auto moving = velocity_dependent_part(gd);
        NT_CHECK(check_conservative(rest).verdict == C);
        const ConservativeReport m = check_conservative(moving);
        NT_CHECK(m.verdict == V && m.velocity_power == VelocityPower::Dissipative);
        const Vec3 v{3.0, -1.0, 2.0}, r{0.1, 0.2, 0.3};
        NT_CHECK(rel_err(rest(0.0, r, v) + moving(0.0, r, v), gd(0.0, r, v)) < 1e-14);  // the split is exact bookkeeping
        NT_CHECK(dot(moving(0.0, r, v), v) < 0.0);
        // Lorentz: the magnetic part does no work; the rest part is a uniform E field with U = -qE·r
        const forces::LorentzForce lf{2.0, {0, 0, 5}, {0, 0, 1}};
        const ConservativeReport lm = check_conservative(velocity_dependent_part(lf));
        NT_CHECK(lm.verdict == V && lm.velocity_power == VelocityPower::NoWork);
        NT_CLOSE(dot(velocity_dependent_part(lf)(0.0, r, v), v), 0.0, 1e-13);
        const auto lr = find_potential(velocity_independent_part(lf));
        NT_CHECK(lr.has_potential());
        NT_CLOSE(lr.U(Vec3{0, 0, 1}) - lr.U(origin), -10.0, 1e-10);
        // shapes other than Force work too: a (t, r) field is split into an equal rest part and a zero moving part
        const auto z = velocity_dependent_part(pulsing_spring);
        NT_CHECK(z(1.0, r, v) == Vec3{});
    }


    // --- velocity power labels (audit 16) -----------------------------------------------------
    {
        auto power = [](const auto& f) { return check_conservative(f).velocity_power; };
        NT_CHECK(power(forces::LinearDrag{0.5}) == VelocityPower::Dissipative);
        NT_CHECK(power(forces::MagneticForce{}) == VelocityPower::NoWork);
        // negative damping +b v: F_v·v > 0 everywhere, it FEEDS energy
        NT_CHECK(power([](double, const Vec3&, const Vec3& v) { return 0.3 * v; }) == VelocityPower::Feeding);
        NT_CHECK(power(forces::sum(forces::Spring{4.0}, [](double, const Vec3&, const Vec3& v) { return 0.3 * v; })) ==
                 VelocityPower::Feeding);
        // Friction written -mu v/|v| is NaN at v = 0, so F_v = F(v) - F(0) has no usable cosine: DELIBERATELY
        // Unknown (the label says "could not tell"); the zero-guarded form -mu normalized(v) is Dissipative.
        NT_CHECK(power([](double, const Vec3&, const Vec3& v) { return -2.0 * v / norm(v); }) == VelocityPower::Unknown);
        NT_CHECK(power([](double, const Vec3&, const Vec3& v) { return -2.0 * normalized(v); }) == VelocityPower::Dissipative);
        // both signs genuinely occur: -x v_x-direction drag flips sign with x
        NT_CHECK(power([](double, const Vec3& r, const Vec3& v) { return r.x * v; }) == VelocityPower::Mixed);
        // A magnetic part tiny against a huge E field: the subtraction's rounding noise hides the angle: Unknown, never a wrong label.
        NT_CHECK(power(forces::LorentzForce{1.0, {1e9, 0, 0}, {0, 0, 1}}) == VelocityPower::Unknown);
        NT_CHECK(power(forces::LorentzForce{1.0, {1e7, 3e6, 0}, {0, 0, 1}}) == VelocityPower::NoWork);
        // wind drag -b(v - w): the label describes the velocity-dependent PART (-b v): Dissipative, although the
        // true drag power -b(v - w)·v can be positive (research-physics.md §8(b)).
        const Vec3 wind{3, 0, 0};
        auto wind_drag = [wind](double, const Vec3&, const Vec3& v) { return -0.5 * (v - wind); };
        NT_CHECK(power(wind_drag) == VelocityPower::Dissipative);
        NT_CHECK(dot(wind_drag(0.0, origin, Vec3{1, 0, 0}), Vec3{1, 0, 0}) > 0.0);
        const std::string d = describe(check_conservative(wind_drag));
        NT_CHECK(d.find("velocity-dependent PART") != std::string::npos && d.find("wind") != std::string::npos);
        NT_CHECK(describe(check_conservative(forces::MagneticForce{})).find("never changes the kinetic energy") != std::string::npos);
        NT_CHECK(describe(check_conservative([](double, const Vec3&, const Vec3& v) { return 0.3 * v; })).find("supplies energy") != std::string::npos);
        NT_CHECK(std::string(to_string(VelocityPower::Feeding)) == "feeding");
        NT_CHECK(std::string(to_string(VelocityPower::Unknown)) == "unknown");
    }

    // --- hole coverage (audit 12 B1): oblique lines, wire pairs, non-cubic boxes ---------------------
    {
        const Vec3 s2{1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0), 0.0};
        const Vec3 s3{1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
        NT_CHECK(check_conservative(wire_field({0.6, -0.6, 0.5}, s2)).verdict == P);  // oblique line along (1,1,0)
        NT_CHECK(check_conservative(wire_field({0.7, 0.7, -0.7}, s3)).verdict == P);  // along (1,1,1)
        auto pair = [](const Vec3& r) {  // antiparallel wires along z at x = ±0.5: the big loop sees 0 total
            return wire_field({0.5, 0, 0}, {0, 0, 1})(r) - wire_field({-0.5, 0, 0}, {0, 0, 1})(r);
        };
        NT_CHECK(check_conservative(pair).verdict == P);
        // box [-3,3]x[-1,1]² with a wire along z at (2,0): outside the old central circles
        NT_CHECK(check_conservative(wire_field({2, 0, 0}, {0, 0, 1}), box({-3, -1, -1}, {3, 1, 1})).verdict == P);
        // lab box [0,10]x[0,2]x[-1,1], wire at (7,1)
        NT_CHECK(check_conservative(wire_field({7, 1, 0}, {0, 0, 1}), box({0, 0, -1}, {10, 2, 1})).verdict == P);
        // and the same box without the wire in it is conservative
        NT_CHECK(check_conservative(wire_field({12, 1, 0}, {0, 0, 1}), box({0, 0, -1}, {10, 2, 1})).verdict == C);
        // the report says where: the worst loop is a cell around the wire (centre within a cell of it)
        const ConservativeReport w = check_conservative(wire_field({2, 0, 0}, {0, 0, 1}), box({-3, -1, -1}, {3, 1, 1}));
        NT_CHECK(w.max_loop_ratio > 0.5);
        NT_CHECK(std::fabs(w.worst_loop_center.x - 2.0) < 3.1 && std::fabs(w.worst_loop_center.y) < 1.1);

        // Seeded random-line sweep: lines through a random point of the inner box [-0.8,0.8]³ (at least
        // 0.2 from every face), random direction. Count the lines classified Conservative (misses).
        detail::Rng rng(987654321);
        int misses = 0;
        const int n_lines = 300;
        for (int k = 0; k < n_lines; ++k) {
            const Vec3 p = rng.in_box(Box{{-0.8, -0.8, -0.8}, {0.8, 0.8, 0.8}});
            Vec3 d;
            do {
                d = rng.in_cube(1.0);
            } while (norm2(d) > 1.0 || norm2(d) < 1e-3);
            d = normalized(d);
            const Verdict v = check_conservative(wire_field(p, d)).verdict;
            if (v == C) ++misses;
        }
        const double miss_rate = static_cast<double>(misses) / n_lines;
        if (print_table) std::printf("random-line sweep: %d of %d lines missed (%.2f%%)\n", misses, n_lines, 100.0 * miss_rate);
        NT_CHECK(miss_rate < 0.01);
    }

    // --- coherent rounding noise is not a curl (audit 12 S4) ---------------------------------------------
    {
        for (const double Cn : {1e8, 1e10}) {
            auto noisy = [Cn](const Vec3& r) {  // n(u) = ((u + C) - C) - u is 0 in exact arithmetic
                auto n = [Cn](double u) { volatile double a = u + Cn; const double b = a - Cn; return b - u; };
                return Vec3{-4.0 * r.x + n(r.y), -4.0 * r.y + n(r.z), -4.0 * r.z + n(r.x)};
            };
            const ConservativeReport r = check_conservative(noisy);
            NT_CHECK(r.verdict != N);
            if (print_table) std::printf("noise C=%g: %s, unconfirmed %d\n", Cn, to_string(r.verdict), r.unconfirmed_curl_candidates);
        }
        // the true-curl rows still confirm (they are in the battery); a weak true curl too:
        NT_CHECK(check_conservative(ab_field(1.0, 1.0 + 1e-4)).verdict == N);
    }

    // --- the step follows the box, not the origin (audit 13 S2, audit 12 S5) -------------------------------
    {
        for (const double X : {1e2, 1e3, 1e4}) {
            const Vec3 c{X + 0.37, -0.21, 0.29};
            CheckOptions o = box({X - 1, -1, -1}, {X + 1, 1, 1});
            auto asym = [c](const Vec3& r) {  // (1 + cos²θ) r̂/r²: central but not spherically symmetric
                const Vec3 d = r - c;
                const double s2 = norm2(d), s = std::sqrt(s2);
                return ((1.0 + d.z * d.z / s2) / (s2 * s)) * d;
            };
            auto coul = [c](const Vec3& r) {
                const Vec3 d = r - c;
                const double s = norm(d);
                return (2.0 / (s * s * s)) * d;
            };
            const ConservativeReport a = check_conservative(asym, o);
            NT_CHECK(a.verdict == N);
            NT_CHECK(check_conservative(coul, o).verdict == C);
        }
        // a (2 + cosθ) r̂ force (central, not spherically symmetric) in a box centred at (1e4,0,0)
        auto asym2 = [](const Vec3& r) {
            const Vec3 c{1e4 + 0.37, -0.21, 0.29};
            const Vec3 d = r - c;
            const double s = norm(d);
            return ((2.0 + d.z / s) / s) * d;
        };
        NT_CHECK(check_conservative(asym2, box({1e4 - 1, -1, -1}, {1e4 + 1, 1, 1})).verdict == N);
        // a molecule-sized swirl (Gaussian width 1e-10 m) in a 1e-9 m box: at the origin, at (1,0,0) and at (1e4,0,0)
        for (const double X : {0.0, 1.0, 1e4}) {
            const Vec3 m{X, 0, 0};
            auto swirl = [m](const Vec3& r) {
                const Vec3 q = r - m;
                return 10.0 * std::exp(-norm2(q) / 1e-20) * Vec3{-q.y, q.x, 0.0};
            };
            const ConservativeReport sw = check_conservative(swirl, box(m - Vec3{5e-10, 5e-10, 5e-10}, m + Vec3{5e-10, 5e-10, 5e-10}));
            if (print_table) std::printf("molecular swirl at x=%g: %s\n", X, to_string(sw.verdict));
            NT_CHECK(sw.verdict == N);
        }
    }

    // --- degenerate regions and notes (audit 12 N4) ------------------------------------------------------------
    {
        const ConservativeReport flat = check_conservative(forces::Spring{4.0}, box({-1, -1, 0}, {1, 1, 0}));
        NT_CHECK(flat.verdict == Verdict::Inconclusive && flat.note[0] != '\0');
        NT_CHECK(describe(flat).find("degenerate") != std::string::npos);
        NT_CHECK(check_conservative(forces::Spring{4.0}, box({1, 1, 1}, {-1, -1, -1})).verdict == Verdict::Inconclusive);
        NT_CHECK(check_conservative(forces::Spring{4.0}, box({kNaN, -1, -1}, {1, 1, 1})).verdict == Verdict::Inconclusive);
        const ConservativeReport nanr = check_conservative([](const Vec3&) { return Vec3{kNaN, 0.0, 0.0}; });
        NT_CHECK(nanr.verdict == Verdict::Inconclusive && nanr.note[0] != '\0');
    }

    // --- find_potential: path-independence validation, boxes, reference rules (audit 15) --------------------------
    {
        // B1: U -> F -> U with U = -6/r singular at the default reference (origin): the numerical force is
        // finite near the singularity but wrong there; the reference must be moved and ΔU correct.
        auto U0 = [](const Vec3& r) { return -6.0 / norm(r); };
        const auto F0 = force_from_potential(U0);
        const auto res = find_potential(F0);
        NT_CHECK(res.report.verdict == C);
        NT_CHECK(res.reference_ok && res.reference_moved && res.r_ref_used != origin);
        const Vec3 p{0.7, 0.2, -0.4}, q{-0.3, 0.9, 0.5};
        NT_CLOSE(res.U(p) - res.U(q), U0(p) - U0(q), 1e-8);
        const Vec3 p2{1.0, 0.0, 0.0}, q2{0.0, 2.0, 0.0};
        NT_CLOSE(res.U(p2) - res.U(q2), -3.0, 1e-8);  // research value ΔU((1,0,0),(0,2,0)) = -3 J... for GMm = 6

        // B2: detours never leave the box. A point mass at (0,0.95,0) ON the straight path, plus a singular
        // line {(x, 1.1, 0)} 0.1 m OUTSIDE the box: an unclamped detour waypoint (0,1.4,0) wraps the line (2π).
        const Vec3 c{0, 0.95, 0};
        auto Fm = [c](const Vec3& r) {
            const Vec3 d = r - c;
            const double n = norm(d);
            const double u = r.y - 1.1, w = r.z;
            const double s = u * u + w * w;
            return (-1.0 / (n * n * n)) * d + Vec3{0.0, -w / s, u / s};
        };
        auto Uex = [c](const Vec3& r) {
            const double u = r.y - 1.1, w = r.z;
            return -1.0 / norm(r - c) - (std::atan2(-w, -u) + pi);
        };
        const Vec3 a{0, 0.95, -0.9};
        const auto m = find_potential(Fm, CheckOptions{}, a, Uex(a));
        NT_CHECK(m.report.verdict == C && m.has_potential());
        NT_CHECK(m.U.region.has_value());
        for (const Vec3& b : {Vec3{0, 0.95, 0.9}, Vec3{0.3, 0.95, 0.6}, Vec3{0, 0.9, 0.9}, Vec3{0.5, -0.5, 0.5}}) {
            NT_CLOSE(m.U(b), Uex(b), 1e-9);
        }

        // B3: a caller's r_ref OUTSIDE the box is rejected (and never moved)
        const auto out = find_potential(hole, hole_b, Vec3{-2, 0, 0});
        NT_CHECK(out.report.verdict == C);
        NT_CHECK(!out.reference_ok && !out.reference_moved && !out.has_potential());
        NT_CHECK(out.r_ref_used == Vec3{-2, 0, 0});
        NT_CHECK(out.note[0] != '\0' && describe(out).find("outside") != std::string::npos);

        // B4: a time-dependent force singular at the default reference only at t_min = 0 (a charge at (0.3 t, 0, 0))
        auto moving_charge = [](double t, const Vec3& r) {
            const Vec3 d = r - Vec3{0.3 * t, 0, 0};
            const double n = norm(d);
            return (1.0 / (n * n * n)) * d;
        };
        const auto mc = find_potential(moving_charge);
        NT_CHECK(mc.report.verdict == T);
        NT_CHECK(mc.reference_ok);
        const Vec3 pp{-0.5, 0.4, 0.2};
        NT_CHECK(std::isfinite(mc.U(0.0, pp)));
        NT_CHECK(std::isfinite(mc.U(pp)));
        NT_CHECK(std::isfinite(mc.U(1.0, pp)));
    }

    // --- force_from_potential: the centre (audit 13 S1, audit 15 S4) and the time-dependent round trip -------------
    {
        const Vec3 c{1000.0, 0, 0};
        auto Uc = [c](const Vec3& r) { return -1.0 / norm(r - c); };
        const Vec3 p = c + Vec3{1.0, 0.0, 0.0};  // 1 m from the centre, 1000 m from the origin
        const Vec3 exact = -(1.0 / 1.0) * normalized(p - c);   // -dU/dr r̂ = -(1/r²) r̂ ... magnitude 1/r² = 1, pointing inward
        const auto centered = force_from_potential(Uc, 1.0, c);
        // the limit is rounding of the coordinates themselves: ulp(1000 m)/h = 1e-13/1e-4 ~ 1e-9 (measured 2.5e-10)
        NT_CHECK(rel_err(centered(0.0, p, Vec3{}), exact) < 1e-8);
        NT_CHECK(rel_err(force_from_potential(Uc)(0.0, p, Vec3{}), exact) > 1e-6);  // the origin-based step is worse (documented)
        // an off-origin Morse well
        auto Um = [c](const Vec3& r) {
            const double e = std::exp(-(norm(r - c) - 1.0) / 0.5);
            return e * e - 2.0 * e;
        };
        const Vec3 q = c + Vec3{1.3, 0.2, -0.4};
        const double e = std::exp(-(norm(q - c) - 1.0) / 0.5);
        const Vec3 fm = (4.0 * (e * e - e)) * normalized(q - c);
        NT_CHECK(rel_err(force_from_potential(Um, 0.5, c)(0.0, q, Vec3{}), fm) < 1e-9);
        // F -> U(t, r) -> F for a time-dependent force: the rotating field
        const auto rot = find_potential(rotating_field, box({-10, -10, -10}, {10, 10, 10}));
        const auto back = force_from_potential(rot.U);
        for (const double t : {0.0, 1.3, 4.0}) {
            const Vec3 pt{0.5, -0.7, 0.9};
            NT_CHECK(rel_err(back(t, pt, Vec3{}), rotating_field(t, pt)) < 1e-9);
        }
        NT_CHECK(check_conservative(back).verdict == T);
    }

    return newton_test::report("potentials");
}
