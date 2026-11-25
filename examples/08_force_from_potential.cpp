// examples/08_force_from_potential.cpp — "From U to F to motion" (Taylor §4.3-4.5, §5.1, §5.4).
//
// Physical situation: example 07 went from a force to its potential energy. This
// one goes the other way, which is how most of physics is actually written down:
// someone hands you U(r) and you need the motion. F = -grad U (Taylor §4.3).
//   (1) A Morse potential in one line (a bond-like well, scaled to everyday size).
//   (2) force_from_potential(U) against the hand-derived -dU/dx.
//   (3) check_conservative on that force: grad U can never have curl (§4.4).
//   (4) A small oscillation in the well: energy_drift, and the period against
//       2 pi / omega with omega = sqrt(U''(0)/m). This is the bridge to chapter 5
//       (§5.1): every smooth minimum looks like a spring. A bigger swing shows the
//       period LENGTHENING, because the well is not a parabola.
//   (5) A time-dependent U(t, r): a spring whose support moves (§4.5). energy_drift
//       is large, energy_balance (which includes dU/dt) is ~0.
//   (6) Add LinearDrag: the energy decays, but energy_balance(traj, U, drag) ~ 0 and
//       work_done(drag) is negative and equals the energy CHANGE E(t) - E(0);
//       the energy lost is -W_drag (§4.1-4.2, §5.4).
//   (7) A CSV of (t, x, E, W_drag) for the damped run: build/morse_damped.csv.
//
// The potential. U(x) = A (1 - e^{-x/S})^2, zero at its minimum x = 0, flattening out to
// A = 1 J far to the right (the "dissociation energy"). S = 0.3 m is the width of the
// well. We use A = 1 J, S = 0.3 m, m = 1 kg so the numbers are easy to read. For a real
// diatomic molecule A ~ 1e-19 J and S ~ 1e-10 m; then pass length_scale ~ 1e-10 to
// force_from_potential, because its finite-difference step must be small compared with the
// distance over which U changes (default 1.0 m would be a catastrophic step for an atom).
//
// Analytic expectation:
//   F(x) = -dU/dx = -(2A/S) e^{-x/S} (1 - e^{-x/S}).
//   U''(0) = 2A/S^2, so omega = sqrt(2A/(m S^2)) = sqrt(2)/0.3 = 4.714 rad/s,
//   period 2 pi/omega = 1.3329 s for small amplitude.
//   The Morse well even has an exact period at any energy E < A:
//       T(E) = T_0 / sqrt(1 - E/A),   E = U(x_start) for a release from rest.
//   (5) Moving support: x(t) = 0.4 (cos 2t - cos 3t), a closed form for the driven spring.
//   (6) Energy lost = work done by drag, to integrator error.
//
// What to look at: the measured period matches 2 pi/omega for the tiny swing and the exact
// T(E) for the big one, and that T(E) is clearly longer than 2 pi/omega (anharmonic).
// The drift in (5) is order 1 and the balance is ~1e-10; in (6) the drift is order 1 and the
// balance is again tiny. All of the "drift" is physics; only the balance measures numerics.
//
// Error column: "rel err" is |numeric - expected| / |expected| (absolute where expected = 0).
//
// The new headers below are also reachable through <newton/newton.hpp>; they're
// listed explicitly so you can see where each name comes from (#pragma once
// makes the double include harmless).

#include <algorithm>
#include <cmath>
#include <cstdio>
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

void compare(const char* label, double numeric, double expected) {
    std::printf("  %-34s numeric % .9g   expected % .9g   rel err %.1e\n", label, numeric, expected,
                rel_err(numeric, expected));
}

void heading(const char* text) { std::printf("\n=== %s ===\n", text); }

// NaN-aware: std::max(m, NaN) silently drops the NaN, so a ledger full of NaN would
// print as "0". Any NaN here makes the result NaN, which prints as such.
double max_abs(const std::vector<double>& v) {
    double m = 0.0;
    for (double x : v) {
        if (std::isnan(x)) return std::nan("");
        m = std::max(m, std::fabs(x));
    }
    return m;
}

// --- (1) the potential, in one line -----------------------------------------------------
constexpr double kA = 1.0;    // well depth [J]
constexpr double kS = 0.3;    // well width [m]
constexpr double kMass = 1.0;  // [kg]
// U(r) = A (1 - e^{-x/S})^2 : a function of r.x only, so the force points along x.
const auto morse_U = [](const Vec3& r) {
    const double e = 1.0 - std::exp(-r.x / kS);
    return kA * e * e;
};
// The same thing by hand: F_x = -dU/dx.
double morse_F(double x) {
    const double e = std::exp(-x / kS);
    return -(2.0 * kA / kS) * e * (1.0 - e);
}
// Exact Morse period at energy E (0 < E < A): T_0 / sqrt(1 - E/A).
double morse_period(double E) {
    const double omega0 = std::sqrt(2.0 * kA / (kMass * kS * kS));
    return 2.0 * constants::pi / omega0 / std::sqrt(1.0 - E / kA);
}

}  // namespace

int main() {
    const double pi = constants::pi;
    const double omega = std::sqrt(2.0 * kA / (kMass * kS * kS));  // sqrt(U''(0)/m)
    const double T_harmonic = 2.0 * pi / omega;

    // ------------------------------------------------------------------------------------
    // (2) force_from_potential vs the analytic derivative (Taylor §4.3).
    //
    // GradientForce is a library Force: it has the (t, r, v) signature, so it drops into
    // simulate() and forces::sum() like Spring does. Inside, F = -grad U by a 4th-order
    // finite difference with step ~1e-4 * max(|r|, length_scale).
    // ------------------------------------------------------------------------------------
    heading("(2) force_from_potential(U) against -dU/dx");
    const auto morse = potentials::force_from_potential(morse_U, kS);  // length_scale = well width
    for (double x : {-0.2, -0.05, 0.1, 0.3, 0.9}) {
        char label[48];
        std::snprintf(label, sizeof label, "F_x at x = % .2f m [N]", x);
        compare(label, morse(0.0, Vec3{x, 0.0, 0.0}, Vec3{}).x, morse_F(x));
    }
    const Vec3 off_axis = morse(0.0, Vec3{0.3, 0.2, -0.1}, Vec3{});
    std::printf("  F_y, F_z at (0.3, 0.2, -0.1): % .1e, % .1e  (U does not depend on y or z, so both are 0)\n",
                off_axis.y, off_axis.z);
    // The stiffness at the bottom: k_eff = -dF/dx = U''(0) = 2A/S^2, from the force itself.
    const double h = 1e-4;  // m; the central difference itself has error ~ (h/S)^2
    const double k_eff = -(morse(0.0, Vec3{h, 0.0, 0.0}, Vec3{}).x - morse(0.0, Vec3{-h, 0.0, 0.0}, Vec3{}).x) / (2.0 * h);
    compare("U''(0) from the force [N/m]", k_eff, 2.0 * kA / (kS * kS));

    // ------------------------------------------------------------------------------------
    // (3) A potential always gives a conservative force (Taylor §4.4: curl grad U = 0).
    // ------------------------------------------------------------------------------------
    heading("(3) classify the derived force");
    const potentials::ConservativeReport rep = potentials::check_conservative(morse);
    std::printf("  verdict: %s, max relative curl = %.1e\n", potentials::to_string(rep.verdict), rep.max_relative_curl);
    std::printf("    %s\n", potentials::describe(rep).c_str());
    std::printf("  -> curl(grad U) = 0 because mixed partial derivatives commute. (In general the finite\n"
                "     differences leave a tiny numerical curl; here it is exactly 0 since U depends on x only.)\n");

    // ------------------------------------------------------------------------------------
    // (4) Small oscillations: the chapter-4 to chapter-5 bridge (Taylor §5.1).
    //
    // Near the minimum, U(x) = 1/2 U''(0) x^2 + ..., so the particle behaves like a spring of
    // stiffness k = U''(0). The release point is on +x so analysis::measure_period (which
    // averages the gaps between UPWARD zero crossings of r.x) sees the oscillation about x = 0.
    // ------------------------------------------------------------------------------------
    heading("(4) oscillation in the well: energy drift and period");
    std::printf("  omega = sqrt(U''(0)/m) = sqrt(2A/(m S^2)) = %.6f rad/s, 2 pi/omega = %.6f s\n", omega, T_harmonic);

    struct Swing {
        const char* name;
        double x0;
    };
    for (const Swing& swing : {Swing{"small swing, x0 = 0.01 m", 0.01}, Swing{"large swing, x0 = S = 0.3 m", kS}}) {
        const double E = morse_U(Vec3{swing.x0, 0.0, 0.0});  // released from rest
        SimOptions opt;
        opt.t_end = 8.0 * morse_period(E);
        opt.dt = morse_period(E) / 1000.0;
        opt.method = Method::RK4;
        const Trajectory traj = simulate(Particle{kMass, Vec3{swing.x0, 0.0, 0.0}, Vec3{}}, morse, opt);

        const std::vector<double> drift = analysis::energy_drift(traj, morse_U);
        const double T = analysis::measure_period(traj);
        std::printf("  -- %s: E = %.6f J --\n", swing.name, E);
        std::printf("  max |energy_drift| = %.1e   (same U as the force came from; RK4, dt = T/1000)\n", max_abs(drift));
        compare("period vs 2 pi/omega [s]", T, T_harmonic);
        compare("period vs exact T(E) [s]", T, morse_period(E));
    }
    std::printf("  -> for a 0.01 m swing the period is only 0.05%% above 2 pi/omega: the well is nearly a\n"
                "     parabola, so the spring picture works. At x0 = S\n"
                "     the particle climbs the soft right wall, the effective stiffness drops, and the\n"
                "     period is about 29%% longer than 2 pi/omega. That is anharmonicity (Taylor §5.1);\n"
                "     the exact Morse T(E) = T_0/sqrt(1 - E/A) predicts it.\n");

    // ------------------------------------------------------------------------------------
    // (5) A time-dependent potential: a spring whose support moves (Taylor §4.5).
    //
    // U(t, r) = 1/2 k |r - s(t)|^2 with s(t) = A cos(w t) x-hat. force_from_potential accepts
    // U(t, r) (it prefers that form when U can take a time argument), and total_energy,
    // energy_drift and energy_balance accept it too. T + U is no longer constant: whoever
    // moves the support does work at the rate dU/dt, so
    //      E(t) - E(0) = integral of dU/dt dt.
    // energy_balance subtracts that integral; energy_drift does not.
    // k = 4 N/m, A = 0.5 m, w = 3 rad/s, m = 1 kg, from rest at the origin.
    // ------------------------------------------------------------------------------------
    heading("(5) a moving spring support: U(t, r)");
    constexpr double k = 4.0, amp = 0.5, w = 3.0;
    const auto moving_U = [](double t, const Vec3& r) {
        return 0.5 * k * norm2(r - Vec3{amp * std::cos(w * t), 0.0, 0.0});
    };
    const auto moving_F = potentials::force_from_potential(moving_U);
    std::printf("  verdict: %s\n", potentials::to_string(potentials::check_conservative(moving_F).verdict));

    SimOptions sopt;
    sopt.t_end = 10.0;
    sopt.dt = 1e-3;
    sopt.method = Method::RK4;
    const Trajectory moving = simulate(Particle{1.0, Vec3{}, Vec3{}}, moving_F, sopt);

    // The closed form: x(t) = 0.4 (cos 2t - cos 3t), v(t) = 1.2 sin 3t - 0.8 sin 2t, so
    // E = 1/2 v^2 + 1/2 k (x - A cos 3t)^2 can be written down independently of the run.
    for (double t : {1.0, 2.5, 5.0, 10.0}) {
        const std::size_t i = static_cast<std::size_t>(std::llround(t / sopt.dt));
        const double x = 0.4 * (std::cos(2.0 * t) - std::cos(3.0 * t));
        const double v = 1.2 * std::sin(3.0 * t) - 0.8 * std::sin(2.0 * t);
        const double E_exact = 0.5 * v * v + 0.5 * k * (x - amp * std::cos(w * t)) * (x - amp * std::cos(w * t));
        char label[48];
        std::snprintf(label, sizeof label, "E(t = %.1f s) [J]", t);
        compare(label, analysis::total_energy(1.0, moving.states[i], moving_U), E_exact);
    }
    const std::vector<double> m_drift = analysis::energy_drift(moving, moving_U);
    const std::vector<double> m_bal = analysis::energy_balance(moving, moving_U, [](double, const Vec3&, const Vec3&) { return Vec3{}; });
    std::printf("  max |energy_drift|   = %.1e   (E really changes: dE/dt = dU/dt)\n", max_abs(m_drift));
    std::printf("  max |energy_balance| = %.1e   (E(t) - E(0) - integral of dU/dt: the numerical error only)\n", max_abs(m_bal));

    // ------------------------------------------------------------------------------------
    // (6) Add drag (Taylor §4.1-4.2, §5.4).
    //
    // Total force = -grad U (the Morse well) + F_drag, with F_drag = -b v, b = 0.3 kg/s. The
    // potential is still the same U, but now E = T + U decays; the ledger is
    //      T + U - (T + U)_0 = W_drag = integral of F_drag . v dt,
    // and energy_balance(traj, U, drag) is the residual of that equation.
    // ------------------------------------------------------------------------------------
    heading("(6) Morse well plus linear drag");
    const forces::LinearDrag drag{0.3};
    const auto damped_force = forces::sum(morse, drag);

    SimOptions dopt;
    dopt.t_end = 10.0;
    dopt.dt = 1e-3;
    dopt.method = Method::RK4;
    const Trajectory damped = simulate(Particle{kMass, Vec3{kS, 0.0, 0.0}, Vec3{}}, damped_force, dopt);

    const std::vector<double> d_drift = analysis::energy_drift(damped, morse_U);
    const std::vector<double> d_bal = analysis::energy_balance(damped, morse_U, drag);
    const std::vector<double> W_drag = analysis::work_done(damped, drag);
    const double E_start = analysis::total_energy(kMass, damped.states.front(), morse_U);
    const double E_end = analysis::total_energy(kMass, damped.states.back(), morse_U);

    std::printf("  E(0) = %.6f J, E(10 s) = %.6f J\n", E_start, E_end);
    std::printf("  energy_drift at t = 10 s: % .6f   (relative loss; the energy really did go away)\n", d_drift.back());
    compare("E(10) - E(0) vs W_drag [J]", E_end - E_start, W_drag.back());
    std::printf("  max |energy_balance(traj, U, drag)| = %.1e   (balance closes: all the loss is drag's work)\n", max_abs(d_bal));
    // For a rough cross-check with chapter 5: light damping makes E fall by e^{-b t/m} on average.
    std::printf("  rough check: E(10) / E(0) = %.4f, e^{-b t/m} = %.4f (energy decays by e^{-bt/m} on average,\n"
                "     modulated by the anharmonic well; not an exact comparison)\n",
                E_end / E_start, std::exp(-0.3 * 10.0 / kMass));

    // ------------------------------------------------------------------------------------
    // (7) CSV: t, x, E, W_drag. Column E + (-W_drag) is flat: E(t) - W_drag(t) = E(0).
    // ------------------------------------------------------------------------------------
    std::vector<std::vector<double>> rows;
    for (std::size_t i = 0; i < damped.size(); i += 10) {
        const State& s = damped.states[i];
        rows.push_back({s.t, s.r.x, analysis::total_energy(kMass, s, morse_U), W_drag[i]});
    }
    const bool wrote = io::write_csv("build/morse_damped.csv", {"t", "x", "E", "W_drag"}, rows);
    std::printf("\n  wrote build/morse_damped.csv (%zu rows): %s\n", rows.size(), wrote ? "ok" : "FAILED");
    std::printf("  plot E and W_drag + E(0) together: they lie on top of each other.\n");

    return 0;
}
