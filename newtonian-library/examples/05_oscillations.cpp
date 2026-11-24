// examples/05_oscillations.cpp — SHM, damping regimes, and driven resonance
// (Taylor ch. 5).
//
// Physical situation, three demos sharing one spring:
//   (1) undamped SHM (forces::Spring alone) — period & amplitude vs
//       analytic::SHM (Taylor §5.2).
//   (2) all three damping regimes (forces::Spring + forces::LinearDrag) —
//       under-, critically-, and over-damped ring-down, vs
//       analytic::DampedOscillator (Taylor §5.4).
//   (3) a driven, lightly-damped oscillator (+ forces::SinusoidalDrive)
//       swept in driving frequency, vs analytic::DrivenOscillator's
//       steady-state amplitude (Taylor §5.5-5.6) — the resonance curve.
//
// Analytic expectation: SHM period/amplitude should match to ~1e-6 or
// better at the dt used here (RK4, dt = period/1000). The damped x(t)
// curves should match the closed forms in each regime. The resonance sweep
// should peak near analytic::DrivenOscillator::resonant_omega() and match
// amplitude() closely once the transient (which decays like e^{-beta t})
// has had time to die out.
//
// What to look at: build/damping_regimes.csv overlays all three regimes on
// one time axis — underdamped should oscillate through zero several times
// while decaying, critically damped should return to zero fastest WITHOUT
// oscillating, and overdamped should crawl back slowly, also without
// oscillating. build/resonance.csv is the amplitude-vs-omega curve; it
// should peak near (slightly below) omega0, at resonant_omega(), and track
// the analytic curve closely across the sweep.

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
    const double omega0 = 2.0;  // natural frequency [rad/s]
    const double k = mass * omega0 * omega0;

    // --- (1) undamped SHM -------------------------------------------------------
    std::printf("=== SHM: period and amplitude ===\n");
    const double x0 = 1.0, v0 = 0.0;
    const Particle shm_particle{mass, Vec3{x0, 0.0, 0.0}, Vec3{v0, 0.0, 0.0}};
    const forces::Spring shm_spring{k};

    const analytic::SHM shm{omega0, x0, v0};

    SimOptions shm_opt;
    shm_opt.t0 = 0.0;
    shm_opt.t_end = 5.0 * shm.period();
    shm_opt.dt = shm.period() / 1000.0;
    shm_opt.method = Method::RK4;

    const Trajectory shm_traj = simulate(shm_particle, shm_spring, shm_opt);
    const double period_numeric = analysis::measure_period(shm_traj);
    const double amplitude_numeric = analysis::measure_amplitude(shm_traj, 0.0);

    row("numeric period [s]", period_numeric);
    row("analytic period [s]", shm.period());
    row("numeric amplitude [m]", amplitude_numeric);
    row("analytic amplitude [m]", shm.amplitude());
    row("energy (const) [J]", shm.energy(mass));
    std::printf("\n");

    // --- (2) three damping regimes -----------------------------------------------
    std::printf("=== damped ring-down: three regimes ===\n");
    struct DampingCase {
        const char* name;
        double beta;
    };
    const DampingCase cases[3] = {
        {"underdamped", 0.3 * omega0},
        {"critically damped", 1.0 * omega0},
        {"overdamped", 3.0 * omega0},
    };

    SimOptions damp_opt;
    damp_opt.t0 = 0.0;
    damp_opt.t_end = 8.0;
    damp_opt.dt = 0.001;
    damp_opt.method = Method::RK4;

    std::vector<Trajectory> damp_trajs;
    for (const DampingCase& c : cases) {
        const forces::LinearDrag drag{2.0 * mass * c.beta};  // b = 2 m beta
        const auto force = forces::sum(shm_spring, drag);
        damp_trajs.push_back(simulate(Particle{mass, Vec3{x0, 0.0, 0.0}, Vec3{}}, force, damp_opt));

        const analytic::DampedOscillator osc{omega0, c.beta, x0, 0.0};
        std::printf("-- %s (beta = % .4g) --\n", c.name, c.beta);
        for (double t : {1.0, 2.0, 4.0}) {
            const std::size_t i = static_cast<std::size_t>(std::llround(t / damp_opt.dt));
            row("  numeric x [m]", damp_trajs.back().states[i].r.x);
            row("  analytic x [m]", osc.x(t));
        }
    }
    std::printf("\n");

    // One CSV, three columns — all three runs share the same t0/dt/t_end, so
    // their sample indices line up one-to-one.
    std::vector<std::vector<double>> damp_rows;
    for (std::size_t i = 0; i < damp_trajs[0].size(); ++i) {
        damp_rows.push_back({damp_trajs[0].states[i].t, damp_trajs[0].states[i].r.x,
                              damp_trajs[1].states[i].r.x, damp_trajs[2].states[i].r.x});
    }
    const bool wrote_damp_csv = io::write_csv(
        "build/damping_regimes.csv", {"t", "x_underdamped", "x_critical", "x_overdamped"}, damp_rows);
    std::printf("%-28s %s\n", "wrote damping_regimes.csv", wrote_damp_csv ? "ok" : "FAILED");
    std::printf("\n");

    // --- (3) driven resonance sweep -----------------------------------------------
    std::printf("=== driven oscillator: resonance sweep ===\n");
    const double beta_driven = 0.1 * omega0;  // light damping -> a sharp, visible peak
    const double f0 = 1.0;                    // driving acceleration amplitude, F0/m [m/s^2]

    // omega is mutated per sweep point below; resonant_omega()/quality_factor()
    // don't depend on it, so its initial value here is arbitrary.
    analytic::DrivenOscillator driven{omega0, beta_driven, f0, omega0};
    row("resonant omega [rad/s]", driven.resonant_omega());
    row("quality factor Q", driven.quality_factor());
    std::printf("\n");

    const int n_sweep = 40;
    const double omega_lo = 0.2 * omega0, omega_hi = 2.0 * omega0;
    std::vector<std::vector<double>> res_rows;
    double best_omega = 0.0, best_amp_numeric = -1.0;

    for (int i = 0; i <= n_sweep; ++i) {
        const double omega = omega_lo + (omega_hi - omega_lo) * static_cast<double>(i) / n_sweep;

        const forces::LinearDrag drag{2.0 * mass * beta_driven};
        const forces::SinusoidalDrive drive{f0 * mass, omega, 0.0, Vec3{1.0, 0.0, 0.0}};
        const auto force = forces::sum(shm_spring, drag, drive);

        SimOptions sweep_opt;
        sweep_opt.t0 = 0.0;
        sweep_opt.t_end = 50.0;  // >> 1/beta_driven, so the transient has decayed
        sweep_opt.dt = 0.01;
        sweep_opt.method = Method::RK4;

        const Trajectory traj = simulate(Particle{mass, Vec3{}, Vec3{}}, force, sweep_opt);
        const double amp_numeric = analysis::measure_amplitude(traj, 30.0);  // steady state only

        driven.omega = omega;
        const double amp_analytic = driven.amplitude();
        res_rows.push_back({omega, amp_numeric, amp_analytic});

        if (amp_numeric > best_amp_numeric) {
            best_amp_numeric = amp_numeric;
            best_omega = omega;
        }
    }

    const bool wrote_res_csv =
        io::write_csv("build/resonance.csv", {"omega", "amplitude_numeric", "amplitude_analytic"}, res_rows);
    std::printf("%-28s %s\n", "wrote resonance.csv", wrote_res_csv ? "ok" : "FAILED");
    row("numeric peak omega [rad/s]", best_omega);
    row("numeric peak amplitude [m]", best_amp_numeric);
    driven.omega = driven.resonant_omega();
    row("analytic peak amplitude [m]", driven.amplitude());

    return 0;
}
