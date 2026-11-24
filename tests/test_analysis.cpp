// tests/test_analysis.cpp — newton::analysis checked against SYNTHETIC data
// built directly from closed-form formulas (not simulate() output): this
// file asks "if the analysis functions are handed exact physics, do they
// report exact physics back?" — isolating the analysis module's own
// correctness from any integrator error. (test_numeric_vs_analytic.cpp is
// the file that mixes analysis functions with real simulate() trajectories.)
#include <cmath>
#include <vector>

#include "test_util.hpp"
#include <newton/newton.hpp>

using namespace newton;

int main() {
    // ======================================================================
    // energy_drift() on exact SHM data ~= 0 everywhere: sample
    // analytic::SHM::x(t)/v(t) directly into a Trajectory (no integrator
    // involved at all), then check the spring-potential energy computed from
    // those exact samples is constant to floating-point precision — SHM is
    // conservative, so this is really testing kinetic_energy()/
    // total_energy()/energy_drift()'s own arithmetic, not any integrator.
    // ======================================================================
    {
        const double mass = 2.0;
        const double omega0 = 1.7;
        const double k = mass * omega0 * omega0;  // F = -k x = m x'' = -m omega0^2 x
        analytic::SHM shm{omega0, 1.5, -0.8};

        Trajectory traj;
        traj.mass = mass;
        const double T = shm.period();
        const int n_samples = 500;
        for (int i = 0; i <= n_samples; ++i) {
            // 3 periods of exact analytic samples, evenly spaced.
            const double t = 3.0 * T * static_cast<double>(i) / n_samples;
            traj.states.push_back(State{Vec3{shm.x(t), 0.0, 0.0}, Vec3{shm.v(t), 0.0, 0.0}, t});
        }

        auto U = analysis::spring_U(k);
        const std::vector<double> drift = analysis::energy_drift(traj, U);
        NT_CHECK(drift.size() == traj.size());
        for (double d : drift) {
            NT_CLOSE(d, 0.0, 1e-9);  // exact analytic samples: only trig-evaluation roundoff
        }

        // kinetic_energy() and total_energy() agree with the closed-form
        // energy() at t=0, where U(x0) = 1/2 k x0^2 = 1/2 m omega0^2 x0^2.
        const double T0 = analysis::kinetic_energy(mass, traj.front().v);
        const double E0 = analysis::total_energy(mass, traj.front(), U);
        NT_CLOSE(T0, 0.5 * mass * shm.v0 * shm.v0, 1e-12);
        NT_CLOSE(E0, shm.energy(mass), 1e-9);
    }

    // ======================================================================
    // find_peaks() and measure_period() on synthetic x(t) = cos(omega t)
    // data, sampled finely — no simulate() involved, so this isolates the
    // peak-finding / zero-crossing logic itself.
    // ======================================================================
    {
        const double omega = 2.3;
        const double T = 2.0 * constants::pi / omega;
        const int per_period = 2000;   // dense sampling: interpolation error negligible
        const int n_periods = 5;

        Trajectory traj;
        traj.mass = 1.0;
        for (int i = 0; i <= n_periods * per_period; ++i) {
            const double t = n_periods * T * static_cast<double>(i) / (n_periods * per_period);
            traj.states.push_back(State{Vec3{std::cos(omega * t), 0.0, 0.0}, Vec3{}, t});
        }

        // Peaks of cos(omega t) occur at t = n*T, x = 1.0. The very first
        // (t=0) and very last (t=5T) samples sit at the array's edges, where
        // find_peaks()'s 3-point comparison can't see both neighbours (see
        // oscillation.hpp), so only the interior peaks at T, 2T, 3T, 4T are
        // guaranteed to be found.
        const analysis::PeakList peaks = analysis::find_peaks(traj);
        NT_CHECK(peaks.t.size() >= 3);  // at least T, 2T, 3T should all be interior
        for (std::size_t i = 0; i < peaks.t.size(); ++i) {
            NT_CLOSE(peaks.x[i], 1.0, 1e-6);
            // Every found peak time should land close to an integer multiple of T.
            const double n_est = std::round(peaks.t[i] / T);
            NT_CLOSE(peaks.t[i], n_est * T, 1e-3 * T);
        }
        // Consecutive peaks are one period apart.
        for (std::size_t i = 1; i < peaks.t.size(); ++i) {
            NT_CLOSE(peaks.t[i] - peaks.t[i - 1], T, 1e-3 * T);
        }

        // measure_period(): mean spacing of upward zero crossings should
        // recover T to interpolation accuracy, far tighter than the
        // find_peaks() tolerance above since it averages over many crossings.
        const double period_measured = analysis::measure_period(traj);
        NT_CLOSE(period_measured, T, 1e-4 * T);

        // measure_amplitude(): with x(t) = cos(omega t), the max |x| after
        // any t_start is exactly 1.0 (there's always another peak ahead
        // within the sampled window).
        const double amp = analysis::measure_amplitude(traj, 0.0);
        NT_CLOSE(amp, 1.0, 1e-6);
    }

    // ======================================================================
    // angular_momentum() of uniform circular motion is constant in time:
    // r(t) = R(cos wt, sin wt, 0), v(t) = R w(-sin wt, cos wt, 0), built
    // directly from the parametric formula (no central-force simulation
    // needed to make the point) — L = m r x v = m R^2 w z_hat at every t,
    // Taylor §3.4's textbook example of angular momentum for circular orbits.
    // ======================================================================
    {
        const double mass = 3.0;
        const double R = 2.0;
        const double w = 1.4;
        const Vec3 L_expected{0.0, 0.0, mass * R * R * w};

        for (double t : {0.0, 0.3, 1.1, 2.7, 5.0}) {
            const Vec3 r{R * std::cos(w * t), R * std::sin(w * t), 0.0};
            const Vec3 v{-R * w * std::sin(w * t), R * w * std::cos(w * t), 0.0};
            State s{r, v, t};
            const Vec3 L = analysis::angular_momentum(mass, s);
            NT_CLOSE_VEC3(L, L_expected, 1e-9);
        }

        // angular_momentum() about a point other than the origin is NOT
        // conserved for this (non-central-about-that-point) motion in
        // general — sanity check that "about" actually shifts the reference
        // point rather than being silently ignored: about = the circle's own
        // center (origin) reproduces L_expected exactly (explicit vs default
        // argument agree), while about = a point ON the circle changes it.
        const State s0{Vec3{R, 0.0, 0.0}, Vec3{0.0, R * w, 0.0}, 0.0};
        const Vec3 L_about_origin = analysis::angular_momentum(mass, s0, Vec3{});
        NT_CLOSE_VEC3(L_about_origin, L_expected, 1e-9);
        const Vec3 L_about_edge = analysis::angular_momentum(mass, s0, Vec3{R, 0.0, 0.0});
        NT_CHECK(norm(L_about_edge) < norm(L_expected));  // r - about = 0 here: L collapses toward zero
    }

    // ======================================================================
    // center_of_mass() and total_momentum(): a two-particle sanity check,
    // independent of the oscillator/circular-motion fixtures above.
    // ======================================================================
    {
        std::vector<Particle> particles;
        particles.push_back(Particle{1.0, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}});
        particles.push_back(Particle{3.0, Vec3{4.0, 0.0, 0.0}, Vec3{-1.0, 0.0, 0.0}});
        // R = (1*0 + 3*4)/(1+3) = 3.0 along x.
        const Vec3 com = analysis::center_of_mass(particles);
        const Vec3 com_expected{3.0, 0.0, 0.0};
        NT_CLOSE_VEC3(com, com_expected, 1e-12);
        // P = 1*(1,0,0) + 3*(-1,0,0) = (-2,0,0).
        const Vec3 P = analysis::total_momentum(particles);
        const Vec3 P_expected{-2.0, 0.0, 0.0};
        NT_CLOSE_VEC3(P, P_expected, 1e-12);

        // Empty / massless input is zero-guarded, not NaN (mirrors
        // Vec3::normalized()'s zero-guard convention).
        const Vec3 com_empty = analysis::center_of_mass({});
        NT_CHECK(com_empty == Vec3{});
        NT_CHECK(!std::isnan(com_empty.x));
    }

    return newton_test::report("analysis");
}
