// tests/test_analysis.cpp — newton::analysis checked against SYNTHETIC data
// built directly from closed-form formulas (not simulate() output): this
// file asks "if the analysis functions are handed exact physics, do they
// report exact physics back?" — isolating the analysis module's own
// correctness from any integrator error. (test_numeric_vs_analytic.cpp is
// the file that mixes analysis functions with real simulate() trajectories.)
#include <algorithm>
#include <cmath>
#include <cstddef>
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

    // ======================================================================
    // Energy bookkeeping with non-conservative and time-dependent pieces
    // (Taylor §4.1 work and power, §4.2 dE = W_nc, §4.5 dU/dt): U(t, r)
    // dispatch, power(), work_done(), energy_balance(). Expected values are
    // closed forms from the damped oscillator (§5.4), the linear-drag
    // projectile (§2.2) and the rotating uniform field, not from this code.
    // ======================================================================
    {
        // NaN-aware: std::max(m, NaN) silently DROPS the NaN, so an all-NaN
        // ledger would read 0 ("balanced"). Any NaN makes the result NaN, and
        // NaN < tol is false, so such a ledger fails every `max_abs(R) < tol`.
        const auto max_abs = [](const std::vector<double>& v) {
            double m = 0.0;
            for (double x : v) {
                if (std::isnan(x)) return std::nan("");
                m = std::max(m, std::fabs(x));
            }
            return m;
        };

        // --- U(t, r) dispatch in total_energy / energy_drift ---------------
        {
            auto U_static = [](const Vec3& r) { return 3.0 * r.x; };
            auto U_timed = [](double t, const Vec3& r) { return 3.0 * r.x + 5.0 * t; };
            const State s{Vec3{2.0, 0.0, 0.0}, Vec3{0.0, 4.0, 0.0}, 0.5};
            // T = 1/2 * 2 * 16 = 16; U = 6 (+ 2.5 when the time term is used).
            NT_CLOSE(analysis::total_energy(2.0, s, U_static), 22.0, 1e-14);
            NT_CLOSE(analysis::total_energy(2.0, s, U_timed), 24.5, 1e-14);

            Trajectory traj;
            traj.mass = 1.0;
            for (int i = 0; i <= 4; ++i) {
                traj.states.push_back(State{Vec3{}, Vec3{}, static_cast<double>(i)});
            }
            // E0 = U(0, 0) = 0, so the 1 J fallback scale applies: drift = 5 t.
            const std::vector<double> dr = analysis::energy_drift(traj, U_timed);
            for (int i = 0; i <= 4; ++i) {
                NT_CLOSE(dr[static_cast<std::size_t>(i)], 5.0 * i, 1e-13);
            }
        }

        // --- power(): F . v -------------------------------------------------
        {
            const State s{Vec3{1.0, 2.0, 3.0}, Vec3{0.3, -1.2, 0.7}, 0.0};
            const forces::MagneticForce mag{2.0, Vec3{0.4, -0.5, 1.5}};
            NT_CLOSE(analysis::power(mag, s), 0.0, 1e-15);  // v x B is perpendicular to v
            const double b = 0.4;
            NT_CLOSE(analysis::power(forces::LinearDrag{b}, s), -b * norm2(s.v), 1e-14);
            // A force along v: P = |F||v|.
            const auto push = [](double, const Vec3&, const Vec3& v) { return 2.0 * v; };
            NT_CLOSE(analysis::power(push, s), 2.0 * norm2(s.v), 1e-14);
        }

        // --- cumulative integrator: exactness, edge cases, odd/even counts --
        {
            // Power P(t) = t^2 via F = (t^2,0,0), v = (1,0,0). The Simpson
            // pair and the parabola-based odd interval are both exact for a
            // quadratic, even on NON-uniform times, so W = t^3/3 for every
            // sample count from 3 up (even and odd alike).
            const auto F2 = [](double t, const Vec3&, const Vec3&) { return Vec3{t * t, 0.0, 0.0}; };
            const double times[] = {0.0, 0.3, 0.5, 1.1, 1.2, 2.0, 2.9, 3.0};
            for (std::size_t n = 3; n <= 8; ++n) {
                Trajectory traj;
                for (std::size_t i = 0; i < n; ++i) {
                    traj.states.push_back(State{Vec3{}, Vec3{1.0, 0.0, 0.0}, times[i]});
                }
                const std::vector<double> W = analysis::work_done(traj, F2);
                NT_CHECK(W.size() == n);
                NT_CHECK(W[0] == 0.0);
                for (std::size_t i = 0; i < n; ++i) {
                    NT_CLOSE(W[i], times[i] * times[i] * times[i] / 3.0, 1e-13);
                }
            }
            // Uniform spacing, P = t^3: Simpson's rule is exact for cubics at
            // the even indices (W = t^4/4), for 6 (even) and 7 (odd) samples.
            const auto F3 = [](double t, const Vec3&, const Vec3&) { return Vec3{t * t * t, 0.0, 0.0}; };
            for (std::size_t n : {std::size_t{6}, std::size_t{7}}) {
                Trajectory traj;
                for (std::size_t i = 0; i < n; ++i) {
                    traj.states.push_back(State{Vec3{}, Vec3{1.0, 0.0, 0.0}, 0.5 * static_cast<double>(i)});
                }
                const std::vector<double> W = analysis::work_done(traj, F3);
                for (std::size_t i = 0; i < n; i += 2) {
                    const double t = traj.states[i].t;
                    NT_CLOSE(W[i], t * t * t * t / 4.0, 1e-13);
                }
            }
            // 0, 1 and 2 samples.
            Trajectory none;
            NT_CHECK(analysis::work_done(none, F2).empty());
            Trajectory one;
            one.states.push_back(State{Vec3{}, Vec3{1.0, 0.0, 0.0}, 1.0});
            const std::vector<double> W1 = analysis::work_done(one, F2);
            NT_CHECK(W1.size() == 1 && W1[0] == 0.0);
            Trajectory two;  // trapezoid: 1/2 (2.0 - 1.0)(1 + 4) = 2.5
            two.states.push_back(State{Vec3{}, Vec3{1.0, 0.0, 0.0}, 1.0});
            two.states.push_back(State{Vec3{}, Vec3{1.0, 0.0, 0.0}, 2.0});
            const std::vector<double> W2 = analysis::work_done(two, F2);
            NT_CHECK(W2.size() == 2 && W2[0] == 0.0);
            NT_CLOSE(W2[1], 2.5, 1e-14);
            // energy_balance on an empty trajectory is empty, not a crash.
            NT_CHECK(analysis::energy_balance(none, analysis::spring_U(1.0), F2).empty());
        }

        // --- damped oscillator, research-physics §8(c) item 1 ---------------
        // m = 1, k = 4, b = 0.4 (beta = 0.2, omega0 = 2), x0 = 1, v0 = 0,
        // E0 = 2 J, W_drag(t) = E(t) - E0 exactly. RK4 at dt = 1e-3, every
        // step recorded. Measured max error (even and odd samples) 2e-12;
        // tolerances below leave ~50x margin.
        const analytic::DampedOscillator damped{2.0, 0.2, 1.0, 0.0};
        const double k_spring = 4.0;
        const double b_drag = 0.4;
        const auto F_damped = forces::sum(forces::Spring{k_spring}, forces::LinearDrag{b_drag});
        const forces::LinearDrag drag{b_drag};
        const auto damped_E = [&](double t) {
            const double x = damped.x(t);
            const double v = damped.v(t);
            return 0.5 * v * v + 0.5 * k_spring * x * x;
        };
        const auto damped_W_err = [&](double dt, int stride) {
            SimOptions opt;
            opt.t_end = 10.0;
            opt.dt = dt;
            opt.record_every = stride;
            const Trajectory traj = simulate(Particle{1.0, Vec3{1.0, 0.0, 0.0}, Vec3{}}, F_damped, opt);
            const std::vector<double> W = analysis::work_done(traj, drag);
            double worst = 0.0;
            for (std::size_t i = 0; i < traj.size(); ++i) {
                worst = std::max(worst, std::fabs(W[i] - (damped_E(traj.states[i].t) - 2.0)));
            }
            return worst;
        };
        {
            SimOptions opt;
            opt.t_end = 10.0;
            opt.dt = 1e-3;
            const Trajectory traj = simulate(Particle{1.0, Vec3{1.0, 0.0, 0.0}, Vec3{}}, F_damped, opt);
            const std::vector<double> W = analysis::work_done(traj, drag);
            NT_CHECK(W.size() == traj.size());
            NT_CHECK(W[0] == 0.0);
            // The three tabulated values (samples 1000, 5000, 10000).
            NT_CHECK(std::fabs(W[1000] - (-0.7369477522)) < 1e-9);
            NT_CHECK(std::fabs(W[5000] - (-1.7043558284)) < 1e-9);
            NT_CHECK(std::fabs(W[10000] - (-1.9596345276)) < 1e-9);
            // Whole-run agreement with the closed form, odd samples included.
            double worst = 0.0;
            for (std::size_t i = 0; i < traj.size(); ++i) {
                worst = std::max(worst, std::fabs(W[i] - (damped_E(traj.states[i].t) - 2.0)));
            }
            NT_CHECK(worst < 1e-10);
            // Drag only removes energy: W never increases.
            for (std::size_t i = 1; i < W.size(); ++i) NT_CHECK(W[i] <= W[i - 1] + 1e-12);

            // Energy ledger: E_spring_U + drag work balances to ~1e-12 of E0.
            const std::vector<double> R = analysis::energy_balance(traj, analysis::spring_U(k_spring), drag);
            NT_CHECK(R.size() == traj.size());
            NT_CHECK(R[0] == 0.0);
            NT_CHECK(max_abs(R) < 1e-10);
            // ... and forgetting the drag in F_nc is caught loudly: energy
            // fell by ~1.96 J of E0 = 2 J, so the residual is ~1.
            const auto no_force = [](double, const Vec3&, const Vec3&) { return Vec3{}; };
            NT_CHECK(max_abs(analysis::energy_balance(traj, analysis::spring_U(k_spring), no_force)) > 0.9);
            // Plain drift sees the same loss as "drift" (it cannot tell
            // physics from error): ~ -0.98 at t = 10.
            NT_CLOSE(analysis::energy_drift(traj, analysis::spring_U(k_spring)).back(),
                     (damped_E(10.0) - 2.0) / 2.0, 1e-9);
        }

        // --- 4th-order convergence: halving dt divides the W error by ~16 ---
        {
            const double e1 = damped_W_err(8e-3, 1);
            const double e2 = damped_W_err(4e-3, 1);
            const double e3 = damped_W_err(2e-3, 1);
            const double r12 = e1 / e2;
            const double r23 = e2 / e3;
            // Measured 16.0 and 16.0 (a trapezoid would give ~4).
            NT_CHECK(r12 > 12.0 && r12 < 20.0);
            NT_CHECK(r23 > 12.0 && r23 < 20.0);
        }

        // --- record_every that does not divide n ---------------------------
        // dt = 0.01, n = 1000, stride 7: samples 0, 7, ..., 994, then the
        // forced last sample at step 1000 (a gap of 6 steps, not 7).
        {
            SimOptions opt;
            opt.t_end = 10.0;
            opt.dt = 0.01;
            opt.record_every = 7;
            const Trajectory traj = simulate(Particle{1.0, Vec3{1.0, 0.0, 0.0}, Vec3{}}, F_damped, opt);
            NT_CHECK(traj.size() == 144);  // 1 + floor(1000/7) = 143, plus the final sample
            NT_CLOSE(traj.states.back().t, 10.0, 1e-12);
            NT_CHECK(traj.states[142].t - traj.states[141].t > traj.states.back().t - traj.states[142].t);  // short last gap
            const std::vector<double> W = analysis::work_done(traj, drag);
            double worst = 0.0;
            for (std::size_t i = 0; i < traj.size(); ++i) {
                worst = std::max(worst, std::fabs(W[i] - (damped_E(traj.states[i].t) - 2.0)));
            }
            NT_CHECK(worst < 2e-4);  // measured 5.0e-5: O(h^4) at h = 0.07 s (~90 samples per period)
            NT_CHECK(worst > 1e-12);  // and it really is the coarse grid being measured
            NT_CHECK(max_abs(analysis::energy_balance(traj, analysis::spring_U(k_spring), drag)) < 1e-4);  // measured 2.5e-5
        }

        // --- linear-drag projectile, research-physics §8(c) item 2 ----------
        // m = 1, b = 0.5, 20 m/s at 45 deg, U = m g y, y0 = 0, E0 = 200 J.
        {
            const double c45 = 20.0 / std::sqrt(2.0);
            const Particle ball{1.0, Vec3{}, Vec3{c45, c45, 0.0}};
            const forces::LinearDrag air{0.5};
            const auto F_proj = forces::sum(forces::UniformGravity{1.0}, air);

            SimOptions opt;
            opt.t_end = 2.0;
            opt.dt = 1e-3;
            const Trajectory traj = simulate(ball, F_proj, opt);
            const std::vector<double> W = analysis::work_done(traj, air);
            NT_CHECK(std::fabs(W[1000] - (-94.6841394160)) < 1e-8);
            NT_CHECK(std::fabs(W[2000] - (-126.7624560775)) < 1e-8);
            const auto Ugrav = analysis::uniform_gravity_U(1.0, constants::g_earth);
            NT_CHECK(max_abs(analysis::energy_balance(traj, Ugrav, air)) < 1e-9);

            // Impact via simulate_until: the last interval is short and the
            // event time is bisected, so this also exercises uneven gaps.
            SimOptions opt2;
            opt2.t_end = 10.0;
            opt2.dt = 1e-3;
            const Trajectory hit = simulate_until(ball, F_proj, opt2, events::ground_impact());
            const std::vector<double> Wh = analysis::work_done(hit, air);
            NT_CHECK(std::fabs(hit.states.back().t - 2.4110895210) < 1e-8);
            NT_CHECK(std::fabs(Wh.back() - (-145.8787809024)) < 1e-8);
            // KE at impact = E0 + W = 54.1212190976 J (U = 0 at y = 0).
            NT_CHECK(std::fabs(analysis::kinetic_energy(1.0, hit.states.back().v) - 54.1212190976) < 1e-8);
            NT_CHECK(max_abs(analysis::energy_balance(hit, Ugrav, air)) < 1e-9);
        }

        // --- time-dependent U: rotating uniform field ----------------------
        // F = (cos t, sin t, 0), m = 1, from rest at the origin, gauge
        // U(t, r) = -F(t) . r. Then E(t) = 2 - 2 cos t - t sin t (§3.3 A),
        // and E(t) - E(0) = integral of dU/dt, with NO non-conservative force.
        {
            const auto F_rot = [](double t, const Vec3&, const Vec3&) {
                return Vec3{std::cos(t), std::sin(t), 0.0};
            };
            const auto U_rot = [](double t, const Vec3& r) {
                return -(std::cos(t) * r.x + std::sin(t) * r.y);
            };
            const auto zero_force = [](double, const Vec3&, const Vec3&) { return Vec3{}; };
            SimOptions opt;
            opt.t_end = 3.5 * constants::pi;
            opt.dt = 1e-3;
            const Trajectory traj = simulate(Particle{1.0, Vec3{}, Vec3{}}, F_rot, opt);

            // E(t) against the closed form, at every sample's own time.
            double worst_E = 0.0;
            for (const State& s : traj.states) {
                const double E_exact = 2.0 - 2.0 * std::cos(s.t) - s.t * std::sin(s.t);
                worst_E = std::max(worst_E, std::fabs(analysis::total_energy(1.0, s, U_rot) - E_exact));
            }
            NT_CHECK(worst_E < 1e-10);

            // The ledger balances only because the dU/dt term is included.
            const std::vector<double> R = analysis::energy_balance(traj, U_rot, zero_force);
            NT_CHECK(R.size() == traj.size());
            NT_CHECK(R[0] == 0.0);
            // E0 = 0 and W_nc = 0 here, so the scale is max|E_i - E_0| ~ 13 J (the largest
            // energy that moves); measured residual ~2e-12.
            NT_CHECK(max_abs(R) < 1e-9);
            // Plain drift is huge here: it is real physics (the field does work on the particle as
            // it rotates), not integrator error.
            const std::vector<double> drift = analysis::energy_drift(traj, U_rot);
            NT_CHECK(max_abs(drift) > 1.0);
            // A static U(r) cannot represent this force, and the ledger knows.
            // (Not -x: for this field T = 1 - cos t = x, so U = -x would balance by coincidence.)
            const auto U_frozen = [](const Vec3& r) { return -r.y; };
            NT_CHECK(max_abs(analysis::energy_balance(traj, U_frozen, zero_force)) > 0.1);
        }

        // --- B1: dU/dt step must follow U's time scale, not the run ---------
        // A driven spring written as U(t, r) = 1/2 k x^2 - F0 cos(w t) x,
        // drive period 1 s, simulated with exactly -dU/dx, so the true
        // ledger is 0. A step tied to |t| or the run length (1e-4 * 1000 s =
        // 0.1 s = a tenth of the period) gave 6e-3 to 4e-2 here.
        {
            const double k = 4.0, F0 = 1.0, w = 2.0 * constants::pi;
            const auto U = [=](double t, const Vec3& r) { return 0.5 * k * r.x * r.x - F0 * std::cos(w * t) * r.x; };
            const auto F = [=](double t, const Vec3& r, const Vec3&) {
                return Vec3{-k * r.x + F0 * std::cos(w * t), 0.0, 0.0};
            };
            const auto zero = [](double, const Vec3&, const Vec3&) { return Vec3{}; };
            const double cases[][2] = {{0.0, 1000.0}, {1000.0, 10.0}, {1000.0, 1000.0}};
            for (const auto& c : cases) {
                SimOptions opt;
                opt.t0 = c[0];
                opt.t_end = c[0] + c[1];
                opt.dt = 1e-3;
                const Trajectory traj = simulate(Particle{1.0, Vec3{0.5, 0.0, 0.0}, Vec3{}}, F, opt);
                // Measured <= 4e-10 for all three (was 5.8e-3, 5.9e-3, 4.3e-2).
                NT_CHECK(max_abs(analysis::energy_balance(traj, U, zero)) < 1e-8);
            }
        }

        // --- B2: simulate_until landing on a step boundary ------------------
        // Vacuum projectile whose flight time is exactly N dt: the bisected
        // last sample lands a hair past the previous one (gap ~ 1e-12 or
        // less), which made the Simpson pair lose ~eps*h1/h2. Gravity is the
        // net force, so W_gravity must equal delta T (RK4 is exact for a
        // constant force). Measured worst |W - dT| over N = 900..1100:
        // about 1e-14 J (was 7.7e-7 J).
        {
            const double g = constants::g_earth, dt = 1e-3;
            const forces::UniformGravity grav{1.0};
            double worst = 0.0;
            int n_tiny = 0;
            for (int N = 900; N <= 1100; ++N) {
                const double vy0 = 0.5 * g * N * dt;
                SimOptions opt;
                opt.t_end = 5.0;
                opt.dt = dt;
                const Trajectory tr = simulate_until(Particle{1.0, Vec3{}, Vec3{3.0, vy0, 0.0}}, grav, opt,
                                                     events::ground_impact());
                const std::vector<double> W = analysis::work_done(tr, grav);
                const double dT = 0.5 * (norm2(tr.states.back().v) - norm2(tr.states.front().v));
                if (tr.states.back().t - tr.states[tr.size() - 2].t < 1e-9 * dt) ++n_tiny;
                worst = std::max(worst, std::fabs(W.back() - dT));
            }
            NT_CHECK(n_tiny > 100);  // the sweep really does produce lopsided last gaps
            NT_CHECK(worst < 1e-11);
            // Direct synthetic version: f = cos t at dt = 1e-3, then one last
            // sample 1e-12 s later (even last index). Integral is sin t.
            std::vector<double> ts;
            for (int i = 0; i <= 999; ++i) ts.push_back(1e-3 * i);
            ts.push_back(0.999 + 1e-12);
            Trajectory traj;
            for (double t : ts) traj.states.push_back(State{Vec3{}, Vec3{1.0, 0.0, 0.0}, t});
            const auto Fc = [](double t, const Vec3&, const Vec3&) { return Vec3{std::cos(t), 0.0, 0.0}; };
            const std::vector<double> W = analysis::work_done(traj, Fc);
            NT_CHECK(traj.size() % 2 == 1);  // last index is even
            NT_CHECK(std::fabs(W.back() - std::sin(ts.back())) < 1e-13);  // was 8e-12
        }

        // --- B3: residual scale includes the dU/dt term ----------------------
        // Rotating field scaled down to F0 = 1e-6 N: E0 = 0, W_nc = 0. The
        // right U balances (~1e-12 of the energy that moved); the wrong static
        // U = -F0 y must read ~1, not ~1e-11 "joules".
        {
            const double F0 = 1e-6;
            const auto F_rot = [=](double t, const Vec3&, const Vec3&) {
                return Vec3{F0 * std::cos(t), F0 * std::sin(t), 0.0};
            };
            const auto U_rot = [=](double t, const Vec3& r) {
                return -F0 * (std::cos(t) * r.x + std::sin(t) * r.y);
            };
            const auto U_wrong = [=](const Vec3& r) { return -F0 * r.y; };
            const auto zero_force = [](double, const Vec3&, const Vec3&) { return Vec3{}; };
            SimOptions opt;
            opt.t_end = 3.5 * constants::pi;
            opt.dt = 1e-3;
            const Trajectory traj = simulate(Particle{1.0, Vec3{}, Vec3{}}, F_rot, opt);
            NT_CHECK(max_abs(analysis::energy_balance(traj, U_rot, zero_force)) < 1e-9);
            const double wrong = max_abs(analysis::energy_balance(traj, U_wrong, zero_force));
            NT_CHECK(wrong > 0.5 && wrong < 2.0);  // measured ~1.0 at F0 = 1 N and 1e-6 N alike
        }

        // --- NaN-aware max_abs, and a NaN force is never "balanced" ---------
        {
            NT_CHECK(std::isnan(max_abs({0.0, std::nan(""), 1.0})));
            const auto nan_force = [](double, const Vec3&, const Vec3&) { return Vec3{std::nan(""), 0.0, 0.0}; };
            SimOptions opt;
            opt.t_end = 1.0;
            opt.dt = 1e-2;
            const Trajectory traj = simulate(Particle{1.0, Vec3{}, Vec3{1.0, 0.0, 0.0}},
                                             [](double, const Vec3&, const Vec3&) { return Vec3{}; }, opt);
            NT_CHECK(!(max_abs(analysis::energy_balance(traj, analysis::spring_U(1.0), nan_force)) < 1e-3));
        }

        // --- correct ledgers with E0 = 0 must read ~0 -----------------------
        {
            // Free fall from rest at y = 0, U = m g y, F_nc = 0. E0 = 0, nothing
            // crosses the ledger, and T_end ~ 4800 J. (A scale containing
            // |E_i - E_0| would cancel the tiny error and read 1.0.)
            const double g = constants::g_earth;
            const auto U = analysis::uniform_gravity_U(1.0, g);
            const auto none = [](double, const Vec3&, const Vec3&) { return Vec3{}; };
            SimOptions opt;
            opt.t_end = 10.0;
            opt.dt = 1e-3;
            const Trajectory traj = simulate(Particle{1.0, Vec3{}, Vec3{}}, forces::UniformGravity{1.0}, opt);
            NT_CHECK(analysis::kinetic_energy(1.0, traj.states.back().v) > 4000.0);
            NT_CHECK(max_abs(analysis::energy_balance(traj, U, none)) < 1e-10);

            // E x B drift from rest: q = 1, E = (0, 0.5, 0), B = z, m = 1,
            // U = -q E.r, F_nc = the magnetic part. E0 = 0, T <= 0.5 J.
            const forces::LorentzForce lor{1.0, Vec3{0.0, 0.5, 0.0}, Vec3{0.0, 0.0, 1.0}};
            const forces::MagneticForce mag{1.0, Vec3{0.0, 0.0, 1.0}};
            const auto U_E = [](const Vec3& r) { return -0.5 * r.y; };
            SimOptions o2;
            o2.t_end = 4.0 * constants::pi;
            o2.dt = 1e-3;
            const Trajectory drift = simulate(Particle{1.0, Vec3{}, Vec3{}}, lor, o2);
            const double drift_speed = drift.states.back().r.x / drift.states.back().t;
            NT_CLOSE(drift_speed, 0.5, 1e-3);  // E x B / B^2 = 0.5 m/s on average
            // W_mag ~ 0 to roundoff, and the ledger balances.
            NT_CHECK(max_abs(analysis::work_done(drift, mag)) < 1e-12);
            NT_CHECK(max_abs(analysis::energy_balance(drift, U_E, mag)) < 1e-9);
            // The whole Lorentz force as F_nc with U = 0 balances too.
            NT_CHECK(max_abs(analysis::energy_balance(drift, [](const Vec3&) { return 0.0; }, lor)) < 1e-9);
            // Forgetting the magnetic part is NOT a mistake the ledger sees:
            // it does no work. Forgetting the electric part is.
            NT_CHECK(max_abs(analysis::energy_balance(drift, U_E, none)) < 1e-9);
            NT_CHECK(max_abs(analysis::energy_balance(drift, [](const Vec3&) { return 0.0; }, none)) > 0.5);
        }

        // --- quadratic drag fall vs closed form (research §8(c) item 3) -----
        // m = 1, c = 0.1, from rest: E - E0 = m v_ter^2 [tanh^2(u)/2 - ln cosh u],
        // u = g t / v_ter, v_ter = sqrt(m g / c); equals W_drag = -c int |v|^3 dt.
        {
            const double g = constants::g_earth, c = 0.1;
            const double v_ter = std::sqrt(g / c);
            const forces::QuadraticDrag qd{c};
            SimOptions opt;
            opt.t_end = 3.0;
            opt.dt = 1e-3;
            const Trajectory traj = simulate(Particle{1.0, Vec3{}, Vec3{}},
                                             forces::sum(forces::UniformGravity{1.0}, qd), opt);
            const std::vector<double> W = analysis::work_done(traj, qd);
            double worst = 0.0;
            for (std::size_t i = 0; i < traj.size(); ++i) {
                const double u = g * traj.states[i].t / v_ter;
                const double exact = v_ter * v_ter * (0.5 * std::tanh(u) * std::tanh(u) - std::log(std::cosh(u)));
                worst = std::max(worst, std::fabs(W[i] - exact));
            }
            NT_CHECK(worst < 1e-8);
            NT_CHECK(std::fabs(W[1000] - (-13.6813677492)) < 1e-8);
            NT_CHECK(std::fabs(W[3000] - (-175.1036576962)) < 1e-7);
            NT_CHECK(max_abs(analysis::energy_balance(traj, analysis::uniform_gravity_U(1.0, g), qd)) < 1e-9);
        }

        // --- driven damped steady state: W_drive = -W_drag = pi b w A^2 -----
        // Over one drive period in the steady state E repeats, so the drive
        // pumps in exactly what the drag takes out (Taylor §5.5-5.6).
        {
            const double m = 1.0, k = 4.0, b = 0.5, F0 = 1.0, w = 1.8;
            const analytic::DrivenOscillator osc{2.0, b / (2.0 * m), F0 / m, w};
            const double A = osc.amplitude();
            const double delta = osc.phase();
            const double T = 2.0 * constants::pi / w;
            // Start ON the steady state (no transient): x = A cos(wt - delta).
            const Vec3 r0{A * std::cos(-delta), 0.0, 0.0};
            const Vec3 v_start{A * w * std::sin(delta), 0.0, 0.0};  // d/dt A cos(wt - delta) at 0
            const forces::SinusoidalDrive drive{F0, w, 0.0, Vec3{1.0, 0.0, 0.0}};
            const forces::LinearDrag damping{b};  // not "drag": that name is taken in an outer scope
            SimOptions opt;
            opt.t_end = T;
            opt.dt = T / 2000.0;  // exactly 2000 steps
            const Trajectory traj = simulate(Particle{m, r0, v_start},
                                             forces::sum(forces::Spring{k}, damping, drive), opt);
            const double expected = constants::pi * b * w * A * A;
            const std::vector<double> Wd = analysis::work_done(traj, drive);
            const std::vector<double> Wg = analysis::work_done(traj, damping);
            NT_CHECK(std::fabs(Wd.back() - expected) < 1e-9 * expected);
            NT_CHECK(std::fabs(Wg.back() + expected) < 1e-9 * expected);
            // Ledger with the drive and the drag both booked as F_nc, U = spring.
            const auto F_nc = forces::sum(damping, drive);
            NT_CHECK(max_abs(analysis::energy_balance(traj, analysis::spring_U(k), F_nc)) < 1e-9);
            // E is periodic: E(T) = E(0).
            NT_CLOSE(analysis::total_energy(m, traj.states.back(), analysis::spring_U(k)),
                     analysis::total_energy(m, traj.states.front(), analysis::spring_U(k)), 1e-9);
        }

        // --- an energy-FEEDING force (+b v): W > 0, E grows, ledger closes --
        {
            const forces::LinearDrag feed{-0.2};  // F = +0.2 v
            SimOptions opt;
            opt.t_end = 10.0;
            opt.dt = 1e-3;
            const Trajectory traj = simulate(Particle{1.0, Vec3{1.0, 0.0, 0.0}, Vec3{}},
                                             forces::sum(forces::Spring{4.0}, feed), opt);
            const std::vector<double> W = analysis::work_done(traj, feed);
            NT_CHECK(W.back() > 1.0);
            const auto U = analysis::spring_U(4.0);
            NT_CHECK(analysis::total_energy(1.0, traj.states.back(), U) > 2.0);
            NT_CHECK(max_abs(analysis::energy_balance(traj, U, feed)) < 1e-9);
            // Booking it with the wrong sign is caught (residual ~ +-1).
            NT_CHECK(max_abs(analysis::energy_balance(traj, U, forces::LinearDrag{0.2})) > 0.5);
        }

        // --- gauge invariance: U + 1e9 J changes nothing ---------------------
        {
            const auto F_damped2 = forces::sum(forces::Spring{4.0}, forces::LinearDrag{0.4});
            const forces::LinearDrag drag2{0.4};
            const auto none = [](double, const Vec3&, const Vec3&) { return Vec3{}; };
            SimOptions opt;
            opt.t_end = 10.0;
            opt.dt = 1e-3;
            const Trajectory traj = simulate(Particle{1.0, Vec3{1.0, 0.0, 0.0}, Vec3{}}, F_damped2, opt);
            const auto U = analysis::spring_U(4.0);
            const auto U_big = [](const Vec3& r) { return 0.5 * 4.0 * norm2(r) + 1e9; };
            const std::vector<double> R = analysis::energy_balance(traj, U, drag2);
            const std::vector<double> Rb = analysis::energy_balance(traj, U_big, drag2);
            double diff = 0.0;
            for (std::size_t i = 0; i < R.size(); ++i) diff = std::max(diff, std::fabs(R[i] - Rb[i]));
            NT_CHECK(diff < 1e-6);  // only the 1e9 J roundoff (~1e-7 J absolute) differs
            NT_CHECK(max_abs(Rb) < 1e-6);
            // A forgotten drag is still detected in the shifted gauge.
            NT_CHECK(max_abs(analysis::energy_balance(traj, U_big, none)) > 0.5);
        }

        // --- scale mutants (review S2): the |W_nc| term of the scale --------
        {
            // (i) At rest at U = 0, pushed by a constant F_nc: E0 = 0, the
            // ledger is exact, residual ~ roundoff relative to max|W| = |dT|.
            const auto push = [](double, const Vec3&, const Vec3&) { return Vec3{3.0, 0.0, 0.0}; };
            SimOptions opt;
            opt.t_end = 2.0;
            opt.dt = 1e-2;
            const Trajectory traj = simulate(Particle{1.0, Vec3{}, Vec3{}}, push, opt);
            const auto U0 = [](const Vec3&) { return 0.0; };
            NT_CHECK(max_abs(analysis::energy_balance(traj, U0, push)) < 1e-12);

            // (ii) A particle coasting at 1 m/s with NO force, but the caller
            // claims F_nc = 5 N along v. W_nc grows to 10 J while E never
            // moves (E0 = 0.5 J), so the right answer is "10 J of work
            // missing from the ledger": residual = -W/max|W| = -1 at the end.
            // If max|W_nc| were left out of the scale it would read -20.
            const auto none = [](double, const Vec3&, const Vec3&) { return Vec3{}; };
            const auto claimed = [](double, const Vec3&, const Vec3&) { return Vec3{5.0, 0.0, 0.0}; };
            const Trajectory coast = simulate(Particle{1.0, Vec3{}, Vec3{1.0, 0.0, 0.0}}, none, opt);
            const std::vector<double> R = analysis::energy_balance(coast, U0, claimed);
            NT_CLOSE(R.back(), -1.0, 1e-9);
            NT_CHECK(max_abs(R) < 1.0 + 1e-9);
        }
    }

    return newton_test::report("analysis");
}
