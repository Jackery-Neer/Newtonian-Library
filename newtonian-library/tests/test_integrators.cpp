// tests/test_integrators.cpp — the numerics contract every integrator and
// the simulate() driver must satisfy, independent of any one force law:
//   (a) constant-force exactness: Verlet & RK4 reproduce ½at² to ~1e-12 at
//       any dt; Euler and semi-implicit Euler do NOT (their O(dt) offset,
//       research numerical-integration.md §4.1(a));
//   (b) SHM convergence order: halving dt should shrink the error by
//       2^p, p = 1 (Euler, semi-implicit), 2 (Verlet), 4 (RK4);
//   (c) SHM energy behaviour over 100 periods: Euler grows without bound,
//       semi-implicit/Verlet stay bounded, RK4 decays negligibly;
//   (d) simulate()'s time bookkeeping: the final sample's t equals t_end to
//       1e-12 (research §7.2 — t = t0 + i*dt every iteration, never
//       accumulated).
#include <cmath>
#include <vector>

#include "test_util.hpp"
#include <newton/newton.hpp>

using namespace newton;

namespace {

/// Runs simulate() on the pure SHM force (Spring with anchor at the origin,
/// natural_length 0, k = m*omega0^2 collapses exactly to F = -k r — see
/// forces/spring.hpp) and returns the max |x_numeric(t) - x_analytic(t)|
/// over every recorded sample, per research §4.2's "error norm over the
/// whole trajectory, not a single endpoint."
double shm_max_error(Method method, double dt, double omega0, double T_end) {
    Particle p{1.0, Vec3{1.0, 0.0, 0.0}, Vec3{}};  // x0=1, v0=0, m=1
    forces::Spring spring{omega0 * omega0, Vec3{}, 0.0};
    SimOptions opt;
    opt.t0 = 0.0;
    opt.t_end = T_end;
    opt.dt = dt;
    opt.method = method;
    Trajectory traj = simulate(p, spring, opt);

    analytic::SHM shm{omega0, 1.0, 0.0};
    double e = 0.0;
    for (const State& s : traj.states) {
        e = std::fmax(e, std::fabs(s.r.x - shm.x(s.t)));
    }
    return e;
}

/// Observed order p = log2(e(dt) / e(dt/2)) between the two smallest dt in a
/// 5-level halving sequence starting at dt0 — using the finest pair keeps
/// the estimate deep in the asymptotic (small-dt) regime, per research
/// §4.2's warning that coarse dt hasn't converged to the leading-order
/// behaviour yet (verified empirically for this exact setup: e.g. Euler's
/// order estimate is 1.27 at the coarsest pair but 1.02 by the finest).
double observed_order(Method method, double omega0, double T_end, double dt0) {
    double e_prev = shm_max_error(method, dt0 / 8.0, omega0, T_end);
    double e_last = shm_max_error(method, dt0 / 16.0, omega0, T_end);
    return std::log2(e_prev / e_last);
}

}  // namespace

int main() {
    // ======================================================================
    // (a) Constant-force exactness (research §4.1(a)).
    // ======================================================================
    // a = F/m is the same vector at every step here, so the true motion is
    // exactly r(t) = r0 + v0 t + 0.5 a t^2, v(t) = v0 + a t. Verlet and RK4
    // must reproduce this to roundoff at ANY dt; Euler and semi-implicit
    // Euler have a genuine O(dt) truncation error even for constant a (the
    // ½dt² term is simply missing from their update equations), so they
    // must NOT match to tight tolerance.
    {
        const double mass = 2.0;
        const Vec3 F{3.0, -5.0, 7.0};
        forces::ConstantForce force{F};
        const Vec3 a = F / mass;
        const Vec3 r0{0.0, 0.0, 0.0};
        const Vec3 v0{1.0, 2.0, 3.0};
        const double dt = 0.1;
        const int n = 10;  // T = 1.0 s

        State s_euler{r0, v0, 0.0};
        State s_semi{r0, v0, 0.0};
        State s_verlet{r0, v0, 0.0};
        State s_rk4{r0, v0, 0.0};
        for (int i = 0; i < n; ++i) {
            s_euler = integrators::euler_step(s_euler, mass, dt, force);
            s_semi = integrators::semi_implicit_euler_step(s_semi, mass, dt, force);
            s_verlet = integrators::velocity_verlet_step(s_verlet, mass, dt, force);
            s_rk4 = integrators::rk4_step(s_rk4, mass, dt, force);
        }

        const double T = n * dt;
        const Vec3 r_exact = r0 + v0 * T + 0.5 * T * T * a;
        const Vec3 v_exact = v0 + a * T;

        // Verlet and RK4: exact to roundoff (1e-12) at dt = 0.1.
        NT_CLOSE_VEC3(s_verlet.r, r_exact, 1e-12);
        NT_CLOSE_VEC3(s_verlet.v, v_exact, 1e-12);
        NT_CLOSE_VEC3(s_rk4.r, r_exact, 1e-12);
        NT_CLOSE_VEC3(s_rk4.v, v_exact, 1e-12);

        // Forward Euler and semi-implicit Euler: NOT exact at this dt — the
        // position error must be well above roundoff (it's O(dt) globally,
        // and n=10 steps at dt=0.1 makes it easily visible).
        NT_CHECK(norm(s_euler.r - r_exact) > 1e-4);
        NT_CHECK(norm(s_semi.r - r_exact) > 1e-4);
    }

    // ======================================================================
    // (b) SHM convergence order (research §4.2 — "the most important test
    // you will write"). Fix T_end = one period, halve dt, and check the
    // observed order matches each method's documented global order to
    // within +-0.2, per DESIGN.md.
    // ======================================================================
    {
        const double omega0 = 1.0;
        const double T = 2.0 * constants::pi / omega0;

        const double p_euler = observed_order(Method::Euler, omega0, T, T / 20.0);
        const double p_semi = observed_order(Method::SemiImplicitEuler, omega0, T, T / 20.0);
        const double p_verlet = observed_order(Method::VelocityVerlet, omega0, T, T / 10.0);
        const double p_rk4 = observed_order(Method::RK4, omega0, T, T / 4.0);

        NT_CLOSE(p_euler, 1.0, 0.2);
        NT_CLOSE(p_semi, 1.0, 0.2);
        NT_CLOSE(p_verlet, 2.0, 0.2);
        NT_CLOSE(p_rk4, 4.0, 0.2);
    }

    // ======================================================================
    // (c) Energy behaviour over 100 periods (research §2.0's summary table
    // and §4.3). dt = T/1000 (comfortably inside the "high precision" band
    // of research §3.1) so that the symplectic methods' bounded wobble
    // (~k/2 and ~k^2/4, k = omega*dt) sits well under 1%, isolating the
    // qualitative behaviour DESIGN.md asks for: Euler grows past 10%,
    // semi-implicit/Verlet stay under 1% for the whole run, RK4's drift
    // stays tiny throughout.
    // ======================================================================
    {
        const double omega0 = 1.0;
        const double T = 2.0 * constants::pi / omega0;
        const double dt = T / 1000.0;
        const double T_end = 100.0 * T;
        const double k = omega0 * omega0;  // spring constant, m=1 -> k=omega0^2

        auto run = [&](Method method) {
            Particle p{1.0, Vec3{1.0, 0.0, 0.0}, Vec3{}};
            forces::Spring spring{k, Vec3{}, 0.0};
            SimOptions opt;
            opt.t0 = 0.0;
            opt.t_end = T_end;
            opt.dt = dt;
            opt.method = method;
            return simulate(p, spring, opt);
        };
        auto U = analysis::spring_U(k);

        const Trajectory traj_euler = run(Method::Euler);
        const Trajectory traj_semi = run(Method::SemiImplicitEuler);
        const Trajectory traj_verlet = run(Method::VelocityVerlet);
        const Trajectory traj_rk4 = run(Method::RK4);

        const std::vector<double> drift_euler = analysis::energy_drift(traj_euler, U);
        const std::vector<double> drift_semi = analysis::energy_drift(traj_semi, U);
        const std::vector<double> drift_verlet = analysis::energy_drift(traj_verlet, U);
        const std::vector<double> drift_rk4 = analysis::energy_drift(traj_rk4, U);

        auto max_abs = [](const std::vector<double>& d) {
            double m = 0.0;
            for (double x : d) m = std::fmax(m, std::fabs(x));
            return m;
        };

        // Forward Euler: unconditionally, exponentially divergent — the
        // final relative energy is many multiples of the initial energy.
        NT_CHECK(std::fabs(drift_euler.back()) > 0.10);

        // Semi-implicit Euler and velocity Verlet: bounded oscillation,
        // never trending — the whole trajectory's energy stays within 1%.
        NT_CHECK(max_abs(drift_semi) < 0.01);
        NT_CHECK(max_abs(drift_verlet) < 0.01);

        // RK4: a real but minuscule monotone decay — "tiny" relative to the
        // other three, several orders of magnitude below the 1% band above.
        NT_CHECK(max_abs(drift_rk4) < 1e-6);
    }

    // ======================================================================
    // (d) simulate()'s time bookkeeping: t = t0 + n*dt exactly (to
    // roundoff), never an accumulated t0 + dt + dt + ... (research §7.2).
    // dt = 0.1 is not exactly representable in binary floating point, which
    // is precisely the case that would expose accumulation drift if the
    // driver got this wrong.
    // ======================================================================
    {
        Particle p{1.0, Vec3{}, Vec3{1.0, 0.0, 0.0}};
        forces::ConstantForce zero_force{};  // F = 0: trivial motion, only the time axis matters here
        SimOptions opt;
        opt.t0 = 0.0;
        opt.t_end = 3.7;
        opt.dt = 0.1;  // t_end/dt = 37 exactly in real arithmetic, inexact in binary

        for (Method m : {Method::Euler, Method::SemiImplicitEuler, Method::VelocityVerlet, Method::RK4}) {
            opt.method = m;
            const Trajectory traj = simulate(p, zero_force, opt);
            NT_CLOSE(traj.back().t, opt.t_end, 1e-12);
        }
    }

    return newton_test::report("integrators");
}
