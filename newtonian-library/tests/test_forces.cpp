// tests/test_forces.cpp — every functor in newton::forces returns the
// hand-computed vector for a chosen (t, r, v); sum() adds; quadratic drag at
// v=0 is exactly the zero vector (no NaN from the singularity-free form).
//
// C++ note: expected results are always bound to a named `expected` local
// before being handed to NT_CLOSE_VEC3, rather than written as an inline
// Vec3{...} literal argument. NT_CLOSE_VEC3 is a 3-parameter macro (a, b,
// tol); the preprocessor splits macro arguments on every top-level comma it
// sees and does not understand braces as grouping, so an inline
// `Vec3{1.0, 2.0, 3.0}` literal passed directly as one argument would be
// misparsed as three. A named local sidesteps the issue entirely (see
// NT_CHECK's own variadic fix in test_util.hpp for the other way to solve
// this).
#include <cmath>

#include "test_util.hpp"
#include <newton/newton.hpp>

using namespace newton;  // test .cpp files may use-namespace; only headers may not (DESIGN.md rule 4)

int main() {
    // --- UniformGravity: F = mass * g -----------------------------------
    {
        forces::UniformGravity grav{2.0, Vec3{0.0, -10.0, 0.0}};
        const Vec3 F = grav(0.0, Vec3{}, Vec3{});
        const Vec3 expected{0.0, -20.0, 0.0};
        NT_CLOSE_VEC3(F, expected, 1e-12);

        // Default constants: F = 1.0 * (0, -g_earth, 0).
        forces::UniformGravity default_grav{};
        const Vec3 F2 = default_grav(0.0, Vec3{}, Vec3{});
        const Vec3 expected_default{0.0, -constants::g_earth, 0.0};
        NT_CLOSE_VEC3(F2, expected_default, 1e-12);
    }

    // --- PointGravity: F = -GM*mass*(r-center)/|r-center|^3 --------------
    {
        // GM=1, mass=1, center at origin, r=(2,0,0): |d|=2, F = -(2,0,0)/8 = (-0.25,0,0).
        forces::PointGravity pg{1.0, 1.0, Vec3{}, 0.0};
        const Vec3 F = pg(0.0, Vec3{2.0, 0.0, 0.0}, Vec3{});
        const Vec3 expected{-0.25, 0.0, 0.0};
        NT_CLOSE_VEC3(F, expected, 1e-12);
        NT_CHECK(F.x < 0.0);  // attractive: force points from r toward center

        // Softening caps the force at short range instead of diverging: at
        // r == center, denom = softening^3 > 0 (not the unsoftened 0/0).
        forces::PointGravity soft{1.0, 1.0, Vec3{}, 0.5};
        const Vec3 Fs = soft(0.0, Vec3{}, Vec3{});
        NT_CHECK(Fs == Vec3{});  // d=0 exactly -> numerator is the zero vector regardless of denom
        NT_CHECK(!std::isnan(Fs.x));

        // Unsoftened, exactly at the center: guarded to zero, not NaN/inf.
        forces::PointGravity unsoft{1.0, 1.0, Vec3{}, 0.0};
        const Vec3 F0 = unsoft(0.0, Vec3{}, Vec3{});
        NT_CHECK(F0 == Vec3{});
        NT_CHECK(!std::isnan(F0.x));
    }

    // --- LinearDrag: F = -b v --------------------------------------------
    {
        forces::LinearDrag drag{0.5};
        const Vec3 v{3.0, -4.0, 0.0};
        const Vec3 F = drag(0.0, Vec3{}, v);
        const Vec3 expected{-1.5, 2.0, 0.0};
        NT_CLOSE_VEC3(F, expected, 1e-12);
    }

    // --- QuadraticDrag: F = -c |v| v, and exactly zero (no NaN) at v=0 ----
    {
        forces::QuadraticDrag drag{0.2};
        const Vec3 v{3.0, 4.0, 0.0};  // |v| = 5
        const Vec3 F = drag(0.0, Vec3{}, v);
        // -0.2 * 5 * (3,4,0) = (-3.0, -4.0, 0.0)
        const Vec3 expected{-3.0, -4.0, 0.0};
        NT_CLOSE_VEC3(F, expected, 1e-12);

        const Vec3 F_at_rest = drag(0.0, Vec3{}, Vec3{});
        NT_CHECK(F_at_rest == Vec3{});  // exact zero: no 0/0 from a v_hat-based formula
        NT_CHECK(!std::isnan(F_at_rest.x) && !std::isnan(F_at_rest.y) && !std::isnan(F_at_rest.z));
    }

    // --- sphere-in-air coefficient factories ------------------------------
    {
        const double D = 0.07;  // 7 cm, Taylor's baseball
        NT_CLOSE(forces::linear_drag_sphere_air(D), constants::sphere_air_beta * D, 1e-12);
        NT_CLOSE(forces::quadratic_drag_sphere_air(D), constants::sphere_air_gamma * D * D, 1e-12);
    }

    // --- Spring: F = -k (|d| - L0) d_hat, d = r - anchor -------------------
    {
        // natural_length = 0 reduces exactly to F = -k r (the SHM case).
        forces::Spring spring{4.0, Vec3{}, 0.0};
        const Vec3 F = spring(0.0, Vec3{1.5, 0.0, 0.0}, Vec3{});
        const Vec3 expected{-6.0, 0.0, 0.0};
        NT_CLOSE_VEC3(F, expected, 1e-12);

        // Nonzero anchor and natural_length: stretched past rest length pulls back.
        forces::Spring anchored{2.0, Vec3{1.0, 0.0, 0.0}, 0.5};
        // r = (2.5, 0, 0): d = (1.5, 0, 0), |d| = 1.5, stretch = 1.5-0.5 = 1.0.
        const Vec3 F2 = anchored(0.0, Vec3{2.5, 0.0, 0.0}, Vec3{});
        const Vec3 expected2{-2.0, 0.0, 0.0};
        NT_CLOSE_VEC3(F2, expected2, 1e-12);

        // Zero-guard: sitting exactly on the anchor with natural_length > 0
        // (an unstable equilibrium) must not produce NaN.
        forces::Spring unstable{3.0, Vec3{5.0, 5.0, 0.0}, 1.0};
        const Vec3 F3 = unstable(0.0, Vec3{5.0, 5.0, 0.0}, Vec3{});
        NT_CHECK(F3 == Vec3{});
        NT_CHECK(!std::isnan(F3.x));
    }

    // --- AnisotropicSpring: F = (-kx x, -ky y, -kz z) ----------------------
    {
        forces::AnisotropicSpring aspring{1.0, 2.0, 3.0};
        const Vec3 F = aspring(0.0, Vec3{2.0, 3.0, 4.0}, Vec3{});
        const Vec3 expected{-2.0, -6.0, -12.0};
        NT_CLOSE_VEC3(F, expected, 1e-12);

        // Equal constants collapse to the isotropic case.
        forces::AnisotropicSpring iso{5.0, 5.0, 5.0};
        forces::Spring iso_ref{5.0, Vec3{}, 0.0};
        const Vec3 r{1.0, -2.0, 0.5};
        const Vec3 iso_F = iso(0.0, r, Vec3{});
        const Vec3 iso_ref_F = iso_ref(0.0, r, Vec3{});
        NT_CLOSE_VEC3(iso_F, iso_ref_F, 1e-12);
    }

    // --- SinusoidalDrive: F(t) = F0 cos(omega t + phase) dir ---------------
    {
        forces::SinusoidalDrive drive{2.0, 3.0, 0.0, Vec3{0.0, 1.0, 0.0}};
        const double t = 0.4;
        const double expected_mag = 2.0 * std::cos(3.0 * t);
        const Vec3 F = drive(t, Vec3{}, Vec3{});
        const Vec3 expected{0.0, expected_mag, 0.0};
        NT_CLOSE_VEC3(F, expected, 1e-12);

        // At t=0 with phase=0, F = F0 * dir exactly.
        const Vec3 F_at_zero = drive(0.0, Vec3{}, Vec3{});
        const Vec3 expected_at_zero{0.0, 2.0, 0.0};
        NT_CLOSE_VEC3(F_at_zero, expected_at_zero, 1e-12);
    }

    // --- ConstantForce: F = F, independent of t, r, v -----------------------
    {
        forces::ConstantForce cf{Vec3{1.0, -2.0, 3.0}};
        const Vec3 F = cf(5.0, Vec3{9.0, 9.0, 9.0}, Vec3{1.0, 1.0, 1.0});
        NT_CLOSE_VEC3(F, cf.F, 1e-12);
    }

    // --- MagneticForce: F = q v x B, always perpendicular to v -------------
    {
        forces::MagneticForce mag{1.0, Vec3{0.0, 0.0, 1.0}};
        const Vec3 v{1.0, 0.0, 0.0};
        const Vec3 F = mag(0.0, Vec3{}, v);
        // x_hat x z_hat = -y_hat.
        const Vec3 expected{0.0, -1.0, 0.0};
        NT_CLOSE_VEC3(F, expected, 1e-12);
        NT_CLOSE(dot(F, v), 0.0, 1e-12);  // does no work: F . v == 0 exactly
    }

    // --- LorentzForce: F = q (E + v x B) ------------------------------------
    {
        forces::LorentzForce lorentz{2.0, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
        const Vec3 v{0.0, 1.0, 0.0};
        // v x B = (0,1,0)x(0,0,1) = (1,0,0). F = 2*((1,0,0)+(1,0,0)) = (4,0,0).
        const Vec3 F = lorentz(0.0, Vec3{}, v);
        const Vec3 expected{4.0, 0.0, 0.0};
        NT_CLOSE_VEC3(F, expected, 1e-12);
    }

    // --- sum(): superposition of forces -------------------------------------
    {
        forces::UniformGravity grav{1.0, Vec3{0.0, -9.8, 0.0}};
        forces::LinearDrag drag{0.5};
        auto combined = forces::sum(grav, drag);
        const Vec3 v{2.0, 0.0, 0.0};
        const Vec3 expected = grav(0.0, Vec3{}, v) + drag(0.0, Vec3{}, v);
        const Vec3 got = combined(0.0, Vec3{}, v);
        NT_CLOSE_VEC3(got, expected, 1e-12);

        // Three forces at once: gravity + drag + a constant thrust.
        forces::ConstantForce thrust{Vec3{0.0, 5.0, 0.0}};
        auto three = forces::sum(grav, drag, thrust);
        const Vec3 expected3 = grav(0.0, Vec3{}, v) + drag(0.0, Vec3{}, v) + thrust(0.0, Vec3{}, v);
        const Vec3 got3 = three(0.0, Vec3{}, v);
        NT_CLOSE_VEC3(got3, expected3, 1e-12);

        // Zero forces: the fold's explicit Vec3{} seed makes this well-formed,
        // returning the zero force (compose.hpp's C++ note).
        auto none = forces::sum();
        NT_CHECK(none(0.0, Vec3{}, Vec3{}) == Vec3{});
    }

    return newton_test::report("forces");
}
