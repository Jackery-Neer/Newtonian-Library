// tests/test_coords.cpp — cylindrical and spherical polar coordinates:
// hand-computed points, the atan2 quadrant fix (Taylor problem 1.42), round
// trips through Cartesian in both directions, agreement between the direct
// cylindrical <-> spherical conversions and the route through Cartesian, and
// the reason theta uses atan2 instead of acos.
#include "test_util.hpp"        // quotes: sits next to this file
#include <newton/newton.hpp>    // angle brackets: reached via -I include

#include <cmath>

using newton::Vec3;
using newton::constants::pi;
using newton::coords::CylindricalPolar;
using newton::coords::SphericalPolar;
using newton::coords::to_cartesian;
using newton::coords::to_cylindrical;
using newton::coords::to_spherical;

int main() {
    const double tol = 1e-12;

    // --- hand-computed point: cylindrical --------------------------------
    // (3, 4, 5): rho is the 3-4-5 hypotenuse; z passes straight through.
    const CylindricalPolar c345 = to_cylindrical(Vec3{3.0, 4.0, 5.0});
    NT_CLOSE(c345.rho, 5.0, tol);
    NT_CLOSE(c345.phi, std::atan2(4.0, 3.0), tol);
    NT_CHECK(c345.z == 5.0);  // exact ==: z is copied, never computed

    // --- the quadrant bug (Taylor problem 1.42) ---------------------------
    // (-1, -1) is in the third quadrant, at phi = -3pi/4. The naive formula
    // atan(y / x) only sees the ratio (-1)/(-1) = +1 and answers +pi/4, the
    // exact opposite direction. atan2 sees both signs.
    const double phi_q3 = to_cylindrical(Vec3{-1.0, -1.0, 0.0}).phi;
    NT_CLOSE(phi_q3, -3.0 * pi / 4.0, tol);
    NT_CLOSE(std::atan(-1.0 / -1.0) - phi_q3, pi, tol);  // naive answer is off by pi

    // --- hand-computed points: spherical, physics convention -------------
    // theta is measured DOWN from +z, phi AROUND z from +x.
    struct Case {
        Vec3 v;
        double r, theta, phi;
    };
    const Case cases[] = {
        {{0.0, 0.0, 2.0}, 2.0, 0.0, 0.0},                       // north pole
        {{0.0, 0.0, -2.0}, 2.0, pi, 0.0},                       // south pole
        {{3.0, 0.0, 0.0}, 3.0, pi / 2.0, 0.0},                  // equator, +x
        {{0.0, 3.0, 0.0}, 3.0, pi / 2.0, pi / 2.0},             // equator, +y
        {{1.0, 1.0, std::sqrt(2.0)}, 2.0, pi / 4.0, pi / 4.0},  // 45 deg down, 45 deg around
    };
    for (const Case& k : cases) {
        const SphericalPolar s = to_spherical(k.v);
        NT_CLOSE(s.r, k.r, tol);
        NT_CLOSE(s.theta, k.theta, tol);
        NT_CLOSE(s.phi, k.phi, tol);
        // Bind the input first: a brace list inside a 3-argument macro call
        // would be split at its commas (see test_util.hpp).
        const SphericalPolar given{k.r, k.theta, k.phi};
        NT_CLOSE_VEC3(to_cartesian(given), k.v, tol);
    }

    // --- round trip: Cartesian -> curvilinear -> Cartesian ----------------
    // Every octant, both poles, the equator, the origin, and sizes from 1e-3
    // to 1e3. This direction works even on the z axis, where phi is
    // undefined, because rho = 0 (or sin(theta) = 0) multiplies phi away.
    const Vec3 points[] = {
        {1.0, 2.0, 3.0},     {-1.0, 2.0, 3.0},    {-1.0, -2.0, 3.0},       {1.0, -2.0, 3.0},
        {1.0, 2.0, -3.0},    {-1.0, 2.0, -3.0},   {-1.0, -2.0, -3.0},      {1.0, -2.0, -3.0},
        {0.0, 0.0, 5.0},     {0.0, 0.0, -5.0},    {4.0, 0.0, 0.0},         {0.0, -4.0, 0.0},
        {-4.0, 0.0, 0.0},    {2e-3, -1e-3, 5e-4}, {700.0, -300.0, 900.0}, {0.0, 0.0, 0.0},
    };
    for (const Vec3& v : points) {
        NT_CLOSE_VEC3(to_cartesian(to_cylindrical(v)), v, tol);
        NT_CLOSE_VEC3(to_cartesian(to_spherical(v)), v, tol);

        // Every angle lands in its documented range.
        const SphericalPolar s = to_spherical(v);
        const CylindricalPolar c = to_cylindrical(v);
        NT_CHECK(s.r >= 0.0 && c.rho >= 0.0);
        NT_CHECK(s.theta >= 0.0 && s.theta <= pi);
        NT_CHECK(s.phi >= -pi && s.phi <= pi);

        // The direct cylindrical <-> spherical conversions must agree with
        // the long way round through Cartesian.
        const SphericalPolar s_direct = to_spherical(c);
        NT_CLOSE(s_direct.r, s.r, tol);
        NT_CLOSE(s_direct.theta, s.theta, tol);
        NT_CLOSE(s_direct.phi, s.phi, tol);
        const CylindricalPolar c_direct = to_cylindrical(s);
        NT_CLOSE(c_direct.rho, c.rho, tol);
        NT_CLOSE(c_direct.phi, c.phi, tol);
        NT_CLOSE(c_direct.z, c.z, tol);
    }

    // --- round trip: curvilinear -> Cartesian -> curvilinear --------------
    // Only meaningful where the angles are defined: theta stays strictly
    // inside (0, pi), off the z axis, and phi strictly inside (-pi, pi), so
    // the seam where -pi and +pi name the same direction never comes up.
    const double radii[] = {0.5, 1.0, 1000.0};
    for (const double r : radii) {
        for (int i = 0; i <= 10; ++i) {
            const double theta = 0.1 + 0.3 * i;  // 0.1 .. 3.1
            for (int j = 0; j <= 12; ++j) {
                const double phi = -3.0 + 0.5 * j;  // -3.0 .. 3.0
                const SphericalPolar s = to_spherical(to_cartesian(SphericalPolar{r, theta, phi}));
                NT_CLOSE(s.r, r, tol);
                NT_CLOSE(s.theta, theta, tol);
                NT_CLOSE(s.phi, phi, tol);

                const double z = (theta - 1.5) * r;  // any z works; vary it too
                const CylindricalPolar c = to_cylindrical(to_cartesian(CylindricalPolar{r, phi, z}));
                NT_CLOSE(c.rho, r, tol);
                NT_CLOSE(c.phi, phi, tol);
                NT_CLOSE(c.z, z, tol);
            }
        }
    }

    // --- why theta uses atan2(rho, z) and not acos(z / r) -----------------
    // A point 1e-10 rad off the +z axis: a 0.1 nm sideways offset at 1 m.
    // z/r rounds to exactly 1.0, so acos snaps the point onto the axis and
    // reports theta = 0. atan2 keeps the angle to full precision.
    const Vec3 near_pole{1e-10, 0.0, 1.0};
    NT_CHECK(std::acos(near_pole.z / newton::norm(near_pole)) == 0.0);  // acos: angle lost
    NT_CLOSE(to_spherical(near_pole).theta / 1e-10, 1.0, tol);         // atan2: angle kept

    // --- the origin: no division anywhere, so no NaN ----------------------
    const SphericalPolar s0 = to_spherical(Vec3{});
    NT_CHECK(s0.r == 0.0 && s0.theta == 0.0 && s0.phi == 0.0);
    const CylindricalPolar c0 = to_cylindrical(Vec3{});
    NT_CHECK(c0.rho == 0.0 && c0.phi == 0.0 && c0.z == 0.0);

    return newton_test::report("coords");
}
