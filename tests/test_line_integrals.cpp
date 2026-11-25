// tests/test_line_integrals.cpp — numerical line integrals and circulation
// (fields/line_integral.hpp). Every expected value is a closed form from
// physics: work = F·Δr for a uniform field, ½k(|a|²-|b|²) for a spring,
// GMm(1/|b| - 1/|a|) for point gravity, and Stokes' theorem for the swirl
// ((-y, x, 0) has curl 2 ẑ, so its circulation is 2 × the enclosed area) and
// for the "hole" field (-y, x, 0)/(x²+y²), whose circulation around the z
// axis is exactly 2π although its curl vanishes everywhere off the axis.
#include "test_util.hpp"        // quotes: sits next to this file
#include <newton/newton.hpp>    // angle brackets: reached via -I include
#include <newton/fields/line_integral.hpp>

#include <cmath>
#include <limits>
#include <vector>

using newton::Vec3;
using newton::constants::pi;
using newton::fields::circulation;
using newton::fields::circulation_detailed;
using newton::fields::line_integral;
using newton::fields::line_integral_detailed;
using newton::fields::LineIntegral;

namespace {

// Relative error of a scalar: |a - b| / max(|b|, tiny).
double rel_err(double a, double b) {
    return std::fabs(a - b) / std::fmax(std::fabs(b), 1e-300);
}

// Evaluation budget: 15 + 30 * (200 - 1) (see detail::kMaxPieces).
constexpr int max_evals = 5985;

// Fields as plain callables of position only (freezing t and v of a Force).
Vec3 swirl(const Vec3& r) { return {-r.y, r.x, 0.0}; }
Vec3 hole(const Vec3& r) {
    const double d = r.x * r.x + r.y * r.y;
    return {-r.y / d, r.x / d, 0.0};
}

}  // namespace

int main() {
    const double tol = 1e-11;

    // --- uniform field: W = F·(b - a), any path ---------------------------
    {
        const newton::forces::UniformGravity g{2.0, {0.0, -9.81, 0.0}};
        const auto F = [&](const Vec3& r) { return g(0.0, r, Vec3{}); };
        const Vec3 a{1.0, 2.0, 3.0}, b{-4.0, 7.5, 0.5};
        const LineIntegral r = line_integral_detailed(F, a, b);
        NT_CHECK(r.converged);
        NT_CLOSE(r.value, dot(g.mass * g.g, b - a), 1e-13);
        NT_CHECK(r.evaluations == 15);  // a constant field needs just one panel
        // Dropping 5 m: W = mgh = 2 * 9.81 * 5.
        NT_CLOSE(line_integral(F, Vec3{0.0, 5.0, 0.0}, Vec3{}), 2.0 * 9.81 * 5.0, tol);
    }

    // --- spring: W = ½k(|a|² - |b|²) --------------------------------------
    {
        const newton::forces::Spring sp{3.5, {}, 0.0};
        const auto F = [&](const Vec3& r) { return sp(0.0, r, Vec3{}); };
        const Vec3 a{1.0, -2.0, 0.5}, b{-0.3, 0.7, 2.0};
        const double exact = 0.5 * 3.5 * (norm2(a) - norm2(b));
        const LineIntegral r = line_integral_detailed(F, a, b);
        NT_CHECK(r.converged);
        NT_CLOSE(r.value, exact, tol);
        NT_CHECK(r.error <= std::fmax(1e-12 * std::fabs(r.value), 1e-14 * r.scale));
        NT_CLOSE(line_integral(F, a, b), exact, tol);
        // Reversing the path flips the sign of the work.
        NT_CLOSE(line_integral(F, b, a), -exact, tol);
    }

    // --- point gravity, lab-ish scale: W = GMm(1/|b| - 1/|a|) -------------
    {
        const newton::forces::PointGravity pg{4.0, 1.5, {}, 0.0};
        const auto F = [&](const Vec3& r) { return pg(0.0, r, Vec3{}); };
        const Vec3 a{3.0, 0.0, 4.0}, b{-1.0, 2.0, 1.0};  // |a| = 5, |b| = sqrt(6)
        const double exact = 4.0 * 1.5 * (1.0 / norm(b) - 1.0 / norm(a));
        const LineIntegral r = line_integral_detailed(F, a, b);
        NT_CHECK(r.converged);
        NT_CLOSE(r.value, exact, tol);
        NT_CHECK(r.evaluations <= max_evals);
        NT_CHECK(r.scale >= std::fabs(r.value));  // S = ∫|F||dr| can never be below |∫F·dr|
    }

    // --- Earth-orbit scale: LEO to far away, GM = 3.986004418e14 ----------
    {
        const double GM = 3.986004418e14, m = 1000.0;
        const newton::forces::PointGravity pg{GM, m, {}, 0.0};
        const auto F = [&](const Vec3& r) { return pg(0.0, r, Vec3{}); };
        const Vec3 a{6.771e6, 0.0, 0.0};                // ~400 km altitude
        const Vec3 b{-2.0e7, 3.5e7, 1.0e7};
        const double exact = GM * m * (1.0 / norm(b) - 1.0 / norm(a));
        const LineIntegral r = line_integral_detailed(F, a, b);
        NT_CHECK(r.converged);
        NT_CLOSE(r.value, exact, tol);
        NT_CHECK(r.evaluations <= max_evals);
        // Work from LEO out to the Moon's distance, along a radial line.
        const Vec3 moon{3.844e8, 0.0, 0.0};
        NT_CLOSE(line_integral(F, a, moon), GM * m * (1.0 / 3.844e8 - 1.0 / 6.771e6), tol);
    }

    // --- polylines: path independence for conservative F ------------------
    {
        const newton::forces::PointGravity pg{1.0, 1.0, {}, 0.0};
        const auto F = [&](const Vec3& r) { return pg(0.0, r, Vec3{}); };
        const Vec3 a{2.0, 0.0, 0.0}, b{0.0, 3.0, 1.0};
        const double direct = line_integral(F, a, b);
        const double detour = line_integral(F, std::vector<Vec3>{a, {2.0, 3.0, 0.0}, {2.0, 3.0, 4.0}, b});
        NT_CLOSE(direct, 1.0 / norm(b) - 1.0 / norm(a), tol);
        NT_CLOSE(detour, direct, tol);
    }

    // --- polylines: the swirl is path dependent, with a known difference --
    // From (1,0,0) to (0,1,0): the direct segment gives 1; the corner route
    // through (1,1,0) gives 2; the route through the origin gives 0. The
    // difference of two routes is the circulation around the closed loop
    // they make (Stokes: curl 2 × area): direct vs corner encloses a
    // triangle of area 1/2 -> 1.
    {
        const Vec3 a{1.0, 0.0, 0.0}, b{0.0, 1.0, 0.0};
        const double direct = line_integral(swirl, a, b);
        const double corner = line_integral(swirl, std::vector<Vec3>{a, {1.0, 1.0, 0.0}, b});
        const double via_origin = line_integral(swirl, std::vector<Vec3>{a, {}, b});
        NT_CLOSE(direct, 1.0, tol);
        NT_CLOSE(corner, 2.0, tol);
        NT_CLOSE(via_origin, 0.0, 1e-13);
        NT_CLOSE(corner - direct, 1.0, tol);
        // A closed polyline (the unit square, counterclockwise) is a circulation: 2 × area.
        const double square = line_integral(
            swirl, std::vector<Vec3>{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0},
                                     {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}});
        NT_CLOSE(square, 2.0, tol);
        // Degenerate polylines: no segments, no work.
        NT_CHECK(line_integral(swirl, std::vector<Vec3>{}) == 0.0);
        NT_CHECK(line_integral(swirl, std::vector<Vec3>{a}) == 0.0);
    }

    // --- circulation of the swirl: 2πR² -----------------------------------
    {
        for (const double R : {0.5, 1.0, 3.0}) {
            const LineIntegral c = circulation_detailed(swirl, Vec3{}, Vec3{0.0, 0.0, 1.0}, R);
            NT_CHECK(c.converged);
            NT_CLOSE(c.value, 2.0 * pi * R * R, tol);
            NT_CHECK(c.evaluations <= max_evals);
        }
        // Off-center loop: circulation = 2 × area, regardless of where the loop sits.
        NT_CLOSE(circulation(swirl, Vec3{5.0, -2.0, 1.0}, Vec3{0.0, 0.0, 1.0}, 2.0), 2.0 * pi * 4.0, tol);
        // A loop in the xz-plane (normal ŷ): the swirl's curl (0,0,2) has no ŷ
        // component, so Stokes gives 0.
        NT_CLOSE(circulation(swirl, Vec3{}, Vec3{0.0, 1.0, 0.0}, 1.0), 0.0, 1e-13);
        // A circle tilted 60 degrees from the z axis: Stokes gives (curl·n̂) × area
        // = 2 cos(60°) × π R² = π R² for R = 1.
        const Vec3 tilt{std::sin(pi / 3.0), 0.0, std::cos(pi / 3.0)};
        NT_CLOSE(circulation(swirl, Vec3{}, tilt, 1.0), pi, tol);
    }

    // --- the hole field: 2π around the axis, 0 beside it ------------------
    {
        const LineIntegral around = circulation_detailed(hole, Vec3{}, Vec3{0.0, 0.0, 1.0}, 1.0);
        NT_CHECK(around.converged);
        NT_CLOSE(around.value, 2.0 * pi, tol);
        // The answer does not depend on the radius or on being centered on the axis
        // as long as the loop encloses it (Ampère's law for a wire).
        NT_CLOSE(circulation(hole, Vec3{}, Vec3{0.0, 0.0, 1.0}, 7.0), 2.0 * pi, tol);
        NT_CLOSE(circulation(hole, Vec3{0.3, 0.2, 5.0}, Vec3{0.0, 0.0, 1.0}, 1.0), 2.0 * pi, tol);
        // Beside the axis: curl is zero inside the loop, so the circulation is 0.
        const LineIntegral beside = circulation_detailed(hole, Vec3{3.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}, 1.0);
        NT_CHECK(beside.converged);
        NT_CLOSE(beside.value, 0.0, 1e-13);
        NT_CHECK(std::fabs(beside.value) <= 1e-12 * beside.scale);
        // A tilted circle through which the axis passes also gives +2π (its normal
        // has positive z component); flipping the normal gives -2π.
        const Vec3 tilt{1.0, 0.0, 1.0};
        NT_CLOSE(circulation(hole, Vec3{}, tilt, 1.0), 2.0 * pi, tol);
        NT_CLOSE(circulation(hole, Vec3{}, -1.0 * tilt, 1.0), -2.0 * pi, tol);
    }

    // --- conservative fields: circulation is zero to scale_tol ------------
    {
        const newton::forces::PointGravity pg{1.0, 1.0, {}, 0.0};
        const auto Fg = [&](const Vec3& r) { return pg(0.0, r, Vec3{}); };
        const newton::forces::Spring sp{2.0, {0.5, 0.0, 0.0}, 0.0};
        const auto Fs = [&](const Vec3& r) { return sp(0.0, r, Vec3{}); };
        // Loops that encircle the origin are fine for point gravity: it is curl-free
        // AND its singularity is a point, so no loop is linked with it.
        const Vec3 n{0.2, -0.5, 1.0};
        for (const double R : {0.5, 2.0}) {
            const LineIntegral cg = circulation_detailed(Fg, Vec3{0.1, 0.0, 0.0}, n, R);
            NT_CHECK(cg.converged);
            NT_CHECK(std::fabs(cg.value) <= 1e-12 * cg.scale);
            const LineIntegral cs = circulation_detailed(Fs, Vec3{1.0, 1.0, 1.0}, n, R);
            NT_CHECK(cs.converged);
            NT_CHECK(std::fabs(cs.value) <= 1e-12 * cs.scale);
        }
    }

    // --- orientation follows the right-hand rule --------------------------
    {
        const Vec3 up{0.0, 0.0, 1.0}, down{0.0, 0.0, -1.0};
        const double ccw = circulation(swirl, Vec3{}, up, 1.0);
        const double cw = circulation(swirl, Vec3{}, down, 1.0);
        NT_CLOSE(ccw, 2.0 * pi, tol);
        NT_CLOSE(cw, -2.0 * pi, tol);
        NT_CLOSE(ccw + cw, 0.0, 1e-12);
        // Only the direction of the normal matters, not its length.
        NT_CLOSE(circulation(swirl, Vec3{}, Vec3{0.0, 0.0, 40.0}, 1.0), ccw, 1e-13);
        // Rotation about +x is the field (0, -z, y); right-hand rule about +x is +.
        const auto swirl_x = [](const Vec3& r) { return Vec3{0.0, -r.z, r.y}; };
        NT_CLOSE(circulation(swirl_x, Vec3{}, Vec3{1.0, 0.0, 0.0}, 1.0), 2.0 * pi, tol);
        NT_CLOSE(circulation(swirl_x, Vec3{}, Vec3{-1.0, 0.0, 0.0}, 1.0), -2.0 * pi, tol);
        // Axis-aligned and general normals: the frame construction never degenerates.
        // F = n × r is a rotation about n with curl 2n, so around a circle with unit
        // normal n̂ = n/|n| and R = 1, Stokes gives (2n·n̂) π = 2π|n|.
        for (const Vec3& n : {Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0},
                              Vec3{1.0, 1.0, 1.0}, Vec3{-0.3, 0.2, -0.9}}) {
            const auto rot_n = [&](const Vec3& r) { return cross(n, r); };
            NT_CLOSE(circulation(rot_n, Vec3{}, n, 1.0), 2.0 * pi * norm(n), tol);
        }
    }

    // --- failure flags ----------------------------------------------------
    {
        const newton::forces::PointGravity pg{1.0, 1.0, {}, 0.0};
        const auto F = [&](const Vec3& r) { return pg(0.0, r, Vec3{}); };
        const double nan_ = std::numeric_limits<double>::quiet_NaN();
        // A path straight through the center: the integral is divergent.
        const LineIntegral through = line_integral_detailed(F, Vec3{-1.0, 0.0, 0.0}, Vec3{2.0, 0.0, 0.0});
        NT_CHECK(!through.converged);
        NT_CHECK(through.evaluations <= max_evals);
        // The minimum-piece-width guard (1e-12 of the range) ends the hunt after about
        // 40 halvings at the singular point: 15 + 30 x 40 = 1215 evaluations, far
        // fewer than the 5985 the piece limit alone would allow.
        NT_CHECK(through.evaluations < 2000);
        NT_CHECK(std::isnan(line_integral(F, Vec3{-1.0, 0.0, 0.0}, Vec3{2.0, 0.0, 0.0})));
        NT_CHECK(std::isnan(line_integral(F, Vec3{0.3, 0.0, 0.0}, Vec3{-0.7, 0.0, 0.0})));
        // The SYMMETRIC path, center hit at the exact midpoint: F·d is odd about
        // the center and cancels in K15 and G7 (work error exactly 0); only the
        // scale integral ∫|F||dr| (even, so it cannot cancel) exposes the
        // divergence. PointGravity returns {} exactly at its center.
        const LineIntegral sym = line_integral_detailed(F, Vec3{-1.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
        NT_CHECK(!sym.converged);
        NT_CHECK(sym.evaluations <= max_evals);
        NT_CHECK(std::isnan(line_integral(F, Vec3{-1.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0})));
        // Same through a circle's diameter-like symmetric polyline.
        NT_CHECK(std::isnan(line_integral(F, std::vector<Vec3>{{-2.0, 0.0, 0.0}, {2.0, 0.0, 0.0}})));
        // And a symmetric path through a NaN-at-center 1/r² field still fails.
        const auto inv_sq = [&](const Vec3& r) { return (-1.0 / (norm2(r) * norm(r))) * r; };  // NaN at 0
        NT_CHECK(!line_integral_detailed(inv_sq, Vec3{-1.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}).converged);
        NT_CHECK(std::isnan(line_integral(inv_sq, Vec3{-1.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0})));
        // Ending exactly on the center (F infinite there) is also divergent.
        NT_CHECK(std::isnan(line_integral(F, Vec3{1.0, 0.0, 0.0}, Vec3{})));
        // A polyline with one bad segment is NaN overall.
        NT_CHECK(std::isnan(line_integral(F, std::vector<Vec3>{{1.0, 1.0, 0.0}, {2.0, 0.0, 0.0},
                                                                {-3.0, 0.0, 0.0}})));
        // A circle through the center of the force is likewise a failure.
        NT_CHECK(!circulation_detailed(F, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}, 1.0).converged);
        NT_CHECK(std::isnan(circulation(F, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0}, 1.0)));
        // Non-finite field values: converged = false, value NaN.
        const auto bad = [&](const Vec3& r) { return r.x > 0.3 ? Vec3{nan_, 0.0, 0.0} : Vec3{1.0, 0.0, 0.0}; };
        const LineIntegral nf = line_integral_detailed(bad, Vec3{}, Vec3{1.0, 0.0, 0.0});
        NT_CHECK(!nf.converged);
        NT_CHECK(std::isnan(nf.value));
        NT_CHECK(std::isnan(line_integral(bad, Vec3{}, Vec3{1.0, 0.0, 0.0})));
        // Bad circles: zero normal, radius <= 0, NaN radius.
        const LineIntegral zn = circulation_detailed(swirl, Vec3{}, Vec3{}, 1.0);
        NT_CHECK(!zn.converged);
        NT_CHECK(zn.evaluations == 0);
        NT_CHECK(!circulation_detailed(swirl, Vec3{}, Vec3{0.0, 0.0, 1.0}, 0.0).converged);
        NT_CHECK(!circulation_detailed(swirl, Vec3{}, Vec3{0.0, 0.0, 1.0}, -1.0).converged);
        NT_CHECK(!circulation_detailed(swirl, Vec3{}, Vec3{0.0, 0.0, 1.0}, nan_).converged);
        NT_CHECK(std::isnan(circulation(swirl, Vec3{}, Vec3{}, 1.0)));
        NT_CHECK(std::isnan(circulation(swirl, Vec3{}, Vec3{0.0, 0.0, 1.0}, -2.0)));
    }

    // --- circle_frame normalizes its argument -----------------------------
    {
        const newton::fields::detail::CircleFrame fr =
            newton::fields::detail::circle_frame(Vec3{0.0, 0.0, 5.0});
        NT_CLOSE(norm(fr.u), 1.0, 1e-14);
        NT_CLOSE(norm(fr.w), 1.0, 1e-14);
        NT_CLOSE(dot(fr.u, fr.w), 0.0, 1e-14);
        NT_CLOSE_VEC3(cross(fr.u, fr.w), (Vec3{0.0, 0.0, 1.0}), 1e-14);
    }

    // --- bookkeeping ------------------------------------------------------
    {
        // Evaluation counts are 15 + 30 × (pieces - 1) and never exceed 5985,
        // even for a hopelessly singular integrand.
        const newton::forces::PointGravity pg{1.0, 1.0, {}, 0.0};
        const auto F = [&](const Vec3& r) { return pg(0.0, r, Vec3{}); };
        const LineIntegral hard = line_integral_detailed(F, Vec3{-1.0, 1e-9, 0.0}, Vec3{1.0, 1e-9, 0.0});
        NT_CHECK(hard.evaluations <= max_evals);
        NT_CHECK(hard.evaluations >= 15);
        NT_CHECK((hard.evaluations - 15) % 30 == 0);
        // A loose tolerance is no more expensive than a tight one.
        const LineIntegral loose = line_integral_detailed(F, Vec3{0.3, 0.0, 0.0}, Vec3{3.0, 4.0, 0.0}, 1e-4);
        const LineIntegral tight = line_integral_detailed(F, Vec3{0.3, 0.0, 0.0}, Vec3{3.0, 4.0, 0.0}, 1e-13);
        NT_CHECK(loose.converged && tight.converged);
        NT_CHECK(loose.evaluations <= tight.evaluations);
        // The scale S = ∮|F||dr| for the unit-circle swirl is |F| = 1 times length 2π.
        const LineIntegral c = circulation_detailed(swirl, Vec3{}, Vec3{0.0, 0.0, 1.0}, 1.0);
        NT_CHECK(c.evaluations <= max_evals);
        NT_CHECK(rel_err(c.scale, 2.0 * pi) < 1e-12);
    }

    return newton_test::report("line_integrals");
}
