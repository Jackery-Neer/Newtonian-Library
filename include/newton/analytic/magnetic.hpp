#pragma once
// newton/analytic/magnetic.hpp — charged particle in a uniform magnetic
// field (Taylor §2.5, §2.7).
//
// Physics: the magnetic force F = q v x B does no work (it's always
// perpendicular to v), so |v| is exactly conserved and the motion is a
// helix: uniform circular motion in the plane perpendicular to B, uniform
// drift along B. That combination of "exactly solvable" and "one of the
// clearest invariants in the book" makes this the standard oracle for
// testing that an integrator preserves |v| (research map §8, row 10) —
// see forces/electromagnetic.hpp for the numeric side.
//
// C++ note: cyclotron_omega() is pure arithmetic and constexpr; the other
// two call std::abs(double), which is not guaranteed constexpr until C++23,
// so they are plain inline.

#include <cmath>

#include <newton/core/constants.hpp>

namespace newton::analytic {

/// Cyclotron (angular) frequency ω = q B / m [rad/s] (Taylor §2.5). Signed:
/// its sign encodes the sense of rotation (which way the particle circles
/// for a given sign of q and B), while cyclotron_radius()/cyclotron_period()
/// below report unsigned geometric quantities.
/// @param q charge [C]
/// @param B magnetic field magnitude (component along the rotation axis) [T]
/// @param m mass [kg]
[[nodiscard]] constexpr double cyclotron_omega(double q, double B, double m) { return q * B / m; }

/// Gyroradius (radius of the circular motion transverse to B):
/// r = m v_perp / (|q| B) [m] (Taylor §2.5-2.7) — equivalently v_perp/|ω|.
/// @param m      mass [kg]
/// @param v_perp speed component perpendicular to B [m/s]
/// @param q      charge [C]
/// @param B      magnetic field magnitude [T]
[[nodiscard]] inline double cyclotron_radius(double m, double v_perp, double q, double B) {
    return m * v_perp / std::abs(q * B);
}

/// Cyclotron period T = 2π m / (|q| B) [s] (Taylor §2.5-2.7) — independent
/// of speed, the property that makes cyclotrons work: every particle takes
/// the same time per revolution regardless of its (transverse) energy.
/// @param m mass [kg]
/// @param q charge [C]
/// @param B magnetic field magnitude [T]
[[nodiscard]] inline double cyclotron_period(double m, double q, double B) {
    return 2.0 * constants::pi * m / std::abs(q * B);
}

}  // namespace newton::analytic
