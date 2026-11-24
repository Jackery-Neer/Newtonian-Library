#pragma once
// newton/analysis/oscillation.hpp — measurements on an oscillating trajectory.
//
// Physics: everything here measures the x-component of position, r.x — the
// convention this library uses throughout Ch. 5 (Taylor's 1D oscillator
// lives along x; see analytic/oscillator.hpp). Feeding a 3D trajectory
// through these functions is harmless — they simply read r.x and ignore
// r.y, r.z — but they exist to turn a simulate()'d Trajectory into the same
// numbers analytic::SHM / DampedOscillator / DrivenOscillator predict in
// closed form, for a numeric-vs-analytic comparison.
//
// C++ note: all three functions take a plain `const Trajectory&` rather than
// a generic component-accessor template. That keeps their signatures the
// same shape as the rest of the analysis module (energy_drift,
// angular_momentum) and avoids a lambda-template parameter for a feature
// (measuring y or z instead of x) nothing in this library currently needs —
// add it if a future example demands it, per DESIGN.md's "correctness and
// readability beat generality."

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <newton/core/trajectory.hpp>

namespace newton::analysis {

/// Times and r.x values of the local maxima found in a Trajectory, in time
/// order; peaks.t[i] and peaks.x[i] describe the same peak. Used for
/// log-decrement extraction from a damped oscillator's envelope
/// (Taylor §5.4: successive peak heights shrink by a constant ratio).
struct PeakList {
    std::vector<double> t;  ///< times of the peaks [s]
    std::vector<double> x;  ///< r.x at each peak [m]
};

/// Local maxima of r.x by 3-point comparison: sample i is a peak iff
/// x[i-1] < x[i] > x[i+1] (strict on both sides). Simple and robust for the
/// smooth, densely-sampled trajectories this library produces; it can miss a
/// peak that lands exactly on a two-sample plateau or falls between grid
/// points, which is why measure_period() below uses interpolated zero
/// crossings rather than peak spacing for anything precision-sensitive.
[[nodiscard]] inline PeakList find_peaks(const Trajectory& traj) {
    PeakList peaks;
    for (std::size_t i = 1; i + 1 < traj.size(); ++i) {
        const double x_prev = traj.states[i - 1].r.x;
        const double x_cur  = traj.states[i].r.x;
        const double x_next = traj.states[i + 1].r.x;
        if (x_cur > x_prev && x_cur > x_next) {
            peaks.t.push_back(traj.states[i].t);
            peaks.x.push_back(x_cur);
        }
    }
    return peaks;
}

/// Period estimated from the mean spacing of upward zero crossings of r.x
/// (a crossing is a sample pair where x goes from negative to non-negative).
/// Each crossing time is linearly interpolated within its bracketing pair,
/// which is far more accurate than snapping to the nearest sample time at
/// typical simulation dt. Returns 0.0 if fewer than two crossings are found
/// (not enough data to measure a spacing at all).
///
/// Physics: for x(t) oscillating about 0 — SHM and its damped/driven
/// relatives (Taylor ch. 5) — consecutive upward crossings are exactly one
/// period apart, so averaging the spacing over many crossings (rather than
/// using just one) also averages out per-crossing interpolation error.
[[nodiscard]] inline double measure_period(const Trajectory& traj) {
    std::vector<double> crossings;
    for (std::size_t i = 1; i < traj.size(); ++i) {
        const double x0 = traj.states[i - 1].r.x;
        const double x1 = traj.states[i].r.x;
        if (x0 < 0.0 && x1 >= 0.0) {
            const double t0 = traj.states[i - 1].t;
            const double t1 = traj.states[i].t;
            // Root of the line through (t0,x0)-(t1,x1): t0 + (0-x0)/(x1-x0) *
            // (t1-t0). x1 - x0 > 0 always holds here because x0 < 0 <= x1, so
            // this never divides by zero.
            crossings.push_back(t0 + (-x0) / (x1 - x0) * (t1 - t0));
        }
    }
    if (crossings.size() < 2) return 0.0;
    // Mean spacing = (last - first) / (count - 1): the sum of consecutive
    // differences between crossings telescopes to exactly this.
    return (crossings.back() - crossings.front()) /
           static_cast<double>(crossings.size() - 1);
}

/// Steady-state amplitude: max |r.x| over samples with t >= t_start. Pick
/// t_start several damping times (roughly 5/beta) after the drive turns on
/// so the transient (Taylor §5.6) has decayed and only the steady-state
/// response remains.
/// @param traj    trajectory to measure
/// @param t_start earliest time [s] to include in the max
[[nodiscard]] inline double measure_amplitude(const Trajectory& traj, double t_start) {
    double amplitude = 0.0;
    for (const auto& s : traj.states) {
        if (s.t >= t_start) amplitude = std::max(amplitude, std::fabs(s.r.x));
    }
    return amplitude;
}

}  // namespace newton::analysis
