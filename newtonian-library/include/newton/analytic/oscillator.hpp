#pragma once
// newton/analytic/oscillator.hpp — closed-form oscillator solutions
// (Taylor ch. 5). The centerpiece of newton::analytic: nearly everything in
// this chapter is exactly solvable, which makes it the richest source of
// validation cases in the whole library (research map §8, rows 23-30).
//
// Physics: three linear ODEs, in increasing generality:
//   SHM              ẍ + ω0² x = 0                  (§5.2)
//   DampedOscillator ẍ + 2β ẋ + ω0² x = 0            (§5.4)
//   DrivenOscillator ẍ + 2β ẋ + ω0² x = f0 cos(ω t)  (§5.5-5.6, steady state)
// All three describe a scalar coordinate x(t) — unlike projectile/drag, this
// is intrinsically 1D (Taylor's spring-mass-damper, or the electrical LRC
// analogue), so no Vec3 appears in this header.
//
// C++ note: std::sin/cos/exp/sqrt/atan2 are not constexpr in C++17, so most
// member functions here are plain (non-constexpr) inline — implicitly so,
// since they are defined inside the class body. The exceptions are noted
// where the formula happens to need no transcendental call.

#include <cmath>

#include <newton/core/constants.hpp>

namespace newton::analytic {

/// Simple harmonic motion x(t) = A cos(ω0 t - δ), fitted to initial
/// conditions x(0) = x0, v(0) = v0 (Taylor §5.2).
struct SHM {
    double omega0{1.0};  ///< natural frequency, sqrt(k/m) [rad/s]
    double x0{0.0};      ///< initial position [m]
    double v0{0.0};      ///< initial velocity [m/s]

    /// Amplitude A = sqrt(x0² + (v0/ω0)²) [m].
    [[nodiscard]] double amplitude() const {
        const double c2 = v0 / omega0;
        return std::sqrt(x0 * x0 + c2 * c2);
    }

    /// Phase δ (rad) such that x0 = A cos(δ), v0/ω0 = A sin(δ).
    [[nodiscard]] double phase() const { return std::atan2(v0 / omega0, x0); }

    /// Position at time t.
    ///
    /// C++ note: implemented as x0 cos(ω0 t) + (v0/ω0) sin(ω0 t) rather than
    /// amplitude()*cos(ω0 t - phase()) — algebraically identical (this is
    /// the trig-identity expansion of the same cosine), but it makes
    /// x(0) == x0 exactly to the last bit (cos(0)=1, sin(0)=0 exactly),
    /// instead of merely to floating-point tolerance through an
    /// atan2-then-cos round trip. That exactness is exactly what
    /// tests/test_analytic.cpp checks.
    [[nodiscard]] double x(double t) const {
        return x0 * std::cos(omega0 * t) + (v0 / omega0) * std::sin(omega0 * t);
    }

    /// Velocity at time t (derivative of x(t) above; same exactness note).
    [[nodiscard]] double v(double t) const {
        return -x0 * omega0 * std::sin(omega0 * t) + v0 * std::cos(omega0 * t);
    }

    /// Period T = 2π/ω0 [s] — independent of amplitude, the defining
    /// feature of SHM (Taylor §5.2). Pure arithmetic: constexpr.
    [[nodiscard]] constexpr double period() const { return 2.0 * constants::pi / omega0; }

    /// Total mechanical energy E = ½ m v0² + ½ m ω0² x0² = ½ m ω0² A²
    /// (Taylor §5.2) — constant for all t since SHM is conservative.
    ///
    /// C++ note: computed as ½m(ω0²x0² + v0²) rather than ½m·ω0²·amplitude()²
    /// — the same identity used in x()/v() above, expanded so no sqrt is
    /// needed (amplitude() calls std::sqrt, which would make this function
    /// non-constexpr for no benefit). This lets energy() be constexpr.
    /// @param m mass [kg]
    [[nodiscard]] constexpr double energy(double m) const {
        return 0.5 * m * (omega0 * omega0 * x0 * x0 + v0 * v0);
    }
};

/// Which of the four textbook regimes a damped oscillator is in
/// (Taylor §5.4), classified by β vs ω0.
enum class Regime { Undamped, Underdamped, CriticallyDamped, Overdamped };

/// Damped oscillator ẍ + 2β ẋ + ω0² x = 0 (Taylor §5.4), fitted to initial
/// conditions x(0) = x0, v(0) = v0. Covers all four regimes; x(t)/v(t)
/// dispatch on regime() so the caller never needs to branch.
struct DampedOscillator {
    double omega0;   ///< undamped natural frequency, sqrt(k/m) [rad/s]
    double beta;     ///< damping rate, b/(2m) [rad/s]
    double x0{1.0};  ///< initial position [m]
    double v0{0.0};  ///< initial velocity [m/s]

    /// Classifies the regime from β vs ω0.
    ///
    /// Physics/numerics: β exactly 0 is Undamped (a meaningful physical
    /// distinction from "very lightly damped", not just a numerical
    /// shortcut). Near β == ω0, both the underdamped formula (divides by
    /// ω1 = sqrt(ω0²-β²)) and the overdamped formula (divides by
    /// s = sqrt(β²-ω0²)) divide by a quantity going to zero — the repeated
    /// root β = ω0 is a genuine mathematical singularity of those two
    /// branches, not of the physics (the critically-damped formula stays
    /// perfectly well-behaved there). So |β - ω0| below a small relative
    /// epsilon is classified CriticallyDamped, routing x(t)/v(t) to the
    /// stable (C1 + C2 t) e^(-βt) formula instead of blowing up 0/0.
    [[nodiscard]] Regime regime() const {
        constexpr double kRelEps = 1e-9;  // relative to omega0
        if (beta == 0.0) return Regime::Undamped;
        if (std::abs(beta - omega0) <= kRelEps * omega0) return Regime::CriticallyDamped;
        return beta < omega0 ? Regime::Underdamped : Regime::Overdamped;
    }

    /// Damped angular frequency ω1 = sqrt(ω0² - β²) [rad/s] (Taylor §5.4).
    /// Only real (and only physically meaningful) for Undamped/Underdamped;
    /// for β > ω0 the argument is negative and this returns NaN — a correct
    /// reflection that overdamped motion doesn't oscillate, so it has no
    /// "damped frequency". See regime() for why callers of x()/v() don't
    /// need to worry about this directly.
    [[nodiscard]] double omega1() const { return std::sqrt(omega0 * omega0 - beta * beta); }

    /// Position at time t, correct in all four regimes (dispatches on
    /// regime()).
    [[nodiscard]] double x(double t) const {
        switch (regime()) {
            case Regime::Overdamped: {
                const double s = std::sqrt(beta * beta - omega0 * omega0);
                const double r_plus = -beta + s;
                const double r_minus = -beta - s;
                const double c_plus = (v0 - r_minus * x0) / (r_plus - r_minus);
                const double c_minus = -(v0 - r_plus * x0) / (r_plus - r_minus);
                return c_plus * std::exp(r_plus * t) + c_minus * std::exp(r_minus * t);
            }
            case Regime::CriticallyDamped: {
                // x = (C1 + C2 t) e^(-βt), C1 = x0, C2 = v0 + β x0.
                const double c1 = x0;
                const double c2 = v0 + beta * x0;
                return (c1 + c2 * t) * std::exp(-beta * t);
            }
            case Regime::Undamped:
            case Regime::Underdamped:
            default: {
                // x = e^(-βt)[x0 cos(ω1 t) + ((v0+βx0)/ω1) sin(ω1 t)].
                // At β=0 this is exactly ω1=ω0, reducing to plain SHM — no
                // special case needed for Undamped.
                const double w1 = omega1();
                const double b_coeff = (v0 + beta * x0) / w1;
                return std::exp(-beta * t) * (x0 * std::cos(w1 * t) + b_coeff * std::sin(w1 * t));
            }
        }
    }

    /// Velocity at time t (derivative of x(t) above), correct in all four
    /// regimes.
    [[nodiscard]] double v(double t) const {
        switch (regime()) {
            case Regime::Overdamped: {
                const double s = std::sqrt(beta * beta - omega0 * omega0);
                const double r_plus = -beta + s;
                const double r_minus = -beta - s;
                const double c_plus = (v0 - r_minus * x0) / (r_plus - r_minus);
                const double c_minus = -(v0 - r_plus * x0) / (r_plus - r_minus);
                return c_plus * r_plus * std::exp(r_plus * t) + c_minus * r_minus * std::exp(r_minus * t);
            }
            case Regime::CriticallyDamped: {
                const double c1 = x0;
                const double c2 = v0 + beta * x0;
                return (c2 - beta * (c1 + c2 * t)) * std::exp(-beta * t);
            }
            case Regime::Undamped:
            case Regime::Underdamped:
            default: {
                // Differentiating x(t) above and simplifying with
                // ω1² = ω0² - β² collapses the cos coefficient to exactly
                // v0 — a clean identity worth trusting rather than
                // re-deriving by finite differences.
                const double w1 = omega1();
                const double sin_coeff = (beta * v0 + omega0 * omega0 * x0) / w1;
                return std::exp(-beta * t) * (v0 * std::cos(w1 * t) - sin_coeff * std::sin(w1 * t));
            }
        }
    }

    /// Logarithmic decrement ln(x_n / x_{n+1}) = 2πβ/ω1 (Taylor §5.4) — the
    /// ratio of successive oscillation peaks, a standard experimental way to
    /// extract β from a ring-down trace. Only meaningful for Underdamped
    /// (peaks require oscillation); NaN propagates naturally from omega1()
    /// otherwise.
    [[nodiscard]] double log_decrement() const { return 2.0 * constants::pi * beta / omega1(); }
};

/// Steady-state response of a damped, sinusoidally driven oscillator
/// ẍ + 2β ẋ + ω0² x = f0 cos(ω t) (Taylor §5.5-5.6), i.e. x_p(t) =
/// amplitude() * cos(ω t - phase()). The transient (which depends on
/// initial conditions and decays like e^(-βt)) is deliberately not modeled
/// here — see DampedOscillator for that piece; superposing the two is the
/// student's exercise (Taylor §5.5).
struct DrivenOscillator {
    double omega0{1.0};  ///< undamped natural frequency, sqrt(k/m) [rad/s]
    double beta{0.0};    ///< damping rate, b/(2m) [rad/s]
    double f0{0.0};      ///< driving acceleration amplitude, F0/m [m/s^2]
    double omega{1.0};   ///< driving angular frequency [rad/s]

    /// Steady-state amplitude A(ω) = f0 / sqrt((ω0²-ω²)² + 4β²ω²)
    /// (Taylor §5.6).
    [[nodiscard]] double amplitude() const {
        const double detune = omega0 * omega0 - omega * omega;
        return f0 / std::sqrt(detune * detune + 4.0 * beta * beta * omega * omega);
    }

    /// Steady-state phase lag δ ∈ (0, π), tan δ = 2βω/(ω0²-ω²) (Taylor
    /// §5.6): δ -> 0 far below resonance, δ = π/2 exactly at ω = ω0, δ -> π
    /// far above resonance.
    ///
    /// C++ note: atan2's numerator 2βω is >= 0 for physical β, ω >= 0, which
    /// is exactly what keeps atan2 in [0, π] rather than the full (-π, π]
    /// range — no extra branching needed to enforce the documented domain.
    [[nodiscard]] double phase() const {
        return std::atan2(2.0 * beta * omega, omega0 * omega0 - omega * omega);
    }

    /// Driving frequency at which amplitude() peaks: ω_res = sqrt(ω0²-2β²)
    /// (Taylor §5.6). Only real for β < ω0/√2 (weak damping); for stronger
    /// damping A(ω) decreases monotonically from ω=0 and has no interior
    /// maximum, so this correctly returns NaN — there is no resonance peak
    /// to report.
    [[nodiscard]] double resonant_omega() const { return std::sqrt(omega0 * omega0 - 2.0 * beta * beta); }

    /// Quality factor Q = ω0/(2β) (Taylor §5.6) — dimensionless measure of
    /// how sharply peaked the resonance is. β = 0 gives Q = +infinity, the
    /// correct limit for an undamped (lossless) oscillator.
    [[nodiscard]] double quality_factor() const { return omega0 / (2.0 * beta); }

    /// Full width at half maximum of A²(ω): FWHM = 2β (Taylor §5.6) — the
    /// width in ω over which the power response A² stays above half its
    /// peak value, for weak damping (β ≪ ω0) where A²(ω) is approximately
    /// Lorentzian.
    [[nodiscard]] constexpr double fwhm() const { return 2.0 * beta; }
};

}  // namespace newton::analytic
