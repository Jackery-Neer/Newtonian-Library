#pragma once
// newton/potentials/force_from_potential.hpp — from a potential energy to the
// force it produces: F = -∇U (Taylor §4.3).
//
// Physics: a conservative force points "downhill" on its potential-energy
// surface, F = -∇U = -(∂U/∂x, ∂U/∂y, ∂U/∂z). Give this file any U you can
// write as a function and it returns a force you can hand straight to
// simulate(), with the gradient taken numerically (fields::gradient, 5-point
// central differences). That makes "invent a potential, watch the motion" a
// one-liner, with no algebra and no chance of a sign error in a hand-derived
// force.
//
// C++ note: the result is a *function object* with the library's Force
// signature, so it composes with forces::sum, check_conservative,
// velocity_dependent_part and everything else that takes a Force.

#include <type_traits>
#include <utility>

#include <newton/core/vec3.hpp>
#include <newton/fields/derivatives.hpp>

namespace newton::potentials {

/// The force F = -∇U of a potential energy U, differentiated numerically.
///
/// Physics: F = -∇U (Taylor §4.3). Any single-valued, smooth U gives a
/// conservative force, because ∇×(∇U) = 0 identically (mixed partial derivatives
/// commute): a potential energy guarantees a conservative force
/// (check_conservative says so for the result). If U also depends on t, the
/// force is curl-free at each instant and E = T + U obeys dE/dt = ∂U/∂t
/// (§4.5); the verdict is then TimeDependentPotential. Only the position
/// gradient is taken: ∂U/∂t is not part of the force.
///
/// Pitfalls (the numerical gradient cannot see them coming):
///  - Single-valuedness. U = -atan2(y, x) looks like a potential but jumps by
///    2π across the negative x axis. Away from that cut its gradient is the
///    hole field, curl-free but not conservative around the axis (the force
///    is not a gradient of any single-valued U), and a stencil straddling
///    the cut returns a spike of order 2π/h.
///  - Kinks. U = k|x| gives a force discontinuous at x = 0. The finite-difference
///    gradient is wrong within 2h of the kink (the 5-point stencil overshoots by
///    7/6, a 17% error, then settles).
///  - Singular points. For U = -GMm/r the gradient is garbage within about 2h
///    of r = 0. GradientForce cannot know where the singularity is.
///  - Length scale. The step must be small against the distance over which U
///    changes: r0 for Lennard-Jones and Morse, the range a for Yukawa. The
///    default 1 m suits everyday mechanics; for a molecule pass
///    length_scale = 1e-10. With the default, an exponentially decaying U
///    (Morse, Yukawa) underflows to exactly 0 at the stencil points and the
///    force comes out exactly 0; a power law (Lennard-Jones) gives a wrong,
///    tiny number (3e-54 N where 1.3e-12 N is right).
///  - The step follows the distance from `center` (default: the ORIGIN), not
///    from the structure: h = 1e-4·max(|r - center|, length_scale). A well
///    centred 1000 m from the origin gets h = 0.1 m and loses accuracy as
///    (|r|/|r - c|)^4 (a force error of 1.7e-4 at 1 m from the well).
///    Pass the centre of the structure as `center`, or write U about the origin.
///  - Never put a reference point (find_potential) within ~2h of a singularity of
///    U: inside that core the numerical force is finite but not a gradient, and
///    integrals from there converge to a wrong, direction-dependent value.
///  - The numerical force is conservative only up to truncation and rounding
///    error, far below the classifier's threshold: check_conservative reports
///    a relative curl of 0 to ~1e-13 for the cases measured.
///
/// C++ note: the potential's callable shape is detected at compile time with
/// `if constexpr` + std::is_invocable_r_v. U(t, r) (two arguments) is preferred
/// when it is callable; otherwise U(r). Anything else is a static_assert with
/// a readable message. A struct with operator() and a factory function
/// (force_from_potential) that deduces the type for you is the library's
/// pattern for functors; compare forces::sum.
template <class PotentialFn>
struct GradientForce {
    PotentialFn U{};            // potential energy: double U(const Vec3& r) or double U(double t, const Vec3& r) [J]
    double length_scale{1.0};   // [m] the step is fields::default_step(r - center, length_scale): a FLOOR on |r - center|
    Vec3   center{};            // [m] where the structure is: the step grows with the distance from here (default: the origin)

    /// F(t, r, v) = -∇U at (r, t) [N]. The velocity is ignored: a potential
    /// force cannot depend on it (Taylor §4.2 condition (i)).
    /// @param t time [s]; used only if U depends on t
    /// @param r position [m]
    [[nodiscard]] Vec3 operator()(double t, const Vec3& r, const Vec3& /*v*/) const {
        constexpr bool time_dependent = std::is_invocable_r_v<double, const PotentialFn&, double, const Vec3&>;
        constexpr bool static_potential = std::is_invocable_r_v<double, const PotentialFn&, const Vec3&>;
        static_assert(time_dependent || static_potential,
                      "force_from_potential() needs a potential energy: a callable double U(const Vec3& r) "
                      "or double U(double t, const Vec3& r), callable on a const object");
        const double h = fields::default_step(r - center, length_scale);
        if constexpr (time_dependent) {
            auto frozen = [this, t](const Vec3& x) { return U(t, x); };
            return -fields::gradient(frozen, r, h);
        } else {
            (void)t;  // a static potential has no time to freeze
            auto field = [this](const Vec3& x) { return U(x); };
            return -fields::gradient(field, r, h);
        }
    }
};

/// Builds the force F = -∇U from a potential energy, so that `simulate` and
/// friends can use it directly:
///     auto F = force_from_potential([](const Vec3& r) { return 2.0 * norm2(r); });  // U = 2 r², F = -4 r
/// See GradientForce for the physics and the pitfalls.
///
/// C++ note: a factory function exists because C++17 cannot deduce the template
/// argument of an aggregate (`GradientForce{lambda, 1.0}` would not compile:
/// class template argument deduction for aggregates arrived in C++20); a
/// function template deduces PotentialFn from the lambda for you.
///
/// @param U            potential energy U(r) or U(t, r) [J], copied
/// @param length_scale smallest distance over which U changes appreciably [m]; default 1 m
/// @param center       where U's structure is centred [m] (the nucleus of a central potential); the
///                     finite-difference step is 1e-4·max(|r - center|, length_scale). Default: the origin.
template <class PotentialFn>
[[nodiscard]] GradientForce<PotentialFn> force_from_potential(PotentialFn U, double length_scale = 1.0,
                                                              const Vec3& center = {}) {
    return GradientForce<PotentialFn>{std::move(U), length_scale, center};
}

}  // namespace newton::potentials
