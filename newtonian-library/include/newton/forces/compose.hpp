#pragma once
// newton/forces/compose.hpp — combine several forces into one.
//
// Physics: Newton's second law only cares about the *net* force,
// F_net = sum_i F_i (superposition — each force acts as if the others
// weren't there). sum(...) makes that superposition explicit in code: a
// projectile with drag is just gravity + drag, evaluated and added.

#include <newton/core/vec3.hpp>

namespace newton::forces {

/// sum(f1, f2, ...) -> a callable satisfying the Force concept that adds up
/// every force it was given: F_net(t, r, v) = f1(t, r, v) + f2(t, r, v) + ...
///
/// Physics: e.g. a projectile under gravity and quadratic drag is
/// `sum(UniformGravity{m}, QuadraticDrag{c})` — no bespoke "gravity + drag"
/// struct needed; superposition of forces (Newton's second law) is exactly
/// vector addition.
///
/// C++ note: a variadic template over the callable types Fs..., each stored
/// by value in the returned lambda's capture — so `sum` owns copies of the
/// forces it was handed and the result stays valid after the temporaries
/// passed in go out of scope. The body is a C++17 *fold expression*:
/// `(Vec3{} + ... + fs(t, r, v))` expands to
/// `(((Vec3{} + f1(t,r,v)) + f2(t,r,v)) + ...)`, i.e. it left-folds vector
/// addition over the pack — and the explicit `Vec3{}` seed means `sum()`
/// with zero forces is well-formed too, returning the zero force. The return
/// type is `auto` because the lambda's closure type is unnameable — see
/// cpp-library-design.md's guidance on when `auto` returns are appropriate.
template <class... Fs>
[[nodiscard]] auto sum(Fs... fs) {
    return [fs...](double t, const Vec3& r, const Vec3& v) -> Vec3 {
        return (Vec3{} + ... + fs(t, r, v));
    };
}

}  // namespace newton::forces
