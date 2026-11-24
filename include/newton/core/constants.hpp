#pragma once
// newton/core/constants.hpp — physical constants and unit helpers.
//
// Policy for the whole library: SI units everywhere, angles in radians,
// `double` everywhere. Textbook constants then work with no conversion layer
// (and no conversion bugs).
//
// C++ note: `inline constexpr` (C++17) gives each constant exactly ONE object
// across the whole program even though this header is pasted into many
// translation units — this is the modern replacement for `#define G 6.67e-11`,
// with a real type, a namespace, and debugger visibility. The digit separator
// (') is C++14 and purely cosmetic.

namespace newton::constants {

inline constexpr double pi = 3.14159265358979323846;

/// Multiply degrees by this to get radians: `45.0 * deg`.
inline constexpr double deg = pi / 180.0;

/// Standard gravitational acceleration [m/s^2] (CODATA conventional value).
inline constexpr double g_earth = 9.80665;

/// Newtonian gravitational constant [m^3 kg^-1 s^-2] (CODATA 2018).
inline constexpr double G = 6.67430e-11;

// --- air drag on a smooth sphere at STP (Taylor §2.1) ----------------------
// Taylor writes the drag magnitude f(v) = b v + c v^2 with b = beta * D and
// c = gamma * D^2, where D is the sphere's diameter. These are his quoted
// coefficients for air at standard conditions.

/// Linear drag geometry coefficient: b = sphere_air_beta * D  [N·s/m^2].
inline constexpr double sphere_air_beta = 1.6e-4;

/// Quadratic drag geometry coefficient: c = sphere_air_gamma * D^2  [N·s^2/m^4].
inline constexpr double sphere_air_gamma = 0.25;

}  // namespace newton::constants
