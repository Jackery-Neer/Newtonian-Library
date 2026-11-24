#pragma once
// newton/newton.hpp — umbrella header. Include this to get the whole library:
//
//     #include <newton/newton.hpp>
//
// In a large project you would include individual module headers to keep
// compile times down; at this library's size the umbrella is always fine.
//
// C++ note: the library's own headers use angle brackets and full paths
// (<newton/core/vec3.hpp>) — the same spelling a user writes. That works
// because `-I include` puts our include root on the search path, and it means
// internal includes never break when files move relative to each other.
//
// This file is also the library's manifest: every public header is listed
// here, grouped by module.

// core — vocabulary types used everywhere
#include <newton/core/constants.hpp>
#include <newton/core/vec3.hpp>
#include <newton/core/particle.hpp>
#include <newton/core/trajectory.hpp>

// coords — cylindrical and spherical polar coordinates; simulate in Cartesian,
// convert at the edges (Taylor §1.7, §4.8)
#include <newton/coords/cylindrical_polar.hpp>
#include <newton/coords/spherical_polar.hpp>
#include <newton/coords/cartesian.hpp>

// forces — force-law functors, composable with newton::forces::sum(...)
#include <newton/forces/forces.hpp>

// integrators — one-step methods State -> State
#include <newton/integrators/euler.hpp>
#include <newton/integrators/semi_implicit_euler.hpp>
#include <newton/integrators/velocity_verlet.hpp>
#include <newton/integrators/rk4.hpp>

// simulation driver — the time-stepping loop, recording, stop events
#include <newton/simulate.hpp>

// analytic — closed-form textbook solutions (the validation oracle)
#include <newton/analytic/projectile.hpp>
#include <newton/analytic/drag.hpp>
#include <newton/analytic/oscillator.hpp>
#include <newton/analytic/magnetic.hpp>
#include <newton/analytic/rocket.hpp>

// analysis — conserved-quantity diagnostics and trajectory measurements
#include <newton/analysis/energy.hpp>
#include <newton/analysis/momentum.hpp>
#include <newton/analysis/oscillation.hpp>

// systems — self-contained special-purpose simulators
#include <newton/systems/rocket.hpp>

// io — trajectory export for plotting
#include <newton/io/csv.hpp>
