#pragma once
// newton/forces/forces.hpp — module umbrella for newton::forces.
//
// Includes every force-law functor and the compose() helper. Include this
// (or newton/newton.hpp) rather than the individual headers unless you're
// trying to minimize compile time in a large project.

#include <newton/forces/compose.hpp>
#include <newton/forces/drag.hpp>
#include <newton/forces/drives.hpp>
#include <newton/forces/electromagnetic.hpp>
#include <newton/forces/gravity.hpp>
#include <newton/forces/spring.hpp>
