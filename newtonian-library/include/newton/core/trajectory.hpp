#pragma once
// newton/core/trajectory.hpp — the recorded output of a simulation.
//
// A Trajectory is just the sequence of States a simulation passed through,
// plus the particle's mass (kept here so analysis functions — energy,
// momentum — can be computed from the trajectory alone). Producers:
// newton/simulate.hpp. Consumers: newton/analysis/*, newton/io/csv.hpp.

#include <cstddef>
#include <vector>

#include <newton/core/particle.hpp>

namespace newton {

struct Trajectory {
    double mass{1.0};           ///< mass of the simulated particle [kg]
    std::vector<State> states;  ///< samples in time order

    [[nodiscard]] std::size_t size() const { return states.size(); }
    [[nodiscard]] bool empty() const { return states.empty(); }
    [[nodiscard]] const State& front() const { return states.front(); }
    [[nodiscard]] const State& back() const { return states.back(); }
};

}  // namespace newton
