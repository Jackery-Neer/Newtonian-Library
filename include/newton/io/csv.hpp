#pragma once
// newton/io/csv.hpp — plain-text CSV export for trajectories and generic data.
//
// C++ note: implemented with <cstdio> (std::FILE*, std::fprintf) rather than
// <fstream> + <iomanip>. The library's precision contract is literally the C
// printf format "%.12g" (12 significant digits, shortest of fixed/scientific
// notation, no trailing zeros) — that is std::fprintf's native behavior, one
// argument, no manipulator soup. Reproducing %g exactly with iostreams needs
// std::defaultfloat + std::setprecision(12) *and* care around trailing zeros
// that %g strips automatically; <cstdio> gets the promised format for free.
// Both approaches are standard-library-only, satisfying DESIGN.md's "no
// dependencies beyond the C++ standard library."

#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include <newton/core/particle.hpp>
#include <newton/core/trajectory.hpp>

namespace newton::io {

/// Write a Trajectory to CSV: header "t,x,y,z,vx,vy,vz", one row per sample
/// in traj.states, each value formatted with %.12g (enough significant
/// digits to round-trip a double without full %.17g's visual noise).
/// Plot with `python3 scripts/plot_csv.py path.csv` or any spreadsheet/tool.
/// @param path path to the file to create/overwrite
/// @param traj trajectory to write
/// @returns false if the file could not be opened (or not fully flushed and
///          closed) for writing; true on success
[[nodiscard]] inline bool write_csv(const std::string& path, const Trajectory& traj) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (f == nullptr) return false;

    std::fprintf(f, "t,x,y,z,vx,vy,vz\n");
    for (const State& s : traj.states) {
        std::fprintf(f, "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                     s.t, s.r.x, s.r.y, s.r.z, s.v.x, s.v.y, s.v.z);
    }

    return std::fclose(f) == 0;
}

/// Generic column-oriented CSV: header row from column_names, then one row
/// per entry of `rows` (each inner vector should have the same length as
/// column_names — a short row just leaves trailing columns blank on that
/// line, a long one keeps writing past the header, so mismatches are visible
/// in the output rather than silently truncated). Values use the same
/// %.12g formatting as the Trajectory overload. Use this for sweep results:
/// write_csv(path, {"omega", "amplitude"}, rows).
/// @param path          path to the file to create/overwrite
/// @param column_names  CSV header, one entry per column
/// @param rows          one inner vector per output row
/// @returns false if the file could not be opened (or not fully flushed and
///          closed) for writing; true on success
[[nodiscard]] inline bool write_csv(const std::string& path,
                                    const std::vector<std::string>& column_names,
                                    const std::vector<std::vector<double>>& rows) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (f == nullptr) return false;

    for (std::size_t i = 0; i < column_names.size(); ++i) {
        std::fprintf(f, "%s%s", column_names[i].c_str(),
                     (i + 1 < column_names.size()) ? "," : "\n");
    }
    for (const std::vector<double>& row : rows) {
        for (std::size_t i = 0; i < row.size(); ++i) {
            std::fprintf(f, "%.12g%s", row[i], (i + 1 < row.size()) ? "," : "\n");
        }
    }

    return std::fclose(f) == 0;
}

}  // namespace newton::io
