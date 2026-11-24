#!/usr/bin/env python3
"""scripts/plot_csv.py — quick-look plotting for newton's CSV output.

Usage:
    python3 scripts/plot_csv.py build/foo.csv [xcol ycol...]

With no column names given, the first CSV column is used as x and every
other column is plotted against it. With column names given, the first
extra argument is the x column and every argument after it is a y column,
e.g.:

    python3 scripts/plot_csv.py build/projectile_drag.csv x y
    python3 scripts/plot_csv.py build/energy_drift.csv t E_euler E_rk4

Degrades gracefully if matplotlib isn't installed: prints the available
column names (so you at least know what you *could* have plotted) and exits
cleanly rather than raising an ImportError traceback. This is deliberately
stdlib (csv, sys) + matplotlib and nothing else, matching the library's own
"no dependencies beyond what's necessary" policy for build/tooling code.
"""

import csv
import sys


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: python3 scripts/plot_csv.py <path.csv> [xcol ycol...]", file=sys.stderr)
        return 1

    path = sys.argv[1]
    extra_cols = sys.argv[2:]

    try:
        with open(path, newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = [[float(v) for v in row] for row in reader if row]
    except FileNotFoundError:
        print(f"error: no such file: {path}", file=sys.stderr)
        return 1
    except (StopIteration, ValueError) as exc:
        print(f"error: could not parse {path} as a newton CSV: {exc}", file=sys.stderr)
        return 1

    if not header or not rows:
        print(f"{path}: no data to plot", file=sys.stderr)
        return 1

    # One list of values per named column, in file order.
    columns = {name: [row[i] for row in rows] for i, name in enumerate(header)}

    if extra_cols:
        xcol, *ycols = extra_cols
    else:
        # Default: first column is x, every remaining column is a y series
        # (for a trajectory CSV that's t vs x,y,z,vx,vy,vz all at once — busy,
        # but a reasonable "just show me something" default).
        xcol, ycols = header[0], header[1:]

    for name in [xcol, *ycols]:
        if name not in columns:
            print(f"error: no column named '{name}' in {path} (have: {', '.join(header)})", file=sys.stderr)
            return 1
    if not ycols:
        print(f"error: need at least one y column (have: {', '.join(header)})", file=sys.stderr)
        return 1

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed -- skipping the plot.")
        print("install it with: pip install matplotlib")
        print(f"columns available in {path}: {', '.join(header)}")
        return 0

    x = columns[xcol]
    fig, ax = plt.subplots()
    for name in ycols:
        ax.plot(x, columns[name], label=name)
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycols[0] if len(ycols) == 1 else "value")
    ax.set_title(path)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
