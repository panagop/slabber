"""Compare Czerny's nine four-edge plates with PyNite Rect and an independent FD grid."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np

from czerny_dataset.czerny import load_table, lookup
from czerny_dataset.pynite_czerny_benchmark import edge_codes
from scripts.four_edge_fd import solve_four_edge_fd
from scripts.run_pynite_ssss import (
    COMBO, E, LOAD, LX, NU, RATIOS, ROOT, THICKNESS, center_displacement,
    moment_extrema, relative_difference_percent, solve_rect_model, support_reactions,
)


CASES = ("SSSS", "F_y0", "F_x0", "FF_y", "FF_x", "FF_adj", "FFF_x", "FFF_y", "FFFF")
FOUR_EDGE_RATIOS = tuple(n / 100 for n in range(100, 201, 5))
FEM_REFINEMENT = (21, 31)
FD_REFINEMENT = (64, 128)
EDGE_NAMES = ("x=0", "y=0", "x=lx", "y=ly")
CORNERS = ((0, 0, 0, 1), (1, 0, 2, 1), (1, 1, 2, 3), (0, 1, 0, 3))
MOMENT_SCALE = LOAD * LX**2
DEFLECTION_SCALE = LOAD * LX**4 / (E * THICKNESS**3)
SCREEN_PERCENT = {"reaction_residual": 1e-6, "reference_difference": 2.0, "last_mesh_change": 1.0}
GLOBAL_QUANTITIES = ("mx_global_max", "my_global_max")


def sample(value, x, y, **details):
    return {"value": float(value), "x_m": float(x), "y_m": float(y), **details}


def czerny_samples(case_id, ratio, table):
    selected = {}
    for case, csv_quantity in table:
        if case != case_id or csv_quantity.startswith(("qx", "qy")):
            continue
        factor, source = lookup(case, ratio, csv_quantity, table=table)
        quantity = csv_quantity
        if case == "FF_adj" and csv_quantity in ("mx_fixed_mid", "my_fixed_mid"):
            # Printed p. 286 labels these m_xermin and m_yermin.
            quantity = csv_quantity.replace("_mid", "_min")
            source["printed_symbol_from_scan"] = f"m_{quantity[1]}ermin"
        scale = DEFLECTION_SCALE if quantity == "deflection_mid" else MOMENT_SCALE
        selected[quantity] = {"value": factor * scale, "coefficient": source["table_value"], "source": source}
    return selected


def rect_fixed_edges(model, nx, ny, edges):
    fixed = {}
    # Normal component index, boundary element index, along-edge count and local face.
    for edge, code, (axis, across, count, face) in zip(
        EDGE_NAMES, edges, ((0, 0, ny, 0), (1, 0, nx, 0), (0, nx - 1, ny, 1), (1, ny - 1, nx, 1)),
    ):
        if code != "F":
            continue

        def at(k, t):
            i, j = (across, k) if axis == 0 else (k, across)
            plate = model.plates[f"P_{i}_{j}"]
            u, v = (face, t) if axis == 0 else (t, face)
            x, y = u * plate.width(), v * plate.height()
            value = plate.moment(x, y, combo_name=COMBO)[axis, 0]
            return sample(value, plate.i_node.X + x, plate.i_node.Y + y, edge=edge, element=plate.name)

        fixed[edge] = {
            "midpoint": at(count // 2, 0.5),
            "minimum": min((at(k, t) for k in range(count) for t in (0, 1)), key=lambda point: point["value"]),
        }
    return fixed


def fd_fixed_edges(grid, edges):
    fixed = {}
    lines = (grid.mx[:, 0], grid.my[0, :], grid.mx[:, -1], grid.my[-1, :])
    for index, (edge, code, line) in enumerate(zip(EDGE_NAMES, edges, lines)):
        if code != "F":
            continue

        def at(k):
            i, j = ((0 if index == 0 else -1), k) if index % 2 == 0 else (k, (0 if index == 1 else -1))
            return sample(line[k], grid.x[i], grid.y[j], edge=edge)

        fixed[edge] = {"midpoint": at(len(line) // 2), "minimum": at(int(line.argmin()))}
    return fixed


def select_source_station_samples(points, fixed):
    # If both opposite edges are fixed, they are symmetric in this benchmark.
    for edge, row in fixed.items():
        axis = edge[0]
        points.setdefault(f"m{axis}_fixed_mid", row["midpoint"])
        points.setdefault(f"m{axis}_fixed_min", row["minimum"])


def corner_samples(edges):
    return [(u, v) for u, v, i, j in CORNERS if edges[i] == edges[j] == "S"]


def corner_magnitude(value, x, y):
    return sample(abs(value), x, y, signed_value=float(value))


def rect_mesh(case, ratio, nx, edges):
    model, ny = solve_rect_model(case, ratio, nx)
    return rect_mesh_results(model, ratio, nx, ny, edges)


def rect_mesh_results(model, ratio, nx, ny, edges):
    ly = LX * ratio
    center = model.plates[f"P_{nx // 2}_{ny // 2}"]
    mx, my, mxy = center.moment(center.width() / 2, center.height() / 2, combo_name=COMBO).ravel()
    w = center_displacement(center)
    points = {
        "mx_mid": sample(mx, LX / 2, ly / 2),
        "my_mid": sample(my, LX / 2, ly / 2),
        "mxy_mid": sample(mxy, LX / 2, ly / 2),
        "deflection_mid": sample(-w, LX / 2, ly / 2, signed_global_z_m=w),
    }
    extrema = moment_extrema(model.plates.values(), LX / 2)
    for field, component, key in (
        ("Mx_global", "Mx", "mx_global_max"),
        ("My_global", "My", "my_global_max"),
        ("My_mid_section", "My", "my_max"),
    ):
        peak = extrema[field]
        points[key] = sample(peak[component], peak["x_m"], peak["y_m"], element=peak["element"])

    fixed = rect_fixed_edges(model, nx, ny, edges)
    select_source_station_samples(points, fixed)
    reactions = support_reactions(model, nx, ny, ly)
    twisting = []
    corner_forces = []
    for u, v in corner_samples(edges):
        plate = model.plates[f"P_{u * (nx - 1)}_{v * (ny - 1)}"]
        mxy = plate.moment(u * plate.width(), v * plate.height(), combo_name=COMBO)[2, 0]
        twisting.append(corner_magnitude(mxy, u * LX, v * ly))
        force = model.nodes[f"N_{u * nx}_{v * ny}"].RxnFZ[COMBO]
        corner_forces.append(corner_magnitude(force, u * LX, v * ly))
    if twisting:
        points["mxy_corner"] = max(twisting, key=lambda point: point["value"])
        points["corner_force"] = max(corner_forces, key=lambda point: point["value"])

    return {
        "nx": nx, "ny": ny, "points": points, "fixed_edge_moments": fixed,
        "reactions": reactions, "simply_supported_corner_twisting": twisting,
    }


def fd_mesh(ratio, nx, edges):
    # Keep a grid node at the geometric midpoint on every printed ratio.
    ny, ly = 2 * round(nx * ratio / 2), LX * ratio
    grid = solve_four_edge_fd(edges, LX, ly, LOAD, E, THICKNESS, nx, ny)
    i, j = nx // 2, ny // 2
    points = {
        "mx_mid": sample(grid.mx[j, i], grid.x[i], grid.y[j]),
        "my_mid": sample(grid.my[j, i], grid.x[i], grid.y[j]),
        "deflection_mid": sample(grid.downward_deflection[j, i], grid.x[i], grid.y[j]),
    }
    for key, field in (("mx_global_max", grid.mx), ("my_global_max", grid.my)):
        peak_j, peak_i = np.unravel_index(field.argmax(), field.shape)
        points[key] = sample(field[peak_j, peak_i], grid.x[peak_i], grid.y[peak_j])
    peak_j = int(grid.my[:, i].argmax())
    points["my_max"] = sample(grid.my[peak_j, i], grid.x[i], grid.y[peak_j])
    fixed = fd_fixed_edges(grid, edges)
    select_source_station_samples(points, fixed)

    twisting = []
    for u, v in corner_samples(edges):
        # Odd reflection in both directions gives the mixed central derivative at an S/S corner.
        inward = grid.downward_deflection[1 if v == 0 else -2, 1 if u == 0 else -2]
        mxy = -grid.rigidity * (1 - 2 * u) * (1 - 2 * v) * inward / (LX / nx * ly / ny)
        twisting.append(corner_magnitude(mxy, u * LX, v * ly))
    if twisting:
        corner = max(twisting, key=lambda point: point["value"])
        points["mxy_corner"] = corner
        points["corner_force"] = corner_magnitude(-2 * corner["value"], corner["x_m"], corner["y_m"])

    return {"nx": nx, "ny": ny, "points": points, "fixed_edge_moments": fixed}


def coefficient(quantity, value):
    return value / DEFLECTION_SCALE if quantity == "deflection_mid" else MOMENT_SCALE / abs(value)


def errors(quantity, actual, reference):
    return {
        "physical_percent": relative_difference_percent(actual, reference),
        "coefficient_percent": relative_difference_percent(coefficient(quantity, actual), coefficient(quantity, reference)),
    }


def compare(source, fem, fd):
    comparisons = {}
    for quantity in dict.fromkeys((*source, "my_mid", *GLOBAL_QUANTITIES)):
        fem_value = fem[-1]["points"][quantity]["value"]
        fd_value = fd[-1]["points"][quantity]["value"]
        czerny = source.get(quantity)
        comparisons[quantity] = {
            "pynite": {**fem[-1]["points"][quantity], "coefficient": coefficient(quantity, fem_value)},
            "finite_difference": {**fd[-1]["points"][quantity], "coefficient": coefficient(quantity, fd_value)},
            "czerny": czerny,
            "pynite_vs_fd": errors(quantity, fem_value, fd_value),
            "pynite_last_refinement": errors(quantity, fem_value, fem[-2]["points"][quantity]["value"]),
            "fd_last_refinement": errors(quantity, fd_value, fd[-2]["points"][quantity]["value"]),
        }
        if czerny is not None:
            comparisons[quantity]["pynite_vs_czerny"] = errors(quantity, fem_value, czerny["value"])
            comparisons[quantity]["fd_vs_czerny"] = errors(quantity, fd_value, czerny["value"])
    return comparisons


def comparison_run(case, ratio, edges, source, fem, fd):
    return {
        "case_id": case, "ratio_ly_lx": ratio, "edges_x0_y0_xlx_yly": edges,
        "pynite_meshes": fem, "fd_meshes": fd, "comparisons": compare(source, fem, fd),
    }


def run_ratio(case, ratio, table, fem_meshes, fd_meshes):
    edges = edge_codes(case)
    source = czerny_samples(case, ratio, table)
    print(f"{case} ly/lx={ratio:.2f} edges={','.join(edges)}", flush=True)
    fem = []
    for nx in fem_meshes:
        mesh = rect_mesh(case, ratio, nx, edges)
        fem.append(mesh)
        print(f"  PyNite {nx}x{mesh['ny']}, signed force residual {mesh['reactions']['residual_percent']:+.3e}%", flush=True)
    fd = [fd_mesh(ratio, nx, edges) for nx in fd_meshes]
    print(f"  FD grids: {', '.join(f'{grid['nx']}x{grid['ny']}' for grid in fd)}", flush=True)
    return comparison_run(case, ratio, edges, source, fem, fd)


def comparison_rows(runs):
    for run in runs:
        for quantity, row in run["comparisons"].items():
            czerny = row["czerny"] or {}
            source = czerny.get("source", {})
            units = "m" if quantity == "deflection_mid" else "kN" if quantity == "corner_force" else "kN m/m"
            yield {
                "case_id": run["case_id"], "ratio_ly_lx": run["ratio_ly_lx"],
                "edges_x0_y0_xlx_yly": "/".join(run["edges_x0_y0_xlx_yly"]),
                "quantity": quantity, "units": units,
                "source_page": source.get("source_page"), "source_csv_quantity": source.get("quantity"),
                "source_symbol": source.get("printed_symbol_from_scan", source.get("original_symbol")),
                "source_value": czerny.get("value"), "source_coefficient": czerny.get("coefficient"),
                **{
                    f"{name}_{field}": row[method][field]
                    for name, method in (("pynite", "pynite"), ("fd", "finite_difference"))
                    for field in ("value", "coefficient", "x_m", "y_m")
                },
                **{
                    f"{comparison}_{field}": row.get(comparison, {}).get(field)
                    for comparison in ("pynite_vs_czerny", "fd_vs_czerny", "pynite_vs_fd", "pynite_last_refinement", "fd_last_refinement")
                    for field in ("physical_percent", "coefficient_percent")
                },
            }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    ratio_options = parser.add_mutually_exclusive_group()
    ratio_options.add_argument("--ratios", type=float, nargs="+", choices=FOUR_EDGE_RATIOS, default=RATIOS)
    ratio_options.add_argument("--all-ratios", action="store_true", help="Use the full 1.00–2.00 printed grid.")
    parser.add_argument("--fem-meshes", type=int, nargs="+", default=[7, 11, 15, *FEM_REFINEMENT])
    parser.add_argument("--fd-meshes", type=int, nargs="+", default=[32, *FD_REFINEMENT])
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "four_edge_benchmark.json")
    args = parser.parse_args()
    for name, meshes, parity, minimum in (("fem", args.fem_meshes, 1, 3), ("fd", args.fd_meshes, 0, 4)):
        if len(meshes) < 2 or meshes != sorted(set(meshes)) or any(n < minimum or n % 2 != parity for n in meshes):
            parser.error(f"--{name}-meshes requires at least two strictly increasing counts of the correct parity")
    table = load_table()
    runtime = {"python": platform.python_version(), "platform": platform.platform(),
               **{package: version(package) for package in ("PyNiteFEA", "numpy", "scipy")}}
    print(runtime, flush=True)
    ratios = FOUR_EDGE_RATIOS if args.all_ratios else args.ratios
    results = [run_ratio(case, ratio, table, args.fem_meshes, args.fd_meshes) for case in args.cases for ratio in ratios]
    report = {
        "runtime": runtime,
        "argv": ["python", "-m", "scripts.run_four_edge_benchmark", *sys.argv[1:]],
        "solver_models": {"pynite": "PyNiteFEA Rect linear static", "fd": "SciPy Kirchhoff finite differences"},
        "inputs": {"lx_m": LX, "load_kN_per_m2": LOAD, "E_kN_per_m2": E, "h_m": THICKNESS, "nu": NU},
        "screening_thresholds_percent": SCREEN_PERCENT,
        "conventions": {
            "moments": "Sagging positive in kN m/m. Fixed-edge values are the normal bending moment.",
            "deflection": "Positive down for comparisons, in m.",
            "reactions": "Global Z up, in kN. Signed boundary nodes summed once. Corner comparisons use magnitude.",
            "twisting": "Corner comparisons use magnitude, with the signed moment retained separately.",
            "fd_corner_force": "Negative twice the FD corner twisting magnitude from the Kirchhoff corner relation.",
        },
        "results": results,
    }
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    csv_path = args.output.with_suffix(".csv")
    rows = list(comparison_rows(results))
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(results)} case/ratio runs and {len(rows)} comparison rows to {args.output} and {csv_path}")


if __name__ == "__main__":
    main()
