"""Compare all six Hahn three-edge cases with PyNite Rect and Kirchhoff Rayleigh–Ritz."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from czerny_dataset.czerny import THREE_EDGE_CSV, load_table, lookup
from czerny_dataset.pynite_czerny_benchmark import edge_codes
from scripts.kirchhoff_ritz import solve_ritz
from scripts.run_four_edge_benchmark import sample
from scripts.run_pynite_ssss import (
    COMBO, E, LOAD, LX, NU, ROOT, THICKNESS, center_displacement, reaction_balance,
    relative_difference_percent, solve_rect_model, support_reactions,
)


CASES = (
    "H17_SSS", "H18_REAR_FIXED", "H19_ONE_SIDE_FIXED", "H20_BOTH_SIDES_FIXED",
    "H21_ADJ_FIXED", "H22_ALL_FIXED",
)
RATIOS = (0.5, 1.0, 1.4)
STATION_GRID = 241
LIMITS = {"reaction_residual": 1e-6, "reference_difference": 2.0, "last_refinement": 1.0}
DIAGNOSTICS = ("mx_global_max", "my_global_max", "my_center", "deflection_ref")


def hahn_reference(case, ratio, table):
    references = {}
    for table_case, quantity in table:
        if case != table_case:
            continue
        factor, metadata = lookup(case, ratio, quantity, table=table)
        if metadata["interpolated"]:
            raise ValueError(f"{case} ratio {ratio} is between printed Hahn knots")
        load = LOAD * LX**2 * ratio
        scales = {"K": load, "K_lx2_over_Ed3": load * LX**2 / (E * THICKNESS**3), "1": 1.0}
        scale = scales[metadata["scale_basis"]]
        references[quantity] = {"value": factor * scale, "scale": scale, "source": metadata}
    return references


def free_mid_deflection(plate):
    # The Rect displacement along its bottom edge is a cubic Hermite trace.
    d = plate.d(COMBO).ravel()
    return (d[2] + d[8]) / 2 + plate.width() * (d[10] - d[4]) / 8


def negative_reach(section, ly):
    """Locate the first zero reached from the rear on the centerline."""
    for (y0, m0), (y1, m1) in zip(section[-2::-1], section[:0:-1]):
        if m0 >= 0 > m1:
            zero = y0 + (y1 - y0) * m0 / (m0 - m1)
            return sample((ly - zero) / ly, LX / 2, zero)
    raise ValueError("No negative-to-positive centerline transition from the rear")


def uses_rear_peak(case, ratio):
    return (case == "H19_ONE_SIDE_FIXED" and ratio > 0.8) or (
        case == "H20_BOTH_SIDES_FIXED" and ratio > 0.6
    )


def rect_moment(plate, x, y, component, **details):
    value = plate.moment(x, y, combo_name=COMBO)[component, 0]
    return sample(value, plate.i_node.X + x, plate.i_node.Y + y, element=plate.name, **details)


def rect_peaks(model, ly):
    mx_peak = my_peak = rear_peak = None
    for plate in model.plates.values():
        width, height = plate.width(), plate.height()
        y0 = plate.i_node.Y
        # Bilinear element moments attain extrema at corners. Include the cut
        # at the mid-depth of the panel for the rear half of the field.
        ys = (0.0, height)
        if y0 < ly / 2 < y0 + height:
            ys += (ly / 2 - y0,)
        for y in ys:
            for x in (0.0, width):
                mx, my, _ = plate.moment(x, y, combo_name=COMBO).ravel()
                px, py = plate.i_node.X + x, y0 + y
                if mx_peak is None or mx > mx_peak["value"]:
                    mx_peak = sample(mx, px, py, element=plate.name)
                if my_peak is None or my > my_peak["value"]:
                    my_peak = sample(my, px, py, element=plate.name)
                if py >= ly / 2 and (rear_peak is None or my > rear_peak["value"]):
                    rear_peak = sample(my, px, py, element=plate.name)
    return mx_peak, my_peak, rear_peak


def rect_run(case, ratio, nx):
    ly = LX * ratio
    model, ny = solve_rect_model(case, ratio, nx)
    i, j = nx // 2, ny // 2
    center = model.plates[f"P_{i}_{j}"]
    bottom = model.plates[f"P_{i}_0"]
    rear = model.plates[f"P_{i}_{ny - 1}"]
    left_mid = model.plates[f"P_0_{j}"]
    lower_left = model.plates["P_0_0"]
    upper_left = model.plates[f"P_0_{ny - 1}"]
    cx, cy = center.width() / 2, center.height() / 2
    mx_peak, my_peak, rear_peak = rect_peaks(model, ly)
    points = {
        "mx_center": rect_moment(center, cx, cy, 0),
        "my_center": rect_moment(center, cx, cy, 1),
        "mx_free_mid": rect_moment(bottom, bottom.width() / 2, 0, 0),
        "my_free_mid_diagnostic": rect_moment(bottom, bottom.width() / 2, 0, 1),
        "deflection_ref": sample(-free_mid_deflection(bottom), LX / 2, 0),
        "deflection_center_diagnostic": sample(-center_displacement(center), LX / 2, ly / 2),
        "mx_global_max": mx_peak,
        "my_global_max": my_peak,
    }
    if case.startswith(("H19", "H20")):
        points["my_center_or_max"] = rear_peak if uses_rear_peak(case, ratio) else points["my_center"]
    if case.startswith(("H19", "H20", "H21", "H22")):
        points["mx_fixed_free_corner"] = rect_moment(lower_left, 0, 0, 0)
        points["mx_fixed_side_mid"] = rect_moment(left_mid, 0, left_mid.height() / 2, 0)
    if case.startswith(("H18", "H21", "H22")):
        points["my_rear_fixed"] = rect_moment(rear, rear.width() / 2, rear.height(), 1)
    if case.startswith(("H17", "H18")):
        points["mxy_free_corner"] = rect_moment(lower_left, 0, 0, 2)
    if case == "H18_REAR_FIXED":
        section = []
        for row in range(ny):
            plate = model.plates[f"P_{i}_{row}"]
            for y in (0.0, plate.height()):
                section.append(rect_moment(plate, plate.width() / 2, y, 1))
        points["my_positive_adjacent"] = max(section, key=lambda point: point["value"])
        points["negative_moment_reach"] = negative_reach([(p["y_m"], p["value"]) for p in section], ly)

    reactions = support_reactions(model, nx, ny, ly)
    open_edges = reactions["open_edge_resultants_kN"]
    free_nodes = reactions["support_nodes_by_location"]["y=0"]
    if any(model.nodes[f"N_{column}_0"].support_DZ for column in range(1, nx)):
        raise ValueError("An interior free-edge node has a vertical constraint")
    reactions["free_edge_max_abs_nodal_fz_kN"] = max(abs(row["fz_kN"]) for row in free_nodes)
    if case == "H17_SSS":
        points["mxy_rear_corner"] = rect_moment(upper_left, 0, upper_left.height(), 2)
        points["reaction_x_total"] = sample((open_edges["x=0"] + open_edges["x=lx"]) / 2, LX / 2, ly / 2)
        points["reaction_y_total"] = sample(open_edges["y=ly"], LX / 2, ly)
    return {"nx": nx, "ny": ny, "points": points, "reactions": reactions}


def ritz_at(field, x, y, component):
    return sample(field.moments(x, y)[component].item(), x, y)


def ritz_reactions(field, ly, degree):
    """Integrate outward effective shear and retain the signed twisting jumps."""
    stations, weights = leggauss(degree + 1)
    gx, gy = (stations + 1) * LX / 2, (stations + 1) * ly / 2
    sides = {
        "x=0": -float(weights @ field.effective_shear(0, gy, "x").ravel()) * ly / 2,
        "x=lx": float(weights @ field.effective_shear(LX, gy, "x").ravel()) * ly / 2,
        "y=ly": float(weights @ field.effective_shear(gx, ly, "y").ravel()) * LX / 2,
    }
    corners = []
    for u, v in ((0, 0), (1, 0), (1, 1), (0, 1)):
        mxy = ritz_at(field, u * LX, v * ly, 2)["value"]
        corners.append({"x_m": u * LX, "y_m": v * ly, "fz_kN": 2 * (1 - 2 * u) * (1 - 2 * v) * mxy})
    free_residual = -float(weights @ field.effective_shear(gx, 0, "y").ravel()) * LX / 2
    return {
        **reaction_balance(-LOAD * LX * ly, [*sides.values(), *(c["fz_kN"] for c in corners)]),
        "open_edge_resultants_kN": sides,
        "signed_corner_twisting_forces": corners,
        "free_edge_effective_shear_integral_kN": free_residual,
    }


def ritz_run(case, ratio, degree, edges):
    ly = LX * ratio
    field = solve_ritz(edges, LX, ly, LOAD, E, THICKNESS, NU, degree)
    x, y = np.linspace(0, LX, STATION_GRID), np.linspace(0, ly, STATION_GRID)
    mx, my, _ = field.moments(x, y)
    mid = STATION_GRID // 2

    def peak(values, first_row=0):
        j, i = np.unravel_index(values.argmax(), values.shape)
        return sample(values[j, i], x[i], y[j + first_row])

    points = {
        "mx_center": sample(mx[mid, mid], LX / 2, ly / 2),
        "my_center": sample(my[mid, mid], LX / 2, ly / 2),
        "mx_free_mid": sample(mx[0, mid], LX / 2, 0),
        "my_free_mid_diagnostic": sample(my[0, mid], LX / 2, 0),
        "deflection_ref": sample(field.partial(LX / 2, 0).item(), LX / 2, 0),
        "deflection_center_diagnostic": sample(field.partial(LX / 2, ly / 2).item(), LX / 2, ly / 2),
        "mx_global_max": peak(mx),
        "my_global_max": peak(my),
    }
    if case.startswith(("H19", "H20")):
        points["my_center_or_max"] = peak(my[mid:], mid) if uses_rear_peak(case, ratio) else points["my_center"]
    if case.startswith(("H19", "H20", "H21", "H22")):
        points["mx_fixed_free_corner"] = sample(mx[0, 0], 0, 0)
        points["mx_fixed_side_mid"] = sample(mx[mid, 0], 0, ly / 2)
    if case.startswith(("H18", "H21", "H22")):
        points["my_rear_fixed"] = sample(my[-1, mid], LX / 2, ly)
    if case.startswith(("H17", "H18")):
        points["mxy_free_corner"] = ritz_at(field, 0, 0, 2)
    if case == "H18_REAR_FIXED":
        j = int(my[:, mid].argmax())
        points["my_positive_adjacent"] = sample(my[j, mid], LX / 2, y[j])
        points["negative_moment_reach"] = negative_reach(list(zip(y, my[:, mid])), ly)

    record = {"degree": degree, "points": points,
              "free_edge_effective_shear_mid_kN_per_m": -field.effective_shear(LX / 2, 0, "y").item()}
    if case == "H17_SSS":
        reactions = ritz_reactions(field, ly, degree)
        sides = reactions["open_edge_resultants_kN"]
        points["mxy_rear_corner"] = ritz_at(field, 0, ly, 2)
        points["reaction_x_total"] = sample((sides["x=0"] + sides["x=lx"]) / 2, LX / 2, ly / 2)
        points["reaction_y_total"] = sample(sides["y=ly"], LX / 2, ly)
        record["reactions_from_effective_shear"] = reactions
    return record


def normalize(value, scale, kind):
    if kind == "denominator":
        return scale / abs(value)
    return value / scale


def comparison_values(quantity, point, scale, kind):
    raw = point["value"]
    value = abs(raw) if quantity.startswith("mxy_") else raw
    return {**point, "value": value, "signed_value": raw, "coefficient": normalize(value, scale, kind)}


def differences(actual, reference):
    return {
        field: relative_difference_percent(actual[field], reference[field])
        for field in ("value", "coefficient")
    }


def compare(ratio, source, fem, ritz):
    rows = {}
    for quantity in dict.fromkeys((*source, *DIAGNOSTICS)):
        reference = source.get(quantity)
        load = LOAD * LX**2 * ratio
        scale = reference["scale"] if reference else (
            load * LX**2 / (E * THICKNESS**3) if quantity == "deflection_ref" else load
        )
        kind = reference["source"]["table_value_kind"] if reference else "denominator"
        px0, px1 = (comparison_values(quantity, mesh["points"][quantity], scale, kind) for mesh in fem[-2:])
        rr0, rr1 = (comparison_values(quantity, mesh["points"][quantity], scale, kind) for mesh in ritz[-2:])
        row = {
            "pynite": px1, "ritz": rr1, "hahn": reference,
            "pynite_vs_ritz": differences(px1, rr1),
            "pynite_last_refinement": differences(px1, px0),
            "ritz_last_refinement": differences(rr1, rr0),
        }
        if reference:
            hahn = {"value": reference["value"], "coefficient": reference["source"]["table_value"]}
            row["pynite_vs_hahn"] = differences(px1, hahn)
            row["ritz_vs_hahn"] = differences(rr1, hahn)
        rows[quantity] = row
    return rows


def run_case(case, ratio, table, meshes, degrees):
    edges = edge_codes(case)
    source = hahn_reference(case, ratio, table)
    print(f"{case} ly/lx={ratio:.2f} edges={','.join(edges)}", flush=True)
    fem = []
    for nx in meshes:
        mesh = rect_run(case, ratio, nx)
        fem.append(mesh)
        reactions = mesh["reactions"]
        print(f"  PyNite {nx}x{mesh['ny']} signed residual {reactions['residual_percent']:+.3e}%", flush=True)
    ritz = [ritz_run(case, ratio, degree, edges) for degree in degrees]
    print(f"  Ritz degrees {degrees}", flush=True)
    return {
        "case_id": case, "ratio_ly_lx": ratio, "edges_x0_y0_xlx_yly": edges,
        "pynite_meshes": fem, "ritz_models": ritz, "comparisons": compare(ratio, source, fem, ritz),
    }


def csv_rows(runs):
    for run in runs:
        for quantity, row in run["comparisons"].items():
            hahn = row["hahn"] or {}
            source = hahn.get("source", {})
            units = (
                "m" if quantity == "deflection_ref" else "kN" if quantity.startswith("reaction_")
                else "fraction" if quantity == "negative_moment_reach" else "kN m/m"
            )
            yield {
                "case_id": run["case_id"], "ratio_ly_lx": run["ratio_ly_lx"], "quantity": quantity, "units": units,
                "source_page": source.get("source_page"), "source_symbol": source.get("original_symbol"),
                "source_value": hahn.get("value"), "source_coefficient": source.get("table_value"),
                **{
                    f"{method}_{field}": row[method][field]
                    for method in ("pynite", "ritz") for field in ("value", "coefficient", "x_m", "y_m")
                },
                **{
                    f"{comparison}_{field}_percent": row.get(comparison, {}).get(field)
                    for comparison in ("pynite_vs_hahn", "ritz_vs_hahn", "pynite_vs_ritz",
                                       "pynite_last_refinement", "ritz_last_refinement")
                    for field in ("value", "coefficient")
                },
            }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    parser.add_argument("--ratios", type=float, nargs="+", default=RATIOS)
    parser.add_argument("--meshes", type=int, nargs="+", default=[7, 11, 15, 21, 31])
    parser.add_argument("--degrees", type=int, nargs="+", default=[8, 14, 22])
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "three_edge_benchmark.json")
    args = parser.parse_args()
    if len(args.meshes) < 2 or args.meshes != sorted(set(args.meshes)) or any(n < 3 or n % 2 == 0 for n in args.meshes):
        parser.error("--meshes needs at least two increasing odd counts of at least 3")
    if len(args.degrees) < 2 or args.degrees != sorted(set(args.degrees)) or args.degrees[0] < 5:
        parser.error("--degrees needs at least two increasing degrees of at least 5")
    if args.output.suffix.lower() != ".json":
        parser.error("--output must have a .json extension")
    table = load_table(THREE_EDGE_CSV)
    for case in args.cases:
        for ratio in args.ratios:
            try:
                hahn_reference(case, ratio, table)
            except ValueError as error:
                parser.error(str(error))

    runtime = {"python": platform.python_version(), "platform": platform.platform(),
               **{package: version(package) for package in ("PyNiteFEA", "numpy", "scipy")}}
    print(runtime, flush=True)
    runs = [run_case(case, ratio, table, args.meshes, args.degrees) for case in args.cases for ratio in args.ratios]
    report = {
        "runtime": runtime,
        "argv": ["python", "-m", "scripts.run_three_edge_benchmark", *sys.argv[1:]],
        "solver_models": {"pynite": "PyNiteFEA Rect linear static", "ritz": "Kirchhoff–Love Rayleigh–Ritz"},
        "inputs": {"lx_m": LX, "load_kN_per_m2": LOAD, "E_kN_per_m2": E, "h_m": THICKNESS, "nu": NU,
                   "pynite_nx": args.meshes, "ritz_degrees_per_axis": args.degrees, "ritz_station_grid": STATION_GRID},
        "screening_thresholds_percent": LIMITS,
        "conventions": {
            "moments": "Sagging positive, kN m/m. For Hahn, M = p*lx*ly/k. Corner twisting comparisons use magnitudes.",
            "deflection": "Positive down, m. Hahn free-edge midspan w = p*lx*ly*lx^2/(omega*E*h^3).",
            "reactions": "Global Z up, kN. PyNite open-edge sums exclude corner nodes. Ritz H17 integrates effective shear.",
            "corners": "PyNite nodal corner reactions and Ritz signed twisting jumps are stored separately.",
        },
        "results": runs,
    }
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    rows = list(csv_rows(runs))
    with args.output.with_suffix(".csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(runs)} case/ratio runs and {len(rows)} comparison rows to {args.output}")


if __name__ == "__main__":
    main()
