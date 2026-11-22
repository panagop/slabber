"""Run the PyNite Rect SSSS pilot against Czerny and the saved SciPy FD reference."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from importlib.metadata import version
from math import fsum
from pathlib import Path

from Pynite import FEModel3D

from czerny_dataset.czerny import load_all_tables, lookup
from czerny_dataset.pynite_czerny_benchmark import edge_codes, nearest_odd


ROOT = Path(__file__).resolve().parents[1]
FD_CSV = ROOT / "czerny_dataset" / "reference_ssss_fd.csv"
DEFAULT_OUTPUT = ROOT / "results" / "pynite_ssss.json"
RATIOS = (1.0, 1.5, 2.0)
LX, LOAD, E, THICKNESS, NU = 4.0, 10.0, 33e6, 0.18, 0.0
COMBO = "P"
LIMITS = {
    "abs_relative_reaction_residual_percent": 1e-6,
    "abs_finest_reference_error_percent": 2.0,
    "abs_last_mesh_change_percent": 1.0,
}


def relative_difference_percent(value, reference):
    return 100 * (value - reference) / reference


def response(value, scale, *, denominator):
    coefficient = scale / value if denominator else value / scale
    return {"value": value, "coefficient": coefficient}


def reaction_balance(applied_fz, support_reactions):
    """Positive global Z is up. Count every signed reaction once."""
    reaction_fz = fsum(support_reactions)
    residual = applied_fz + reaction_fz
    return {
        "applied_fz_kN": applied_fz,
        "support_fz_kN": reaction_fz,
        "residual_fz_kN": residual,
        "residual_percent": 100 * residual / abs(applied_fz),
    }


def center_displacement(plate):
    """Evaluate the Rect displacement polynomial at its center.

    Local bending DOFs are w, RX = dw/dy, and RY = -dw/dx. The center
    formula is exact on the Rect's 12-term basis for these nodal DOFs.
    """
    d = plate.d(COMBO).ravel()
    corner_average = fsum(d[[2, 8, 14, 20]]) / 4
    x_correction = plate.width() * fsum((-d[4], d[10], d[16], -d[22])) / 16
    y_correction = plate.height() * fsum((d[3], d[9], -d[15], -d[21])) / 16
    return corner_average + x_correction + y_correction


def support_reactions(model, nx, ny, ly):
    """Record vertical reactions by open edge and at individual corners."""
    indices = {
        "x=0": [(0, j) for j in range(1, ny)],
        "y=0": [(i, 0) for i in range(1, nx)],
        "x=lx": [(nx, j) for j in range(1, ny)],
        "y=ly": [(i, ny) for i in range(1, nx)],
        "corners": [(0, 0), (nx, 0), (nx, ny), (0, ny)],
    }
    groups = {}
    for location, positions in indices.items():
        nodes = [model.nodes[f"N_{i}_{j}"] for i, j in positions]
        groups[location] = [
            {"x_m": node.X, "y_m": node.Y, "fz_kN": float(node.RxnFZ[COMBO])}
            for node in nodes
        ]
    return {
        **reaction_balance(
            -LOAD * LX * ly,
            (row["fz_kN"] for rows in groups.values() for row in rows),
        ),
        "open_edge_resultants_kN": {
            edge: fsum(row["fz_kN"] for row in rows)
            for edge, rows in groups.items() if edge != "corners"
        },
        "support_nodes_by_location": groups,
    }


def moment_extrema(plates, mid_x):
    """Find maxima of the unsmoothed Rect Mx and My fields.

    Each is bilinear inside an element, so a whole-panel peak is attained
    at an element corner. Evaluate each adjacent element's own limiting
    value. On x=lx/2 the field is linear within each crossed element.
    """
    global_samples = []
    section_samples = []
    for plate in plates:
        width, height = plate.width(), plate.height()
        x0, y0 = plate.i_node.X, plate.i_node.Y
        for y in (0.0, height):
            for x in (0.0, width):
                mx, my, _ = map(float, plate.moment(x, y, combo_name=COMBO).ravel())
                global_samples.append({
                    "x_m": x0 + x, "y_m": y0 + y, "element": plate.name, "Mx": mx, "My": my,
                })
            if x0 < mid_x < x0 + width:
                my = float(plate.moment(mid_x - x0, y, combo_name=COMBO)[1, 0])
                section_samples.append({
                    "x_m": mid_x, "y_m": y0 + y, "element": plate.name, "My": my,
                })
    return {
        "Mx_global": max(global_samples, key=lambda sample: sample["Mx"]),
        "My_global": max(global_samples, key=lambda sample: sample["My"]),
        "My_mid_section": max(section_samples, key=lambda sample: sample["My"]),
    }


def solve_rect_model(case_id, ratio, nx, *, check_statics=False, poisson=NU):
    ly, ny = LX * ratio, nearest_odd(nx * ratio)
    left, bottom, right, top = edge_codes(case_id)
    if not {left, right, top} <= {"S", "F"} or bottom not in {"S", "F", "FREE"}:
        raise ValueError(f"Expected three supported edges and an optional free y=0 edge for {case_id}")
    model = FEModel3D()
    model.add_material("Concrete", E, G=E / (2 * (1 + poisson)), nu=poisson, rho=0.0)
    nodes = {}
    for j in range(ny + 1):
        for i in range(nx + 1):
            name = f"N_{i}_{j}"
            nodes[i, j] = model.add_node(name, LX * i / nx, ly * j / ny, 0.0)
            clamped = (
                (i == 0 and left == "F") or (j == 0 and bottom == "F")
                or (i == nx and right == "F") or (j == ny and top == "F")
            )
            # In-plane and drilling DOFs are fixed for this pure-bending model.
            model.def_support(
                name, support_DX=True, support_DY=True, support_RZ=True,
                support_DZ=i in (0, nx) or j == ny or (j == 0 and bottom != "FREE"),
                support_RX=clamped, support_RY=clamped,
            )

    for j in range(ny):
        for i in range(nx):
            name = model.add_plate(
                f"P_{i}_{j}", nodes[i, j], nodes[i + 1, j],
                nodes[i + 1, j + 1], nodes[i, j + 1], THICKNESS, "Concrete",
            )
            model.add_plate_surface_pressure(name, -LOAD, case=COMBO)
    model.add_load_combo(COMBO, {COMBO: 1.0})
    model.analyze_linear(check_statics=check_statics)
    return model, ny


def solve_mesh(ratio, nx):
    ly = LX * ratio
    model, ny = solve_rect_model("SSSS", ratio, nx, check_statics=True)
    center = model.plates[f"P_{nx // 2}_{ny // 2}"]
    x, y = center.width() / 2, center.height() / 2
    mx, my, mxy = map(float, center.moment(x, y, combo_name=COMBO).ravel())
    qx, qy = map(float, center.shear(x, y, combo_name=COMBO).ravel())
    w = center_displacement(center)
    extrema = moment_extrema(model.plates.values(), LX / 2)
    reactions = support_reactions(model, nx, ny, ly)
    return {
        "nx": nx, "ny": ny, "element_count": nx * ny,
        "center_element_local_axes_in_global": center.T()[:3, :3].tolist(),
        "center": {"x_m": LX / 2, "y_m": ly / 2, "Mx": mx, "My": my, "Mxy": mxy, "w_m": w},
        "center_transverse_shear": {"Qx": qx, "Qy": qy},
        "positive_moment_extrema": extrema,
        "reactions": reactions,
        "responses": {
            "Mx_mid": response(mx, LOAD * LX**2, denominator=True),
            "My_mid_section_max": response(extrema["My_mid_section"]["My"], LOAD * LX**2, denominator=True),
            "My_global_max": response(extrema["My_global"]["My"], LOAD * LX**2, denominator=True),
            "downward_deflection_mid": response(-w, LOAD * LX**4 / (E * THICKNESS**3), denominator=False),
        },
    }


def czerny_reference(ratio, table):
    references = {}
    for key, quantity, scale in (
        ("Mx_mid", "mx_mid", LOAD * LX**2),
        ("My_mid_section_max", "my_max", LOAD * LX**2),
        ("downward_deflection_mid", "deflection_mid", LOAD * LX**4 / (E * THICKNESS**3)),
    ):
        factor, source = lookup("SSSS", ratio, quantity, table=table)
        references[key] = {"value": scale * factor, "coefficient": source["table_value"], "source": source}
    return references


def fd_references():
    with FD_CSV.open(newline="", encoding="utf-8") as stream:
        rows = {float(row["ratio_ly_lx"]): row for row in csv.DictReader(stream) if int(row["nx"]) == 64}
    return {
        ratio: {
            "nx": int(row["nx"]), "ny": int(row["ny"]),
            "My_global_max_location_m": [float(row["My_max_x"]), float(row["My_max_y"])],
            "responses": {
                key: {"value": float(row[field]), "coefficient": float(row[coefficient])}
                for key, field, coefficient in (
                    ("Mx_mid", "Mx_mid", "kx"),
                    ("My_global_max", "My_max", "ky"),
                    ("downward_deflection_mid", "w_mid", "deflection_multiplier"),
                )
            },
        }
        for ratio, row in rows.items()
    }


def comparison(actual, reference):
    return {
        key: {
            field: relative_difference_percent(actual[key][field], expected[field])
            for field in ("value", "coefficient")
        }
        for key, expected in reference.items()
    }


def max_abs_error(differences):
    return max(abs(error) for row in differences.values() for error in row.values())


def run_case(ratio, meshes, table, fd):
    source = czerny_reference(ratio, table)
    results = []
    for nx in meshes:
        print(f"SSSS ly/lx={ratio:.2f}, nx={nx}", flush=True)
        mesh = solve_mesh(ratio, nx)
        actual = mesh["responses"]
        mesh["czerny_error_percent"] = comparison(actual, source)
        mesh["fd_error_percent"] = comparison(actual, fd["responses"])
        if results:
            mesh["change_from_previous_percent"] = comparison(actual, results[-1]["responses"])
        results.append(mesh)
        residual = mesh["reactions"]["residual_percent"]
        print(f"  {nx}x{mesh['ny']} | signed reaction residual {residual:+.3e}%")

    finest = results[-1]
    metrics = {
        "abs_relative_reaction_residual_percent": max(abs(row["reactions"]["residual_percent"]) for row in results),
        "abs_finest_reference_error_percent": max(
            max_abs_error(finest["czerny_error_percent"]), max_abs_error(finest["fd_error_percent"]),
        ),
        "abs_last_mesh_change_percent": max_abs_error(finest["change_from_previous_percent"]),
    }
    passed = all(value <= LIMITS[key] for key, value in metrics.items())
    print(f"  max finest error {metrics['abs_finest_reference_error_percent']:.3f}% | "
          f"max last mesh change {metrics['abs_last_mesh_change_percent']:.3f}% | "
          f"pilot criteria {'pass' if passed else 'not met'}")
    return {
        "ratio_ly_lx": ratio, "ly_m": LX * ratio,
        "czerny_reference": source, "fd_reference": fd,
        "meshes": results, "acceptance": {"passes": passed, "metrics_percent": metrics},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meshes", type=int, nargs="+", default=[7, 11, 15, 21, 31])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    meshes = args.meshes
    if len(meshes) < 2 or meshes != sorted(set(meshes)) or any(n < 3 or n % 2 == 0 for n in meshes):
        parser.error("--meshes requires at least two increasing odd counts >= 3")

    runtime = {
        "python": platform.python_version(), "platform": platform.platform(),
        **{package: version(package) for package in ("PyNiteFEA", "numpy", "scipy")},
    }
    print(runtime)
    table, fd = load_all_tables(), fd_references()
    report = {
        "solver": "PyNiteFEA Rect, linear static",
        "runtime": runtime,
        "argv": ["python", "-m", "scripts.run_pynite_ssss", *sys.argv[1:]],
        "case_id": "SSSS",
        "edges_x0_y0_xlx_yly": edge_codes("SSSS"),
        "inputs": {
            "lx_m": LX, "load_kN_per_m2": LOAD, "E_kN_per_m2": E,
            "thickness_m": THICKNESS, "nu": NU, "nx": meshes,
        },
        "conventions": {
            "axes": "Rect local x, y, z coincide with global X, Y, Z",
            "load": "Pressure is negative Z. Deflection w is signed. Compared deflection is positive down.",
            "moment": "Positive sagging, in kN m/m. Coefficients are p*lx^2/M in this SSSS pilot.",
            "deflection": "Compared multiplier is -w*E*h^3/(p*lx^4).",
            "shear": "Qx and Qy are element transverse shears in kN/m. Reactions are nodal forces in kN.",
        },
        "fd_source_file": FD_CSV.relative_to(ROOT).as_posix(),
        "criteria_percent": LIMITS,
        "cases": [run_case(ratio, meshes, table, fd[ratio]) for ratio in RATIOS],
    }
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved {args.output}")
    return 0 if all(case["acceptance"]["passes"] for case in report["cases"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
