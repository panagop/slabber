"""Shared Streamlit page for the two slab benchmark families."""

from dataclasses import dataclass
import json
from pathlib import Path

import streamlit as st

from czerny_dataset.czerny import THREE_EDGE_CSV, load_table
from czerny_dataset.pynite_czerny_benchmark import edge_codes, nearest_odd
from scripts.run_four_edge_benchmark import (
    FD_REFINEMENT, FEM_REFINEMENT, FOUR_EDGE_RATIOS, comparison_run, czerny_samples, fd_mesh,
    rect_mesh, rect_mesh_results,
)
from scripts.run_pynite_ssss import E, LOAD, LX, NU, ROOT, THICKNESS, solve_rect_model
from slabber_app.charts import EDGE_LABELS, MX_COLOR, MY_COLOR, moment_diagram, section_profiles, support_plan
from slabber_app.field_results import element_center_values, show_field_results


@dataclass(frozen=True)
class SlabFamily:
    slug: str
    title: str
    source: str
    citation: str
    reference: str
    reference_key: str
    reference_error_key: str
    page_offset: int
    cases: dict[str, str]
    default_case: str
    ratios: tuple[float, ...]


FOUR_EDGE = SlabFamily(
    "four_edge", "Four-edge slabs", "Czerny", "Beton-Kalender, 1996",
    "Finite difference", "finite_difference", "fd", 276,
    {
        "SSSS": "All edges simply supported",
        "F_y0": "Fixed at y = 0",
        "F_x0": "Fixed at x = 0",
        "FF_y": "Both y edges fixed",
        "FF_x": "Both x edges fixed",
        "FF_adj": "Two adjacent edges fixed",
        "FFF_x": "Three fixed edges, x = lₓ simply supported",
        "FFF_y": "Three fixed edges, y = lᵧ simply supported",
        "FFFF": "All edges fixed",
    },
    "FFFF", FOUR_EDGE_RATIOS,
)
THREE_EDGE = SlabFamily(
    "three_edge", "Three-edge slabs", "Hahn", "Tafeln 17–22, 1970",
    "Rayleigh–Ritz", "ritz", "ritz", 360,
    {
        "H17_SSS": "Three simply supported edges",
        "H18_REAR_FIXED": "Fixed rear edge",
        "H19_ONE_SIDE_FIXED": "One fixed side",
        "H20_BOTH_SIDES_FIXED": "Both sides fixed",
        "H21_ADJ_FIXED": "Fixed side and rear edge",
        "H22_ALL_FIXED": "Three fixed edges",
    },
    "H17_SSS", (0.25, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0),
)
QUANTITY_LABELS = {
    "mx_mid": "Mₓ at center",
    "mx_center": "Mₓ at center",
    "my_mid": "Mᵧ at center",
    "my_center": "Mᵧ at center",
    "my_max": "Peak Mᵧ on centerline",
    "mx_global_max": "Mₓ peak across slab",
    "my_global_max": "Mᵧ peak across slab",
    "mx_fixed_mid": "Mₓ at fixed edge midpoint",
    "my_fixed_mid": "Mᵧ at fixed edge midpoint",
    "mx_fixed_min": "Minimum Mₓ on fixed edge",
    "my_fixed_min": "Minimum Mᵧ on fixed edge",
    "mx_free_mid": "Mₓ at free edge midpoint",
    "mx_fixed_free_corner": "Mₓ at fixed/free corner",
    "mx_fixed_side_mid": "Mₓ at fixed side midpoint",
    "my_positive_adjacent": "Positive Mᵧ peak on centerline",
    "my_rear_fixed": "Mᵧ at fixed rear midpoint",
    "mxy_corner": "Corner twisting moment magnitude",
    "mxy_free_corner": "Free corner twisting moment magnitude",
    "mxy_rear_corner": "Rear corner twisting moment magnitude",
    "corner_force": "Corner force magnitude",
    "reaction_x_total": "Resultant on one side, excluding corners",
    "reaction_y_total": "Rear resultant, excluding corners",
    "deflection_mid": "Center deflection",
    "deflection_ref": "Free edge midpoint deflection",
    "negative_moment_reach": "Negative moment reach β",
}


@st.cache_data(show_spinner=False)
def read_results(path: Path, modified_ns: int):
    # The modification time invalidates the cache when the benchmark is rerun.
    return json.loads(path.read_text(encoding="utf-8"))


@st.cache_data(show_spinner="Computing PyNite results…", max_entries=48)
def plot_data(case, ratio, nx, poisson, *, compare_four_edge=False):
    model, ny = solve_rect_model(case, ratio, nx, poisson=poisson)
    run = None
    if compare_four_edge:
        edges = edge_codes(case)
        source = czerny_samples(case, ratio, load_table())
        fem = [rect_mesh(case, ratio, FEM_REFINEMENT[0], edges), rect_mesh_results(model, ratio, nx, ny, edges)]
        fd = [fd_mesh(ratio, divisions, edges) for divisions in FD_REFINEMENT]
        run = comparison_run(case, ratio, edges, source, fem, fd)
    return section_profiles(model, nx, ny), element_center_values(model), run


def quantity_label(quantity, case, ratio):
    if quantity == "my_center_or_max":
        uses_peak = (case.startswith("H19") and ratio > 0.8) or (case.startswith("H20") and ratio > 0.6)
        return "Peak Mᵧ near rear edge" if uses_peak else "Mᵧ at center"
    return QUANTITY_LABELS.get(quantity, quantity)


def display_units(quantity):
    if quantity.startswith("deflection_"):
        return "mm", 1000
    if quantity == "corner_force" or quantity.startswith("reaction_"):
        return "kN", 1
    if quantity == "negative_moment_reach":
        return "fraction", 1
    return "kN·m/m", 1


def show_metrics(run, family):
    finest = run["pynite_meshes"][-1]
    points = finest["points"]
    mx_key, my_key, w_key = (
        ("mx_mid", "my_mid", "deflection_mid") if family.slug == "four_edge"
        else ("mx_center", "my_center", "deflection_center_diagnostic")
    )
    reaction = finest["reactions"]
    for column, (label, value, help_text) in zip(st.columns(4), (
        ("Mₓ at center", f"{points[mx_key]['value']:.3f} kN·m/m", "Positive denotes sagging."),
        ("Mᵧ at center", f"{points[my_key]['value']:.3f} kN·m/m", "Positive denotes sagging."),
        ("Center deflection", f"{1000 * points[w_key]['value']:.3f} mm", "Positive downward."),
        ("Vertical reaction sum", f"{reaction['support_fz_kN']:.2f} kN",
         f"Signed PyNite nodal reaction sum. Equilibrium residual {reaction['residual_percent']:.2e}%."),
    ), strict=True):
        column.metric(label, value, help=help_text, border=True)


def show_figures(case, ratio, nx, ny, poisson, key, sections, cells):
    lx, ly = LX, LX * ratio
    edges = edge_codes(case)
    plan, diagrams = st.columns([1, 2.25], gap="medium")
    with plan, st.container(border=True):
        st.subheader("Support layout", anchor=False)
        st.altair_chart(support_plan(lx, ly, edges), width="content",
                        alt=f"Support plan for {case}, with Mx and My section lines")
        for edge_group in ((0, 2), (1, 3)):
            boundaries = ("x = 0", "y = 0", "x = lₓ", "y = lᵧ")
            st.caption(" · ".join(f"{boundaries[i]}  {EDGE_LABELS[edges[i]]}" for i in edge_group))
        st.caption("Hatching marks fixed edges. A dashed perimeter marks a free edge.")

    with diagrams, st.container(border=True):
        st.subheader("Centerline moment diagrams", anchor=False)
        st.caption(f"PyNite Rect  ·  {nx} × {ny} mesh  ·  Positive denotes sagging")
        x_points, y_points = sections
        values = [point["moment"] for point in (*x_points, *y_points)]
        low, high = min(*values, 0), max(*values, 0)
        pad = max(0.08 * (high - low), 0.01)
        domain = [low - pad, high + pad]
        st.markdown("**Mₓ along x** · at y = lᵧ/2")
        st.altair_chart(moment_diagram(x_points, "x", lx, domain, MX_COLOR), width="stretch",
                        alt="Signed Mx moment along the horizontal slab centerline")
        st.markdown("**Mᵧ along y** · at x = lₓ/2")
        st.altair_chart(moment_diagram(y_points, "y", ly, domain, MY_COLOR), width="stretch",
                        alt="Signed My moment along the vertical slab centerline")
        st.caption("The plots retain each element’s moment values at shared boundaries.")
    show_field_results(case, ratio, nx, ny, poisson, cells, key)


def show_comparison(run, family):
    st.subheader("Benchmark comparisons", anchor=False)
    tab_values, tab_coefficients, tab_refinement = st.tabs(["Physical values", "Coefficients", "Mesh refinement"])
    source_key = family.source.lower()
    value_key = "physical_percent" if source_key == "czerny" else "value"
    coefficient_key = "coefficient_percent" if source_key == "czerny" else "coefficient"
    values, coefficients, refinement = [], [], []
    source_metadata = None
    matching_rows = source_rows = 0

    for quantity, record in run["comparisons"].items():
        label = quantity_label(quantity, run["case_id"], run["ratio_ly_lx"])
        units, scale = display_units(quantity)
        source = record[source_key]
        numerical = record[family.reference_key]
        if source:
            source_metadata = source["source"]
            source_rows += 1
            errors = (record[f"pynite_vs_{source_key}"], record[f"{family.reference_error_key}_vs_{source_key}"])
            if all(abs(error[value_key]) <= 2 for error in errors):
                matching_rows += 1
        values.append({
            "Quantity": label, "Unit": units, family.source: source["value"] * scale if source else None,
            "PyNite": record["pynite"]["value"] * scale,
            family.reference: numerical["value"] * scale,
            "PyNite vs source (%)": record[f"pynite_vs_{source_key}"][value_key] if source else None,
            "PyNite vs reference (%)": record[f"pynite_vs_{family.reference_error_key}"][value_key],
        })
        if source:
            coefficients.append({
                "Quantity": label, "Symbol": source["source"].get("printed_symbol_from_scan",
                                                                 source["source"]["original_symbol"]),
                family.source: source["source"]["table_value"], "PyNite": record["pynite"]["coefficient"],
                family.reference: numerical["coefficient"],
                "PyNite vs source (%)": record[f"pynite_vs_{source_key}"][coefficient_key],
            })
        refinement.append({
            "Quantity": label,
            "PyNite change (%)": record["pynite_last_refinement"][value_key],
            f"{family.reference} change (%)": record[f"{family.reference_error_key}_last_refinement"][value_key],
        })

    with tab_values:
        if source_rows:
            st.caption(f"Both solvers are within 2% of {family.source} on {matching_rows} of {source_rows} physical values.")
        show_table(values, "Physical benchmark values", 3)
    with tab_coefficients:
        basis = r"p,l_x^2" if source_key == "czerny" else r"K=p,l_xl_y"
        st.caption(f"Moment denominator scale ${basis}$. Deflection and other coefficients follow the source definitions.")
        show_table(coefficients, "Source and numerical coefficients", 3)
    with tab_refinement:
        meshes = run["pynite_meshes"]
        previous, finest = meshes[-2:]
        mesh_pair = f"{previous['nx']} × {previous['ny']} → {finest['nx']} × {finest['ny']}"
        if source_key == "czerny":
            old_ref, new_ref = run["fd_meshes"][-2:]
            steps = f"{old_ref['nx']} × {old_ref['ny']} → {new_ref['nx']} × {new_ref['ny']}"
        else:
            old_ref, new_ref = run["ritz_models"][-2:]
            steps = f"degree {old_ref['degree']} → {new_ref['degree']}"
        st.caption(f"Change in physical values on the final refinement  ·  PyNite {mesh_pair}  ·  {family.reference} {steps}")
        show_table(refinement, "Final refinement changes", 3)
    if source_metadata:
        printed_page = source_metadata["source_page"]
        page_url = f"{source_metadata['source_url']}#page={printed_page - family.page_offset}"
        st.caption(f"Historical data from [{family.source}, printed p. {printed_page}]({page_url}).")


def show_table(rows, label, precision):
    numeric = {
        column: st.column_config.NumberColumn(format=f"%.{precision}f")
        for column in rows[0]
        if any(isinstance(row[column], (int, float)) for row in rows)
    }
    st.dataframe(rows, width="stretch", height="content", hide_index=True, column_config=numeric, alt=label)


def render_dashboard(family):
    st.title(family.title, anchor=False)
    st.caption(f"{family.source} · {family.citation}  /  {family.reference} reference")

    path = ROOT / "results" / f"{family.slug}_benchmark.json"
    try:
        report = read_results(path, path.stat().st_mtime_ns)
    except (OSError, json.JSONDecodeError):
        report = None
    selector, aspect = st.columns([2.25, 1], vertical_alignment="bottom")
    cases = list(family.cases)
    case = selector.selectbox("Slab case", cases, index=cases.index(family.default_case),
                              format_func=lambda item: f"{item} · {family.cases[item]}",
                              key=f"{family.slug}_case")
    case_runs = {run["ratio_ly_lx"]: run for run in report["results"] if run["case_id"] == case} if report else {}
    four_edge = family.slug == "four_edge"
    ratios = family.ratios
    ratio = aspect.selectbox("Aspect ratio  lᵧ/lₓ", ratios, index=ratios.index(1.0),
                             format_func=lambda item: f"{item:.2f}", key=f"{family.slug}_ratio")
    run = case_runs.get(ratio)
    inputs = report["inputs"] if report and run else {}
    poisson = inputs.get("nu", NU)
    st.caption(f"Uniform pressure {inputs.get('load_kN_per_m2', LOAD):g} kN/m²  ·  "
               f"Thickness {inputs.get('h_m', THICKNESS):g} m  ·  E = {inputs.get('E_kN_per_m2', E) / 1e6:g} GPa  ·  ν = {poisson:g}")
    live_comparison = four_edge and run is None
    if run:
        finest = run["pynite_meshes"][-1]
        nx, ny = finest["nx"], finest["ny"]
    else:
        nx = FEM_REFINEMENT[-1]
        ny = nearest_odd(nx * ratio)

    if live_comparison:
        st.caption("Calculated on demand for this ratio using the final two PyNite and finite-difference meshes.")
    elif run is None:
        hahn_table = load_table(THREE_EDGE_CSV)
        last_ratio = max(rows[-1]["ratio_ly_lx"] for (case_id, _), rows in hahn_table.items() if case_id == case)
        if ratio > last_ratio:
            st.caption(f"PyNite results calculated on demand. Hahn tabulates this case through lᵧ/lₓ = {last_ratio:.2f}.")
        else:
            st.caption("PyNite results calculated on demand. A saved Hahn and Rayleigh–Ritz comparison is unavailable for this selection.")
            if report is None:
                st.code("uv run --locked python -m scripts.run_three_edge_benchmark", language="bash")

    sections, cells, live_run = plot_data(case, ratio, nx, poisson, compare_four_edge=live_comparison)
    run = run or live_run
    if run:
        show_metrics(run, family)
    show_figures(case, ratio, nx, ny, poisson, family.slug, sections, cells)
    if run:
        show_comparison(run, family)
