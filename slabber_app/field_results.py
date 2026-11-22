"""Finite-element field maps and optional PyNite-rendered figures."""

from importlib.util import find_spec
import logging
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import altair as alt
import pandas as pd
import streamlit as st
from Pynite import FEModel3D

from czerny_dataset.pynite_czerny_benchmark import edge_codes
from scripts.run_pynite_ssss import COMBO, LX, ROOT, center_displacement
from slabber_app.charts import EDGE_COLORS, EDGE_LABELS


RESULT_LABELS = {"Mx": "Mₓ", "My": "Mᵧ", "Mxy": "Mₓᵧ", "deflection": "Vertical deflection"}
LOGGER = logging.getLogger(__name__)


def element_center_values(model: FEModel3D) -> list[dict[str, str | float]]:
    cells = []
    for plate in model.plates.values():
        width, height = plate.width(), plate.height()
        origin = plate.i_node
        mx, my, mxy = plate.moment(width / 2, height / 2, combo_name=COMBO).ravel()
        cells.append({
            "element": plate.name,
            "x0": float(origin.X),
            "x1": float(origin.X + width),
            "y0": float(origin.Y),
            "y1": float(origin.Y + height),
            "x": float(origin.X + width / 2),
            "y": float(origin.Y + height / 2),
            "Mx": float(mx),
            "My": float(my),
            "Mxy": float(mxy),
            "deflection": -1000 * center_displacement(plate),
        })
    return cells


def field_map(cells, result, ly, supports):
    label = RESULT_LABELS[result]
    units = "mm" if result == "deflection" else "kN·m/m"
    values = [cell[result] for cell in cells]
    low, high = min(values), max(values)
    extent = max(abs(low), abs(high), 1e-12)
    if result == "deflection" and low >= 0:
        colors = alt.Scale(domain=[0, extent], range=["#F0F7F7", "#086861"])
    else:
        colors = alt.Scale(domain=[-extent, 0, extent], range=["#335E9F", "#F4F5F6", "#B15324"])

    x = alt.X("x0:Q", title="x (m)", scale=alt.Scale(domain=[0, LX], nice=False),
              axis=alt.Axis(grid=False, tickCount=5))
    y = alt.Y("y0:Q", title="y (m)", scale=alt.Scale(domain=[0, ly], nice=False),
              axis=alt.Axis(grid=False, tickCount=5))
    rows = pd.DataFrame(cells, columns=["element", "x0", "x1", "y0", "y1", "x", "y", result])
    heatmap = alt.Chart(rows, mark=alt.MarkDef(type="rect")).encode(
        x=x, x2="x1:Q", y=y, y2="y1:Q",
        color=alt.Color(f"{result}:Q", title=f"{label} ({units})", scale=colors,
                        legend=alt.Legend(orient="bottom", gradientLength=160, format=".3~g")),
        tooltip=[
            alt.Tooltip("element:N", title="Element"),
            alt.Tooltip("x:Q", title="Center x (m)", format=".3f"),
            alt.Tooltip("y:Q", title="Center y (m)", format=".3f"),
            alt.Tooltip(f"{result}:Q", title=f"{label} ({units})", format=".5g"),
        ],
    )
    edge_coordinates = (
        (0, 0, 0, ly, "x = 0"), (0, 0, LX, 0, "y = 0"),
        (LX, 0, LX, ly, "x = lₓ"), (0, ly, LX, ly, "y = lᵧ"),
    )
    edges = []
    for support, (x0, y0, x1, y1, name) in zip(supports, edge_coordinates, strict=True):
        edge = alt.Chart(alt.Data(values=[
            {"x0": x0, "y0": y0, "x1": x1, "y1": y1, "edge": name, "support": EDGE_LABELS[support]},
        ]), mark=alt.MarkDef(type="rule", color=EDGE_COLORS[support], strokeWidth=2.4,
                            strokeDash=[5, 4] if support == "FREE" else [1, 0])).encode(
            x=x, x2="x1:Q", y=y, y2="y1:Q",
            tooltip=[alt.Tooltip("edge:N", title="Edge"), alt.Tooltip("support:N", title="Support")],
        )
        edges.append(edge)
    pixels_per_meter = 280 / max(LX, ly)
    # Axes and the color legend sit outside the physical plotting area.
    return (alt.layer(heatmap, *edges).properties(
        width=round(LX * pixels_per_meter),
        height=round(ly * pixels_per_meter),
        autosize=alt.AutoSizeParams(type="pad", contains="content"),
    ).configure_view(stroke=None))


@st.cache_data(show_spinner="Rendering the PyNite figure…", max_entries=48)
def native_figure(case: str, ratio: float, nx: int, poisson: float, result: str) -> tuple[bytes, float]:
    with TemporaryDirectory(prefix="slabber-render-", dir=ROOT / "results") as folder:
        image_path = Path(folder) / "figure.png"
        completed = subprocess.run(
            [sys.executable, "-m", "slabber_app.render_plate", case, str(ratio), str(nx),
             str(poisson), result, str(image_path)],
            cwd=ROOT, capture_output=True, check=True, text=True, encoding="utf-8", errors="replace",
            timeout=60, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        return image_path.read_bytes(), float(completed.stdout.strip().splitlines()[-1])


def show_field_results(case, ratio, nx, ny, poisson, cells, key):
    with st.container(border=True):
        st.subheader("PyNite field results", anchor=False)
        result = st.selectbox("Result", RESULT_LABELS.keys(), format_func=RESULT_LABELS.__getitem__,
                              key=f"{key}_field", required=True)
        maps, native = st.tabs(["Interactive maps", "PyNite render"], key=f"{key}_field_view", on_change="rerun")
        label = RESULT_LABELS[result]

        if maps.open:
            with maps:
                chart, summary = st.columns([2, 1], gap="medium")
                chart.altair_chart(field_map(cells, result, LX * ratio, edge_codes(case)), width="content",
                                   alt=f"Full-slab {label} map for {case} with one value per element center")
                units = "mm" if result == "deflection" else "kN·m/m"
                values = [cell[result] for cell in cells]
                summary.caption(f"PyNite Rect · {nx} × {ny} mesh · Element centers")
                summary.metric("Minimum", f"{min(values):.4g} {units}")
                summary.metric("Maximum", f"{max(values):.4g} {units}")
                if result == "deflection":
                    summary.caption("Positive values are downward.")
                elif result != "Mxy":
                    summary.caption("Positive moments denote sagging.")

        if native.open:
            with native:
                if find_spec("pyvista") is None:
                    st.info("Launch Slabber with the rendering extra to enable PyNite figures.")
                    st.code("uv run --locked --extra rendering streamlit run app.py", language="bash")
                    return
                try:
                    png, scale = native_figure(case, ratio, nx, poisson, result)
                except (subprocess.SubprocessError, OSError, ValueError, IndexError) as error:
                    detail = error.stderr if isinstance(error, subprocess.CalledProcessError) else str(error)
                    LOGGER.warning("PyNite render failed for %s at %s (%s): %s", case, ratio, result, (detail or "")[-2500:])
                    st.warning("PyNite could not render this figure. The interactive map is available.")
                    return
                caption = (f"PyNite deformed shape · Displacements exaggerated {scale:.1f}×"
                           if result == "deflection" else f"PyNite {label} contour · Nodally averaged · kN·m/m")
                st.image(png, caption=caption, width=740, output_format="PNG", alt=f"{caption} for {case}")
