"""Support schematics and signed, piecewise-linear PyNite section diagrams."""

import altair as alt

from scripts.run_pynite_ssss import COMBO


MX_COLOR = "#087C74"
MY_COLOR = "#AE571E"
EDGE_LABELS = {"F": "Fixed", "S": "Simply supported", "FREE": "Free"}
EDGE_COLORS = {"F": "#18283C", "S": "#4D7199", "FREE": "#718297"}


def section_profiles(model, nx, ny):
    horizontal, vertical = [], []
    # Rect normal moments are linear along a straight cut within an element.
    # Keep each element's own limiting values at shared boundaries.
    for i in range(nx):
        plate = model.plates[f"P_{i}_{ny // 2}"]
        for x in (0.0, plate.width()):
            mx = plate.moment(x, plate.height() / 2, combo_name=COMBO)[0, 0]
            horizontal.append({"segment": i, "position": plate.i_node.X + x, "moment": float(mx)})
    for j in range(ny):
        plate = model.plates[f"P_{nx // 2}_{j}"]
        for y in (0.0, plate.height()):
            my = plate.moment(plate.width() / 2, y, combo_name=COMBO)[1, 0]
            vertical.append({"segment": j, "position": plate.i_node.Y + y, "moment": float(my)})
    return horizontal, vertical


def support_plan(lx, ly, edges):
    span = max(lx, ly)
    pad = span * 0.17
    half = span / 2 + pad
    x = alt.X("x:Q", scale=alt.Scale(domain=[lx / 2 - half, lx / 2 + half], nice=False), axis=None)
    y = alt.Y("y:Q", scale=alt.Scale(domain=[ly / 2 - half, ly / 2 + half], nice=False), axis=None)

    def lines(items):
        return alt.Chart(alt.Data(values=items)).encode(x=x, x2="x2:Q", y=y, y2="y2:Q")

    fill = lines([{"x": 0, "x2": lx, "y": 0, "y2": ly}]).mark_rect(color="#E9F1F6")
    sides = (
        (0, 0, 0, ly, -1, 0, "x = 0"),
        (0, 0, lx, 0, 0, -1, "y = 0"),
        (lx, 0, lx, ly, 1, 0, "x = lₓ"),
        (0, ly, lx, ly, 0, 1, "y = lᵧ"),
    )
    layers = [fill]
    for code, (x0, y0, x1, y1, out_x, out_y, name) in zip(edges, sides, strict=True):
        if code == "F":
            hatch_depth = 0.026 * span
            hatch = []
            for step in range(1, 12):
                t = step / 12
                px, py = x0 + (x1 - x0) * t, y0 + (y1 - y0) * t
                hatch.append({"x": px, "y": py,
                              "x2": px + (out_x + abs(out_y) * 0.55) * hatch_depth,
                              "y2": py + (out_y + abs(out_x) * 0.55) * hatch_depth})
            layers.append(lines(hatch).mark_rule(color=EDGE_COLORS[code], strokeWidth=1))
        edge = [{"x": x0, "y": y0, "x2": x1, "y2": y1, "edge": name, "support": EDGE_LABELS[code]}]
        layers.append(lines(edge).mark_rule(
            color=EDGE_COLORS[code],
            strokeWidth=3 if code == "F" else 2,
            strokeDash=[5, 4] if code == "FREE" else [],
        ).encode(tooltip=[alt.Tooltip("edge:N", title="Edge"), alt.Tooltip("support:N", title="Support")]))

    cuts = [
        {"x": 0, "x2": lx, "y": ly / 2, "y2": ly / 2, "color": MX_COLOR},
        {"x": lx / 2, "x2": lx / 2, "y": 0, "y2": ly, "color": MY_COLOR},
    ]
    layers.append(lines(cuts).mark_rule(strokeDash=[5, 3], strokeWidth=1.6).encode(
        color=alt.Color("color:N", scale=None, legend=None)
    ))
    labels = [
        {"x": lx * 0.76, "y": ly / 2 + pad * 0.38, "text": "Mₓ", "color": MX_COLOR},
        {"x": lx / 2 + pad * 0.38, "y": ly * 0.78, "text": "Mᵧ", "color": MY_COLOR},
        {"x": lx / 2, "y": -pad * 0.58, "text": f"lₓ {lx:.2f} m  →", "color": "#43546B"},
    ]
    label_style = alt.MarkDef(type="text", fontSize=12, fontWeight=600)
    layers.append(alt.Chart(alt.Data(values=labels), mark=label_style).encode(
        x=x, y=y, text="text:N", color=alt.Color("color:N", scale=None, legend=None)
    ))
    ly_label = alt.Chart(alt.Data(values=[
        {"x": -pad * 0.58, "y": ly / 2, "text": f"lᵧ {ly:.2f} m  →"},
    ]), mark=alt.MarkDef(type="text", angle=270, fontSize=12, fontWeight=600, color="#43546B")).encode(
        x=x, y=y, text="text:N"
    )
    layers.append(ly_label)
    return alt.layer(*layers).properties(width=215, height=215).configure_view(stroke=None)


def moment_diagram(points, axis, length, domain, color):
    position = alt.X("position:Q", title=f"{axis} (m)", scale=alt.Scale(domain=[0, length], nice=False),
                     axis=alt.Axis(tickCount=7, grid=False))
    moment = alt.Y("moment:Q", title="kN·m/m", scale=alt.Scale(domain=domain, nice=False),
                   axis=alt.Axis(tickCount=5, gridOpacity=0.55))
    base = alt.Chart(alt.Data(values=points)).encode(
        x=position, y=moment, detail="segment:N", order=alt.Order("position:Q")
    )
    zero = alt.Chart(alt.Data(values=[{"moment": 0}]), mark=alt.MarkDef(
        type="rule", color="#718297", strokeDash=[4, 3]
    )).encode(y=moment)
    filled = base.properties(mark=alt.MarkDef(type="area", color=color, opacity=0.15)).encode(y2=alt.Y2(datum=0))
    curve = base.properties(mark=alt.MarkDef(type="line", color=color, strokeWidth=2.2))
    stations = base.properties(mark=alt.MarkDef(type="circle", color=color, size=27, opacity=0.66)).encode(
        tooltip=[
            alt.Tooltip("position:Q", title=f"{axis} (m)", format=".3f"),
            alt.Tooltip("moment:Q", title="Moment (kN·m/m)", format=".3f"),
        ]
    )
    return alt.layer(filled, zero, curve, stations).properties(height=135).configure_view(stroke=None)
