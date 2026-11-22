"""Off-screen PyNite renderer, run in a separate process by Streamlit."""

import argparse
from pathlib import Path

from scripts.run_pynite_ssss import COMBO, LX, solve_rect_model
from slabber_app.field_colors import FIELD_COLORS, field_extent


def render_plate(case: str, ratio: float, nx: int, poisson: float, result: str, output: Path) -> float:
    import pyvista
    from matplotlib.colors import LinearSegmentedColormap
    from Pynite.Rendering import Renderer

    pyvista.OFF_SCREEN = True
    model, _ = solve_rect_model(case, ratio, nx, poisson=poisson)
    deformed = result == "deflection"
    scale = 1.0
    if deformed:
        max_displacement = max(abs(node.DZ[COMBO]) for node in model.nodes.values())
        if max_displacement > 0:
            scale = 0.1 * min(LX, LX * ratio) / max_displacement

    renderer = Renderer(model)
    renderer.combo_name = COMBO
    renderer.theme = "print"
    renderer.window_width = 1024
    renderer.window_height = 768
    renderer.render_nodes = False
    renderer.render_loads = False
    renderer.show_labels = False
    renderer.deformed_shape = deformed
    renderer.deformed_scale = scale
    renderer.color_map = "dz" if deformed else result
    renderer.scalar_bar = True

    try:
        renderer.update(reset_camera=False, off_screen=True)
        plotter = renderer.plotter
        contours = plotter.mesh.point_data["Contours"]
        if deformed:
            # PyNite provides local z displacement in metres, positive upward.
            contours *= -1000
            plotter.view_isometric()
            title = "Vertical deflection (mm)"
        else:
            plotter.view_xy()
            title = f"{result} (kN m/m)"
        extent = field_extent(float(contours.min()), float(contours.max()))
        cmap = LinearSegmentedColormap.from_list("slabber", FIELD_COLORS, N=257)
        mapper = plotter.mapper
        mapper.lookup_table.apply_cmap(cmap, n_values=cmap.N)
        plotter.update_scalar_bar_range((-extent, extent))
        plotter.scalar_bar.SetTitle(title)
        plotter.scalar_bar.SetPosition(0.2, 0.025)
        for actor in plotter.actors.values():
            if isinstance(actor, pyvista.Actor) and actor.mapper is mapper:
                actor.prop.lighting = False
        plotter.enable_parallel_projection()
        plotter.reset_camera()
        plotter.camera.zoom(0.86)  # Leave room for the horizontal legend below the slab.
        plotter.screenshot(str(output))
    finally:
        renderer.plotter.close()
    return scale


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case")
    parser.add_argument("ratio", type=float)
    parser.add_argument("nx", type=int)
    parser.add_argument("poisson", type=float)
    parser.add_argument("result", choices=["Mx", "My", "Mxy", "deflection"])
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(render_plate(args.case, args.ratio, args.nx, args.poisson, args.result, args.output))


if __name__ == "__main__":
    main()
