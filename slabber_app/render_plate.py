"""Off-screen PyNite renderer, run in a separate process by Streamlit."""

import argparse
from pathlib import Path

from scripts.run_pynite_ssss import COMBO, LX, solve_rect_model


def render_plate(case: str, ratio: float, nx: int, poisson: float, result: str, output: Path) -> float:
    import pyvista
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
    renderer.color_map = None if deformed else result
    renderer.scalar_bar = not deformed
    renderer.scalar_bar_text_size = 16

    def set_view(plotter: pyvista.Plotter) -> None:
        if deformed:
            plotter.view_isometric()
        else:
            plotter.view_xy()
            plotter.scalar_bar.SetTitle(f"{result} (kN m/m)")
        plotter.enable_parallel_projection()

    renderer.post_update_callbacks.append(set_view)
    try:
        renderer.screenshot(str(output), interact=False, reset_camera=True)
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
