"""Independent finite-difference Kirchhoff plate reference for four supported edges."""

from dataclasses import dataclass

import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import spsolve


@dataclass
class PlateGrid:
    x: np.ndarray
    y: np.ndarray
    downward_deflection: np.ndarray
    mx: np.ndarray
    my: np.ndarray
    rigidity: float


def _second_difference(n, step):
    return diags((1.0, -2.0, 1.0), (-1, 0, 1), shape=(n, n), format="csc") / step**2


def solve_four_edge_fd(edges, lx, ly, pressure, elastic_modulus, thickness, nx, ny):
    """Solve D * Δ²(w) = p with nu=0 and odd/even reflection at S/F edges.

    All physical boundary displacements are zero. The squared interior
    Laplacian uses odd ghost displacement across each edge (w=0, w_nn=0).
    A fixed edge instead has an even ghost (w=0, w_n=0), adding 2/step**4
    to the adjacent interior diagonal of the biharmonic operator.
    """
    if len(edges) != 4 or not set(edges) <= {"S", "F"}:
        raise ValueError("Expected four simply supported or fixed edge codes")
    if min(nx, ny) < 4 or nx % 2 or ny % 2:
        raise ValueError("FD mesh counts must be even and at least 4")

    left, bottom, right, top = (edge == "F" for edge in edges)
    x, y = np.linspace(0, lx, nx + 1), np.linspace(0, ly, ny + 1)
    dx, dy = lx / nx, ly / ny
    ix, iy = nx - 1, ny - 1
    laplacian = (
        kron(eye(iy), _second_difference(ix, dx), format="csc")
        + kron(_second_difference(iy, dy), eye(ix), format="csc")
    )
    ghost_correction = np.zeros((iy, ix))
    ghost_correction[:, 0] += left * 2 / dx**4
    ghost_correction[:, -1] += right * 2 / dx**4
    ghost_correction[0, :] += bottom * 2 / dy**4
    ghost_correction[-1, :] += top * 2 / dy**4
    biharmonic = laplacian @ laplacian + diags(ghost_correction.ravel(), format="csc")
    rigidity = elastic_modulus * thickness**3 / 12
    interior = spsolve(biharmonic.tocsc(), np.full(ix * iy, pressure / rigidity))
    w = np.pad(interior.reshape(iy, ix), 1)
    mx, my = np.zeros_like(w), np.zeros_like(w)
    mx[:, 1:-1] = -rigidity * (w[:, :-2] - 2 * w[:, 1:-1] + w[:, 2:]) / dx**2
    my[1:-1, :] = -rigidity * (w[:-2, :] - 2 * w[1:-1, :] + w[2:, :]) / dy**2
    # Use the same even ghost when evaluating the normal curvature at a clamp.
    mx[:, 0] = -rigidity * 2 * left * w[:, 1] / dx**2
    mx[:, -1] = -rigidity * 2 * right * w[:, -2] / dx**2
    my[0, :] = -rigidity * 2 * bottom * w[1, :] / dy**2
    my[-1, :] = -rigidity * 2 * top * w[-2, :] / dy**2
    return PlateGrid(x, y, w, mx, my, rigidity)
