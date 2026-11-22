"""Tensor-product Rayleigh–Ritz reference for a uniformly loaded Kirchhoff plate."""

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import legder, leggauss, legvander
from scipy.linalg import null_space, solve


def _legendre_derivatives(t, degree, derivative):
    if derivative == 0:
        return legvander(t, degree)
    coefficients = legder(np.eye(degree + 1), m=derivative, axis=0)
    return legvander(t, degree - derivative) @ coefficients


@dataclass
class SpanBasis:
    length: float
    degree: int
    transform: np.ndarray

    def evaluate(self, stations, derivative=0):
        t = 2 * np.atleast_1d(stations) / self.length - 1
        return (2 / self.length)**derivative * (
            _legendre_derivatives(t, self.degree, derivative) @ self.transform
        )


def _span_basis(length, degree, low, high):
    constraints = []
    for endpoint, boundary in ((-1.0, low), (1.0, high)):
        if boundary in ("S", "F"):
            constraints.append(_legendre_derivatives([endpoint], degree, 0)[0])
        if boundary == "F":
            constraints.append(_legendre_derivatives([endpoint], degree, 1)[0])
    transform = null_space(np.vstack(constraints)) if constraints else np.eye(degree + 1)
    return SpanBasis(length, degree, transform)


@dataclass
class RitzPlate:
    x: SpanBasis
    y: SpanBasis
    amplitudes: np.ndarray
    rigidity: float
    poisson: float

    def partial(self, x, y, dx=0, dy=0):
        """Evaluate a displacement derivative on the tensor grid, indexed [y, x]."""
        return (self.x.evaluate(x, dx) @ self.amplitudes @ self.y.evaluate(y, dy).T).T

    def moments(self, x, y):
        xx = self.partial(x, y, dx=2)
        yy = self.partial(x, y, dy=2)
        xy = self.partial(x, y, dx=1, dy=1)
        return (
            -self.rigidity * (xx + self.poisson * yy),
            -self.rigidity * (yy + self.poisson * xx),
            -self.rigidity * (1 - self.poisson) * xy,
        )

    def effective_shear(self, x, y, normal_axis):
        """Effective Kirchhoff shear in the positive axis direction, in kN/m."""
        transverse = 2 - self.poisson
        if normal_axis == "x":
            return self.rigidity * (self.partial(x, y, dx=3) + transverse * self.partial(x, y, dx=1, dy=2))
        return self.rigidity * (self.partial(x, y, dy=3) + transverse * self.partial(x, y, dx=2, dy=1))


def solve_ritz(edges, lx, ly, pressure, elastic_modulus, thickness, poisson, degree):
    """Minimize Kirchhoff bending energy over polynomials satisfying essential restraints.

    S imposes w=0. F also imposes the normal derivative w_n=0.
    The free edge has no essential constraint, so its moment and effective
    shear conditions are natural conditions of the variational problem.
    """
    if len(edges) != 4 or not set(edges) <= {"S", "F", "FREE"}:
        raise ValueError("Expected S, F, or FREE edge restraints")
    if degree < 5:
        raise ValueError("Use a polynomial degree of at least five")
    left, bottom, right, top = edges
    x = _span_basis(lx, degree, left, right)
    y = _span_basis(ly, degree, bottom, top)
    points, weights = leggauss(degree + 1)
    wx, wy = weights * lx / 2, weights * ly / 2
    xs = [x.evaluate((points + 1) * lx / 2, order) for order in range(3)]
    ys = [y.evaluate((points + 1) * ly / 2, order) for order in range(3)]

    def tensor_form(dx_test, dy_test, dx_trial, dy_trial):
        gram_x = xs[dx_test].T @ (wx[:, None] * xs[dx_trial])
        gram_y = ys[dy_test].T @ (wy[:, None] * ys[dy_trial])
        return np.kron(gram_x, gram_y)

    cross = tensor_form(2, 0, 0, 2)
    stiffness = (
        tensor_form(2, 0, 2, 0) + tensor_form(0, 2, 0, 2)
        + poisson * (cross + cross.T) + 2 * (1 - poisson) * tensor_form(1, 1, 1, 1)
    )
    rigidity = elastic_modulus * thickness**3 / (12 * (1 - poisson**2))
    rhs = pressure / rigidity * np.kron(wx @ xs[0], wy @ ys[0])
    coefficients = solve(stiffness, rhs, assume_a="pos").reshape(xs[0].shape[1], ys[0].shape[1])
    return RitzPlate(x, y, coefficients, rigidity, poisson)
