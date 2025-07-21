import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def vectorized_version(v, alpha, lambda_val):
    v = int(v)
    alpha = int(alpha) * np.pi / 180
    sx, sy = -0.1, 0
    s = sx + 1j * sy
    r = 1.1
    beta = alpha
    k = -2 * int(r) * int(v) * np.sin(beta)
    tol = 5e-2

    # Mesh or grid generation in the circle or complex z plane
    x, y = np.meshgrid(np.arange(-5, 5.01, 0.01), np.arange(-5, 5.01, 0.01))
    z = x + 1j * y

    w = v * np.exp(1j * alpha)  # free stream rotation to angle of attack, deg

    # Grid tolerance check for flow visualization
    for p in range(len(x)):
        for q in range(len(y)):
            if np.abs(z[p, q] - s) <= r - tol:
                z[p, q] = np.nan

    # Total complex aerodynamic potential function and grid in airfoil plane
    f = np.zeros_like(z, dtype=complex)
    J = np.zeros_like(z, dtype=complex)

    for p in range(len(x)):
        for q in range(len(y)):
            if not np.isnan(z[p, q]):
                f[p, q] = (
                    w * z[p, q]
                    + (v * np.exp(-1j * alpha) * r**2) / (z[p, q] - s)
                    + 1j * k * np.log(z[p, q])
                )
                J[p, q] = z[p, q] + lambda_val**2 / z[p, q]  # grid in the airfoil plane
            else:
                f[p, q] = np.nan
                J[p, q] = np.nan

    # JOUKOWSKI MAPPING FUNCTION
    phi = np.arange(0, 2 * np.pi + 0.05, 0.05)
    zcirc = np.zeros(len(phi), dtype=complex)
    zair = np.zeros(len(phi), dtype=complex)

    for p in range(len(phi)):
        zcirc[p] = r * (np.cos(phi[p]) + 1j * np.sin(phi[p])) + s
        zair[p] = -zcirc[p] + lambda_val**2 / (-zcirc[p])

    return f, J, zcirc, zair
