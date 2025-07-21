import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


def solve_thin_airfoil_theory(x_airfoil, camber_line, v_inf, alpha, camber_slope):
    # Collocation points (quarter-chord rule)
    n_colloc = 150
    c = 1
    alpha = np.radians(alpha)
    theta = np.linspace(
        np.pi / (2 * n_colloc), np.pi - np.pi / (2 * n_colloc), n_colloc
    )
    x_colloc = c / 2 * (1 - np.cos(theta))

    # Interpolate camber slope at collocation points
    camber_slope_interp = np.interp(x_colloc, x_airfoil, camber_slope)

    # Right-hand side: boundary condition
    rhs = alpha - camber_slope_interp

    # Solve for Fourier coefficients
    n_coeffs = 4

    # A_0 coefficient (related to lift)
    A_0 = np.mean(rhs)

    # Higher order coefficients A_n
    A_n = np.zeros(n_coeffs)
    for n in range(1, n_coeffs):
        A_n[n] = 2 * np.mean(rhs * np.cos(n * theta))

    A = np.concatenate([[A_0], A_n])

    # Vortex strength distribution
    gamma = np.zeros_like(x_colloc)
    for i, x in enumerate(x_colloc):
        theta_i = np.arccos(1 - 2 * x / c)
        if np.sin(theta_i) > 1e-10:  # Avoid division by zero
            gamma[i] = (
                2
                * v_inf
                * (
                    A[0] * (1 + np.cos(theta_i)) / np.sin(theta_i)
                    + np.sum([A[n] * np.sin(n * theta_i) for n in range(1, len(A))])
                )
            )

    # Calculate lift coefficient
    cl = 2 * np.pi * (A[0] + A[1] / 2)

    # Calculate moment coefficient about quarter chord
    cm_c4 = -np.pi / 2 * (A[1] - A[2])

    return x_colloc, gamma, cl, cm_c4


def velocity_field(
    x_airfoil, camber_line, v_inf, alpha, camber_slope, x_field, y_field
):
    u = np.zeros_like(x_field)
    v = np.zeros_like(y_field)
    c = 1
    # Free stream contribution (CAMBIA DA STRINGHE A NUMERI)
    u += v_inf * np.cos(int(np.degrees(alpha)))
    v += v_inf * np.sin(int(np.degrees(alpha)))
    x_colloc, gamma, cl, cm_c4 = solve_thin_airfoil_theory(
        x_airfoil, camber_line, v_inf, alpha, camber_slope
    )
    # Vortex contribution
    if gamma is not None:
        # Vortex panel positions
        n_colloc = len(gamma)
        theta = np.linspace(
            np.pi / (2 * n_colloc), np.pi - np.pi / (2 * n_colloc), n_colloc
        )
        x_panels = c / 2 * (1 - np.cos(theta))
        dx_panel = c / n_colloc

        for i, (x_panel, gamma_i) in enumerate(zip(x_panels, gamma)):
            # Distance from panel to field point
            dx = x_field - x_panel
            dy = y_field - np.interp(x_panel, x_airfoil, camber_line)
            r_sq = dx**2 + dy**2

            # Avoid singularities
            r_sq = np.maximum(r_sq, 1e-10)

            # Vortex-induced velocity
            u += -gamma_i * dy / (2 * np.pi * r_sq) * dx_panel
            v += gamma_i * dx / (2 * np.pi * r_sq) * dx_panel

    return u, v, x_colloc, gamma, cl, cm_c4


def calculate_thin_airfoil(
    x_airfoil,
    camber_line,
    v_inf,
    alpha,
    camber_slope,
    x_upper,
    y_upper,
    x_lower,
    y_lower,
):
    # Create computational grid
    x_range = (-1.5, 1.5)
    y_range = (-1.5, 1.5)
    n_points = (100, 100)
    x = np.linspace(x_range[0], x_range[1], n_points[0])
    y = np.linspace(y_range[0], y_range[1], n_points[1])
    X, Y = np.meshgrid(x, y)
    c = 1
    # Calculate velocity field
    U, V, x_colloc, gamma, cl, cm_c4 = velocity_field(
        x_airfoil, camber_line, v_inf, alpha, camber_slope, X, Y
    )

    # Mask points inside airfoil
    mask_upper = np.zeros_like(X, dtype=bool)
    mask_lower = np.zeros_like(X, dtype=bool)

    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            x_point = X[i, j]
            y_point = Y[i, j]

            if 0 <= x_point <= c:
                y_upper_interp = np.interp(x_point, x_upper, y_upper)
                y_lower_interp = np.interp(x_point, x_lower, y_lower)

                if y_lower_interp <= y_point <= y_upper_interp:
                    mask_upper[i, j] = True
                    mask_lower[i, j] = True

    # Set velocity to zero inside airfoil
    U[mask_upper] = 0
    V[mask_upper] = 0
    cp = -2 * gamma / v_inf
    return X, Y, U, V, x_colloc, gamma, cl, cm_c4, cp
