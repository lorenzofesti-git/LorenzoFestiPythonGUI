from tkinter import *
from tkinter import ttk
import sv_ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib import *
from numpy import *
import math
import numpy as np
import simpy as sp
import scipy.spatial as spatial
import scipy.integrate as integrate
import pathlib as path


def naca4_coordinates(naca, n_points, x):
    # Parse NACA number
    m = int(naca[0]) / 100.0  # max camber
    p = int(naca[1]) / 10.0  # location of max camber
    t = int(naca[2:]) / 100.0  # thickness

    # Thickness distribution
    yt = (
        5
        * t
        * (
            0.2969 * np.sqrt(x)
            - 0.1260 * x
            - 0.3516 * x**2
            + 0.2843 * x**3
            - 0.1015 * x**4
        )
    )

    # Camber line and its slope
    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)
    for i, xi in enumerate(x):
        if xi < p and p != 0:
            yc[i] = m / p**2 * (2 * p * xi - xi**2)
            dyc_dx[i] = 2 * m / p**2 * (p - xi)
        elif p != 0:
            yc[i] = m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * xi - xi**2)
            dyc_dx[i] = 2 * m / (1 - p) ** 2 * (p - xi)
        else:
            yc[i] = 0
            dyc_dx[i] = 0

    theta = np.arctan(dyc_dx)

    # Upper and lower surfaces
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)

    # Combine upper and lower surfaces
    xf = np.concatenate((np.flip(xu), xl))
    yf = np.concatenate((np.flip(yu), yl))
    return xu, xl, yu, yl, xf, yf, yc, yt, dyc_dx


def naca5_coordinates(naca, n_points, x):
    # Parse NACA 5-digit code
    naca = str(naca)
    pp = int(naca[1])
    cld = int(naca[0]) * 0.15  # design lift coefficient
    p = int(naca[1]) / 20.0  # position of max camber
    reflex = int(naca[3])  # reflex flag
    t = (int(naca[3]) * 10 + int(naca[4])) / 100.0  # thickness
    rmatstd = [0.0580, 0.1260, 0.2025, 0.2900, 0.3910]
    k1matstd = [361.400, 51.640, 15.957, 6.643, 3.230]
    rmatrfl = [0.1300, 0.2170, 0.3180, 0.4410]
    k1matrfl = [51.990, 15.793, 6.520, 3.191]
    k1k2matrfl = [0.000764, 0.00677, 0.0303, 0.1355]

    # Thickness distribution
    def yt(x):
        yt_v = (
            5
            * t
            * (
                0.2969 * np.sqrt(x)
                - 0.1260 * x
                - 0.3516 * x**2
                + 0.2843 * x**3
                - 0.1015 * x**4
            )
        )
        return yt_v

    # Camber line and its derivative
    def camber(x):
        if reflex == 0:
            # standard camber       da agg:     ([pp-1])
            yc = np.where(
                x < p,
                k1matstd[pp - 1]
                / 6
                * (
                    x**3
                    - 3 * rmatstd[pp - 1] * x**2
                    + rmatstd[pp - 1] ** 2 * (3 - rmatstd[pp - 1]) * x
                ),
                k1matstd[pp - 1] * rmatstd[pp - 1] ** 3 / 6 * (1 - x),
            )
            dyc_dx = np.where(
                x < p,
                k1matstd[pp - 1]
                / 6
                * (
                    3 * x**2
                    - 6 * rmatstd[pp - 1] * x
                    + rmatstd[pp - 1] ** 2 * (3 - rmatstd[pp - 1])
                ),
                -k1matstd[pp - 1] * rmatstd[pp - 1] ** 3 / 6,
            )
        else:
            # reflex camber  da agg: [pp-2]
            yc = np.where(
                x < p,
                k1matrfl[pp - 2]
                / 6
                * (
                    (x - rmatrfl[pp - 2]) ** 3
                    - k1k2matrfl[pp - 2] * x * (1 - rmatrfl[pp - 2]) ** 3
                    - x * rmatrfl[pp - 2] ** 3
                    + rmatrfl[pp - 2] ** 3
                ),
                k1matrfl[pp - 2]
                / 6
                * (
                    k1k2matrfl[pp - 2] * (x - rmatrfl[pp - 2]) ** 3
                    - k1k2matrfl[pp - 2] * x * (1 - rmatrfl[pp - 2]) ** 3
                    - x * rmatrfl[pp - 2] ** 3
                    + rmatrfl[pp - 2] ** 3
                ),
            )
            dyc_dx = np.where(
                x < p,
                k1matrfl[pp - 2]
                / 6
                * (
                    3 * (x - rmatrfl[pp - 2]) ** 2
                    - k1k2matrfl[pp - 2] * (1 - rmatrfl[pp - 2]) ** 3
                    - rmatrfl[pp - 2] ** 3
                ),
                k1matrfl[pp - 2]
                / 6
                * (
                    3 * k1k2matrfl[pp - 2] * (x - rmatrfl[pp - 2]) ** 2
                    + -k1k2matrfl[pp - 2] * (1 - rmatrfl[pp - 2]) ** 3
                    - rmatrfl[pp - 2] ** 3
                ),
            )
        return yc, dyc_dx

    # Cosine spacing (better point distribution)
    # beta = np.linspace(0, np.pi, n_points)
    m = cld / 0.3
    yt_vals = yt(x)
    yc, dyc_dx = camber(x)
    theta = np.arctan(dyc_dx)

    xu = x - yt_vals * np.sin(theta)
    yu = yc + yt_vals * np.cos(theta)
    xl = x + yt_vals * np.sin(theta)
    yl = yc - yt_vals * np.cos(theta)

    # Combine upper and lower surfaces
    xf = np.concatenate((np.flip(xu), xl))
    yf = np.concatenate((np.flip(yu), yl))
    return xu, xl, yu, yl, xf, yf, yc, yt_vals, dyc_dx


def nacanum(naca, x):
    corda = 1
    arr_naca = list(naca)
    naca_subtype = len(arr_naca)
    xn = np.asarray(x / float(corda))
    i = 0
    ymediana = 0 * x
    delty = 0 * x
    if naca_subtype == 4:
        [xu, xl, yu, yl, xf, yf, yc, yt, dyc_dx] = naca4_coordinates(naca, len(xn), x)
        delty = (yu + yl) / 2
        # print(delty)
        return xu, xl, yu, yl, xf, yf, yc, yt, dyc_dx
    elif naca_subtype == 5:
        [xu, xl, yu, yl, xf, yf, yc, yt, dyc_dx] = naca5_coordinates(naca, len(xn), x)
        delty = (yu + yl) / 2
        # print(delty)
        return xu, xl, yu, yl, xf, yf, yc, yt, dyc_dx
    else:
        return None  # Restituisce None se il numero non ha 4 o 5 arr_naca
