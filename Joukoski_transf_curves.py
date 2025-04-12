import numpy as np
import math
from matplotlib import *
import matplotlib.pyplot as plt
import plotly.figure_factory as ff


def conformal_map():

   # Joukowski transform parameters
   c = 1.0                             # transform parameter
   h, k = 0.1, 0.1                     # center of circle in z plane                   
   R_0 =math.sqrt((c - h)**2 + k**2)   # circle radius

   # =================================

   # curve in z plane
   n = 500
   x = np.linspace(-R_0 + h, R_0 + h, n)
   # yu = np.sqrt(R_0**2 - (x - h)**2) + k    # upper curve
   # yl = -np.sqrt(R_0**2 - (x - h)**2) + k   # lower curve

   # # hack to fix NaNs in yu, yl if sqrt(very small number) occurs
   # yu[np.argwhere(np.isnan(yu))] = k
   # yl[np.argwhere(np.isnan(yl))] = k

   # zu = x + yu * 1j   # upper curve
   # zl = x + yl * 1j   # lower curve

   # # zeta plane curve
   # zeta_l = zl + c**2 / zl
   # zeta_u = zu + c**2 / zu


   #===== generating the grid ===========================

   # grid in polar coordinates
   Rlim = 5                          # domain limit in r 
   Nr, Ntheta = 100, 145             # number of grid points in r and theta
   r = np.linspace(R_0, Rlim, Nr)
   theta = np.linspace(0, 2 * np.pi, Ntheta)   
   R, T = np.meshgrid(r, theta)

   # convert polar grid to cartesian 
   X = R * np.cos(T) + h
   Y = R * np.sin(T) + k

   # Joukoski transform on grid
   Z = X + Y*1j
   zeta_grid = Z + c**2 / Z


   # ============ Solving flow over the airfoil ===============
   U = 1.0                           # Uniform flow velocity
   aoa = 20.0 * math.pi / 180         # angle of attack
   Dstr = R_0**2 * 2 * math.pi * U   # doublet strength


   # grid in the zp (z prime) reference frame
   Xp = (X - h) * np.cos(aoa) + (Y - k) * np.sin(aoa)
   Yp = (Y - k) * np.cos(aoa) - (X - h) * np.sin(aoa)

   # Kutta condition (stagnation point at trailing edge)
   Vstr = -Yp[0, 0] * 4 * np.pi * U

   #velocity field in zp plane 

   v_r = np.zeros(np.shape(X))
   v_t = np.zeros(np.shape(X))
   psi = np.zeros(np.shape(X))
   v_r=(U*(1-((R_0**2)/(r**2)))*np.cos(theta))
   v_t=(-U*(1+((R_0**2)/(r**2)))*np.sin(theta))-(Vstr/(2 * np.pi * r))
   u=v_r*np.cos(theta)-v_t*np.sin(theta)
   v=v_r*np.sin(theta)+v_t*np.cos(theta)
   # velocity field in zeta plane
   dzeta_dz = 1 - (c/Z)**2
   V_zeta = (u - v * 1j) / dzeta_dz
   u_zeta = V_zeta.real
   v_zeta = -V_zeta.imag

   x_dd = np.linspace(-3, 3, 100)
   y_dd = np.linspace(-3, 3, 100)
   Y_d, X_d = np.meshgrid(x_dd, y_dd)
   # da mettere nello script principale
   fig = ff.create_streamline(x_dd, y_dd, u_zeta, v_zeta, arrow_scale=.1)
   fig.show() 
  