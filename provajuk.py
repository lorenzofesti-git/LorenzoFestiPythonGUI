import numpy as np
import math
from matplotlib import *
import matplotlib.pyplot as plt

def conformal_map():

   # Joukowski transform parameters
   c = 1.0                             # transform parameter
   h, k = -0.15, 0                     # center of circle in z plane
   R_0 = 1.15                          # circle radius 

   # ===== Joukoski transform curves ============================

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

   # plot z plane and zeta plane grids

   #----------------------------------------------------------

   # ============ Solving flow over the airfoil ===============
   U = 1.0                           # Uniform flow velocity
   aoa = 20.0 * math.pi / 180         # angle of attack
   Dstr = R_0**2 * 2 * math.pi * U   # doublet strength


   # grid in the zp (z prime) reference frame
   Xp = (X - h) * np.cos(aoa) + (Y - k) * np.sin(aoa)
   Yp = (Y - k) * np.cos(aoa) - (X - h) * np.sin(aoa)
   plt.plot(Xp,Yp)
   plt.show