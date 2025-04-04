import numpy as np
import os
import scipy.spatial as spatial
import scipy.integrate as integrate
import sympy as sp
import matplotlib.pyplot as plt




class ThinAirfoilTheory(object):
	def __init__(self):
		self.dat_file, self.txt_file = None, None
		self.x_coords_surf, self.y_coords_surf = np.array([]), np.array([])
		self.coords_surf_mat = np.array([])
		self.x_coords_mean, self.y_coords_mean = np.array([]), np.array([])
		self.coords_mean_mat = np.array([])
		self.poly_order = 5  # Set polynomial order to fit the mean line camber.
		self.aoa = None  # Preallocate the angle of attack.
		self.vor = None  # Object of the Voronoi diagram.
		self.best_pol_dict = dict()  # Dictionary containing the main information of the camber line's polynomial.
		self.best_pol = None
		self.best_pol_der = None  # First derivative of the best polynomial for mean camber line.
		self.sympy_poly_der = None  # Sympy polynomial of best_pol_der.
		self.sympy_der = None  # Sympy equation of sympy_poly_der with the change of variable required.
		self.lambda_der = None  # Lambda function of sympy_der.
		self.header = ''
		self.cl, self.zero_lift_angle = None, None


	def solve_theory(self, angle_of_attack, chord=1, n_coefficients=3, report=True):
		"""
		Solve the thin airfoil theory for a given angle of attack and chord length using the number of A coefficients
		as given in n_coefficients.
		:param angle_of_attack: Angle of attack of the airfoil IN RADIANS.
		:param chord: Chord length with same units as the x coordinates. Optional, default is normalized (c=1).
		:param n_coefficients: Number of coefficients to be solved. Must be greater than 2. Optional, default is 3.
		:param report: Boolean indicating if a report containing all relevant calculated information should be printed.
		:return: Tuple containing the corresponding coefficients.
		"""
		# Perform the change of variable. To do so, first transform the numpy polynomial to sympy.
		self.sympy_poly_der = sp.Poly(self.best_pol_der.coefficients, sp.Symbol('x'))

		# Perform the change of variable: x = c/2 * (1-cos(theta)).
		theta = sp.Symbol('theta')
		self.sympy_der = self.sympy_poly_der.subs({sp.Symbol('x'): 0.5*chord*(1-sp.cos(theta))})
		self.lambda_der = sp.lambdify(theta, self.sympy_der.as_expr(), modules='numpy')

		# Compute the coefficients.
		coefficients = []
		A0 = angle_of_attack - (1/np.pi)*integrate.quad(self.lambda_der, 0, np.pi)[0]
		coefficients.append(A0)
		assert n_coefficients >= 2, 'More than 1 coefficient should be computed in order to derive data from this theory'
		for i in np.arange(1, n_coefficients):
			coefficients.append((2/np.pi)*integrate.quad(lambda angle: self.lambda_der(angle)*np.cos(i*angle), 0, np.pi)[0])

		# Compute data derived from the theory.
		self._compute_relevant_data(coefficients, chord)

		return tuple(coefficients)

	def _compute_relevant_data(self, coefficients, chord):
		"""
		Compute all the information that can be derived from the thin airfoil theory given the necessary coefficients.
		:param coefficients: Array-like of all the computed coefficients.
		:param chord: The length of chord of the airfoil, with same units as the x coordinates.
		:return:
		"""
		# Compute the lift coefficient.
		self.cl = 2*np.pi*(coefficients[0] + 0.5*coefficients[1])

		# Compute the zero lift angle of attack.
		factor = -(1/np.pi)
		self.zero_lift_angle = factor*integrate.quad(lambda angle: self.lambda_der(angle)*(np.cos(angle)-1), 0, np.pi)[0]

		# Compute the moment coefficient about the Leading Edge.
		self.cm_le = -(self.cl/4 + np.pi/4*(coefficients[1] - coefficients[2]))

		# Compute the moment coefficient about the quarter chord.
		self.cm_quarter = np.pi/4 * (coefficients[2] - coefficients[1])

		# Compute the center of pressure.
		self.x_cp = chord/4*(1+np.pi/self.cl*(coefficients[1] - coefficients[2]))


if __name__ == '__main__':
	thin_theory = ThinAirfoilTheory()
	thin_theory.read_file(filename='NACA 2408.dat')
	thin_theory.build_mean_camber_line(max_pol_degree=7)
	A0, A1, A2 = thin_theory.solve_theory(angle_of_attack=5*np.pi/180, n_coefficients=3)