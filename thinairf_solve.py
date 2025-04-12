import numpy as np
import os
import scipy.spatial as spatial
import scipy.integrate as integrate
import sympy as sp
import matplotlib.pyplot as plt



def solve_thin_airf_theory(angle_of_attack, chord=1, n_coefficients=3, report=True):
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
	sympy_poly_der = sp.Poly(best_pol_der.coefficients, sp.Symbol('x'))

	# Perform the change of variable: x = c/2 * (1-cos(theta)).
	theta = sp.Symbol('theta')
	sympy_der = sympy_poly_der.subs({sp.Symbol('x'): 0.5*chord*(1-sp.cos(theta))})
	lambda_der = sp.lambdify(theta, sympy_der.as_expr(), modules='numpy')

	# Compute the coefficients.
	coefficients = []
	A0 = angle_of_attack - (1/np.pi)*integrate.quad(lambda_der, 0, np.pi)[0]
	coefficients.append(A0)
	assert n_coefficients >= 2, 'More than 1 coefficient should be computed in order to derive data from this theory'
	for i in np.arange(1, n_coefficients):
		coefficients.append((2/np.pi)*integrate.quad(lambda angle: lambda_der(angle)*np.cos(i*angle), 0, np.pi)[0])

	# Compute data derived from the theory.
	compute_relevant_data(coefficients, chord)

	return tuple(coefficients)

def compute_relevant_data(coefficients, chord):
	"""
	Compute all the information that can be derived from the thin airfoil theory given the necessary coefficients.
	:param coefficients: Array-like of all the computed coefficients.
	:param chord: The length of chord of the airfoil, with same units as the x coordinates.
	:return:
	"""
	# Compute the lift coefficient.
	cl = 2*np.pi*(coefficients[0] + 0.5*coefficients[1])

	# Compute the zero lift angle of attack.
	factor = -(1/np.pi)
	zero_lift_angle = factor*integrate.quad(lambda angle: lambda_der(angle)*(np.cos(angle)-1), 0, np.pi)[0]

	# Compute the moment coefficient about the Leading Edge.
	cm_le = -(cl/4 + np.pi/4*(coefficients[1] - coefficients[2]))

	# Compute the moment coefficient about the quarter chord.
	cm_quarter = np.pi/4 * (coefficients[2] - coefficients[1])

	# Compute the center of pressure.
	x_cp = chord/4*(1+np.pi/cl*(coefficients[1] - coefficients[2]))

