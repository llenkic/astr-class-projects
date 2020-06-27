'''
Authors: Laura Lenkic

Date: 10/29/2017

Filename: integrators.py

Fourth-Order Runge-Kutta and Second-Order Leapfrog integrators and definition of equations to be integrated.
'''

import numpy as np

#---------------------------Functions for Problem 1---------------------------
def x_dot(cond,coeff):
	return coeff[0]*cond[0] - coeff[1]*cond[0]*cond[1] - coeff[2]*cond[0]

def y_dot(cond,coeff):
	return -1*coeff[3]*cond[1] + coeff[4]*cond[0]*cond[1] - coeff[5]*cond[1]
#-----------------------------------------------------------------------------

#---------------------------Functions for Problem 2---------------------------
def vx(cond,coeff):
	return cond[1]
	
def vx_dot(cond, coeff):
      return (-coeff[0]*cond[0])/(1 + coeff[0]*cond[0]**2 + coeff[1]*cond[2]**2)**(3/2)
 
def vy(cond,coeff):
	return cond[3]

def vy_dot(cond, coeff):
    return (-coeff[1]*cond[2])/(1 + coeff[0]*cond[0]**2 + coeff[1]*cond[2]**2)**(3/2)    
    
#-----------------------------------------------------------------------------


def rk4(functions,ics,coeff,t):
	"""
	Fourth-order Runge Kutta Integrator

	param functions: array of functions that are to be integrated
	param ics: array of initial conditions of problem
	param coeff: array of coefficients for functions
	param t: array of start point, end point, and time step to use in integration
	"""	
	# Vals stores each new calculated step in solution, to be used to calculate next step.
	vals = ics

	# Store estimated slopes, intermediate values of solutions, and final value of solution.
	k1 = np.zeros(len(functions))
	k2 = np.zeros(len(functions))
	k3 = np.zeros(len(functions))
	k4 = np.zeros(len(functions))
	int_vals = np.zeros(len(functions))
	dt = t[2]
	solutions = np.zeros((len(functions),t[1]/dt))
	for i in range(len(functions)):
		solutions[i][0] = ics[i]

	# Integrate from t = 0 to 100, in time steps dt.
	for j in np.arange(0 + t[2], t[1]/t[2]):
		# First estimate of next point in solution.
		for i in range(0,len(functions)):
			k1[i] = functions[i](vals,coeff)
		for i in range(0,len(functions)):
			int_vals[i] = vals[i] + 0.5*dt*k1[i]

		# Second estimate, based on first.
		for i in range(0,len(functions)):
			k2[i] = functions[i](int_vals,coeff)
		for i in range(0,len(functions)):
			int_vals[i] = vals[i] + 0.5*dt*k2[i]

		# Third estimate, based on second.
		for i in range(0,len(functions)):
			k3[i] = functions[i](int_vals,coeff)
		for i in range(0,len(functions)):
			int_vals[i] = vals[i] + dt*k3[i]

		# Final estimate based on four slope estimates; store in solution.
		for i in range(0,len(functions)):
			k4[i] = functions[i](int_vals,coeff)
		# for i in range(0,len(functions)):
			solutions[i][j] = vals[i] + dt*((1./6.)*k1[i] + (1./3.)*k2[i] + (1./3.)*k3[i] + (1./6.)*k4[i])

		# Update to new point to use in integration of next point.
		vals = []
		for i in range(0, len(functions)):
			vals.append(solutions[i][j])

	return solutions
 
   
def leapfrog(functions,ics,coeff,t):
	"""
	Second-order Leap Frog Integrator

	param functions: array of functions that are to be integrated
	param ics: array of initial conditions of problem
	param coeff: array of coefficients for functions
	param t: array of start point, end point, and time step to use in integration
	"""	
	# Vals stores each new calculated step in solution, to be used to calculate next step.
	vals = ics
	
	# Store intermediate values of solutions, and final value of solution.
	int_vals = np.zeros(len(functions))
	dt = t[2]
	solutions = np.zeros((len(functions),t[1]/dt))
	# setting initial conditions
	for i in range(len(functions)):
		solutions[i][0] = ics[i]
 
 	# Integrate from t = 0 to 100, in time steps dt.
	for j in np.arange(0 + t[2], t[1]/t[2]):
		
		# OPEN DRIFT: Calculate midstep of the in between functions
		for i in np.arange(0,len(functions), 2):
			int_vals[i] = vals[i] + 0.5*dt*vals[i+1]

		# KICK: Calculating the final for the main function
		for i in np.arange(1,len(functions), 2):
			solutions[i][j] = vals[i] + dt*functions[i](int_vals, coeff)
			
		# CLOSING DRIFT: Calculating the final for the in between functions
		for i in np.arange(0,len(functions), 2):
			solutions[i][j] = int_vals[i] + 0.5*dt*solutions[i+1][j]

		# Update to new point to use in integration of next point.
		vals = []
		for i in range(0, len(functions)):
			vals.append(solutions[i][j])


	return solutions
 
