'''
Authors: Laura Lenkic

Date: 10/06/2017

Filename: prob2a.py

Estimates the volume of a d dimensional hypershpere using Monte Carlo integration with 10 to 10^7 samples, in steps of powers of 10.
'''

import random
import numpy as np
from scipy.misc import factorial
from scipy.misc import factorial2

def hsphere_vol(d):
	"""
	This function computes the volume of a d dimensional hypersphere, using Monte Carlo integration. It generates N random points within the smallest cube enclosing the
	sphere, then finds the number of points that lie within the hypersphere. The volume of the hypersphere is then the ratio of points withing the hypersphere and the
	number of random points generated, times the volume of the cube.

	param d: dimension of the hypersphere
	"""
	# Calculate the value of the gamma function based on whether d is even or odd.
	if (d%2 == 0):
		gamma_fnc = factorial(d/2.)
	else:
		gamma_fnc = np.sqrt(2)*factorial2(d-2)/2.**((d-1)/2.)
	
	# Calculate the volume of the hypershpere from the analytic formula and print to screen.
	analytic_vol = np.pi**(d/2.)/(gamma_fnc)
	print "The volume obtained from the analytic formula is: "+str(analytic_vol)+"."

	# Estimate the volume with samples ranging from 10 to 10^7 in steps of powers of 10.
	for i in range(1,8):
		N = 10**i
		cube_pts = []
		sphe_pts = 0

		# Generate N*d random points within the smallest hypercube containing the hypersphere. Hypersphere has radius 1, so the hypercube has sides of length 2.
		for j in range(0,N):
			for k in range(0,d):
				pt = (2*random.random())-1
				cube_pts.append(pt)

		cube_pts = np.reshape(np.asarray(cube_pts), (N,d))

		# Take each random point and calculate it's length, if it is less than 1, then add 1 to the count of points that lie in the hypersphere.
		for j in range(0,N):
			length = 0
			for k in range(0,d):
				length += cube_pts[j][k]**2

			if (length <= 1):
				sphe_pts += 1

		# Caclculate volume and uncertainties of hypersphere by multiplying the ratio of points in the hypersphere to the number of points in the hypercube by the 
		# volume of the hypercube.
		sphe_vol = sphe_pts*(2.0**d)/N
		err = (2.0**d)*np.sqrt((float(sphe_pts)/N - (float(sphe_pts)/N)**2)/N)

		# Print result to screen.
		print "Number of Samples: "+str(N)+", Estimated Volume of "+str(d)+"D sphere: "+str(sphe_vol)+" +/- "+str(err)+", Number of Points in Sphere: "+str(sphe_pts)

