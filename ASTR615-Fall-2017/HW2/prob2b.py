'''
Authors: Laura Lenkic

Date: 10/06/2017

Filename: prob2b.py

Estimates the volume of the surface described in the homework sheet using Monte Carlo integration with 10 to 10^7 samples, in steps of powers of 10.
'''

import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

# Define arrays for storing estimated masses for each given number of samples.
masses = []
Ns = []
errs = []

# For powers of 10 ranging from 1 to 7, the volume element is the volume of the cube divided by the number of samples.
for k in range(1,8):
	N = 10**k
	dV = (6.**3)/N
	mass = 0
	err = 0
	mean = 0
	pts_in = 0

	# Generate N*3 (because 3D problem) random points within the cube, then for each set of (x,y,z) coordinates, check if the point lies within the boundaries defined by the 
	# surface of the shape of interest.
	x = []
	y = []
	z = []
	for j in range(0,N):
		pt = []
		for i in range(0,3):
			pt.append((6*random.random())-3)

		# If the point lies within the shape, calculate the density at that point and multiply by dV to get the mass element. Sum all mass elements to get the total mass.
		if ((pt[0]**2 + pt[1]**2 + pt[2]**2 <= 9) and (pt[0]**2 - pt[1]**2 >= 0)):
			rho = 1 + pt[0]**2 + 2*(pt[1] + pt[2])**2
			mass += rho*dV
			err += (rho)**2
			mean += rho
			pts_in += 1
			if (j%500 == 0):
				x.append(pt[0])
				y.append(pt[1])
				z.append(pt[2])
			
	# Store total mass, number of samples, and 1 sigma uncertainty in arrays for plotting, and print the result to the screen.
	masses.append(mass)
	Ns.append(N)
	err = err/N
	mean = mean/N
	errs.append(6**3*(np.sqrt((err-mean**2)/N)))
	print "Number of Samples: "+str(N)+", Total Mass Estimated: "+str(mass)+" +/- "+str(6**3*(np.sqrt((err-mean**2)/N)))+", Number of Points in Sphere: "+str(pts_in)

# Plot estimated total mass as a function of number of samples.
plt.scatter(np.log10(Ns),masses)
plt.errorbar(np.log10(Ns),masses, yerr=errs, fmt='none')
plt.xlabel('log(Number of samples in integration, N)')
plt.ylabel('Total Mass Estimated')
plt.xlim([0.5,7.5])
plt.ylim([np.min(masses)-100,np.max(masses)+100])
plt.savefig('Masses.png', dpi=300)

# Plot the 3D volume.
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z);
plt.savefig('Volume.png', dpi=300)

