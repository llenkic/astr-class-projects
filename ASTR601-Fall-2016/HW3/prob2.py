'''
Authors: Laura Lenkic

Date: 10/13/2016

Filename: prob2.py

Evolution of high and low scattering cross section photons traveling through a
material.
'''

# Import all necessary packages
import matplotlib.pyplot as plt
import numpy as np
import random
import sys

def photon_scat(L_init):
	"""
    Follows 10,000 photons all initially in the high scattering cross-section mode as they travel through a material. Both photon polarization modes have their own specific 
    probability of interacting. After they interact, they have equal probability of scattering in any direction. Determines the fraction of photons in each mode that cross 
    the surface of the material, and the angular distribution of these photons.
    
	param L_init: starting point below surface for all photons; must be negative
	"""
	# Check that user input for L_init is a negative number
	if L_init > 0:
		print('You must input a negative number!')
		sys.exit()

    # Create array of 10,000 photons
	N = 10000
	photon = np.ones(N)

    # Mean free path of low and high scattering cross-section modes respectively
	l_low = 1000
	l_high = 1

	# Probability that a photon in the low (high) mode will convert to the high (low) mode respectively
	P12 = 0.25
	P21 = 0.001

    # Arrays to store angles of low mode photons, high mode photons and the polarization of escaped photons and initialize the values of distance traveled between scattering events
	angles_low = []
	angles_high = []
	mode = []
	d_low = 0				
	d_high = 0

    # Initialize the count of photons that leave the surface of the material, the number of photons exiting that are in the high mode, and the angle (direction) of the photon travel after a scattering event
	leave = N
	num_high = 0
	angle = 0
	i = 0
    # Track each of the N photons as they travel through the material until it leaves the surface or until 20,000 scattering events take place; assume photon will stay in material and scatter indefinitely after 20,000 events
	while (i < N):
		scat = 1
		L = L_init
        # Generate three random numbers: 1. for calculating the distance traveled by the photon before it scatter, 2. for determining if the photon switches polarization modes, 3. for determining the direction of travel of the photon after scattering
		while (L < 0 and scat < 20000):
			x = random.random()
			y = random.random()
			d_low = -l_low*np.log(y)
			d_high = -l_high*np.log(y)				

			angle = np.random.uniform(0,2*np.pi)	
			if (photon[i] == 1 and x <= P21):
				photon[i] = 0
			elif (photon[i] == 0 and x <= P12):
				photon[i] = 1
			scat = scat + 1

			if (photon[i] == 1):
				L = L + d_high*np.cos(angle)
			if (photon[i] == 0):
				L = L + d_low*np.cos(angle)

			if (scat == 20000):
				leave = leave - 1
                    
        # Store the angle with which each photon exits the surface of the material and count the number of photons in the high mode
		if (scat < 20000):
			if (angle > np.pi):
				angle = 2*np.pi-angle
			else: 
				angle = angle
			angles_high.append(angle) if photon[i] == 1 else angles_low.append(angle)
			mode.append(photon[i])
		if (photon[i] == 1 and scat < 20000):
			num_high = num_high + 1

		i = i + 1

	# Print results to terminal and plot the angular distribution of photons exiting the surface for both polarization modes
	print('Number of escaped photons that are in the high mode: '+str(num_high))
	print('Total number of photons that escape: '+str(leave))
	print('Fraction of escaped high-mode photons: '+str(num_high/float(leave)))
	plt.hist(angles_high,50, alpha=0.5, color='#88CCEE', label='High-Mode Photons')
	plt.hist(angles_low, 50, alpha=0.5, color='#DDCC77', label='Low-Mode Photons')
	plt.legend(loc='upper right')
	plt.xlabel('Degrees [radians]')
	plt.ylabel('Number')
	plt.grid(which='major',axis='both',color='grey',alpha=0.2,linewidth=1, linestyle='--')
	plt.minorticks_on()
	plt.savefig('./'+str(abs(L_init))+'_Below.png', dpi=300, bbox_inches='tight', pad_inches=0.1)
	plt.clf()
