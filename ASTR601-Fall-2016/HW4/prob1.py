'''
Authors: Laura Lenkic

Date: 11/03/2016

Filename: prob1.py

Simplified Comptonization: redistribution of photon energies from scattering off electrons, assuming that initial photon energies are always 0.001 * electron rest mass energy, 
that electrons have a speed of 0.1 * speed of light, that the probability of scattering is independent of photon energy and the direction relative to the photon
'''

# Import all packages
import matplotlib.pyplot as plt
import numpy as np
import random

# Going to compare distribution of photon energies after 10, 100, and 1,000 scatterings
scat = [10,100,1000]

# Defining physical constants: electron mass (g), speed of light (cm/s), electron speed, gamma factor, and electron rest mass energy
m_e = 9.10938356e-28
c = 2.99792458e+10
b = 0.1
v_e = b*c
g = 1/(np.sqrt(1-b**2))
E_rest = m_e*c**2

# Transparency factors and colours for histograms that will be plotted
trans = [0.8,0.5,0.3]
cs = ['#EAC124','#3563EC','#D92405']

# Follow the energies 100,000 photons after 10, 100, and 1,000 scatterings
for k in range(3):
	E_scat = []
	for i in range(100000):
		E_init = 0.001*E_rest
		for j in range(scat[k]):
			#Boost from lab frame to initial rest frame of electron.
			theta_i = np.random.uniform(-1,1)
			E_in = E_init*(g*(1-b*theta_i))
			
			#Angle relative to direction of motion of electron, in electron rest frame.
			aberr = (theta_i-b)*(1-b*theta_i)**(-1)

			#Direction of photon after scattering, in electron rest frame, relative to electron.
			theta_f = np.random.uniform(-1,1)
			phi_f = np.random.uniform(0,2*np.pi)

			compton = theta_f*theta_i+np.sin(np.arccos(theta_i))*np.sin(np.arccos(theta_f))*np.cos(phi_f)
			recoil = E_in/(1.+(E_in/E_rest)*(1.-compton))

			#Boost out of electron rest frame.
			E_init = recoil*(g*(1+b*theta_f))
	
		E_scat.append(E_init)
	# Plot the distribution of photon energies for all three scattering cases and save the output
	plt.hist(np.log10(E_scat),100,range=(-11.0,-7.0),log=True,alpha=trans[k],color=cs[k],label=str(scat[k])+' Scatterings')

plt.legend(loc='upper left')
plt.xlabel('log(Photon Energy) [erg]')
plt.ylabel('Counts')
plt.grid(which='major',axis='both',color='grey',alpha=0.2,linewidth=1, linestyle='--')
plt.minorticks_on()
plt.savefig('./Compton.png', dpi=300, bbox_inches='tight', pad_inches=0.1)
