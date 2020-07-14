'''
Authors: Laura Lenkic

Date: 02/23/2017

Filename: prob2.py

Specific energy release of matter settling on to the equator of a rotating neutron star, as a function of the neutron star rotation rate.
'''

# Import necessary packages
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Define physical constants and neutron star properties in CGS units: gravitational constant, speed of light, neutron star mass, neutron star radius
G = 6.674e-8			#cm^3/g/s^2
c = 2.998e10			#cm/s
M = 1.5*(1.989e33)
R = 1000000.0
const = 2*G*M/(R*c**2)

# Array for specific energy and rotation rate
u_t = []
omega = []
# Calculate the specific energy released by matter falling in for a range of rotation rates
for nu in np.arange(0,c/(R+10000000),0.1):
	w = nu*2*np.pi
	omega.append(nu)
	u_t.append(c**2-np.sqrt((1-const)**2/((1-const)-(R**2)*(w**2)/c**2))*c**2)

plt.plot(omega,u_t)
plt.grid(which='major',axis='both',color='grey',alpha=0.2,linewidth=1,linestyle='--')
plt.minorticks_on()
plt.xlabel('Rotation Frequency [Hz]')
plt.ylabel('Specific Energy Release [erg/g]')
plt.xlim([np.amin(omega),np.amax(omega)+100])
plt.ylim([np.amin(u_t),np.amax(u_t)+0.2E+20])
plt.savefig('SpecificEnergy.png',dpi=300,bbox_inches='tight',pad_inches=0.1)
