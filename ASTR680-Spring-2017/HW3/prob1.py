'''
Authors: Laura Lenkic

Date: 03/09/2017

Filename: prob1.py

Calculating the amount of mass a black hole must accrete to reach different values of black hole spin for two cases:
1. If the black hole has initial spin of 0 and accretes matter at the innermost stable circular orbit in only the prograde direction.
2. If the black hole has initial spin of 0.999 and accretes matter at the innermost stable circular orbit in only the retrograde directon.
'''

# Import necessary packages
import matplotlib.pyplot as plt
import numpy as np

# Set initial black hole parameters: dimensionless spin parameter, mass, spin parameter
j = 0.0
M = 1.0
a = j*M

# Initialize arrays to hold all values of dimensionless spin parameters and black hole masses, the incremental change in dimensionless spin parameter and mass
j_all = []
M_all = []
dJ = 0.0
dM = 0.0
M2 = M
m = 0.000001*M2
J = 0.0
# Accrete matter onto the black hole until its dimensionless sping parameter becomes nearly 1
while (j < 0.9999999999):
	# The radius of the innermost circular orbit depends on the spin paramenter; as the black hole accretes matter and spins up, this radius changes. Calculate this radius 
	# based on current black hole properties
	Z1 = 1.0 + (1.0-j**(2.0))**(1.0/3.0)*((1.0+j)**(1.0/3.0)+(1.0-j)**(1.0/3.0))
	Z2 = (3.0*j**(2.0)+Z1**2.0)**(0.5)

	r_ISCO = M*(3.0+Z2-((3.0-Z1)*(3.0+Z1+2.0*Z2))**(0.5))

	# Calculate constant appears frequently in equations for the specific angular momentum and specific energy of a particle being accreted
	C = np.sqrt(M*r_ISCO)

	# Specific angular momentum and specific energy (i.e., per unit rest mass) of a particle in a circular geodesic around a black hole of mass M and spin parameter jM
	u_phi = (C*(r_ISCO**(2.0)-2.0*a*C+a**(2.0)))/(r_ISCO*(r_ISCO**(2.0)-3.0*M*r_ISCO+2.0*a*C)**(0.5))
	u_t = (r_ISCO**(2.0)-2.0*M*r_ISCO+a*C)/(r_ISCO*(r_ISCO**(2.0)-3.0*M*r_ISCO+2.0*a*C)**(0.5))

	# Save current j and M into array
	j_all.append(j)
	M_all.append(M)

	# Calculate the increment in j and M, then increment these values, and set new spin parameter of black hole
	dJ = m*u_phi
	dM = m*u_t
	J=J+dJ
	M=M+dM
	j = J/(M**2)
	a = j*M

# Mass needed for black hole to accrete to reach various levels of spin
index1 = np.abs(np.subtract.outer(j_all,0.5)).argmin(0)
index2 = np.abs(np.subtract.outer(j_all,0.9)).argmin(0)
index3 = np.abs(np.subtract.outer(j_all,0.9999999999)).argmin(0)

# Print mass accreted to reach different levels of spin
print('Part a:')
print('The black hole accretes '+str(M_all[index1]-1)+'M0 to get to j = 0.5')
print('The black hole accretes '+str(M_all[index2]-1)+'M0 to get to j = 0.9')
print('The black hole accretes '+str(M_all[index3]-1)+'M0 to get to j = 0.9999999999')

# Plot the change in black hole spin as a function of mass as it accretes matter at the innermost stable circular orbit
plt.plot(M_all,j_all)
plt.grid(which='major',axis='both',color='grey',alpha=0.2,linewidth=1,linestyle='--')
plt.minorticks_on()
plt.xlim([np.min(M_all),np.max(M_all)])
plt.ylim([0,1.005])
plt.xlabel('Black Hole Mass')
plt.ylabel('Spin parameter, j')
plt.savefig('Prograde.png',dpi=300,bbox_inches='tight',pad_inches=0.1)

# Repeat above for a black hole with j close to 1 and calculate how much matter it needs to accrete to reach spin of 0. 
j = 0.999
M = 1.0
a = j*M

m = 0.001*M
j_all = []
M_all = []
dJ = 0.0
dM = 0.0
M2 = M
J = 0.999
while (j > 0):
	Z1 = 1.0 + (1.0-j**(2.0))**(1.0/3.0)*((1.0+j)**(1.0/3.0)+(1.0-j)**(1.0/3.0))
	Z2 = (3.0*j**(2.0)+Z1**2.0)**(0.5)

	r_ISCO = M*(3.0+Z2+((3.0-Z1)*(3.0+Z1+2.0*Z2))**(0.5))

	C = np.sqrt(M*r_ISCO)

	u_phi = -(C*(r_ISCO**(2.0)+2.0*a*C+a**(2.0)))/(r_ISCO*(r_ISCO**(2.0)-3.0*M*r_ISCO-2.0*a*C)**(0.5))
	u_t = (r_ISCO**(2.0)-2.0*M*r_ISCO-a*C)/(r_ISCO*(r_ISCO**(2.0)-3.0*M*r_ISCO-2.0*a*C)**(0.5))

	j_all.append(j)
	M_all.append(M)

	dJ = m*u_phi
	dM = m*u_t
	J = J + dJ
	M = M + dM
	j = J/(M**2)
	a = j*M

index1 = np.abs(np.subtract.outer(j_all,0.99)).argmin(0)
index2 = np.abs(np.subtract.outer(j_all,0.9)).argmin(0)
index3 = np.abs(np.subtract.outer(j_all,0.5)).argmin(0)
index4 = np.abs(np.subtract.outer(j_all,0.0)).argmin(0)

print('Part b:')
print('The black hole accretes '+str(M_all[index1]-1)+'M0 to get to j = 0.99')
print('The black hole accretes '+str(M_all[index2]-1)+'M0 to get to j = 0.9')
print('The black hole accretes '+str(M_all[index3]-1)+'M0 to get to j = 0.5')
print('The black hole accretes '+str(M_all[index4]-1)+'M0 to get to j = 0')

plt.clf()
plt.plot(M_all,j_all)
plt.grid(which='major',axis='both',color='grey',alpha=0.2,linewidth=1,linestyle='--')
plt.minorticks_on()
plt.xlabel('Black Hole Mass')
plt.ylabel('Spin parameter, j')
plt.savefig('Retrograde.png',dpi=300,bbox_inches='tight',pad_inches=0.1)
