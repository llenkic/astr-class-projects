'''
Authors: Laura Lenkic

Date: 09/29/2016

Filename: prob1b.py

Evolution of states of 10,000 independent two-level atoms, when the Einstein coefficients A_21, B_21, B_12 are non-zero.
'''

# Import necessary packages
import matplotlib.pyplot as plt
import numpy as np
import random
import sys

# Set the initial parameters evoking the ground state and first excited state of hydrogen, and the number of atoms (nH)
x = 0.001
A_21 = 0.001
E = 10.4
kT = 5.0
g1 = 2.0
g2 = 8.0
nH = 10000.0
t = np.arange(0.001,10.001,0.001)
atom = np.zeros(10000)

# Calculate the probability that an atom will go from the ground state (1) to the excited state (2) and vice-versa. These two equations are derived from the fact that P_12 = B_12 * J and P_21 = A_21 + B_21 * J
P_12 = (g2/g1)*x/(np.exp(E/kT)-1)
P_21 = A_21 + x/(np.exp(E/kT)-1)

# Set the initial number of atoms in the ground state and first excited state -- these will store observed values
n1 = nH
n2 = 0

# List to hold the fraction of atoms in each state as a function of time -- these will store observed values
f1 = []
f2 = []

# Set the initial number of atoms in the ground state and first excited state -- these will store theoretical values for comparison
n1_th = nH
n2_th = 0

# List to hold the fraction of atoms in each state as a function of time -- these will store theoretical values for comparison
f1_th = []
f2_th = []

# Run through all the time states and calculate the number of atoms in each state based on the probabilities calculated on line 28 and 29 (P_12 and P21)
i = 0
while i < (len(t)):
	i2=0
    # Generate a random number between 0 and 1 and check if it is greater than P_12 or P_21. If it is, flip the atom to the opposite state
	while i2 < len(atom):
		if (atom[i2] == 0 and P_12 >= random.random()):
			atom[i2] = 1
		elif (atom[i2] == 1 and P_21 >= random.random()):
			atom[i2] = 0
		i2 = i2 + 1
    # Count the number of atoms in the excited state and calculate the number of atoms in the ground state, then calculate the fraction of atoms in each state and store the results.
	n2 = np.count_nonzero(atom)
	n1 = nH-n2
	print(i,n2,n1)
	f1.append(float(n1)/nH)
	f2.append(float(n2)/nH)

    # Calculate the theoretical expectation for the number of atoms in each state and the fractions; store the results.
	n1_th = n1_th + (nH-n1_th)*P_21 - n1_th*P_12
	f2_th.append((nH-n1_th)/float(nH))

    # Move to the next time step
	i = i + 1

# Plot the resulting fraction of atoms in the excited state and compare it to the theoretical prediction
plt.plot(t,f2,color='red',label='Observed Fraction')
plt.plot(t,f2_th,color='black',label='Theoretical Fraction')
plt.xlim(0,10)
plt.ylim(0,np.amax(f2)*1.1)
plt.xlabel(r't/t$_{0}$')
plt.ylabel('n$_{2}$/n$_{H}$')
plt.grid(which='major',axis='both',color='grey',alpha=0.2,linewidth=1,linestyle='--')
plt.minorticks_on()
plt.legend(loc='lower right')
plt.savefig('./n2-nH-fraction.png',dpi=300,bbox_inches='tight',pad_inches=0.1)
