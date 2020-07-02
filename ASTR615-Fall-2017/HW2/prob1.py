'''
Authors: Laura Lenkic

Date: 10/06/2017

Filename: prob2a.py

This code fits a Lorentzian and Gaussian to the provided data file using
non-linear least squares fitting routines.
'''

# import necessary modules
import numpy as np
import os.path
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import matplotlib.pyplot as plt

######################## FUNCTIONS ########################
# Lorentzian function
def lorentzian(nu, nu0, aL):
	return (1./np.pi)*(aL/( np.power(nu-nu0,2.) + np.power(aL,2.) ))

# Gaussian function	
def gaussian(nu, nu0, aD):
	return (1./aD)*np.sqrt(np.log(2.)/np.pi)*np.exp(-1.*np.log(2.)*np.power(nu-nu0,2.)/np.power(aD,2.))
	
# Chi Square function
def chisq(data, data_err, expected, npar):
	return np.sum(np.power(data-expected,2.)/np.power(data_err,2.))/float(len(data)-npar)
###########################################################

# the input file
# yes, I renamed it...why would you put a space in a file name :( 
file = "./homework_2.txt"

if os.path.isfile(file)==False:
	print
	print "The expected data file is not in the current directory."
	print "Please move homework_2.txt to the same directory as"
	print "hw2_1.py or adjust the path in the file definition in line 29"
	print "of hw2_1.py to the current location of homework_2.txt"
	print
	exit()

# read columns of data file in arrays
nu, phi, phi_err = np.genfromtxt(file, unpack=True)

# calculate the mean of nu to approximate nu_0
# this will be provided to the fitting routines as a first guess for the value of nu_0
nu_mean = np.mean(nu)

# use scipy's curve_fit to fit the Lorentzian
# provide a guess for nu_0, but not for alpha_L
# one-sigma error bars on parameter estimates are the sqrt of the diagonals of the covariance matrix l_pcov
l_popt, l_pcov = curve_fit(lorentzian, nu, phi, p0=[nu_mean, 1.], sigma=phi_err)
l_perr = np.sqrt(np.diag(l_pcov))
l_mod = lorentzian(nu,*l_popt)
l_chi = chisq(phi,phi_err,l_mod, 2.)
print
print "Lorentzian Fit with 1-sigma Errors"
print "   nu_0 =",l_popt[0],"+/-", l_perr[0]
print "   alpha_L =",l_popt[1],"+/-", l_perr[1]
print "   Reduced Chi Square =",l_chi
print 

# use scipy's curve_fit to fot the Gaussian
# provide a guess for nu_0, but not for alpha_D
# one-sigma error bars on parameter estimates are the sqrt of the diagonals of the covariance matrix g_pcov
g_popt, g_pcov = curve_fit(gaussian, nu, phi, p0=[nu_mean, 1.], sigma=phi_err)
g_perr = np.sqrt(np.diag(g_pcov))
g_mod = gaussian(nu,*g_popt)
g_chi = chisq(phi,phi_err,g_mod, 2.)
print
print "Gaussian Fit with 1-sigma Errors"
print "   nu_0 =",g_popt[0],"+/-", g_perr[0]
print "   alpha_D =",g_popt[1],"+/-", g_perr[1]
print "   Reduced Chi Square =",g_chi
print

# make a plot!
plt.plot(nu, lorentzian(nu, *l_popt), color='red', label='Lorentzian')
plt.plot(nu, gaussian(nu, *g_popt), color='lime', label='Gaussian')
plt.errorbar(nu,phi,yerr=phi_err,fmt='.', color='gray', label='data')
plt.ylim([-.005, .045])
plt.xlabel(r'$\nu$')
plt.ylabel(r'$\phi$')
plt.legend(frameon=False, numpoints=1)
plt.show()
