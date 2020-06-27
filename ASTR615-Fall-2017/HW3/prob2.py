#!/usr/bin/env python2
# -*- coding: utf-8 -*-

'''
Authors: Laura Lenkic

Date: 10/30/2017

Filename: prob2.py

Finding and plotting the two-dimensional orbit for a point particle from the given potential 
'''

import integrators
import matplotlib.pyplot as plt
import numpy as np 

# Equations needed for part c
# potential
def pot(x, y, alpha, beta):
	return -1/(np.sqrt(1 + alpha*x**2 + beta*y**2))

# energy	
def E(vx, vy, potential):
	return ((vx**2 + vy**2)/2) + potential

# constants
alpha = 2
beta = 2
x = 1
y = 0
xdot = 0 
ydot = 0.1
start = 0.
end = 100. 

# time steps
tsteps = [1, 0.5, 0.25, 0.1]

# looping through the time steps
for i in range(len(tsteps)):
	t = [start, end, tsteps[i]]
	time = np.arange(t[0], t[1], t[2])
	
	# running integrators
	lf_sol = integrators.leapfrog([integrators.vx,integrators.vx_dot, integrators.vy, integrators.vy_dot],[x, xdot, y, ydot],[alpha, beta],t)
	rk4_sol = integrators.rk4([integrators.vx,integrators.vx_dot, integrators.vy, integrators.vy_dot],[x, xdot, y, ydot],[alpha, beta],t)

	# plotting the position
	plt.figure(1)
	plt.clf()
	plt.plot(lf_sol[0], lf_sol[2], c = 'b', label = 'Leapfrog')
	plt.plot(rk4_sol[0], rk4_sol[2], c = 'g', label = 'RK4')
	plt.title('Orbit of a Point Particle with dt = {}'.format(tsteps[i]))
	plt.xlabel('x')
	plt.ylabel('y')
	plt.legend(loc='best')
	plt.savefig('orbit_tstep{}'.format(i))
	
	# plotting the energy
	plt.figure(2)
	plt.clf()
	plt.plot(time, E(lf_sol[1], lf_sol[3], pot(lf_sol[0], lf_sol[2], alpha, beta)), c='b', label = 'Leapfrog')
	plt.plot(time, E(rk4_sol[1], rk4_sol[3], pot(rk4_sol[0], rk4_sol[2], alpha, beta)), c='G', label = 'RK4')
	plt.title('Energy of a Point Particle with dt = {}'.format(tsteps[i]))
	plt.xlabel('Time')
	plt.ylabel('Energy')
	plt.legend(loc='best')
	plt.savefig('energy_tstep{}'.format(i))
