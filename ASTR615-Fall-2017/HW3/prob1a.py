#!/usr/bin/env python2
# -*- coding: utf-8 -*-

'''
Authors: Laura Lenkic

Date: 10/29/2017

Filename: prob1a.py

Determining how populations of rabbits and foxes changes when we use coefficients of: A = 1, B = 0.1, C = 1.5, D = 0.03, a = b = 0.
'''

import integrators
import matplotlib.pyplot as plt

# Time steps for integrations
time_step = [1,0.5,0.25,0.1]

# Integrate Lotka-Volterra Predator-Prey equations and plot the results for each time given time step.
for i in range(0,len(time_step)):
	solutions = integrators.rk4([integrators.x_dot,integrators.y_dot],[30,3],[1,0.1,0,1.5,0.03,0],[0.,100., time_step[i]])
	plt.plot(solutions[0],solutions[1],label='h = '+str(time_step[i]))
	plt.xlabel('Population Density of Rabbits')
	plt.ylabel('Population Density of Foxes')
	plt.legend(loc='best')
	print("Finished the integration for time step: "+str(time_step[i])+".")

plt.savefig('Pop_Evolution.png',dpi=300)
