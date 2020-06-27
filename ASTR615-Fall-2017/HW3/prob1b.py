#!/usr/bin/env python2
# -*- coding: utf-8 -*-

'''
Authors: Laura Lenkic and Elizabeth Tarantino

Date: 10/29/2017

Filename: prob1b.py

Integration of Lotka-Volterra Predator-Prey equations with a = b ~ 1.23.
'''

import integrators
import matplotlib.pyplot as plt

# Integrate equations, plot result, and print final populations to screen.
solutions = integrators.rk4([integrators.x_dot,integrators.y_dot],[30,3],[1,0.1,1.217,1.5,0.03,1.217],[0,100, 0.1])
plt.plot(solutions[0],solutions[1])
plt.xlabel('Population Density of Rabbits')
plt.ylabel('Population Density of Foxes')

print '\n'
print 'The rabbit population is: '+str(solutions[0][-1])+'.\nThe fox populatino is: '+str(solutions[1][-1])+'.'
print '\n'

plt.savefig('Population_0.png',dpi=300)
