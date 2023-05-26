#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 05:59:00 2023

@author: slw
"""

import euler
import heun
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from timeit import timeit
import analytic_solution

#%% Federauslenkung
value_range = [0,70]
stepps = 500
m = 10.0
c = 1.0
k = 3.0
t0 = 0.0
y0 = 10
y1 = 0
rate_of_change = lambda t, y , y1: (-c * y1 - k * y) / m
spring_x, spring_y0 = euler.euler_ode2(rate_of_change, [t0, y0,y1],
                                       value_range, stepps)

y = []
for i in spring_x:
    y.append(analytic_solution.function_federposition(i))
    
spring_x_heun, spring_y0_heun = heun.heun_ode2(rate_of_change, [t0, y0,y1],
                                       value_range, stepps)

#rate_of_change_scipy = lambda y , y1, t: (-c * y1 - k * y) / m
def spring(y: list, t):
    m = 10.0
    c = 1.0
    k = 3.0
    y0 = y[0]
    y1 = y[1]
    dydt = np.array([float(y1), float((-c * y1 - k * y0) / m)])
    return dydt
y_scipy = odeint(spring, [y1, y0], spring_x)
print(y_scipy)
fig, ax = plt.subplots()
ax.plot(spring_x, y, color='y')
ax.plot(spring_x_heun, spring_y0_heun, '.', color='g')
ax.plot(spring_x, spring_y0, '.', color = 'b')
ax.plot(spring_x, y_scipy[:, 1], '.', color = 'r')
ax.set_xlabel('t [sec]')
ax.set_ylabel('s [cm]')
plt.show()
