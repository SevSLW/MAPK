#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 02:05:36 2023

@author: slw
"""
import euler
import heun
import matplotlib.pyplot as plt
import analytic_solution
import numpy as np
from scipy.integrate import odeint
#import cProfile
#import pstats
from timeit import timeit

#%% Flüssigkeit abkühlen T'(t) = - 0.1 (T(t) -20) T'(t) + 0.1 (T(t) -20) = 0
value_range = [0,60]
stepps = 20
rate_of_cange = lambda x, y: - 0.1 * (y - 20)
cooling_points = euler.euler_method(rate_of_cange, 
                                    [0, 95], 
                                    value_range, 
                                    stepps)
#with cProfile.Profile() as pr:
    #euler.euler_method(rate_of_cange, [0, 95], value_range, stepps)
#cooling_points = euler_method(lambda x, y: -sin(x), [0, 1], [0, 20], 100)
time = timeit(lambda: euler.euler_method(rate_of_cange, [0, 95], value_range, stepps), 
              number=1000)
print(f'time euler: {time}')


expr = analytic_solution.function_cooling
x = np.linspace(value_range[0], value_range[1], stepps)
y = []
for i in x:
    y.append(expr(i))



distance_sum = 0
distance_list = []
for i, value in enumerate(x):
    distance = np.absolute(y[i] - cooling_points['y'][i])
    distance_list.append(distance)
    distance_sum += distance
print(distance_sum)
print(np.mean(distance_list))
#stats = pstats.Stats(pr)
#stats.sort_stats(pstats.SortKey.TIME)
#stats.print_stats()


cooling_points_heun = heun.heun_method(rate_of_cange, [0, 95], value_range, stepps)
time = timeit(lambda: heun.heun_method(rate_of_cange, [0, 95], value_range, stepps), 
              number=1000)
print(f'time heun: {time}')
distance_sum = 0
distance_list = []
for i, value in enumerate(x):
    distance = np.absolute(y[i] - cooling_points_heun['y'][i])
    distance_list.append(distance)
    distance_sum += distance
print(distance_sum)
print(np.mean(distance_list))

rate_of_cange = lambda y, x: - 0.1 * (y - 20)
y_scipy = odeint(rate_of_cange, 95, x)
time = timeit(lambda: odeint(rate_of_cange, 95, x), 
              number=1000)
print(f'time odeint: {time}')
distance_sum = 0
distance_list = []
for i, value in enumerate(x):
    distance = np.absolute(y[i] - y_scipy[i])
    distance_list.append(distance)
    distance_sum += distance
print(distance_sum)
print(np.mean(distance_list))



fig, ax  = plt.subplots()
ax.plot(cooling_points['x'], cooling_points['y'], '>', color='b')
ax.plot(x, y, color='y')
ax.plot(cooling_points_heun['x'], cooling_points_heun['y'], '<', color='g')
ax.plot(x, y_scipy, '.', color='r')
ax.set_xlabel('t [min]')
ax.set_ylabel('T [°C]')
plt.show()

time_euler = []
time_heun = []
time_odeint = []
disttot_euler = []
disttot_heun = []
disttot_odeint = []
distav_euler = []
distav_heun = []
distav_odeint = []

for i in np.arange(5, 400, 5):
    alue_range = [0,60]
    stepps = i
    rate_of_cange = lambda x, y: - 0.1 * (y - 20)
    cooling_points = euler.euler_method(rate_of_cange, 
                                        [0, 95], 
                                        value_range, 
                                        stepps)
    #with cProfile.Profile() as pr:
        #euler.euler_method(rate_of_cange, [0, 95], value_range, stepps)
    #cooling_points = euler_method(lambda x, y: -sin(x), [0, 1], [0, 20], 100)
    time = timeit(lambda: euler.euler_method(rate_of_cange, [0, 95], value_range, stepps), 
                  number=1000)
    time_euler.append(time)


    expr = analytic_solution.function_cooling
    x = np.linspace(value_range[0], value_range[1], stepps)
    y = []
    for i in x:
        y.append(expr(i))



    distance_sum = 0
    distance_list = []
    for i, value in enumerate(x):
        distance = np.absolute(y[i] - cooling_points['y'][i])
        distance_list.append(distance)
        distance_sum += distance
    disttot_euler.append(distance_sum)
    distav_euler.append(np.mean(distance_list))
    #stats = pstats.Stats(pr)
    #stats.sort_stats(pstats.SortKey.TIME)
    #stats.print_stats()


    cooling_points_heun = heun.heun_method(rate_of_cange, [0, 95], value_range, stepps)
    time = timeit(lambda: heun.heun_method(rate_of_cange, [0, 95], value_range, stepps), 
                  number=1000)
    time_heun.append(time)
    distance_sum = 0
    distance_list = []
    for i, value in enumerate(x):
        distance = np.absolute(y[i] - cooling_points_heun['y'][i])
        distance_list.append(distance)
        distance_sum += distance
    disttot_heun.append(distance_sum)
    distav_heun.append(np.mean(distance_list))


    rate_of_cange = lambda y, x: - 0.1 * (y - 20)
    y_scipy = odeint(rate_of_cange, 95, x)
    time = timeit(lambda: odeint(rate_of_cange, 95, x), 
                  number=1000)
    time_odeint.append(time)
    distance_sum = 0
    distance_list = []
    for i, value in enumerate(x):
        distance = np.absolute(y[i] - y_scipy[i])
        distance_list.append(distance)
        distance_sum += distance
    disttot_odeint.append(distance_sum)
    distav_odeint.append(np.mean(distance_list))

x = np.arange(5, 400, 5)
fig, axs = plt.subplots(2)
axs[0].plot(x[0:12], disttot_euler[0:12], color='b')
axs[0].plot(x[0:12], disttot_heun[0:12], color='g')
axs[0].plot(x[0:12], disttot_odeint[0:12], color='r')
axs[1].plot(x[0:12], distav_euler[0:12], color='b')
axs[1].plot(x[0:12], distav_heun[0:12], color='g')
axs[1].plot(x[0:12], distav_odeint[0:12], color='r')
plt.show()

fig, ax = plt.subplots()
ax.plot(x, time_euler, color='b')
ax.plot(x, time_heun, color='g')
ax.plot(x, time_odeint, color='r')
plt.show()




