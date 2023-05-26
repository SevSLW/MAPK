#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 23:31:00 2023

@author: slw
"""
import analytic_solution
import numpy as np
import matplotlib.pyplot as plt

def euler_method(ode, start_condition: list, value_range: list, stepps: int):
    x_list = [start_condition[0]]
    y_list = [start_condition[1]]
    x = start_condition[0]
    x = 0
    y = start_condition[1]
    for i in np.linspace(value_range[0], value_range[1], stepps):
        if i == 0:
            continue
        x = i
        y += value_range[1] / stepps * ode(x, y)
        x_list.append(x)
        y_list.append(y)
    return {'x': x_list,
                'y': y_list}

def euler_ode2(ode, start_condition: list, value_range: list, stepps: int):
    x_list = [start_condition[0]]
    y_list = [np.array([start_condition[1], start_condition[2]])]
    x = start_condition[0]
    y = y_list[-1]
    y0_list = [start_condition[1]]
    y1_list = [start_condition[2]]
    for i in np.linspace(value_range[0], value_range[1], stepps):
        if i == 0:
            continue
        x = i
        y0 = y_list[-1][0]
        y1 = y_list[-1][1]
        y = y_list[-1] + value_range[1]/ stepps * np.array([y1_list[-1], ode(x, y0, y1)])
        x_list.append(x)
        y_list.append(y)
        y0_list.append(y[0])
        y1_list.append(y[1])
    return x_list, y0_list
        

# #%% Flüssigkeit abkühlen T'(t) = - 0.1 (T(t) -20) T'(t) + 0.1 (T(t) -20) = 0
# cooling_points = euler_method(lambda x, y: - 0.1 * (y - 20), [0, 95], [0,60], 61)
# #cooling_points = euler_method(lambda x, y: -sin(x), [0, 1], [0, 20], 100)
# print(cooling_points)
#
# fig, ax  = plt.subplots()
# ax.plot(cooling_points['x'], cooling_points['y'], 'x')
# expr = analytic_solution.function_cooling
# x = np.linspace(0, 60, 61)
# y = []
# for i in x:
#     y.append(expr(i))
# ax.plot(x, y)
