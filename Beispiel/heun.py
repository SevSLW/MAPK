#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 03:44:50 2023

@author: slw
"""

import analytic_solution
import numpy as np

def heun_method(ode, start_condition: list, value_range: list, stepps: int):
    x_list = [start_condition[0]]
    y_list = [start_condition[1]]
    x = start_condition[0]
    x = 0
    y = start_condition[1]
    step_width = value_range[1] / stepps
    for i in np.linspace(value_range[0], value_range[1], stepps):
        if i == 0:
            continue
        x = i
        changerate_current = ode(x, y)
        y_new = y + step_width * changerate_current
        x_new = x + step_width
        changerate_new = ode(x_new, y_new)
        changerate_mean = (changerate_current + changerate_new) / 2
        y += step_width * changerate_mean
        x_list.append(x)
        y_list.append(y)
    return {'x': x_list,
            'y': y_list}

def heun_ode2(ode, start_condition: list, value_range: list, stepps: int):
    x_list = [start_condition[0]]
    y_list = [np.array([start_condition[1], start_condition[2]])]
    y0_list = [y_list[-1][0]]
    x = start_condition[0]
    y = y_list[-1]
    step_width = value_range[1] / stepps
    for i in np.linspace(value_range[0], value_range[1], stepps):
        if i == 0:
            continue   
        x = i
        changerate_current = ode(x, y_list[-1][0], y_list[-1][1])
        y_new = y_list[-1] + step_width * np.array([y_list[-1][1], changerate_current])
        x_new = x + step_width
        changerate_new = ode(x_new, y_new[0], y_new[1])
        changerate_mean = (changerate_current + changerate_new) / 2
        y = y_list[-1] + step_width * np.array([y_list[-1][1], changerate_mean])
        x_list.append(x)
        y_list.append(y)
        y0_list.append(y[0])
    return x_list, y0_list


        
        