#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:42:51 2023

@author: slw
"""

import euler
import heun
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from timeit import timeit
import analytic_solution

#%% cosschwingung y'' = - y 
value_range = [0,60]
stepps = 100
rate_of_change = lambda t, y, y1: - y

euler.euler_method(rate_of_change, start_condition, value_range, stepps)