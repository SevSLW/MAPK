#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 23:43:15 2023

@author: slw
"""

from sympy import *

#%% Flüssigkeit abkühlen T'(t) = - 0.1 (T(t) -20) T'(t) + 0.1 (T(t) -20) = 0

t = symbols('t')
T = Function('T')
cooling_constant = 0.1
room_temp = 20
start_condition = 95

function_cooling = dsolve(diff(T(t), t) + cooling_constant * (T(t) - room_temp))
print(function_cooling)
constant = solve(function_cooling.rhs.subs(t,0) - start_condition)
print(constant)
function_cooling_c = function_cooling.rhs.subs("C1", constant[0])
print(function_cooling_c)
function_cooling = lambdify(t, function_cooling_c)
print(function_cooling)
print(function_cooling(5))


#%% Feder gedämpft myii + cyi + ky = 0 m = 10, c = 1, k = 3 Anfang y(0) = 10, yi(0) = 0 
m = 10
c = 1
k = 3
x = symbols('x')
y = Function('y')
print()
print('Federauslenkun')
function_federposition = dsolve(m*diff(y(x), x, 2) + c * diff(y(x), x) + k * y(x))
print(function_federposition)
constants = solve([function_federposition.rhs.subs(x, 0) - 10,
                  diff(function_federposition.rhs, x).subs(x, 0)])
function_federposition = function_federposition.rhs.subs(constants)
print(constants)
print(function_federposition)
function_federposition = lambdify(x, function_federposition)



