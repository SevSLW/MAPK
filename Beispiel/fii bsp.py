#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:45:19 2023

@author: slw
"""

from pylab import*
import matplotlib.pyplot as plt
import numpy as np

# we are approximating the solution of y'' = f(x,y,y') for x in [x_0, x_1] satisfying the Cauchy condition of order 2:
# y(x_0) = y0 and y'(x_0) = y1

def f(x, y_d_0, y_d_1):
    print(y_d_1)
    return -y_d_0

# here f defines the equation y'' = -y

def explicit_euler(x0, x1, y0, y1, N,):
    # The following formula relates h and N
    h = (x1 - x0)/(N+1)

    xd = list()
    yd = list()

    xd.append(x0)
    # to allow group operations in R^2, we use the numpy library
    yd.append(np.array([y0, y1]))

    for i in range (1,N+1) :
        # We use the explicite Euler scheme y_{i+1} = y_i + h * f(x_i, y_i)
        # remember that now, yd is a list of vectors

        # the equivalent order 1 equation is [y, y']' = [y', f(x,y,y')]
        y = yd[-1] + h * np.array([yd[-1][1], f(xd[-1], yd[-1][0], yd[-1][1])])  # vector of dimension 2
        # print(y)
        # you can replace the above scheme by any other (R-K 4 for example !)

        x = xd[-1] + h  # vector of dimension 1

        yd.append(y)
        xd.append(x)

    return xd, yd

x0 = 0
x1 = 30
y0 = 1
y1 = 0

# the only function satisfying y(0) = 1, y'(0) = 0 and y'' = -y is y(x) = cos(x)
N = 5000
xd, yd =explicit_euler(x0, x1, y0, y1, N)
# I only want the first variable of yd
yd_1 = list(map(lambda y: y[0], yd))

plt.plot(xd,yd_1)
plt.show()

from pylab import*
import matplotlib.pyplot as plt

# we are approximating the solution of y' = f(x,y) for x in [x_0, x_1] satisfying the Cauchy condition y(x_0) = y0

def f(x, y0):
    return y0

# here f defines the equation y' = y

def explicit_euler(x0, x1, y0, N,):
    # The following formula relates h and N
    h = (x1 - x0)/(N+1)

    xd = list()
    yd = list()

    xd.append(x0)
    yd.append(y0)

    for i in range (1,N+1) :
        # We use the explicite Euler scheme y_{i+1} = y_i + h * f(x_i, y_i)
        y = yd[-1] + h * f(xd[-1], yd[-1])
        # you can replace the above scheme by any other (R-K 4 for example !)
        x = xd[-1] + h

        yd.append(y)
        xd.append(x)

    return xd, yd



N = 250
x1 = 5
x0 = 0
y0 = 1
# the only function which satisfies y(0) = 1 and y'=y is y(x)=exp(x).
xd, yd =explicit_euler(x0, x1, y0, N)
plt.plot(xd,yd)
plt.show()
# this plot has the right shape !



from pylab import*
import matplotlib.pyplot as plt
import numpy as np

# we are approximating the solution of y'' = f(x,y,y') for x in [x_0, x_1] satisfying the Cauchy condition of order 2:
# y(x_0) = y0 and y'(x_0) = y1

def f(x, y_d_0, y_d_1):
    m = 10.0
    c = 1.0
    k = 3.0
    return (-c * y_d_1 - k * y_d_0) / m

# here f defines the equation y'' = -y

def explicit_euler(x0, x1, y0, y1, N,):
    # The following formula relates h and N
    h = (x1 - x0)/(N+1)

    xd = list()
    yd = list()

    xd.append(x0)
    # to allow group operations in R^2, we use the numpy library
    yd.append(np.array([y0, y1]))

    for i in range (1,N+1) :
        # We use the explicite Euler scheme y_{i+1} = y_i + h * f(x_i, y_i)
        # remember that now, yd is a list of vectors

        # the equivalent order 1 equation is [y, y']' = [y', f(x,y,y')]
        y = yd[-1] + h * np.array([yd[-1][1], f(xd[-1], yd[-1][0], yd[-1][1])])  # vector of dimension 2
        # print(y)
        # you can replace the above scheme by any other (R-K 4 for example !)

        x = xd[-1] + h  # vector of dimension 1

        yd.append(y)
        xd.append(x)

    return xd, yd

x0 = 0
x1 = 30
y0 = 1
y1 = 0

# the only function satisfying y(0) = 1, y'(0) = 0 and y'' = -y is y(x) = cos(x)
N = 5000
xd, yd =explicit_euler(x0, x1, y0, y1, N)
# I only want the first variable of yd
yd_1 = list(map(lambda y: y[0], yd))

plt.plot(xd,yd_1)
plt.show()
