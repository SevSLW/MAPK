#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 00:55:10 2023

@author: slw
"""

import numpy as np

space = np.linspace(0, 3, 50)

for i in space:
    print(i)

print(range(10))
for i in range(10):
    print(i)
    
a = np.arange(5, 400, 5)
print(a[:3])