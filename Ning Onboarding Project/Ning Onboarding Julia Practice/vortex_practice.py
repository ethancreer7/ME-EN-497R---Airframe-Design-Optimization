#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 17:38:57 2025

@author: ethancreer
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

d = 1.0
tao = np.array([0,0,1])
pi = np.pi
P1 = np.array([0.0,-0.5,0.0])
P2 = np.array([0.0,0.5,0.0])
P3 = np.array([1.0,0.5,0.0])
P4 = np.array([1.0,-0.5,0.0])


def radius_calc(p1, p2):
    return np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

def velocity_influence(tao, p1, p2):
    rad = p2 - p1
    cross_vector = np.cross(tao, rad)
    r = radius_calc(p1, p2)
    return cross_vector/(2*pi*r*r)

test1 = velocity_influence(tao,P1,P2)
test2 = velocity_influence(tao,P3,P2)
test3 = velocity_influence(tao,P4,P2)
final = test1 + test2 + test3
p1 = final*0.01 + P2
print(p1)