#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 17:38:57 2025

@author: ethancreer
"""

import numpy as np
import matplotlib.pyplot as plt

dt = 0.01
time = 0
tao1 = np.array([0,0,1])
tao2 = np.array([0,0,-1])
pi = np.pi


P1 = np.array([0.0,-0.5,0.0])
P2 = np.array([0.0,0.5,0.0])
P3 = np.array([1.0,0.5,0.0])
P4 = np.array([1.0,-0.5,0.0])

V1 = np.array([0,0,0])
V2 = np.array([0,0,0])
V3 = np.array([0,0,0])
V4 = np.array([0,0,0])


vortex1 = np.array([0.0,-0.5,0.0])
vortex2 = np.array([0.0,0.5,0.0])
vortex3 = np.array([1.0,0.5,0.0])
vortex4 = np.array([1.0,-0.5,0.0])

def radius_calc(p1, p2):
    return np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

def velocity_influence(p1, p2):
    rad = p2 - p1
    if p1[1] < 0:
        cross_vector = np.cross(tao2, rad)
    elif p1[1] > 0:
        cross_vector = np.cross(tao1, rad)
    else:
        raise Exception("Sorry, your tao logic is incorrect.")
    r = radius_calc(p1, p2)
    return cross_vector/(2*pi*r*r)

def vortex_1():
    inf2 = velocity_influence(P2,P1)
    inf3 = velocity_influence(P3,P1)
    inf4 = velocity_influence(P4,P1)
    final_inf = inf2 + inf3 + inf4
    return final_inf * dt + P1

def vortex_2():
    inf1 = velocity_influence(P1,P2)
    inf3 = velocity_influence(P3,P2)
    inf4 = velocity_influence(P4,P2)
    final_inf = inf1 + inf3 + inf4
    return final_inf * dt + P2
    

def vortex_3():
    inf1 = velocity_influence(P1,P3)
    inf2 = velocity_influence(P2,P3)
    inf4 = velocity_influence(P4,P3)
    final_inf = inf1 + inf2 + inf4
    return final_inf * dt + P3

def vortex_4():
    inf1 = velocity_influence(P1,P4)
    inf2 = velocity_influence(P2,P4)
    inf3 = velocity_influence(P3,P4)
    final_inf = inf1 + inf2 + inf3
    return final_inf * dt + P4


def calc_final_positions():
    final1 = vortex_1()
    final2 = vortex_2()
    final3 = vortex_3()
    final4 = vortex_4()
    return final1, final2, final3, final4
    

while time < 40:
    P1, P2, P3, P4 = calc_final_positions()
    vortex1 = np.vstack([vortex1,P1])
    vortex2 = np.vstack([vortex2,P2])
    vortex3 = np.vstack([vortex3,P3])
    vortex4 = np.vstack([vortex4,P4])
    time += 0.01
    
plt.plot(vortex1[:,0],vortex1[:,1],color='orange')
plt.plot(vortex2[:,0],vortex2[:,1],color='orange')
plt.plot(vortex3[:,0],vortex3[:,1],color='purple')
plt.plot(vortex4[:,0],vortex4[:,1], color='purple')
plt.yticks([0])
plt.xticks(np.arange(0.0,12.5,2.5))
plt.tight_layout()
plt.show()


xdata1 = vortex1[:,0]
xdata2 = vortex2[:,0]
xdata3 = vortex3[:,0]
xdata4 = vortex4[:,0]
ydata1 = vortex1[:,1]
ydata2 = vortex2[:,1]
ydata3 = vortex3[:,1]
ydata4 = vortex4[:,1]




# "animation that will take forever below"
# for i in range(4000):
#     plt.plot(xdata1[:i],ydata1[:i], color='orange')
#     plt.plot(xdata1[i],ydata1[i], color='orange', markersize=10, marker='o')
#     plt.plot(xdata2[:i],ydata2[:i],color='orange')
#     plt.plot(xdata2[i],ydata2[i], color = 'orange', markersize=10, marker='o')
#     plt.plot(xdata3[:i],ydata3[:i],color='purple')
#     plt.plot(xdata3[i],ydata3[i], color = 'purple', markersize=10,marker='o')
#     plt.plot(xdata4[:i],ydata4[:i],color='purple')
#     plt.plot(xdata4[i],ydata4[i], color = 'purple', markersize=10,marker='o')
#     plt.pause(0.01)
   