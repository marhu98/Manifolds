#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 18:28:31 2020

@author: Paco
"""


from __future__ import division
from sympy import *
from sympy.plotting import *
import math

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

def plotSolution(f,tvars,x0,interval,steps):

    (t,p)=rk4(f,tvars,x0,interval, steps)
    p1 = [i[0] for i in p]
    p2 = [i[1] for i in p]
    plt.plot(p1,p2)

def listProd(l1,l2):
    l = []
    for i in range(len(l1)):
        l.append(l1[i]*l2[i])
    return l
        
def listSum(l1,l2):
    l = []
    for i in range(len(l1)):
        l.append(l1[i]+l2[i])
    return l

def listSum(*args):
    l = []
    for i in range(len(args[0])):
      l.append(sum([x[i] for x in args]))
    return l

def rk4(f,variables,x0,interval, steps):
    n = len(variables)
    t = [interval[0]]
    p = [x0]
    h = (interval[1]-interval[0])/steps
    
    for i in range(steps):
        t.append(t[-1]+h)
        
        currentPoint = p[-1]
        
        k1=f(currentPoint)
        
        k2=f(listSum(p[-1],listProd(n*[h/2],k1)))
        
        k3=f(listSum(p[-1],listProd(n*[h/2],k2)))
        
        k4=f(listSum(p[-1],listProd(n*[h],k3)))
        
        p.append(listSum(currentPoint,listProd(n*[1/6*h],listSum(k1,listProd(n*[2],k2),listProd(n*[2],k3),k4))))   
    return (t,p)


# x,y = symbols("x y")

# r=0.5

# f = lambda z: [r*z[0]*(1-z[0])]

# (t,p)=rk4(f,[x],[0.2],(0,10), 1000)
# x = [i[0] for i in p]
# plt.plot(t,x)


