#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 19:32:02 2020

@author: Paco
"""

from __future__ import division
from sympy import *
from sympy.plotting import *
import math

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

from rk4 import *

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




def geodesics_eq_function(Gamma):
    dim = len(Gamma)
    eqs = list(symbols("eq1:{}".format(2*dim+1)))
    #xs = symbols("x1:{}".format(dim+1))
    dxs = list(symbols("dx1:{}".format(dim+1)))
    ddxs = list(symbols("ddx1:{}".format(dim+1)))
    
    for k in range(dim):
        eqs[k]=dvars[k]
    
    #print(eqs)
    for k in range(dim):
        el = 0
        for i in range(dim):
            for j in range(dim):
                el += Gamma[i][j][k]*dvars[i]*dvars[j]
                #print(el,i,j,k)
        eqs[dim+k]=(-el).simplify()
        
    return eqs
    


def geodesics_eq(Gamma):
    dim = len(Gamma)
    eqs = list(symbols("eq1:{}".format(dim+1)))
    #xs = symbols("x1:{}".format(dim+1))
    dxs = list(symbols("dx1:{}".format(dim+1)))
    ddxs = list(symbols("ddx1:{}".format(dim+1)))
    
    for k in range(dim):
        el = 0
        for i in range(dim):
            for j in range(dim):
                el += Gamma[i][j][k]*dvars[i]*dvars[j]
                #print(el,i,j,k)
        eqs[k]=(ddvars[k]+el).simplify()
    return eqs


def dotprod(l1,l2):
    Sum = 0
    
    for i in range(len(l1)):
        Sum = Sum + l1[i]*l2[i]
    
    return Sum

def putInTheValues(f,x,y):
    for i in range(len(f)):
        for j in range(len(x)):
            f[i] = f[i].subs(x[j],y[j])
    return f

    

dim = 2
g=[]

x,y,z,r,phi,theta,p = symbols("x y z r phi theta p")
dx,dy,ddx,ddy= symbols("dx dy ddx ddy")
variables = [x,y]

dvars = [dx,dy]
ddvars = [ddx,ddy]




for i in range(dim):
    g.append([])
    for j in range(dim):
        g[i].append(Symbol("g{}{}".format(i+1,j+1)))
        
    
        
g[0][0]=1/y**2
g[0][1]=0
g[1][0]=0
g[1][1]=1/y**2


invg = Matrix(g).inv()
invg = invg.tolist()

Gamma = []

for i in range(dim):
    Gamma.append([])
    for j in range(dim):
        Gamma[i].append([])
        for m in range(dim):
            el = 0
            for k in range(dim):
                el = invg[k][m]*(diff(g[j][k],variables[i]) + diff(g[k][i],variables[j]) - diff(g[i][j],variables[k]) ) + el
            Gamma[i][j].append((1/2*el).simplify())


# func = geodesics_eq_function(Gamma)

# f = lambda z: [i.subs({x:z[0],y:z[1],dx:z[2],dy:z[3]}) for i in func]
# #f = lambda z: putInTheValues(func.copy(),[x,y,dx,dy],z)

# (t,p)=rk4(f,[x,y,dx,dy],[0,1,1,10],(0,1), 100)
# p1 = [i[0] for i in p]
# p2 = [i[1] for i in p]
# plt.plot(p1,p2)

func = geodesics_eq_function(Gamma)

phi = np.linspace(-math.pi/2,math.pi/2,10)

x0s0=np.cos(phi)
x0s1=np.sin(phi)

tvars = variables+dvars
f = lambdify([tvars],func)
# f = lambda z: [i.subs({r:z[0],theta:z[1],dr:z[2],dtheta:z[3]}) for i in func]

# x0=[0,1,1,10]

# plotSolution(f, x0, (0,1), 1000)

for i in range(1,len(phi)-1):
    plotSolution(f,tvars, [0,1,x0s0[i],x0s1[i]], (0,10), 10**3)

for i in range(1,len(phi)-1):
    plotSolution(f,tvars, [6,5,x0s0[i],x0s1[i]], (0,10), 10**3)
    
        