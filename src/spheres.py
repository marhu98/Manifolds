#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 01:06:52 2020

@author: Paco
"""

from __future__ import division
from sympy import *
from sympy.plotting import *
import math


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
    

dim = 2
g=[]

t = Symbol("t")

x,y,z,r,phi,theta,p = symbols("x y z r phi theta p")
dphi,dtheta,ddphi,ddtheta = symbols("dphi dtheta ddphi ddtheta")
variables = [phi,theta]

phi0 = Symbol("phi0")

phi0 = 3*pi/4

r=1

x = r*cos(phi)*cos(theta)
y = r*cos(phi)*sin(theta)
z = r*sin(phi)

dvars = [dphi,dtheta]
ddvars = [ddphi,ddtheta]




p=Matrix([x,y,z])




#plot3d_parametric_surface(x, y, z*cte(theta), (phi, 0, 2*pi),(theta,-pi,pi))
#plot3d_parametric_line(x.subs(phi,phi0), y.subs(phi,phi0),z.subs(phi,phi0)+sin(theta)/10**6, (theta, -pi, pi))


for i in range(dim):
    g.append([])
    for j in range(dim):
        g[i].append(dotprod(diff(p,variables[i]),diff(p,variables[j])).simplify())
        

        
#g[0][0] = dotprod(diff(p,phi),diff(p,phi))
#g[0][1] = dotprod(diff(p,phi),diff(p,theta))
#g[1][0] = dotprod(diff(p,theta),diff(p,phi))
#g[1][1] = dotprod(diff(p,theta),diff(p,theta))

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


    
        