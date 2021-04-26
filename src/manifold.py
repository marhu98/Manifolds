#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 23:43:39 2020

@author: Paco
"""

from __future__ import division
from sympy import *
from sympy.plotting import *
from sympy.abc import *
import math



def mLatex(self,a):
    result = ""
    
    for eq in a:
        result += "$$"+latex(eq).replace("\\","~Ñ")+"$$"
    
    return result

def plotSolution(f,tvars,x0,interval,steps,**kwargs):    
    (t,p)=rk4(f,tvars,x0,interval, steps)
    p1 = [i[0] for i in p]
    p2 = [i[1] for i in p]
    plt.plot(p1,p2)
    
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



def listProd(l1,l2):
    l = []
    for i in range(len(l1)):
        l.append(l1[i]*l2[i])
    return l
        
def listDiff(l1,l2):
    l = []
    for i in range(len(l1)):
        l.append(l1[i]-l2[i])
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


class ManifoldPatch:
    
    def __init__(self,coordinates,**kwargs):
        """
    

        Parameters
        ----------
        f : TYPE
            DESCRIPTION.
        tvars : TYPE
            DESCRIPTION.
        x0 : TYPE
            DESCRIPTION.
        interval : TYPE
            DESCRIPTION.
        steps : TYPE
            DESCRIPTION.
        **kwargs : TYPE
            metric:Boolean.
                True if metric is given
    
        Returns
        -------
        None.
    
        """
        
        self.name = ""
        self.coords = coordinates
        self.vars=[]
        
        for i in coordinates:
            for j in i.free_symbols:
                if j not in self.vars:
                    self.vars.append(j)
        
        self.dim = len(self.vars)
        
        self.e=[]
        for i in range(self.dim):
            ei = self.dim*[0]
            ei[i]=1
            self.e.append(ei)
        
        self.xs = symbols("x1:{}".format(self.dim+1))
        
        self.vs = symbols("v1:{}".format(self.dim+1))
        self.dvs = symbols("dv1:{}".format(self.dim+1))
        
        self.toSubstitute = {}
        
        for i in range(self.dim):
            self.toSubstitute[self.vars[i]]=self.xs[i]
        
        self.dxs = []
        for i in range(self.dim):
            self.dxs.append(symbols("d{}".format(self.xs[i],self.dim+1)))
        self.ddxs = []
        for i in range(self.dim):
            self.ddxs.append(symbols("dd{}".format(self.xs[i],self.dim+1)))
            
        self.ddxs = symbols("ddx1:{}".format(self.dim+1))
        # self.xs = coordinates
        
        if "metric" not in kwargs:
            self.setFirstFundamentalForm()
        else:
            g = kwargs["metric"]
            self.invg = Matrix(self.g).inv()
            self.invg = self.invg.tolist()
        self.setChrystoffelSymbols()
        
        self.geodesic = self.geodesics_eq()
        self.geodesic_func = self.geodesics_eq_function()
        
        self.parallelTransportEq = self.parallelTransportEquations()
        
        self.R = self.setCurvature()
        self.Ricc = self.setRicciCurvature()
        self.S = self.setScalarCurvature()
        
    def setCurvature(self):
        R = []
        for i in range(self.dim):
            R.append([])
            for j in range(self.dim):
                R[i].append([])
                for k in range(self.dim):
                    R[i][j].append([])
                    for s in range(self.dim):
                        diffTerm = 0
                        # print(i,j,k,s)
                        if self.Gamma[i][k][s]!=0:
                            # print(self.Gamma[i][k][s])
                            diffTerm = diff(self.Gamma[i][k][s],self.vars[j])
                        if self.Gamma[j][k][s]!=0:
                            # print(self.Gamma[i][k][s])
                            diffTerm -= diff(self.Gamma[j][k][s],self.vars[i])
                        
                        buff = 0
                        
                        for l in range(self.dim):
                            buff += self.Gamma[i][k][l]*self.Gamma[j][l][s]
                            buff -= self.Gamma[j][k][l]*self.Gamma[i][l][s]
                        
                        R[i][j][k].append((buff+diffTerm).simplify())
                        
        
        return R
    
    def setRicciCurvature(self):
        Ricc = []
        for i in range(self.dim):
            Ricc.append([])
            for j in range(self.dim):
                buff = 0
                
                for k in range(self.dim):
                    buff += self.R[k][i][j][k]
                
                Ricc[i].append(buff.simplify())
        return Ricc
    
    def setScalarCurvature(self):
        S = 0
        
        for j in range(self.dim):
            for m in range(self.dim):
                S += self.invg[j][m]*self.Ricc[m][j]
        
        S = trigsimp(expand_trig(S.simplify()))
        
        return S
                
    
    def connection(self,X,Y):
        """
        Takes the coordinates of two vector fields
        in the form of two lists
        """
        Z = self.dim*[0]
        
        for k in range(self.dim):
            for i in range(self.dim):
                for j in range(self.dim):
                    
                    # print(Y[j])
                    # diffTerm = Y[j]
                    
                    Z[k] += X[i](self.vars[0])*Y[j](self.vars[0])*self.Gamma[i][j][k]
                    
                    Z[k]=Z[k].simplify()
                    
        for i in range(self.dim):
            for j in range(self.dim):
                diffTerm = Function("diffTerm")
                diffTerm = diff(Y[k](self.vars[i]),self.vars[i])
                Z[j] += X[i](self.vars[0])*diffTerm
                Z[j]=Z[j].simplify()
        return Z
        
    def bracket(self,X,Y):
        return listDiff(self.connection(X,Y),self.connection(Y,X))
    
    def setName(self,name):
        self.name = name
        
        
    def __str__(self):
        result = ""
        if self.name != "":
            result += "Manifold name: {}\n\n\n".format(self.name)
        result += "Coordinates\n"
        
        for i in range(len(self.xs)):
            result += "x_{} : {}\n".format(i+1,self.xs[i])   
        
        result += "\n"
            
        result += "Metric\n"
        for i in range(self.dim):
            for j in range(self.dim):
                result += "g_{}{} : {}\n".format(i+1,j+1,self.g[i][j])  
                
        result += "\n"
                
        result += "Chrystoffel Symbols\n"
        for k in range(self.dim):
            for i in range(self.dim):
                for j in range(self.dim):
                    result += "Gamma_{}{}^{} : {}\n".format(i+1,j+1,k+1,self.Gamma[i][j][k]) 
                    
        result += "\n"
        
        result += "Geodesics equations\n"
        for eq in self.geodesic:
            result += "{} \n".format(eq)
        
        return result
    
    def setFirstFundamentalForm(self):
        g = []
        # print(self.vars)
        # print(self.coords)
        p = Matrix(self.coords)
        for i in range(self.dim):
            g.append([])
            for j in range(self.dim):
                g[i].append(Symbol("g{}{}".format(i+1,j+1)))
        for i in range(self.dim):
            for j in range(self.dim):
                g[i][j]=self.dotprod(diff(p,self.vars[i]),diff(p,self.vars[j])).simplify()
        
        self.g = g
        
        self.invg = Matrix(self.g).inv()
        self.invg = self.invg.tolist()
        
    def norm(self,X,variable):
        """
        Works out the norm of the tangent vector
        field along a curve

        """
        Y = [z.diff(variable) for z in X]
        
        res =sqrt(self.dotProduct(Y,Y)).simplify()
        
        s = {}
        
        for i in range(self.dim):
            s["{}".format(self.vars[i])]=X[i]
        res = res.subs(s).simplify()
        res = trigsimp(expand_trig(res))
        
        return res
        
    def dotProduct(self,X,Y):
        res = 0
        
        for i in range(self.dim):
            for j in range(self.dim):
                # res += self.g[i][j]*X[i](self.vars[0])*Y[j](self.vars[0])
                res += self.g[i][j]*X[i]*Y[j]
        return res
        
    def setChrystoffelSymbols(self):
        self.Gamma = self.getChrystoffelSymbols()
        
        
    def getFirstFundamentalForm(self):   
        return self.g
            
        
    def getChrystoffelSymbols(self):
        Gamma = []

        for i in range(self.dim):
            Gamma.append([])
            for j in range(self.dim):
                Gamma[i].append([])
                for m in range(self.dim):
                    el = 0
                    for k in range(self.dim):
                        el = self.invg[k][m]*(diff(self.g[j][k],self.vars[i]) + diff(self.g[k][i],self.vars[j]) - diff(self.g[i][j],self.vars[k]) ) + el
                    Gamma[i][j].append((1/2*el).simplify())
        return Gamma
        
    def geodesics_eq_function(self):
        dim = len(self.Gamma)
        eqs = list(symbols("eq1:{}".format(2*dim+1)))
        
        
        for k in range(self.dim):
            eqs[k]=self.dxs[k]
        
        for k in range(self.dim):
            el = 0
            for i in range(self.dim):
                for j in range(self.dim):
                    el += self.Gamma[i][j][k].subs(self.toSubstitute)*self.dxs[i]*self.dxs[j]
                    #print(el,i,j,k)
            eqs[dim+k]=(-el).simplify()
            
        return lambdify([self.vars+self.dxs],eqs)
        
    
    
    def geodesics_eq(self):
        dim = len(self.Gamma)
        eqs = list(symbols("eq1:{}".format(dim+1)))
        
        for k in range(dim):
            el = 0
            for i in range(dim):
                for j in range(dim):
                    el += self.Gamma[i][j][k].subs(self.toSubstitute)*self.dxs[i]*self.dxs[j]
            eqs[k]=(self.ddxs[k]+el).simplify()
        return eqs
    
    def parallelTransportEquations(self):
        eqs = self.dim*[0]
        
        for k in range(self.dim):
            for i in range(self.dim):
                for j in range(self.dim):
                    eqs[k] += self.Gamma[i][j][k].subs(self.toSubstitute)*self.vs[j]*self.dxs[i]
            eqs[k] += self.dvs[k]
            eqs[k]=eqs[k].simplify()
        return eqs
    def getParallelTransportEquationsFunction(self,xs):
        eqs = self.dim*[0]
        
        for k in range(self.dim):
            for i in range(self.dim):
                for j in range(self.dim):
                    eqs[k] += -self.Gamma[i][j][k].subs(self.toSubstitute)*self.vs[j]*diff(xs[i])
            eqs[k]=eqs[k].simplify()
        return eqs
    
    def dotprod(self,l1,l2):
        Sum = 0
        
        for i in range(len(l1)):
            Sum = Sum + l1[i]*l2[i]
        
        return Sum
      
        
class PlanePolarCoordinates(ManifoldPatch):
    def __init__(self):
        r,theta = symbols("r theta")
        super().__init__([r*cos(theta),r*sin(theta)])
        self.setName("Plane with Polar Coordinates")
        
class S2(ManifoldPatch):
    def __init__(self,r,**kwargs):
        if "stereo" in kwargs and kwargs["stereo"]==True:
            u,v = symbols("u v")
            den = 1+u**2+v**2
            super().__init__([2*u/den,2*v/den,(-1+u**2+v**2)/den])            
        else:
            theta,phi = symbols("theta phi")
            super().__init__([r*cos(phi)*cos(theta),r*sin(phi)*cos(theta),r*sin(theta)])
        self.setName("2-Sphere")
        
class S3(ManifoldPatch):
    def __init__(self,r):
        theta,phi,tau = symbols("theta phi tau")
        super().__init__([r*cos(phi)*cos(theta)*cos(tau),r*sin(phi)*cos(theta)*cos(tau),r*sin(theta)*cos(tau),r*sin(tau)])
        self.setName("3-Sphere")
 


class PlaneWithMetric(ManifoldPatch):
    
    def __init__(self,g):
        theta,phi,x,y = symbols("theta phi x y")
        
        self.g = g
        
        
        self.coords = [x,y]
        super(PlaneWithMetric,self).__init__(self.coords,metric = g)
        
        # self.vars=[]
        
        # for i in self.coords:
        #     for j in i.free_symbols:
        #     # self.vars.update(i.free_symbols)
        #         if j not in self.vars:
        #             self.vars.append(j)
        
        # self.dim = len(self.vars)
        
        # self.xs = symbols("x1:{}".format(self.dim+1))
        
        # self.dxs = []
        # for i in range(self.dim):
        #     self.dxs.append(symbols("d{}".format(self.vars[i],self.dim+1)))
        # self.ddxs = []
        # for i in range(self.dim):
        #     self.ddxs.append(symbols("dd{}".format(self.vars[i],self.dim+1)))
            
        # self.ddxs = symbols("ddx1:{}".format(self.dim+1))
        # self.xs = self.coords
        
        # self.g = g
        
        # self.invg = Matrix(self.g).inv()
        # self.invg = self.invg.tolist()
        
        # self.setChrystoffelSymbols()
        # self.geodesic = self.geodesics_eq()

class Plane(PlaneWithMetric):
    
    def __init__(self):
        theta,phi,x,y = symbols("theta phi x y")
        
        g = [2*[0],2*[0]]
        g[0][0]=1
        g[0][1]=0
        g[1][0]=0
        g[1][1]=1
        
        
        super().__init__(g)
        
        self.setName("Normal Plane")
       
class Lobatchevsky(PlaneWithMetric):
    
    def __init__(self):
        theta,phi,x,y = symbols("theta phi x y")
        
        g = [2*[0],2*[0]]
        g[0][0]=1/y**2
        g[0][1]=0
        g[1][0]=0
        g[1][1]=1/y**2
        
        
        super().__init__(g)
        
        self.setName("Lobatchevsky Plane")
    def simpConn(self, i, j ):
        res = []
        for k in range(2):
            res.append(self.Gamma[i-1][j-1][k])
        return res
# myPlane = PlanePolarCoordinates()
# print(myPlane)

a = [Function("a{}".format(i)) for i in range(2)]
b = [Function("b{}".format(i)) for i in range(2)]

# mSphere = S2(1)
# print(mSphere)

# # my3Sphere = S3(1)
# # print(my3Sphere)

plane = Plane()

myLobat = Lobatchevsky()

l = [[myLobat.simpConn(i+1,j+1)]for j in range(2) for i in range(2)]
#print(myLobat.Gamma)
# print("")
# print(myLobat.geodesic)
# print("")
#print(myLobat.parallelTransporteq)
#print(myLobat)