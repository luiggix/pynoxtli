#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Luis M. de la Cruz [Updated on jue abr  2 12:49:41 CST 2020]
"""
#-----------------------------------------------------------
# PARA DEFINIR EL PATH ABSOLUTO DE LOS MÓDULOS DE PYNOXTLI
#
import os, sys
if not("pynoxtli/base" in sys.path[0][-13:]):
    sys.path.insert(0, os.path.abspath('..'))
#-----------------------------------------------------------

import numpy as np
from fvm.numScheme import NumericalScheme

class sAdvDiff1D(NumericalScheme):
    
    def __init__(self, mesh, Su, Gamma = 1.0, rho = 1.0):
        super().__init__(mesh, Su)
        self.__Sp = 0            
        self.__Gamma = Gamma
        self.__rho = rho
        
    @property
    def Gamma(self):
        return self.__Gamma
    
    @Gamma.setter
    def Gamma(self, gamma):
        self.__Gamma = gamma
    
    def setVelocity(self, u):
        self.__u = u
        
    def Sp(self, Sp):
        self.__Sp = Sp * self.dx
    
    def source(self, phi):
        return self.Su * self.dx # You need to create a new array!! the mult create the new array
    
    def calc(self, i):
        dE = self.__Gamma / self.dx
        dW = self.__Gamma / self.dx
        dP = dE + dW - self.__Sp
#
# Average
#        
#        cE = -self.__rho * self.__u[i] * 0.5
#        cW =  self.__rho * self.__u[i-1] * 0.5
#
# Upwind
# 
        cE = max((-self.__u[i],0)) 
        cW = max((self.__u[i-1],0))
        cP = cE + cW + self.__rho * (self.__u[i] - self.__u[i-1])
         
        aE = dE + cE 
        aW = dW + cW
        aP = dP + cP

        return np.array([-aW, aP, -aE])

class sAdvDiff2D(NumericalScheme):
    
    def __init__(self, mesh, Su, Gamma = 1.0, rho = 1.0):
        super().__init__(mesh, Su)
        self.__Sp = 0            
        self.__Gamma = Gamma
        self.__rho = rho
        
    @property
    def Gamma(self):
        return self.__Gamma
    
    @Gamma.setter
    def Gamma(self, gamma):
        self.__Gamma = gamma

    def setVelocity(self, u):
        self.__u = u
        
    def Sp(self, Sp):
        self.__Sp = Sp * self.dx * self.dy
    
    def source(self, phi):
        return self.Su * self.dx * self.dy
    
    def calc(self, i, j):
        dE = self.__Gamma * self.dy / self.dx
        dW = self.__Gamma * self.dy / self.dx
        dN = self.__Gamma * self.dx / self.dy
        dS = self.__Gamma * self.dx / self.dy                
        dP = dE + dW + dN + dS #- self.__Sp

        cE = - self.__rho * self.__u[i  ,j  ] * 0.5
        cW =   self.__rho * self.__u[i-1,j  ] * 0.5
        cN = - self.__rho * self.__u[i  ,j  ] * 0.5
        cS =   self.__rho * self.__u[i  ,j-1] * 0.5
            # Upwind
 #           CE = max((-self.__u[i],0)) 
 #           CW = max((self.__u[i-1],0))
        cP = cE + cW + cN + cS + self.__rho * (self.__u[i,j] - self.__u[i-1,j] + self.__u[i,j] - self.__u[i,j-1])
         
        aE = dE + cE 
        aW = dW + dE
        aN = dN + cN
        aS = dS + cS
        aP = dP + cP

        return np.array([-aS, -aW, aP, -aE, -aN])
    
        
if __name__ == '__main__':

    from utils.displayInfo import printInfo   
    from geo.line import Line 
    
    N = 6

    rod = Line(1.0)
    malla = rod.constructMesh(N)
    ivx, _, _ = malla.bounds(bi = 1, ei = N-1)
    su = np.ones(ivx)
    vel = np.ones(ivx)
    laplace = sAdvDiff1D(malla, su)
    laplace.setVelocity(vel)
    printInfo(Descr = 'Testing Diffusion1D', dx = laplace.dx, vx = malla.vx)    
    print(laplace.calc(3))
    print(laplace.source(su))
    
#    from geo.Rectangle import Rectangle
#    cuadro = Rectangle(1,1)
#    malla2 = cuadro.constructMesh(N,N)
#    ivx, ivy, _ = malla2.bounds(bi = 1, ei = N-1,bj = 1, ej = N-1)
#    su = np.ones((ivy, ivx))
#    laplace2 = Diffusion2D(malla2, su) 
#    printTest(Descr = 'Testing Diffusion2D', hx = laplace2.hx, hy = laplace2.hy,
#              vx = malla.vx, vy = malla.vy)    
#    print(laplace2.calc(3,4))
#    print(laplace2.source(su))    
    
