#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Engine for Numerical Schrodinger Equation Solver

    # Copyright (C) 2017 Zhang Chang-kai #
    # Contact via: phy.zhangck@gmail.com #
    # General Public License version 3.0 #

'''Numerical Schrodinger Equation Solve Engine'''

# Principle of Calculation
# ========================
#
# This program solves 1D Schrodinger equation through brute-force integrat-
# ion by Runge-Kutta forth order numetrical integration method. The Schrod-
# inger equation reads
#                    iħ ∂ψ/∂t = - ħ²/2m ∂²ψ/∂x² + Vψ
# In this program, physical constants are neglected. Thus, the actual equa-
# tion solved is
#                     i ∂ψ/∂t = - ∂²ψ/∂x² + Vψ
# To solve this equation, Runge-Kutta method is used to estimate the value
# of the wave function at each time step.
#


# import dependencies
import numpy as np
from config import *


# class of methods for Schrodinger equation Solver
# read the configuration from config.py

class Schdger(object):

    # initialization according to configuration from config.py
    def __init__(self):

        # spacial step size
        self.dx = spc_size
        # temporal step size
        self.dt = tmp_size
        # total spacial steps
        self.nx = spc_steps
        # total temporal steps
        self.nt = int(rec_num)

        # integration parameter
        self.r = 1j * self.dt / self.dx**2

        # set coordinates
        self.x = np.linspace(spc_range[0], spc_range[1], self.nx)

        # initialize wave function ψ(t,x)
        self.psi = np.empty([self.nt, self.nx], dtype=complex)
        # initialize probability density
        self.pbd = np.empty([self.nt, self.nx], dtype=float)

        # set potential
        self.potential = potential(self.x)

        # initial condition
        self.psi[0] = initial(self.x)
        self.pbd[0] = np.power(self.psi[0].real, 2) \
                    + np.power(self.psi[0].imag, 2)

        # boundary conditions
        self.periodic = True if bound == "periodic" else False
        self.boundary = bound if not self.periodic else None
        
    
    # set time step size
    def _setTimeSize(self, newnt):
        self.dt = newnt


    # generate Laplacian
    def _genLaplace(self, psi):
        
        # initialize derivative
        laplace = np.empty_like(psi)
        # central difference
        laplace[1:-1] = (psi[2:] - 2*psi[1:-1] + psi[0:-2])
        # dirichlet boundary condition
        laplace[0] = laplace[-1] = 0

        # return Laplacian
        return laplace / self.dx ** 2


    # generate boundary conditions
    def _genBoundary(self, i):

        if not self.periodic:
            self.psi[i, 0] = self.boundary[0]
            self.psi[i, -1] = self.boundary[1]


    # solve the equation
    def solveEq(self):
        
        # integration
        for i in range(self.nt-1):
            
            # set boundary conditions
            self._genBoundary(i)
            
            # generate Hamiltonian
            hamiltonian = lambda x: -self._genLaplace(x)
            
            # Runge-Kutta methods
            first = -1j * hamiltonian(self.psi[i])
            secnd = -1j * (hamiltonian(self.psi[i] + self.dt * first / 2))
            third = -1j * (hamiltonian(self.psi[i] + self.dt * secnd / 2))
            forth = -1j * (hamiltonian(self.psi[i] + self.dt * third))

            # estimation
            self.psi[i+1] = self.psi[i] + \
                self.dt * (first + 2 * secnd + 2 * third + forth) / 6
            self.pbd[i+1] = abs(self.psi[i+1])**2


# End of Engine for Numerical Schrodinger Equation Solver
