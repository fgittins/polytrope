# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import ode, simps
import matplotlib.pyplot as plt

import warnings
warnings.simplefilter('ignore', (RuntimeWarning, UserWarning))

class Polytrope(object):
    """
    A polytrope class. Largely used formulae from 
    https://ocw.mit.edu/courses/physics/8-282j-introduction-to-astronomy-spring-2006/assignments/compset2_soln.pdf
    and some of my own calculations.

    Parameters
    ----------

    n : float
        Polytropic index.

    Attributes
    -------

    xi_1 : float
        Dimensionless radius.
    m : float
        Dimensionless mass.
    i : float
        Dimensionless moment of inertia.
    Omega : float
        Dimensionless potential energy.
    """
    def __init__(self, n):
        self.n = n
        # to ensure the calculation is not done multiple times
        self._calculation = False
    
    @property
    def xi_1(self):
        """Dimensionless radius."""
        if not self._calculation:
            self._calculate()
        xi_1 = self._xis[-1]
        return xi_1
    
    @property
    def m(self):
        """Dimensionless mass."""
        if not self._calculation:
            self._calculate()
        m = 4*np.pi * self.xi_1**(-1) * np.abs(self._dtheta_dxis[-1])
        return m
    
    @property
    def i(self):
        """Dimensionless moment of inertia."""
        if not self._calculation:
            self._calculate()
        i = ((8*np.pi/3) * simps(self._thetas**n*self._xis**4, self._xis) 
                  / (self.m * self.xi_1**5))
        return i
    
    @property
    def Omega(self):
        """Dimensionless potential."""
        if not self._calculation:
            self._calculate()
        Omega = (-16*np.pi**2 
                      * simps(self._xis**3*self._thetas**self.n*
                              self._dtheta_dxis, 
                              self._xis)
                      / (self.m**2 * self.xi_1**5))
        return Omega

    def _calculate(self):
        """Run integration."""
        self._xis, self._thetas, self._dtheta_dxis = self._integrate()
        self._calculation = True

    def _LaneEmden(self, xi, y):
        """Defining the Lane-Emden equation."""
        theta, dtheta_dxi = y
        dy_dxi = np.array([dtheta_dxi, - theta**self.n - 2*dtheta_dxi/xi])
        return dy_dxi
    
    def _start(self, xi):
        """Starting point for integration."""
        theta = 1 - xi**2/6
        dtheta_dxi = - xi/3
        y = np.array([theta, dtheta_dxi])
        return y
    
    def _integrate(self):
        """Integrate the Lane-Emden equation using `ode`."""
        xi_0 = 1e-6
        y_0 = self._start(xi_0)

        xis, thetas, dtheta_dxis = [], [], []

        r = (ode(self._LaneEmden).set_integrator('dopri5')
                                 .set_initial_value(y_0, xi_0))
        dxi = 1e-3
        while r.successful() and r.y[0] > 0:
            r.integrate(r.t + dxi)
            xis.append(r.t)
            thetas.append(r.y[0])
            dtheta_dxis.append(r.y[1])
        return np.array(xis), np.array(thetas), np.array(dtheta_dxis)
    
    def plot_temperature(self):
        """Plot the temperature against radius."""
        plt.figure()
        plt.plot(self._xis/self.xi_1, self._thetas)
        plt.xlabel(r'Dimensionless radius, $\xi / \xi_1$')
        plt.ylabel(r'Dimensionless temperature, $\theta(\xi)$')
        plt.show()
    
    def plot_density(self):
        """Plot the density against radius."""
        plt.figure()
        plt.plot(self._xis/self.xi_1, self._thetas**self.n)
        plt.xlabel(r'Dimensionless radius, $\xi / \xi_1$')
        plt.ylabel(r'Dimensionless density, $\theta(\xi)^n$')
        plt.show()
