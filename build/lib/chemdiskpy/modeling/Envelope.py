#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_____________________________________________________________________________________________________________
file name: Envelope
last update: may 2021
language: PYTHON 3.8
short description:  geometrical model of a class 0/I envelope. 
                    Every input and output in au, except output for densities returned in g.cm-3. 
_____________________________________________________________________________________________________________
"""


import numpy as np
from .. constants.constants import mu, autocm, amu, Ggram, kb, M_sun

import parameters as p



#___________________________________________
#   _________   ___    __  ___        ___
#  |   ______| |   \  |  | \  \      /  /
#  |  |______  |    \ |  |  \  \    /  /
#  |   ______| |  |\ \|  |   \  \  /  /
#  |  |______  |  | \    |    \  \/  /
#  |_________| |__|  \___|     \____/
#___________________________________________
class Envelope:
    def __init__(self, rmin=p.rin, rmax=p.rmax, r_centri=p.r_centri, \
                       acc_rate=p.acc_rate, star_mass=p.star_mass, env_mass=p.env_mass, \
                       cavpl=p.cavpl, cav_fact=p.cav_fact, cavz0=p.cavz0):
        self.rmin=rmin
        self.rmax=rmax
        self.r_centri=r_centri
        self.acc_rate=acc_rate
        self.star_mass=star_mass
        self.env_mass=env_mass
        self.cavpl=cavpl
        self.cav_fact=cav_fact
        self.cavz0=cavz0









#---------------_________EXAMPLE OF HEADER FROM NUMPY
"""
Return the standard deviation of the array elements along the given axis.

Refer to `numpy.std` for full documentation.

See Also
--------
numpy.std

Notes
-----
This is the same as `ndarray.std`, except that where an `ndarray` would
be returned, a `matrix` object is returned instead.

"""





#_______________________________________________#
#                     CLASS                     #
#_______________________________________________#

#_____________Outer class: ENVELOPE_____________
#   _________   ___    __  ___        ___
#  |   ______| |   \  |  | \  \      /  /
#  |  |______  |    \ |  |  \  \    /  /
#  |   ______| |  |\ \|  |   \  \  /  /
#  |  |______  |  | \    |    \  \/  /
#  |_________| |__|  \___|     \____/
#_______________________________________________
class Envelope:
    def __init__(self, dist, mu, R_centri=p.R_centri, acc_rate=p.acc_rate, star_mass=p.star_mass):
        self.R_centri=R_centri
        self.acc_rate=acc_rate
        self.star_mass=star_mass



    """
	The following methods give the physical parameters related to the gas phase in the envelope.
	The gas is given by Hydrogren nuclei. Equation references are related to the documentation.
	Instance methods:
        A) solution (parameter for the analytical solution)
        B) mu_0()
        C) env_density(self, acc_rate, star_mass, dist, R_centri, mu, mu_0)
    """

    """ A)
    Equation (1). Initializer. Returns the parameter for the solution of the initial streamline of the infalling material 
    Arguments:
        -mu:            cos(theta) where theta is the angle above the midplane 
        -dist:          distance to the central core in the spherical coordinates (dist, theta, phi) [au]
        -R_centri:      centrifugal radius [au]
    """
    def solution(self, dist, mu, R_centri):
        sol = (27*mu*dist*R_centri**2+np.sqrt(729*mu**2*dist**2*R_centri**4+108*R_centri**3*(dist-R_centri)**3))**(1/3)/(3*R_centri*2**(1/3))
        #solution = 729*mu**2*dist**2*R_centri**4 + 108*R_centri**3*(dist-R_centri)**3
        return sol

    """ B)
    Equation (2). Returns the solution of the initial streamline of the infalling material (cos(theta_0))
    Arguments:
        -dist:          distance to the central core in the spherical coordinates (dist, theta, phi) [au]
        -R_centri:      centrifugal radius [au]
    """	
    def mu_0(self, dist, R_centri):
        solution = self.solution()
        return solution - ((1)/(3*solution))

    """ C)
    Equation (2). Returns the density of gas of the envelope [g.cm-3].
    Arguments:
        -acc_rate:      accretion rate in solar mass per year (converted in g per second)
        -dist:          distance to the central core in the spherical coordinates (dist, theta, phi) [au]
        -R_centri:      centrifugal radius [au]
        -mu:            cos(theta) where theta is the angle above the midplane
        -mu0:           initial streamline of the infalling material
    """	
    def env_density(self, acc_rate, star_mass, dist, R_centri, mu, mu_0):
        return (((acc_rate*M_sun)/3.154e7)/(4*np.pi))*(Ggram*star_mass*M_sun*(dist*autocm)**3)**(-1/2)*(1 + mu/mu_0)**(-1/2)*(mu/mu_0 + 2*mu_0**2*(R_centri/dist))**(-1)

    # #class Dust:
    #     """
    #     This Class gives the physical quantities related to the grains in the envelope.
    #     Equation references are those of the documentation.
    #     Instance methods:
	#         A) initializer ...
    #     """
        
