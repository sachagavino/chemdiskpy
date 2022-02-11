#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_____________________________________________________________________________________________________________
file name: geometrical_model
last update: may 2021
language: PYTHON 3.8
short description:  geometrical model of a class 0/I object consisting of the disk and the envelope. 
                    Every input and output in au, except output for densities returned in g.cm-3. 
_____________________________________________________________________________________________________________
"""


import numpy as np
from .. constants.constants import mu, autocm, amu, Ggram, kb, M_sun

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
    """Outer Class"""	

    class Gas:
        """
	    This Class gives the physical parameters related to the gas phase in the envelope.
	    The gas density is given by number of Hydrogen nuclei. Equation references are related to the documentation.
	    Instance methods:
		    A) initializer (parameter for the analytical solution)
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
        def __init__(self, dist, mu, R_centri):
            solution = (27*mu*dist*R_centri**2+np.sqrt(729*mu**2*dist**2*R_centri**4+108*R_centri**3*(dist-R_centri)**3))**(1/3)/(3*R_centri*2**(1/3))
            #solution = 729*mu**2*dist**2*R_centri**4 + 108*R_centri**3*(dist-R_centri)**3
            self.solution = solution

        """ B)
        Equation (2). Returns the solution of the initial streamline of the infalling material (cos(theta_0))
        Arguments:
            -dist:          distance to the central core in the spherical coordinates (dist, theta, phi) [au]
            -R_centri:      centrifugal radius [au]
        """	
        def mu_0(self, dist, R_centri):
            return self.solution - ((1)/(3*self.solution))

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
        
