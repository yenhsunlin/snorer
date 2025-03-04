# Created by Yen-Hsun Lin (Academia Sinica) in 04/2024.
# Copyright (c) 2024 Yen-Hsun Lin.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version (see <http://www.gnu.org/licenses/>).
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

__all__ = ['Geometry',
           'Propagation',]


#---------- Import required utilities ----------#

from numpy import sqrt,cos,sin,clip,finfo
from constants import constant



##########################################################################
#                                                                        #
#   General Classes and Functions for Numerics                           #
#                                                                        #
##########################################################################


"""
This script contains following classes and functions

Classes
------
1. Geometry

Functions
------
None

The docstrings should be sufficient for their self-explanations
"""


class Geometry:
    """
    The class constructs the static geomatrical relations for d, rprime and cos(psi)
    when (l,theta,phi) and (Rs,Re,beta) are specified. 

    /*-----------------------------------------------------------------------------*/

    ********************
    *                  *
    *   Class Inputs   *
    *                  *
    ********************
    
        l: The line-of-sight distance, kpc
    theta: The zenith angle at Earth, centers SN, rad
      phi: Azimuthal angle at Earth, centers SN, rad
       Rs: Distance from Earth to SN, kpc
       Re: Distance from Earth to GC, kpc
     beta: Off-center angle, rad


    ********************
    *                  *
    *    Attributes    *
    *                  *
    ********************

    When a Geometry instance is initialized, the following attributes will be assigned,

          d: Distance from SN to boost point, kpc
     rprime: Distance from GC to boost point, kpc
    cos_psi: cos(psi) at boost point where psi is the direction for BDM at B pointing
             Earth
    
    Suppose,
    
    >>> l,theta,phi,Rs,Re,beta = 5.160e-9,1e-4,0,11.6,8.5,0.71
    >>> bdmGeo = Geometry(l,theta,phi,Rs,Re,beta)
    
    Then the lengths and cos(psi) in the propagation geometry are specified by 
    
    >>> print(bdmGeo.d)
    11.59999999483988
    >>> print(bdmGeo.rprime)
    8.49999999608676
    >>> print(bdmGeo.cos_psi)
    0.9999999721078604

    See API/Propagation/Positioning for detail
    """

    def __init__(self,l,theta,phi,Rs,Re,beta):
        self.l = l 
        self.theta = theta
        self.phi = phi
        self.Rs = Rs
        self.Re = Re
        self.beta = beta

        # Evaluated quantities
        self._h = self.l*sin(theta)
        self._b = self.l*cos(theta)
        self._rho2 = self.get_rho2(self._b)
        self._cos_delta = self.get_cos_delta()
        self._a2 = self.get_a2()
        self._d2 = self.get_d2()
        
        # Following will be assigned as class attributes
        self._rprime = sqrt(self.get_rprime2())
        self._d = sqrt(self._d2)
        self._cos_psi = self.get_cos_psi()
        
    @property
    def rprime(self):
        return self._rprime
    
    @property
    def d(self):
        return self._d
    
    @property
    def cos_psi(self):
        return self._cos_psi

    # --- The following geometrical quantities are evaluated in squared terms
    def get_rho2(self,b):
        """
        Get rho-squared, see Eq. (5) in document
        """
        Re, beta = self.Re, self.beta
        rho2 = b**2 + Re**2 - 2 * b * cos(beta)
        rho2 = clip(rho2,0,None) # restricts rho-squared to be positive
        return rho2
    
    def get_cos_delta(self):
        """
        Get cos(delta), see Eq. (6) in document
        """
        Re, rho2, b = self.Re, self._rho2, self._b
        rho = sqrt(rho2)
        cos_delta = (Re**2 - rho2 - b**2)/2/rho/b
        return clip(cos_delta,-1,1) # restrics sin value in [-1,1]
    
    def get_a2(self):
        """
        Get a-squared, see Eq. (4) in document
        """
        rho2 = self._rho2
        rho = sqrt(rho2)
        h = self._h
        cos_delta = self._cos_delta
        sin_delta = sqrt(1 - cos_delta**2)
        sin_phi = sin(self.phi)
        a2 = rho2 + h**2 * sin_phi**2 - 2 * rho * h * sin_delta * sin_phi
        return clip(a2,0,None) # restricts rho-squared to be positive
    
    def get_rprime2(self):
        """
        Get r'-squared, see Eq. (3) in document
        """
        a2, h2 = self._a2, self._h**2
        cos_phi = cos(self.phi)
        rprime2 = a2 + h2 * cos_phi**2
        return clip(rprime2,0,None) # restricts rho-squared to be positive
    
    def get_d2(self):
        """
        Get d-squared, see Eq. (8) in document
        """
        l, Rs, theta = self.l, self.Rs, self.theta
        d2 = l**2 + Rs**2 - 2 * l * Rs * cos(theta)
        return clip(d2,0,None) # restricts rho-squared to be positive
    
    def get_cos_psi(self):
        """
        Get cos(psi), see Eq. (7) in document
        """
        Rs, d2, l = self.Rs, self._d2, self.l
        d = self.d
        cos_psi = (Rs**2 - d2 - l**2)/2/d/l
        return clip(cos_psi,-1,1)


class Propagation(Geometry):
    """
    Superclass: Geometry

    The class constructs the dynamical geomatrical relations for d, rprime and cos(psi)
    when (t,vx,theta,phi) and (Rs,Re,beta) are specified.

    Unlike its superclass Geometry, the class parameter l is now replaced by a specific
    time t and dimensionless BDM velocity vx. This allows it to incorporate time-dependent
    feature when evaluating the geometrical quantities during propagation.  

    This class is also not exclusively for SN in MW or LMC, it can be generalized to SN
    in arbitrary distant galaxy as long as the aforementioned inputs are determined.
    The BDM emissivity along the line-of-sight then can be determined when calculate
    the BDM flux and event at Earth associated to that particular SN.

    /*-----------------------------------------------------------------------------*/

    ********************
    *                  *
    *   Class Inputs   *
    *                  *
    ********************
    
        t: The BDM at specific time, seconds
           Time-zero is set to be the arrival of SNnu at Earth
       vx: BDM dimesionless velocity, in the unit of c
    theta: The zenith angle at Earth, rad
      phi: Azimuthal angle at Earth, rad
       Rs: Distance from Earth to SN, kpc
       Re: Distance from Earth to GC, kpc
     beta: Off-center angle, rad


    ********************
    *                  *
    *    Attributes    *
    *                  *
    ********************

    When a Geometry instance is initialized, the following attributes will be assigned,

          l: The line-of-sight distance, kpc
          d: Distance from SN to upscattered point, kpc
     rprime: Distance from GC to upscattered point, kpc
    cos_psi: cos(psi) at boost point where psi is the direction for BDM at B pointing
             Earth
    
    Suppose,
    
    >>> bdmProp = Propagation(t=59,vx=0.9,theta=1e-4,phi=0,Rs=11.6,Re=8.5,beta=0.71)
    
    Then the lengths in the propagation geometry are specified by 
    
    >>> print(bdmProp.l)  
    5.160120751743069e-09
    >>> print(bdmProp.d)  
    11.59999999483988
    >>> print(bdmProp.rprime)
    8.49999999608676
    >>> print(bdmProp.cos_psi) 
    1.0

    See API/Propagation/Positioning for detail
    """

    def __init__(self,t,vx,theta,phi,Rs,Re,beta):
        self.t = t
        self.vx = vx
        self.Rs = Rs
        self.theta = theta

        # Evaluated quantities
        self._zeta = self.get_zeta()
        self._l = self.get_ell()

        # Now initialize superclass
        super().__init__(self._l,theta,phi,Rs,Re,beta)
    
    def get_zeta(self):
        """
        Get zeta, see Eq. (10) in the document, kpc
        """
        Rs, t = self.Rs, self.t
        zeta = Rs + constant.c * t / constant.kpc2cm
        return zeta
    
    def get_ell(self):
        """
        Get ell, see Eq. (11) in the document, kpc
        """
        vx, Rs, theta, zeta = self.vx, self.Rs, self.theta, self._zeta
        cos_theta = cos(theta)
        alpha = sqrt((Rs**2 - zeta**2) * (1 - vx**2) + (Rs * vx * cos_theta - zeta)**2)
        gamma = Rs * vx * cos_theta
        ell = - vx * (alpha + gamma - zeta) / (1 - vx**2)
        return ell

