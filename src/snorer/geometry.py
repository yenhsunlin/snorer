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

__all__ = ['Geometry',]


#---------- Import required utilities ----------#

from numpy import sqrt,isclose,cos,sin
from .constants import constant
from .sysmsg import FlagError



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
    This class constructs the propagation geometry when the loactions of SN, GC and
    Earth, ToF and BDM velocity are specified.

    This class is not exclusively for SN in MW or LMC, it can be generalized to SN
    in arbitrary distant galaxy as long as the aforementioned inputs are determined.
    The BDM emissivity along the line-of-sight then can be determined when calculate
    the BDM flux and event at Earth associated to that particular SN.

    /*-----------------------------------------------------------------------------*/

    ********************
    *                  *
    *   Class Inputs   *
    *                  *
    ********************
    
        t: The BDM ToF, seconds
    theta: The zenith angle at Earth, rad
      phi: Azimuthal angle at Earth, rad
       vx: BDM velocity, in the unit of c
    Rstar: Distance from Earth to SN, kpc
       Re: Distance from Earth to GC, kpc
     beta: Off-center angle, rad


    ********************
    *                  *
    *    Attributes    *
    *                  *
    ********************

    When a Geometry instance is initialized, the following attributes will be assigned,

         d: Distance from Earth to upscattered point, kpc
         D: Distance from SN to upscattered point, kpc
    rprime: Distance from GC to upscattered point, kpc
    cosPsi: The cosine value of the required scattering angle that points toward Earth
            at the upscattered point
    
    Suppose,
    
    >>> bdmGeo = Geometry(t=59,theta=1e-4,vx=0.9,Rstar=11.6,Re=8.5,beta=0.71,phi=0)
    
    Then the lengths in the propagation geometry are specified by 
    
    >>> bdmGeo.d
    5.160120751743069e-09
    >>> bdmGeo.D
    11.59999999483988
    >>> bdmGeo.rprime
    8.49999999608676

    See Phys. Rev. D 108, 083013 (2023), arXiv:2307.03522 for theoretical foundation
    """

    def __init__(self,t,theta,phi,vx,Rstar,Re,beta):
        self.t = t
        self.theta = theta
        self.phi = phi
        self.vx = vx
        self.Rstar = Rstar
        self.Re = Re
        self.beta = beta
        
    @property
    def d(self):
        return self.__class__.get_d(self.t,self.vx,self.Rstar,self.theta)

    @property
    def D(self):
        return self.__class__.get_D(self.d,self.Rstar,self.theta)

    @property
    def rprime(self):
        return self.__class__.get_rPrime(self.d,self.Re,self.theta,self.phi,self.beta)

    @property
    def cosPsi(self):
        return self.__class__.get_cosPsi(self.d,self.Rstar,self.theta)
    
    @classmethod
    def get_cosPsi(cls,d,Rstar,theta) -> float:
        """
        Get the cosine value of scattering angle cos(psi).
        If we did it with law of cosine, then for the case of psi > pi/2,
        it will always return pi - psi which cannot reflect the pratical
        situation
        
        Input
        ------
        d: the l.o.s distance, kpc
        Rstar: the distance between Earth and SN, kpc
        theta: the open-angle, rad
        
        Output
        ------
        psi: scattering angle, rad
        """
        
        # Get D^2
        D2 = cls.get_D(d,Rstar,theta,is_squared = True)
        D = sqrt(D2)
        # Get cos(psi)
        denominator = 2*D*d # check if the denominator in the law of cosine is not 0.0
        if ~isclose(0,denominator,atol=1e-100):
            numerator = Rstar**2 - D2 - d**2
            cosPsi = numerator/denominator
            # Dealing with round-off error
            if cosPsi > 1: cosPsi = 1
            elif cosPsi < -1: cosPsi = -1
            else: pass
        else:
            # the denominator is 0.0, which means d = 0, applying L'Hospital's rule to get cos(psi)
            cosPsi = 0
        return cosPsi
    
    @classmethod
    def get_D(cls,d,Rstar,theta,is_squared = False) -> float:
        """
        Calculate the distance between SN and boosted point D
        
        Input
        ------
        d: the l.o.s distance, kpc
        Rstar: the distance between Earth and SN, kpc
        theta: the open-angle, rad
        is_squared: return the square of such distance, default is False
        
        Output
        ------
        D: the distance D
        """
        # Calculate D^2 via law of cosine
        D2 = d**2 + Rstar**2 - 2*d*Rstar*cos(theta)
        # D2 might turn minus due to round-off error, it shoud truncate at 0
        if D2 < 0: D2 = 0
    
        if is_squared is True:
            return D2
        elif is_squared is False:
            return sqrt(D2)
        else:
            raise FlagError('Keyword argument \'is_squared\' must be a boolean.')
        
    @classmethod
    def get_ell(cls,d,Re,theta,beta,is_squared = False) -> float:
        """
        Calculate the distance ell
        
        Input
        ------
        d: the l.o.s distance, kpc
        Re: the distance between Earth and GC, kpc
        theta: the open-angle, rad
        beta: the off-center angle, rad
        is_square: return the squared ell, default is False
        
        Output
        ------
        ell: the distance ell
        """
        # Calculate ell^2 via law of cosine
        ell2 = Re**2 + (d*cos(theta))**2 - 2*Re*d*cos(theta)*cos(beta)
        # ell2 might turn minus due to round-off error, it should truncate at 0
        if ell2 < 0: ell2 = 0.0
    
        if is_squared is True:
            return ell2
        elif is_squared is False:
            return sqrt(ell2)
        else:
            raise FlagError('Flag \'is_squared\' must be a boolean.')
    
    @classmethod
    def get_rPrime(cls,d,Re,theta,phi,beta) -> float:
        """
        Calculate the distance from boosted point to GC r'
        
        Input
        ------
        d: the l.o.s distance, kpc
        Re: the distance between Earth and GC, kpc
        theta: the open-angle, rad
        phi: the azimuthal angle, rad
        beta: the off-center angle, rad
        
        Output
        ------
        r': kpc
        """
        # ell^2
        ell2 = cls.get_ell(d,Re,theta,beta,is_squared=True)
        # h
        h = d*sin(theta)
        
        # Calculate cos(iota) and iota
        denominator = 2*cos(theta)*sqrt(ell2)*d # check if the denomator in the law of cosine is not 0.0
        if ~isclose(0,denominator,atol=1e-100):
            numerator = Re**2 - ell2 - (d*cos(theta))**2
            cosIota = numerator/denominator
            # Dealing with round-off error
            if cosIota > 1: cosIota = 1
            elif cosIota < -1: cosIota = -1
            else: pass
        else:
            # the denominator is 0, which means d = 0, applying L'Hospital to get cos(iota)
            cosIota = 0
        # Using sin(arccos(x)) = sqrt(1-x^2)
        sinIota = sqrt(1 - cosIota**2)
        
        # Calculate r'^2
        rp2 = ell2*cosIota**2 + (sqrt(ell2)*sinIota - h*sin(phi))**2 + h**2*cos(phi)**2
        return sqrt(rp2)
    
    @classmethod
    def get_d(cls,t,vx,Rstar,theta) -> float:
        """
        Calculate the distance line-of-sight d
        
        Input
        ------
        t: the arrival time of BDM at Earth calibrated by SN neutrino, second
        vx: BDM velocity in the unit of light speed
        Rstar: the distance between Earth and SN, kpc
        theta: the open-angle, rad
        
        Output
        ------
        d: the l.o.s, kpc
        """
        zeta = Rstar + constant.c*t/constant.kpc2cm
        cos_theta = cos(theta)
        prefactor = vx/(1 - vx**2)
        root_squared = (Rstar**2 - zeta**2)*(1 - vx**2) + (Rstar*vx*cos_theta - zeta)**2
        if root_squared < 0:
            root_squared = 0
        d = prefactor*(zeta - Rstar*vx*cos_theta - sqrt(root_squared))
        return abs(d)