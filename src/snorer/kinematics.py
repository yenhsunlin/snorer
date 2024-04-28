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

__all__ = ['Neutrino',
           'get_vx',
           'get_maxPsi',
           'fx_lab',
           'get_thetaRange',
           'get_tof',]


#---------- Import required utilities ----------#

from numpy import sin,cos,tan,arccos,sqrt,pi
from scipy.optimize import root_scalar
from .constants import constant


##########################################################################
#                                                                        #
#   General Classes and Functions for Numerics                           #
#                                                                        #
##########################################################################


"""
This script contains following classes and functions

Classes
------
1. Neutrino

Functions
------
1. get_vx
2. get_maxPsi
3. fx_lab
4. get_thetaRange
5. get_tof

The docstrings should be sufficient for their self-explanations
"""


class Neutrino:
    """
    This class constructs the required neturino energy to have BDM compound (Tx,mx,psi) 

    /*-----------------------------------------------------------------------------*/

    ********************
    *                  *
    *   Class Inputs   *
    *                  *
    ********************
    
     Tx: BDM kinetic energy, MeV
     mx: DM mass, MeV
    psi: Lab frame scattering angle,rad


    ********************
    *                  *
    *    Attributes    *
    *                  *
    ********************

    When a Neutrino instance is established, the following attributes will be assigned,

           Ev: The required SNv energy to produce the given (Tx,mx,psi), MeV
          dEv: Jacobian dEv/dTx
    is_sanity: Does the given (Tx,mx,psi) satisfy energy conservation, bool
    
    Lets check if we want to upscatter mx = 1 keV DM to have kinetic energy Tx = 15 MeV with
    lab frame scattering angle psi = 0.05 rad
    
    >>> snv = Neutrino(Tx=15,mx=1e-3,psi=0.05)
    
    Then,
    
    >>> snv.Ev
    -0.8451953159963379
    >>> snv.dEv
    0.04756415761959332
    >>> snv.is_sanity
    False

    The last attribute is_sanity yields False as the required SNv energy is negative, which
    is clearly insane. This class assists the user to check whether such BDM kinematics can
    be accomplished or not.
    
    See Phys. Rev. D 108, 083013 (2023), arXiv:2307.03522 for details.
    """
    
    def __init__(self,Tx,mx,psi):
        self.Tx = Tx
        self.mx = mx
        self.psi = psi

    @property
    def is_sanity(self):
        """
        If returns False, it means the combination (Tx,mx,psi)
        violates energy conservation and is unphysical
        """
        return self.__class__.sanity_check(self.Tx,self.mx,self.psi)
    
    @property
    def Ev(self):
        return self.__class__.get_Ev(self.Tx,self.mx,self.psi)

    @property
    def dEv(self):
        return self.__class__.get_dEv(self.Tx,self.mx,self.psi)
    
    @classmethod
    def sanity_check(cls,Tx,mx,psi) -> bool:
        """
        Check if the combination (Tx,mx,psi) not violating
        energy conservation
        """
        px = sqrt(Tx*(Tx + 2*mx))
        if Tx < px*cos(psi):
            return True
        else:
            return False
    
    @classmethod
    def get_Ev(cls,Tx,mx,psi) -> float:
        """
        Get the required neutrino energy to boost DM up with kinetic
        energy Tx at angle psi
        """
        px = sqrt(Tx*(Tx + 2*mx))
        return - mx*Tx/(Tx - px*cos(psi))
    
    @classmethod
    def get_dEv(cls,Tx,mx,psi) -> float:
        """
        Get the dEv/dTx
        """
        px = sqrt(Tx*(Tx + 2*mx))
        x = cos(psi)
        return mx**2*Tx*x/(Tx - px*x)**2


def get_vx(Tx,mx) -> float:
    """
    Get the BDM velocity
    
    Input
    ------
    Tx: Kinetic energy of the particle, MeV
    mx: Mass of the particle, MeV
    
    Output
    ------
    vx: BDM velocity, dimensionless
    """
    return sqrt(Tx*(Tx + 2*mx))/(Tx + mx)
    

def fx_lab(Ev,mx,psi) -> float:
    """
    Calculate the angular distribution for cross section
    in lab frame with scattering angle psi. This is for
    model-independent case where the total cross section
    is independent of energy. If a particular model is
    introduced, one may obtain the angular distribution
    via the scattering amplitude instead of this.

    Input
    ------
    Tx: BDM kinetic energy, MeV
    mx: DM mass, MeV
    psi: Lab frame scattering angle, rad
    
    Output
    ------
    prob. dist.: differential distribution at psi
    """
    if 0 <= psi <= pi/2 and Ev > 0:
        # evaluate Lorentz factor in CM frame
        s = mx**2 + 2*Ev*mx
        Ecm = 0.5*(s + mx**2)/sqrt(s)
        gamma = Ecm/mx
        # evaluate the angular distribution
        sec = 1/cos(psi)
        fx = gamma**2*sec**3/pi/(1 + gamma**2*tan(psi)**2)**2
    else:
        fx = 0
    return fx


def get_maxPsi(Tx,mx) -> float:
    """
    Get the maximumly allowed scattering angle psi
    
    Input
    ------
    Tx: BDM kinetic energy, MeV
    mx: DM mass, MeV
    
    Output
    ------
    psi_max: rad
    """
    maxCosValue = sqrt(Tx/(Tx + 2*mx))
    return arccos(maxCosValue)


def get_thetaRange(t,Tx,mx,Rstar) -> tuple:
    """
    Find the range of zenith angle theta that yields non-zero SNv BDM flux 
    
    Input
    ------
    t: the BDM ToF, if t > t_van, the result is unphysical
    Tx: DM kinetic energy, MeV
    mx: DM mass, MeV
    Rstar: Distance between Earth and SN, in kpc
    
    Output
    ------
    tup: (theta_min,thata_max), rad
    """
    vx = get_vx(Tx,mx)
    t0 = Rstar*constant.kpc2cm/constant.c
    # find a_min
    psiM = get_maxPsi(Tx,mx)
    
    # We will use Newton-Raphson method to find the root that corresponds to theta_max
    def f(theta):
        """The target function to be solve"""
        return sin(theta) + sin(psiM - theta)/vx - (t/t0 + 1)*sin(psiM)
    def fp(theta):
        """The derivative of the target function, f_prime"""
        return cos(theta) - cos(psiM - theta)/vx
    
    # find solutions to theta bound using Newton-Raphson method
    sol_theta_min = root_scalar(f,method='newton',x0=0,fprime=fp).root
    sol_theta_max = root_scalar(f,method='newton',x0=pi/2,fprime=fp).root
    # If thetaMin < 0, it's not physical. Set it to 0
    if sol_theta_min < 0:
        sol_theta_min = 0
    return sol_theta_min,sol_theta_max


def get_tof(Tx,mx,Rstar) -> tuple:
    """
    Get the ToF: peak time and vanishing time
    
    Input
    ------
    Tx: BDM kinetic energy, MeV
    mx: DM mass, MeV
    Rstar: Distance from SN to the Earth, kpc
    
    Output
    ------
    tup: (t_peak,t_van), seconds
    """
    # Get maximum psi and BDM velocity
    psiM = get_maxPsi(Tx,mx)
    vx = get_vx(Tx,mx)
    
    # Solving the corresponding theta that maximizes t
    def theta(theta):
        """ Target function """
        return cos(psiM - theta)/cos(theta) - vx
    theta = root_scalar(theta, method='brentq', bracket=[0,pi/2]).root

    # Evaluating the vanishing time and peak time
    t0 = Rstar*constant.kpc2cm/constant.c
    t_van = ((sin(theta) + sin(psiM - theta)/vx)/sin(psiM) - 1)*t0
    t_peak = Rstar*constant.kpc2cm/vx/constant.c - t0
    return t_peak,t_van