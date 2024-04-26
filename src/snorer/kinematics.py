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

from np import sin,cos,tan,arccos,sqrt,isclose,pi
from scipy.optimize import root_scalar as _root_scalar
from .constant import Constants



##########################################################################
#                                                                        #
#   General Classes and Functions for Numerics                           #
#                                                                        #
##########################################################################


constant = Constants()


class Kinematics:
    """
    This class constructs the BDM kinematics after it gets boosted by SNv

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

    When a Kinematics instance is initialized, the following attributes will be assigned,

           Ev: The required SNv energy to produce the given (Tx,mx,psi), MeV
          dEv: Jacobian dEv/dTx
    is_sanity: Does the given (Tx,mx,psi) satisfy energy conservation, bool
    
    Lets check if we want to upscatter mx = 1 keV DM to have kinetic energy Tx = 15 MeV with
    lab frame scattering angle psi = 0.05 rad
    
    >>> bdmKinetic = Kinematics(Tx=15,mx=1e-3,psi=0.05)
    
    Then,
    
    >>> bdmKinetic.Ev
    -0.8451953159963379
    >>> bdmKinetic.dEv
    0.04756415761959332
    >>> bdmKinetic.is_sanity
    False

    The last attribute is_sanity yields False as the required SNv energy is negative, which
    is clearly insane. This class assists the user to check whether such BDM kinematics can
    be accomplished or not.

    The docstring in each class method should be suffcient for self-explanation.
    Mathematical details are documented in the Appendix of 
    
    See Phys. Rev. D 108, 083013 (2023), arXiv:2307.03522 for theoretical foundation
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
        return self.sanity_check()
    
    # @property
    # def vx(self):
    #     return self.__class__.get_vx(self.Tx,self.mx)
    
    @property
    def Ev(self):
        return self.get_Ev()

    @property
    def dEv(self):
        return self.get_dEv()
    
    def sanity_check(self) -> bool:
        """
        Check if the combination (Tx,mx,psi) not violating
        energy conservation
        """
        Tx = self.Tx
        mx = self.mx
        psi = self.psi
        tan_psi2 = tan(psi)**2
        if isclose(0,tan_psi2,atol=1e-100):
            return False
        elif Tx < 2*mx/tan_psi2:
            return True
        else:
            return False
    
    def get_Ev(self) -> float:
        """
        Get the required neutrino energy to boost DM up with kinetic
        energy Tx at angle psi
        """
        Tx = self.Tx
        mx = self.mx
        psi = self.psi
        px = sqrt(Tx*(Tx + 2*mx))
        return - mx*Tx/(Tx - px*cos(psi))
    
    def get_dEv(self) -> float:
        """
        Get the dEv/dTx
        """
        Tx = self.Tx
        mx = self.mx
        psi = self.psi
        px = sqrt(Tx*(Tx + 2*mx))
        x = cos(psi)
        return mx**2*Tx*x/(Tx - px*x)**2

    @staticmethod
    def get_vx(Tx,mx) -> float:
        """
        Get the velocity of particle with mass mx and kinetic energy Tx
        
        Input
        ------
        Tx: Kinetic energy of the particle, MeV
        mx: Mass of the particle, MeV
        
        Output
        ------
        vx: Particle velocity in the unit of light speed
        """
        return sqrt(Tx*(Tx + 2*mx))/(Tx + mx)
    
    @staticmethod
    def dmdOmega_lab(Ev,mx,psi) -> float:
        """
        Calculate the angular distribution for cross section
        in lab frame with scattering angle psi. This is for
        model-independent case where the total cross section
        is independent of energy. If a particular model is
        introduced, one may obtain the angular distribution
        via the scattering amplitude instead of this.
        """
        if 0 <= psi <= pi/2 and Ev > 0:
            # evaluate Lorentz factor in CM frame
            s = mx**2+2*Ev*mx
            Ecm = 0.5*(s + mx**2)/sqrt(s)
            gamma = Ecm/mx
            # evaluate the angular distribution
            sec = 1/cos(psi)
            dndOmega = gamma**2*sec**3/pi/(1+gamma**2*tan(psi)**2)**2
        else:
            dndOmega = 0
        return dndOmega

    @staticmethod
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
    t: the BDM ToF, if t > t_van, the result is not trustable
    Tx: DM kinetic energy, MeV
    mx: DM mass, MeV
    Rstar: Distance between Earth and SN, in kpc
    
    Output
    ------
    tup: (theta_min,thata_max), rad
    """
    vx = Kinematics.get_vx(Tx,mx)
    t0 = Rstar*constant.kpc2cm/constant.c
    # find a_min
    psiM = Kinematics.get_maxPsi(Tx,mx)
    
    # We will use Newton-Raphson method to find the root that corresponds to theta_max
    # This method requires the target function _f and its derivative _fp
    def _f(theta):
        return sin(theta) + sin(psiM - theta)/vx - (t/t0 + 1)*sin(psiM)
    
    def _fp(theta):
        return cos(theta) - cos(psiM - theta)/vx
    
    # find solutions to theta bound using Newton-Raphson method
    sol_theta_min = _root_scalar(_f,method='newton',x0=0,fprime=_fp).root
    sol_theta_max = _root_scalar(_f,method='newton',x0=pi/2,fprime=_fp).root
    # If thetaMin < 0, it's not physical. Set it to 0
    if sol_theta_min < 0: sol_theta_min = 0
        
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
    psiM = Kinematics.get_maxPsi(Tx,mx)
    vx = Kinematics.get_vx(Tx,mx)
    
    # Solving the corresponding theta that maximizes t
    def _theta(theta):
        """ Target function """
        return cos(psiM - theta)/cos(theta) - vx
    theta = _root_scalar(_theta, method='brentq', bracket=[0,pi/2]).root

    # Evaluating the vanishing time and peak time
    t0 = Rstar*constant.kpc2cm/constant.c
    t_van = ((sin(theta) + sin(psiM - theta)/vx)/sin(psiM) - 1)*t0
    t_peak = Rstar*constant.kpc2cm/vx/constant.c - t0
    return t_peak,t_van