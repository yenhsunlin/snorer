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

__all__ = ['Kinematics',
           'Mandelstam',
           'Neutrino',
           'get_vx',
           'get_maxPsi',
           'fx_lab',
           'get_thetaRange',
           'get_tof',
           'KallenLambda',
           'get_tBound',]


#---------- Import required utilities ----------#

from numpy import sin,cos,tan,arccos,sqrt,pi,abs
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
1. Kinematics
2. Mandelstam
3. Neutrino

Functions
------
1. get_vx
2. get_maxPsi
3. fx_lab
4. get_thetaRange
5. get_tof
6. KallenLambda
7. get_tBound

The docstrings should be sufficient for their self-explanations
"""


class Kinematics:
    """
    This class constructs the required kinetic energy T1 of incoming particle with
    mass m1 to boost the target with mass m2 to kinetic energy T2 along the direction
    psi. See the following scheme
                                      
                                  ●
      m1          m2             /
      ●----> T1   ◯         ------------
                                 \ )psi
                               T2 ◯
         
         Before                  After                            
        
    /*-----------------------------------------------------------------------------*/

    ********************
    *                  *
    *   Class Inputs   *
    *                  *
    ********************
    
     m1: Incident particle mass, MeV
     m2: Target mass, MeV
     T2: Kinetic energy received by the target, MeV
    psi: Lab frame scattering angle, rad


    ********************
    *                  *
    *     Methods      *
    *                  *
    ********************

    When an associated instance is established, the following methods can be called,

        get_T1(): Get the required kietic energy T1 for the incident particle that
                  can have (m2,T2,psi)
       get_dT1(): Get the Jacobian dT1/dT2
    get_sanity(): Does the given (m1,m2,T2,psi) satisfy energy conservation, bool
     get_dLips(): Get the diff. Lorentz invariant phase space associated with psi

    Note that the total energy of the target E2 = T2 + m2, similarily, E1 = T1 + m1


    ********************
    *                  *
    *     Examples     *
    *                  *
    ********************
    
    In terms of SNv BDM, we would like to have T1 = Ev, m1 = mv = 0 and T2 = Tx,
    m2 = mx. Suppose, mx = 0.001 MeV, Tx = 15 MeV with psi = 0.05 rad, we want to obtain
    the required Ev to make this happens
    
    >>> snv = Kinematics(m1=0,m2=1e-3,T2=15,psi=0.05)
    
    Then,

    >>> snv.get_T1()
    -0.8451953159962898
    >>> snv.get_dT1()
    0.0031707324661873464
    >>> snv.get_sanity()
    False
    >>> snv.get_dLips()
    1583.2490337942002

    The last methid get_sanity() yields False as the required SNv energy is negative, which
    is clearly insane.
    """
    def __init__(self,m1,m2,T2,psi):
        self.m1 = m1
        self.m2 = m2
        self.T2 = T2
        self.psi = psi

    @property
    def _x(self):
        return cos(self.psi)
    
    def get_sanity(self) -> bool:
        return self.get_T1() >= 0 
    
    def get_T1(self) -> float:
        T2,m1,m2,x = self.T2,self.m1,self.m2,self._x
        p2sq = T2*(T2+2*m2)
        # Mathematica generated expression
        numerator = sqrt(m1**2*p2sq**2*x**4 + p2sq*T2**2*x**2*(m2**2-m1**2)) + m2*T2**2
        denominator = p2sq*x**2 - T2**2
        return numerator/denominator - m1

    def get_dT1(self) -> float:
        T2,m1,m2,x = self.T2,self.m1,self.m2,self._x
        # Mathematica generated expression
        numerator = (m2*x**2*(m1**2*T2*(-T2+(2*m2+T2)*x**2)+m2*(m2*T2*(T2+(2*m2+T2)*x**2)+2*sqrt(T2**2*(2*m2+T2)*x**2*(m2**2*T2+m1**2*(-T2+(2*m2+T2)*x**2))))))
        denominator = ((T2-(2*m2+T2)*x**2)**2*sqrt(T2**2*(2*m2+T2)*x**2*(m2**2*T2+m1**2*(-T2+(2*m2+T2)*x**2))))
        return numerator/denominator

    def get_dLips(self) -> float:
        T1,m1,m2 = self.get_T1(),self.m1,self.m2
        s = m1**2 + m2**2 + 2*m2*(T1 + m1)
        psq = (s - (m1+m2)**2)*(s-(m1-m2)**2)/4/s
        return abs(self._du_dx()/64/pi/s/psq/2/pi)

    def _du_dx(self) -> float:
        T2,m1,m2,x = self.T2,self.m1,self.m2,self._x
        E1 = self.get_T1() + m1
        E2 = T2 + m2
        p2sq = T2*(T2+2*m2)
        p2 = sqrt(p2sq)
        # Mathematica generated expression
        num = (p2sq*T2**2*x*((m1-m2)*(m1+m2)*T2**2-(m1**2+m2**2)*p2sq*x**2-2*m2*sqrt((-m1**2+m2**2)*p2sq*T2**2*x**2+m1**2*p2sq**2*x**4)))
        den = ((T2-p2*x)**2*(T2+p2*x)**2*sqrt((-m1**2+m2**2)*p2sq*T2**2*x**2+m1**2*p2sq**2*x**4))
        return -2*(num/den)*(E2 - p2*x) + 2*E1*p2


class Mandelstam(Kinematics):
    """
    Superclass: Kinematics

    This class constructs the associated Mandelstam variables associated with the
    following scattering scheme with (m2,T2,psi) are given. See the docstring of
    class Kinematics for detail.
                                      
                                  ●
      m1          m2             /
      ●----> T1   ◯         ------------
                                 \ )psi
                               T2 ◯
         
         Before                  After                            
        
    /*-----------------------------------------------------------------------------*/

    ********************
    *                  *
    *   Class Inputs   *
    *                  *
    ********************
    
     m1: Incident particle mass, MeV
     m2: Target mass, MeV
     T2: Kinetic energy received by the target, MeV
    psi: Lab frame scattering angle, rad


    ********************
    *                  *
    *    Attributes    *
    *                  *
    ********************

    When an associated instance is established, the following attributes are assigned,

    T1: The required kietic energy for the incident particle that can have
        (m2,T2,psi), MeV
     s: Associated Mandelstam variable s, MeV^2
     t: Associated Mandelstam variable t, MeV^2
     u: Associated Mandelstam variable u, MeV^2
    """
    def __init__(self,m1,m2,T2,psi):
        super().__init__(m1,m2,T2,psi)

    @property
    def T1(self) -> float:
        """Get the required T1 and assigend it as an instance attribute"""
        return self.get_T1()
    
    @property
    def s(self):
        E1 = self.T1 + self.m1
        return 2*E1*self.m2 + self.m1**2 + self.m2**2
    
    @property
    def t(self):
        E2 = self.T2 + self.m2
        return -2*E2*self.m2 + self.m1**2 + self.m2**2

    @property
    def u(self):
        return 2*(self.m1**2 + self.m2**2) - self.s - self.t


class Neutrino(Kinematics):
    """
    Superclass: Kinematics
    
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
    is clearly insane.
    
    See Phys. Rev. D 108, 083013 (2023), arXiv:2307.03522 for details.
    """
    
    def __init__(self,Tx,mx,psi):
        super().__init__(0,mx,Tx,psi)

    @property
    def is_sanity(self):
        """
        If returns False, it means the combination (Tx,mx,psi)
        violates energy conservation and is unphysical
        """
        return self.get_sanity()
    
    @property
    def Ev(self):
        """Get the required SNv energy"""
        return self.get_T1()

    @property
    def dEv(self):
        """Get the Jacobian dEv/dTx"""
        return self.get_dT1()


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


def KallenLambda(x,y,z) -> float: 
    """
    Kallen lambda function

    Input
    ------
    x: scalar variable
    y: scalar variable
    z: scalar variable

    Output
    ------
    scalar
    """
    return x**2 + y**2 + z**2 - 2*(x*y + y*z + z*x)


def get_tBound(m1,m2,T1) -> tuple:
    """
    Get the allowed range for Mandelstam variable t

    Input
    ------
    m1: Incident particle mass, MeV
    m2: Target mass, MeV
    T1: Incident particle kinetic energy, MeV

    Output
    ------
    tuple: (t_min,t_max), MeV^2 

    See Eqs. (3.34,35) in V. Ilisie, Concepts in QFT, Springer (2016)
    """
    E1 = T1 + m1
    s = 2*m2*E1 + m1**2 + m2**2
    factor1 = m1**2 + m2**2 - s/2 - (m1**2 - m2**2)**2/2/s
    factor2 = KallenLambda(s,m1**2,m2**2)/2/s
    return factor1 - factor2,factor1 + factor2