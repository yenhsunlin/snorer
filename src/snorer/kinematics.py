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
           'get_psiMax',
           'get_gx',
           'get_thetaMax',
           'get_tof',
           'KallenLambda',
           'get_tBound',]


#---------- Import required utilities ----------#

from numpy import sin,cos,tan,arccos,sqrt,pi,abs,clip,atleast_1d,broadcast_arrays,zeros_like,nditer
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
    
     T2: Kinetic energy received by the target, MeV
     m1: Incident particle mass, MeV
     m2: Target mass, MeV
    psi: Lab frame scattering angle, rad


    ********************
    *                  *
    *     Methods      *
    *                  *
    ********************

    When an associated instance is established, the following attributes can be retrieved,

        T1: The required kietic energy T1 for the incident particle that can have (T2,m2,psi)
       dT1: The Jacobian dE1/dT2
         x: The cos(psi)
    sanity: Does the given (T2,m1,m2,psi) satisfy energy conservation, bool
     dLips: The differential Lorentz invariant phase space associated with psi

    Note that the total energy of the particle 2 is E2 = T2 + m2, similarily, E1 = T1 + m1


    ********************
    *                  *
    *     Examples     *
    *                  *
    ********************
    
    In terms of SNv BDM, we would like to have T1 = Ev, m1 = mv = 0 and T2 = Tx,
    m2 = mx. Suppose, mx = 0.001 MeV, Tx = 15 MeV with psi = 0.05 rad, we want to obtain
    the required Ev to make this happens
    
    >>> snv = Kinematics(T2=15,m1=0,m2=1e-3,psi=0.05)
    
    Then,

    >>> snv.T1
    -0.8451953159962898
    >>> snv.dT1
    0.0031707324661873464
    >>> snv.sanity
    False
    >>> snv.dLips
    1583.2490337942002

    The last attribute "sanity" yields False as the required SNv energy is negative, which
    is clearly insane.

    See "API/Particle kinematics/2-2 elastic scattering" for details
    """
    def __init__(self,T2,m1,m2,psi):
        self.T2 = T2
        self.m1 = m1
        self.m2 = m2
        self.psi = psi

        # Evaluated quantities
        self._x = cos(self.psi) # cos(psi)
        self._p2squared = T2*(T2 + 2*m2) # momentum squared of particle 2
        self._p2 = sqrt(T2*(T2 + 2*m2)) # momentum of particle 2
        self._sanity = (self.T2 < self._p2 * self.x) # does energy conservation hold?
        self._T1 = self.get_T1() # required kinetic energy of particle 1
        self._dT1 = self.get_dT1() # Jacobian dT1/dT2
        self._dLips = self.get_dLips() # differential Lorentz invariant phase space
    
    @property
    def T1(self) -> float:
        return self._T1
    
    @property
    def dT1(self) -> float:
        return self._dT1
    
    @property
    def x(self) -> float:
        return self._x
    
    @property
    def sanity(self) -> bool:
        return self._sanity
    
    @property
    def dLips(self) -> float:
        return self._dLips
    
    def get_T1(self) -> float:
        """
        Eq. (2) in API/Particle kinematics/2-2 elastic scattering
        """
        T2,m1,m2,p2,p2sq,x = self.T2,self.m1,self.m2,self._p2,self._p2squared,self.x
        numerator = T2**2 * m2 + p2 * x * sqrt(m1**2 * p2sq * x**2 + T2**2 * (m2**2 - m1**2))
        denominator = p2sq * x**2 - T2**2
        E1 = numerator/denominator # total energy of particle 1
        T1 = E1 - m1 # kinetic energy of particle 1
        return T1
        # -- Following are the old expressions in v1
        # # Mathematica generated expression
        # p2sq = T2*(T2+2*m2)
        # # Mathematica generated expression
        # numerator = sqrt(m1**2*p2sq**2*x**4 + p2sq*T2**2*x**2*(m2**2-m1**2)) + m2*T2**2
        # denominator = p2sq*x**2 - T2**2
        # return numerator/denominator - m1
        
        # following expressions may yield invalid value in sqrt
        # numerator = (2*m1*T2+2*m2*T2-4*m1*m2*x**2-2*m1*T2*x**2+sqrt(-4*(-m1**2*T2-2*m1*m2*T2  \
        #              -m2**2*T2)*(-T2+2*m2*x**2+T2*x**2)+(-2*m1*T2-2*m2*T2+4*m1*m2*x**2        \
        #              +2*m1*T2*x**2)**2))
        # denominator = (2*(-T2+2*m2*x**2+T2*x**2))
        # return numerator/denominator

    def get_dT1(self) -> float:
        """
        Eq. (3) in API/Particle kinematics/2-2 elastic scattering
        """
        T2,m1,m2,x = self.T2,self.m1,self.m2,self.x
        delta = -T2 + (T2 + 2 * m2) * x**2
        alpha = m1**2 * delta
        kappa = sqrt((T2 + 2 * m2) * (alpha + T2 * m2**2))
        beta = m2**2 * (2 * T2 + delta)
        gamma = 2 * m2 * x * kappa
        eta = delta**2 * x * kappa
        dT1 = m2 * x**2 * (alpha + beta + gamma)/eta
        return dT1
        # -- Following are the old expressions in v1
        # # Mathematica generated expression
        # numerator = (m2*x**2*(m1**2*T2*(-T2+(2*m2+T2)*x**2)+m2*(m2*T2*(T2+(2*m2+T2)*x**2)+    \
        #              2*sqrt(T2**2*(2*m2+T2)*x**2*(m2**2*T2+m1**2*(-T2+(2*m2+T2)*x**2)))))) 
        # denominator = ((T2-(2*m2+T2)*x**2)**2*sqrt(T2**2*(2*m2+T2)*x**2*(m2**2*T2+m1**2       \
        #               *(-T2+(2*m2+T2)*x**2))))
        # return numerator/denominator

    def get_dLips(self) -> float:
        T1,m1,m2 = self.T1,self.m1,self.m2
        s = m1**2 + m2**2 + 2*m2*(T1 + m1)
        psq = (s - (m1+m2)**2)*(s-(m1-m2)**2)/4/s
        return abs(self._du_dx()/64/pi/s/psq/2/pi)

    def _du_dx(self) -> float:
        T2,m1,m2,x = self.T2,self.m1,self.m2,self.x
        E1 = self.T1 + m1
        E2 = T2 + m2
        p2sq = self._p2squared
        p2 = self._p2
        # Mathematica generated expression
        num = (p2sq*T2**2*x*((m1-m2)*(m1+m2)*T2**2-(m1**2+m2**2)*p2sq*x**2-2*m2*sqrt((-m1**2  \
               +m2**2)*p2sq*T2**2*x**2+m1**2*p2sq**2*x**4)))
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
    
     T2: Kinetic energy received by the target, MeV
     m1: Incident particle mass, MeV
     m2: Target mass, MeV
    psi: Lab frame scattering angle, rad


    ********************
    *                  *
    *    Attributes    *
    *                  *
    ********************

    When an associated instance is established, the following attributes are assigned,

    T1: The required kietic energy for the incident particle that can have (T2,m2,psi), MeV
     s: Associated Mandelstam variable s, MeV^2
     t: Associated Mandelstam variable t, MeV^2
     u: Associated Mandelstam variable u, MeV^2

    See "API/Particle kinematics/2-2 elastic scattering" for details
    """
    def __init__(self,T2,m1,m2,psi):
        super().__init__(T2,m1,m2,psi)
        self._s = self.get_s()
        self._t = self.get_t()
        self._u = self.get_u()
    
    @property
    def s(self):
        return self._s
    
    @property
    def t(self):
        return self._t

    @property
    def u(self):
        return self._u
    
    def get_s(self):
        E1 = self.T1 + self.m1
        return 2 * E1 * self.m2 + self.m1**2 + self.m2**2
    
    def get_t(self):
        E2 = self.T2 + self.m2
        return -2 * self.m2 * (E2 - self.m2)

    def get_u(self):
        return 2 * (self.m1**2 + self.m2**2) - self.s - self.t


class Neutrino(Kinematics):
    """
    Superclass: Kinematics
    
    This class constructs the required neturino energy to have BDM with (Tx,mx,psi) 

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
    sanity: Does the given (Tx,mx,psi) satisfy energy conservation, bool
    
    Lets check if we want to upscatter mx = 1 keV DM to have kinetic energy Tx = 15 MeV with
    lab frame scattering angle psi = 0.05 rad
    
    >>> snv = Neutrino(Tx=15,mx=1e-3,psi=0.05)
    
    Then,
    
    >>> snv.Ev
    -0.8451953159963379
    >>> snv.dEv
    0.00317073246618734
    >>> snv.sanity
    False

    The last attribute "sanity" yields False as the required SNv energy is negative, which
    is clearly insane.
    
    See "User Manual/Physics Overview" for details
    """
    
    def __init__(self,Tx,mx,psi):
        super().__init__(Tx,0,mx,psi)

    @property
    def Ev(self) -> float:
        """
        Get the required SNv energy
        """
        return self.T1
    
    @property
    def dEv(self) -> float:
        """
        Get the dEv/dTx
        """
        return self.dT1


def get_vx(Tx,mx) -> float:
    """
    Get dimensionless BDM velocity
    
    Parameters
    ----------
    Tx : array_like
        BDM kinetic energy, MeV
    mx : array_like
        DM mass, MeV
    
    Returns
    -------
    out : scalar/ndarray
        BDM velocity, dimensionless
    """
    return sqrt(Tx*(Tx + 2*mx))/(Tx + mx)
    

# def fx_lab(Ev,mx,psi) -> float:
#     """
#     Calculate the angular distribution for cross section
#     in lab frame with scattering angle psi. This is for
#     model-independent case where the total cross section
#     is independent of energy. If a particular model is
#     introduced, one may obtain the angular distribution
#     via the scattering amplitude instead of this.

#     Input
#     ------
#     Tx: BDM kinetic energy, MeV
#     mx: DM mass, MeV
#     psi: Lab frame scattering angle, rad
    
#     Output
#     ------
#     prob. dist.: differential distribution at psi
#     """
#     if 0 <= psi <= pi/2 and Ev > 0:
#         # evaluate Lorentz factor in CM frame
#         s = mx**2 + 2*Ev*mx
#         Ecm = 0.5*(s + mx**2)/sqrt(s)
#         gamma = Ecm/mx
#         # evaluate the angular distribution
#         sec = 1/cos(psi)
#         fx = gamma**2*sec**3/pi/(1 + gamma**2*tan(psi)**2)**2
#     else:
#         fx = 0
#     return fx


def _gx_lab(Ev,mx,psi) -> float:
    """
    Unit component of get_gx, see its docstring
    """
    # evaluate Lorentz factor in CM frame
    s = mx**2 + 2*Ev*mx
    Ecm = 0.5*(s + mx**2)/sqrt(s)
    gamma = Ecm/mx
    # evaluate the angular distribution
    sec = 1/cos(psi)
    gx = gamma**2*sec**3/pi/(1 + gamma**2*tan(psi)**2)**2
    return gx


def get_gx(Ev,mx,psi) -> float:
    """
    Calculate the probability density for cross section at scattering
    angle psi and averaged over azimuthal angle in lab frame. This is
    for energy-independent cross section.

    Parameters
    ----------
    Ev : array_like
        Incoming neutrino energy, MeV
    mx : array_like
        DM mass, MeV
    psi : array_like
        Lab frame scattering angle, rad
    
    Returns
    -------
    out : scalar/ndarray
        Probability density for cross section at psi and averaged
        over azimuthal angle 2*pi. The result is a scalar if the three
        inputs are all scalars. The unit is 1/sr.
    """
    psi = atleast_1d(psi) # Let psi be at least 1d array for easy manipulation
    Ev,mx,psi = broadcast_arrays(Ev,mx,psi) # Broadcast input into the same dim
    invalid_region = ~((0 <= psi) & (psi < 0.5 * pi)) # psi that is outside the valid region
    gx = _gx_lab(Ev,mx,psi) # Evaluate gx
    gx[invalid_region] = 0 # mask the invalid region as 0
    return gx if invalid_region.size > 1 else gx.item()


def get_psiMax(Tx,mx) -> float:
    """
    Get the maximumly allowed scattering angle psi
    
    Parameters
    ----------
    Tx : array_like
        BDM kinetic energy, MeV
    mx : array_like
        DM mass, MeV
    
    Returns
    -------
    out : scalar/ndarray
        Maximum allowed scattering angle psi, rad
    """
    maxCosValue = sqrt(Tx/(Tx + 2*mx))
    maxCosValue = clip(maxCosValue,-1,1)
    return arccos(maxCosValue)


def get_thetaMax(t,Tx,mx,Rstar) -> float:
    """
    Find the maximum BDM field-of-view that centers SN at particular time t
    
    Input
    ------
    t : scalar
        The BDM at particular time t, seconds. If t > t_van, the result is
        unphysical
    Tx : scalar 
        BDM kinetic energy, MeV
    mx : scalar
        DM mass, MeV
    Rs : scalar
        Distance to supernova, kpc
    
    Returns
    -------
    out : scalar
        Maximum field-of-view centers supernova, theta_MAX, rad

    See Eq. (24) in "User Manual/Physics Overview" for detail.
    """
    vx = get_vx(Tx,mx)
    tv = Rstar*constant.kpc2cm/constant.c
    # find a_min
    psiM = get_psiMax(Tx,mx)
    
    # We will use Newton-Raphson method to find the root that corresponds to theta_max
    def f(theta):
        """The target function to be solve"""
        return sin(theta) + sin(psiM - theta)/vx - (t/tv + 1)*sin(psiM)
    
    # find solutions to theta bound using Newton-Raphson method
    theta_M = root_scalar(f,method='newton',x0=pi/2).root
    return theta_M


def _get_tof(Tx,mx,Rs) -> tuple:
    """
    Get peak time and vanishing time for BDM
    
    Parameters
    ----------
    Tx : scalar 
        BDM kinetic energy, MeV
    mx : scalar
        DM mass, MeV
    Rs : scalar
        Distance to SN, kpc
    
    Returns
    -------
    out : tuple
        (t_peak,t_van), seconds

    See Eqs. (22-23) in "User Manual/Physics Overview" for detail.
    """
    # Get maximum psi and BDM velocity
    psiM = get_psiMax(Tx,mx)
    vx = get_vx(Tx,mx)
    
    # Solving the corresponding theta that maximizes t
    def theta(theta):
        """ Target function """
        return cos(psiM - theta)/cos(theta) - vx
    theta_MAX = root_scalar(theta, method='brentq', bracket=[0,pi/2]).root

    # Evaluating the vanishing time and peak time
    tv = Rs*constant.kpc2cm/constant.c
    t_van = ((sin(theta_MAX) + sin(psiM - theta_MAX)/vx)/sin(psiM) - 1)*tv
    t_peak = Rs*constant.kpc2cm/vx/constant.c - tv
    return t_peak,t_van


def get_tvan(Tx,mx,Rs) -> float:
    """
    Get the BDM vanishing time. The time-zero is set as
    the arrival of SN neutrinos at Earth.

    Parameters
    ----------
    Tx : array_like
        BDM kinetic energy, MeV
    mx : array_like
        DM mass, MeV
    Rs : array_like
        Distance to supernova, kpc

    Returns
    -------
    out : scalar/ndarray
        BDM vanishing time, seconds
    """
    Rs = atleast_1d(Rs) # Let Rs be at least 1d array for easy manipulation
    Tx,mx,Rs = broadcast_arrays(Tx,mx,Rs)
    # Setup empty array to store t_van values for every Rs
    t_van = zeros_like(Rs)
    # Use nditer to iterate and get the t_van
    with nditer([Tx,mx,Rs,t_van],op_flags=[['readonly'],['readonly'],['readonly'],['writeonly']]) as it:
        for T,m,r,tvan in it:
            tvan[...] = _get_tof(T,m,r)[-1]
    return t_van if t_van.size > 1 else t_van.item()


def KallenLambda(x,y,z) -> float: 
    """
    Kallen lambda function, A useful function for evaluating
    kinetical quantities in particle physics.

    Parameters
    ----------
    x : array_like
    y : array_like
    z : array_like

    Returns
    -------
    out : scalar/ndarray
    """
    return x**2 + y**2 + z**2 - 2 * (x * y + y * z + z * x)


def get_tBound(T1,m1,m2) -> tuple:
    """
    Get the allowed range for Mandelstam variable t

    Parameters
    ----------
    T1 : array_like
        Incident particle kinetic energy, MeV
    m1 : array_like
        Incident particle mass, MeV
    m2 : array_like   
        Target mass, MeV

    Returns
    -------
    out : tuple/array of tuple
        If all inputs are scalar: (t_min,t_max), MeV^2
        One of them is array, after broadcasting: (array_t_min, array_t_max) MeV^2

    See Eqs. (3.34,35) in V. Ilisie, Concepts in QFT, Springer (2016)
    """
    T1,m1,m2 = broadcast_arrays(T1,m1,m2)
    E1 = T1 + m1
    s = 2*m2*E1 + m1**2 + m2**2
    factor1 = m1**2 + m2**2 - s/2 - (m1**2 - m2**2)**2/2/s
    factor2 = KallenLambda(s,m1**2,m2**2)/2/s
    return factor1 - factor2,factor1 + factor2
