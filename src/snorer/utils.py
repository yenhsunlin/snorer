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

__all__ = ['snorerGeneralInterface',]


#---------- Import required utilities ----------#

from numpy import pi,sin,cos,arccos,isclose
import vegas
from astropy import units as u
from astropy.coordinates import SkyCoord
from .snorerMain import snNuSpectrum
from .kinematics import Mandelstam,Neutrino,get_vx,get_thetaRange,get_tof
from .geometry import Geometry
from .halo import rhox
from .constant import Constants


##########################################################################
#                                                                        #
#   General Classes and Functions for Numerics                           #
#                                                                        #
##########################################################################


"""
This script contains following classes and functions

Classes
------
1. snorerGeneralInterface

Functions
------
None

The docstrings should be sufficient for their self-explanations
"""


class snorerGeneralInterface(Constants):
    """
    Superclass: Constants
    
    Class for evaluating SNv BDM coming from SN in arbitrary distant galaxy with DM-v
    and DM-e interaction cross sections descrbied by a sepcific particle model. This
    class has an dependency on Astropy for handling SN/GC coordinates expressed in
    ICRS J2000.0 system.


    /*-----------------------------------------------------------------------------*/

    ********************
    *                  *
    *   Class Inputs   *
    *                  *
    ********************
    
    SN_coord: SN coordinate list, [RA,DEC,dist]: the first two are right ascension
              and declination of the celestial object respectively, both in string
              type. The last one is the distance to the object in kpc
    GC_coord: Identical to SN_coord, but is for the galactic center
     amp2_xv: Func type, amplitude squared for DM-v interaction, 4 positioning
              arguments.
              amp2_xv := some_func(s,t,u,mx): the first 3 are Mandelstam variables
              and the last one is the DM mass.
     amp2_xe: Identical to amp2_xv, but is for DM-e interaction
    **kwargs: Keywaord arguments that will be passed to SkyCoord class in Astropy.
              Leave it(them) blank if you have no idea to hack


    /*-------------------------------IMPORTANT NOTE--------------------------------*/

    The first FOUR inputs are POSITIONING-ONLY arguments, so their ORDER matters.


    ********************
    *                  *
    *    Attributes    *
    *                  *
    ********************

    When an instance is initialized, the following attributes will be assigned,

            Rstar: Distance to SN, kpc
               Re: Distance to GC, kpc
             beta: Off-center angle, rad
    separation_3d: Separation distance between the two celestial objects, kpc
    
    None of these attributes can be reassigned manually as they are completely
    determined when SN_coord and GC_coord are specified. If you want to have different
    attributes, you can only do by updating SN_coord/GC_coord.

    Suppose, we want to know the SN1987a in LMC, we already documented the coordinates
    of these two in constant, thus

    >>> SN_coord = constant.SN1987a_coord
    >>> GC_coord = constant.LMC_coord
    >>> sn1987a = snoreGeneralInterface(SN_coord,GC_coord,amp2_xv,amp2_xe)

    We explain the amplitude-square later. Now we can view the intrinsic properties
    of this instance by

    >>> sn1987a
          Dist to GC: 4.997e+01 kpc
          Dist to SN: 5.170e+01 kpc
     Seperation dist: 2.008e+00 kpc
    Off-center angle: 2.004e-02 rad

    They can be retrieved by, eg.,
    
    >>> sn1987a.beta     # off-center angle
    0.02004146106280973


    /*------------------------SNv Spectrum and DM Profile--------------------------*/

    When an instance is initialized, it also automatically assigned a SNv spectrum
    (in flux unit) and a DM halo profile. The default SNv spectrum is borrowed from
    `snNuSpectrum` and can be viewed

    >>> Ev,D=15,20      # MeV, kpc
    >>> sn1987a.snNuSpectrum(Ev,D)
    837920069.3938683

    You are free to modify it if you have your own SNv model to test, say `foo_snv`.
    The modified function takes two positioning-only arguments and can be accomplished
    by

    >>> sn1987a.snNuSpectrum = lambda Ev,D: foo_snv(Ev,D)   # change to my own SNv flux

    Similarily, we also borrowed `rhox` to generate DM number density. The default
    is MW NFW profile (# per cm^3) without spike and can be viewed

    >>> r,mx = 15,1e-4  # kpc, MeV
    >>> sn1987a.nx(15,1e-4)
    1149556.450912663

    If you want to implement your own halo profile, say `foo_nx`, sure no problem.
    Same as above, this new function takes two positioning-only arguments

    >>> sn1987a.nx = lambda r,mx: foo_nx(r,mx)

    Whatever features you wish to have in your `foo_snv` and `foo_nx`, you have to
    hardcode them in the beginning as the associated attributes only take two
    positioning-only arguments.
    
    If you have more than two inputs in these two functions, you might have to consider
    removing them or apply `partial` in the functools module to assign them in the
    beginning.


    /*---------------------------Constructing Amplitudes---------------------------*/
    
    We take the model discussed in Phys. Rev. D 108, 083013 (2023) for instance. Both
    DM-v and DM-e have the same amplitude square show by Eq. (3) in terms of Mandelstam
    variables (s,t,u),

        M^2 = 2 (Q/t - mV^2)^2 (s^2 + u^2 + 4t (m1^2 + m2^2) - 2 (m1^2 + m2^2)^2)

    where Q is the multiplication of coupling constants, m1 and m2 are the masses of
    incident and target particles respectively and mV the mediator mass.

    Thus we can construct DM-v amplitude square by letting, m1 = 0 and assume mV = mx/3,
    gV = 1e-06, gx = 1e-02,

    def amp2_xv(s,t,u,mx) -> float:
        mV = mx/3
        gV,gx = 1e-06,1e-02
        Q = gV*gx
        return (s**2 + u**2 + 4*t*(mx**2) - 2*(mx**2)**2)*(Q/(t - mV**2))**2

    Similarily, for DM-e scattering, m1 = mx, m2 = me and kinetic mixing eps = 1e-06

    def amp2_xe(s,t,u,mx) -> float:
        mV = mx/3
        me = constant.me
        gx,eps = 1e-02,1e-06
        Q = gx*eps
        return 2*(s**2 + u**2 + 4*t*(me**2 + mx**2) - 2*(me**2 + mx**2)**2)*(Q/(t - mV**2))**2

    These are the desired amplitudes and serve as the inputs in the class.


    ********************
    *                  *
    *     Methods      *
    *                  *
    ********************

    This class has the following method

         sigma_xe(Tx,mx): Yields the total DM-e cross section for a given (Tx,mx), cm^2
    dsigma_xv(Tx,mx,psi): Yields the differential DM-v cross section at given (Tx,mx,psi), cm^2/sr
           flux(t,Tx,mx): Yields the SNv BDM flux at Earth given (t,Tx,mx), #/MeV/cm^2/s
               event(mx): Yields the SNv BDM event in a given period at Earth given mx, # per electron

    The cross section parts are easy to understand as they are generated by the amplitude-squares
    inputed by user.

    >>> snv1987a.sigma_xe(Tx=15,mx=1e-3)             # cm^2
    3.5462696696305866e-36
    >>> snv1987a.dsigma_xv(Tx=15,mx=1e-3,psi=0.01)   # cm^2/sr
    1.541723841734109e-33

    To generate SNv BDM flux at Earth, the `flux` method requires the BDM ToF t in seconds, BDM
    kinetic energy Tx in MeV and DM mass mx in MeV. Other optional keyword arguments can be checked
    via the associated docstring.

    >>>snv1987a.flux(t=100,Tx=15,mx=1e-3)            # 1/cm^2/MeV/s
    2.4216253535745897e-09

    For BDM event
    
    >>>snv1987a.event(mx=1e-3,neval=50000)           # BDM event number per electron
    1.1068108958437556e-32
    
    Note that to convert event into a real case, we have to multiply the total electron number in the
    detector. For Super-Kamionkande, its total electron number is about 7.4e+33. Thus the total event
    before its vanishing, SK can accumulate around 7.34e+33*1.1068e-32 ~ 81.2 events before SNv BDM
    vanished.
    """
    def __init__(self,SN_coord,GC_coord,amp2_xv,amp2_xe,/,**kwargs): 
        self.SN_coord = SN_coord
        self.GC_coord = GC_coord
        self.amp2_xv = amp2_xv
        self.amp2_xe = amp2_xe
        self.kwargs = kwargs                                    # keyword arguments that will be passed to SkyCoord in Astropy
        self.__dict__.update(self.kwargs)
        
        #*------ Default attributes at the inception -----*#
        self.nx = lambda r,mx: rhox(r,184,24.42,2)/mx           # default profile is MW NFW w/o spike
        self.snNuSpectrum = lambda Ev,D: snNuSpectrum(Ev,D)     # default SN nu spectrum, flux unit

    def __repr__(self):
        return '{:>18s}'.format('Dist to GC:') + ' {:>.3e} kpc'.format(self.Re) + '\n' +                 \
               '{:>18s}'.format('Dist to SN:') + ' {:>.3e} kpc'.format(self.Rstar) + '\n' +              \
               '{:>18s}'.format('Seperation dist:') + ' {:<.3e} kpc'.format(self.separation_3d) + '\n' + \
               '{:>18s}'.format('Off-center angle:') + ' {:<.3e} rad'.format(self.beta) 

    @property
    def separation_3d(self):
        return self._get_geometry_3d()[1]

    @property
    def beta(self):
        return self._get_geometry_3d()[0]

    @property
    def Rstar(self):
        return self.SN_coord[2]

    @property
    def Re(self):
        return self.GC_coord[2]

    def _get_geometry_3d(self) -> float:
        """Get the celestial geometry of SN and GC from user-input"""
        snRA,snDEC,snDist = self.SN_coord  # right ascension, declination, distance in kpc
        gcRA,gcDEC,gcDist = self.GC_coord  # right ascension, declination, distance in kpc
        
        # Call SkyCoord in Astropy to tackle the celestial coordinates of these objects
        sn = SkyCoord(snRA,snDEC,distance=snDist*u.kpc,**self.kwargs)
        gc = SkyCoord(gcRA,gcDEC,distance=gcDist*u.kpc,**self.kwargs)
        # Get the 3d seperation of these two stellar objects, kpc
        sepDist = sn.separation_3d(gc).value
        # Get beta, law of cosine
        cos_beta = (snDist**2 + gcDist**2 - sepDist**2)/2/snDist/gcDist
        beta = arccos(cos_beta)
        return beta,sepDist
        
    def sigma_xe(self,Tx,mx) -> float:
        """Obtain total sigma_xe for a given (Tx,mx), cm^2"""
        me = self.me
        Ex = Tx + mx                       # Total BDM energy, kinetic + mass
        s = 2*me*Ex + mx**2 + me**2
        p_squared = (s - (me+mx)**2)*(s-(me-mx)**2)/4/s
        t_min,t_max = get_tBound(mx,me,Tx) # Allowed t range for dsigma/dt integration
        sigma = quad(lambda t: self.amp2_xe(s,t,2*(mx**2+me**2)-s-t,mx),t_min,t_max)[0]
        sigma = sigma/32/s/p_squared       # integration over azimuthal angle 2pi is incoporated
        sigma = sigma*self.perMeV2cm**2    # convert 1/MeV^2 to cm^2
        return sigma

    def dsigma_xv(self,Tx,mx,psi) -> float:
        """Obtain total sigma_xe for a given (Tx,mx), cm^2"""
        varMandelstam = Mandelstam(0,mx,Tx,psi)  # Get the Mandelstam variables and dLips at psi for (Tx,mx,psi)
        s,t,u,dLips = varMandelstam.s,varMandelstam.t,varMandelstam.u,varMandelstam.get_dLips()
        sigma = self.amp2_xv(s,t,u,mx)*dLips     # Evaluate the differential cross section at psi
        sigma = sigma*self.perMeV2cm**2          # convert 1/MeV^2 to cm^2
        return sigma
            
    def _emissivity(self,Ev,dEv,Tx,mx,psi,r,D) -> float:
        dfv = snNuSpectrum(Ev,D)                 # Get the SNv flux at D
        dsigma = self.dsigma_xv(Tx,mx,psi)       # Get the diff. DM-v cross section at psi
        jx = dfv*dsigma*dEv*self.nx(r,mx)        # Compute the emissivity at boost point
        return jx

    def _diff_flux(self,t,Tx,mx,theta,phi,tau) -> float:
        """Get the BDM differential flux"""
        Rstar,Re,beta = self.Rstar,self.Re,self.beta
        vx = get_vx(Tx,mx)
        
        # propagation geometry
        bdmGeometry = Geometry(t,theta,phi,vx,Rstar,Re,beta)
        d = bdmGeometry.d
        D = bdmGeometry.D 
        rprime = bdmGeometry.rprime
        psi = arccos(bdmGeometry.cosPsi)
        
        # Required SNv energy
        snv = Neutrino(Tx,mx,psi)
        Ev = snv.Ev 
        dEv = snv.dEv 
        is_sanity = snv.is_sanity
        
        if is_sanity:
            jx = self._emissivity(Ev,dEv,Tx,mx,psi,rprime,D)
            # Jacobian
            if ~isclose(0,D,atol=1e-100):
                J = self.c/((d - Rstar*cos(theta))/D + 1/vx)
            else:
                J = self.c*vx
            # BDM flux
            return J*jx*vx*sin(theta)*tau
        else:
            return 0
        
    def flux(self,t,Tx,mx,tau=10,nitn=10,neval=30000) -> float:
        """
        The SNv BDM flux at Earth after integrated over zenith angle theta and
        azimuthal angle phi
    
        Input
        ------
        t: The BDM ToF, relative to the first SN neutrino's arrival
        Tx: BDM kinetic energy, MeV
        mx: DM mass, MeV
        tau: The duration of SN explosion, default 10 s
        nitn: Numer of chains in for iterations, vegas variable, unsigned int
        neval: Number of evaluations in each chain, vegas variable, unsigned int
        
        Output
        ------
        scalar: The diff. BDM flux at Earth, # per MeV per cm^2 per second
    
        See Eq. (6) in Phys. Rev. D 108, 083013 (2023) 
        """
        Rstar,Re,beta = self.Rstar,self.Re,self.beta
        _,t_van = get_tof(Tx,mx,Rstar)                           # get the vanishing time t_van 
        if t <= t_van:
            theta_min,theta_max = get_thetaRange(t,Tx,mx,Rstar)  # get the zenith angle range that contains non-zero BDM flux
            integ = vegas.Integrator([[theta_min,theta_max],[0,2*pi]])  # (theta,phi)
            flux = integ(lambda x: self._diff_flux(t=t,Tx=Tx,mx=mx,theta=x[0],phi=x[1],tau=tau),
                         nitn=nitn,neval=neval).mean
            return flux
        else:                                                    # t > t_van will yield zero BDM
            return 0
    
    def event(self,mx,TxRange=[5,30],tRange=[10,35*constant.year2Seconds],tau=10,
              nitn=10,neval=30000) -> float:
        """
        The SNv BDM evnet at Earth after integrated over exposure time t, BDM
        kinetic energy Tx, zenith angle theta and azimuthal angle phi
    
        Input
        ------
        mx: DM mass, MeV
        TxRange: BDM kinetic energy range of interest, [Tx_min,Tx_max], MeV
        tRange: Detector exposure time, [t_min,t_max], seconds
        tau: The duration of SN explosion, default 10 s
        nitn: Numer of chains in for iterations, vegas variable, unsigned int
        neval: Number of evaluations in each chain, vegas variable, unsigned int
        
        Output
        ------
        scalar: Total SNv BDM event in the given period and energy range
    
        See Eq. (16) in Phys. Rev. D 108, 083013 (2023) and the discussion in
        the maintext.
        """
        Rstar,Re,beta = self.Rstar,self.Re,self.beta
        Tx_min,Tx_max = TxRange
        t_min,t_max = tRange
        
        _,t_van = get_tof(Tx_min,mx,Rstar)                             # get ToF
        theta_min,theta_max = get_thetaRange(t_min,Tx_min,mx,Rstar)    # get the theta range with non-zero BDM flux
    
        
        if t_van <= t_max: # check if user-input maximum exposure time t_max is smaller than the vanishing time
            t_max = t_van  # if so, reset t_max as t_van
        
        if t_min < t_max:  # sometimes the user-input beginning time could be larger than the vanishing time if DM mass is very light 
            integ = vegas.Integrator([[t_min,t_max],[Tx_min,Tx_max],[theta_min,theta_max],[0,2*pi]])  #(t,Tx,theta,phi)
            event = integ(lambda x: self._diff_flux(t=x[0],Tx=x[1],mx=mx,theta=x[2],phi=x[3],tau=tau)*self.sigma_xe(Tx=x[1],mx=mx),
                          nitn=nitn,neval=neval).mean
            return event
        else:
            return 0    
        