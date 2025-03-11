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

__all__ = ['BoostedDarkMatter',
           'galactic_to_beta',
           'equatorial_to_beta',]


#---------- Import required utilities ----------#

from numpy import pi,sin,cos,arccos,asarray,clip
from scipy.integrate import quad
import vegas
from astropy.coordinates import SkyCoord
from .snorerMain import sn_nu_spectrum
from .kinematics import Mandelstam,Neutrino,get_vx,get_thetaMax,_get_tof,get_tBound
from .geometry import Propagation
from .halo import rhox,HaloSpike
from .constants import Constants,constant
from .params import params
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
1. BoostedDarkMatter

Functions
------
1. galactic_to_beta
2. equatorial_to_beta

The docstrings should be sufficient for their self-explanations
"""


class BoostedDarkMatter(Constants):
    """
    Superclass: Constants
    
    Class with medoths that evaluate SNv BDM coming from supernova in arbitrary distant
    galaxy with DM-v and DM-e interaction cross sections descrbied by a specific particle
    model.

    /*-----------------------------------------------------------------------------*/

    ********************
    *                  *
    *   Class Inputs   *
    *                  *
    ********************
    
    Rs : float
        Distance from Earth to SN, kpc.
    Rg : float
        Distance from Earth to the center of a distant galaxy, kpc.
    amp2_xv : func
        Amplitude squared for DM-v interaction, 4 positioning arguments.
        amp2_xv = some_func(s,t,u,mx): the first 3 are Mandelstam variables and the last
        one is the DM mass.
    amp2_xe : func
        Identical to amp2_xv, but is for DM-e interaction.
    is_spike : bool
        Is spike feature included? Default is False.
    **kwargs
        Keyword arguments for characteristic parameters of NFW profile and spike halo. If
        'is_spike = False', the parameters for configuring spiky halo will not be used.
        Default values assume Milky Way.

    ********************
    *                  *
    *    Attributes    *
    *                  *
    ********************

    Suppose, we want to know the SN1987a in LMC, we can extract the equatorial coordinates
    by constant.LMC_coord and constant.SN1987a_coord. They documents [RA,DEC,dist] of the
    stellar objects. One can further use 'snorer.equatorial_to_beta' to get the off-cetner
    angle beta. Thus, we have

    >>> Rs, Rg, beta = 49.97, 51.4, 0.02004
    >>> sn1987a = snorer.BoostedDarkMatter(Rs,Rg,beta,amp2_xv,amp2_xe)

    When an instance is initialized, it also automatically assigns DM halo profile. 
    
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

    This class has the following methods

                nx(r,mx): Yields DM number density at place distant r to GC, 1/cm^3
         sigma_xe(Tx,mx): Yields total DM-e cross section for a given (Tx,mx), cm^2
    dsigma_xv(Tx,mx,psi): Yields differential DM-v cross section at given (Tx,mx,psi), cm^2/sr
           flux(t,Tx,mx): Yields SNv BDM flux at Earth given (t,Tx,mx), #/MeV/cm^2/s
               event(mx): Yields SNv BDM event in a given period at Earth given mx, # per electron

    The cross sections are easy to understand as they are generated by the amplitude-squared inputed
    by the user.

    >>> snv1987a.sigma_xe(Tx=15,mx=1e-3)             # cm^2
    3.5462696696305866e-36
    >>> snv1987a.dsigma_xv(Tx=15,mx=1e-3,psi=0.01)   # cm^2/sr
    1.541723841734109e-33

    To generate SNv BDM flux at Earth, the 'flux' method requires the a particular time t in seconds,
    BDM kinetic energy Tx in MeV and DM mass mx in MeV. Other optional keyword arguments can be
    checked via the associated docstring.

    >>>snv1987a.flux(t=100,Tx=15,mx=1e-3)            # 1/cm^2/MeV/s
    2.4216253535745897e-09

    For BDM event
    
    >>>snv1987a.event(mx=1e-3,neval=50000)           # BDM event number per electron
    1.1068108958437556e-32
    
    Note that to convert event into a real case, we have to multiply the total electron number in the
    detector. For Super-Kamionkande, its total electron number is about 7.4e+33. Thus the total event
    before its vanishing, SK can accumulate around 7.34e+33*1.1068e-32 ~ 81 event before SNv BDM
    vanished.
    """
    def __init__(self,Rs,Rg,beta,amp2_xv,amp2_xe,is_spike=False,**kwargs): 
        self.Rs = Rs
        self.Rg = Rg
        self.beta = beta
        self.amp2_xv = amp2_xv
        self.amp2_xe = amp2_xe
        self.is_spike = is_spike

        # Setup internal function _nx for DM number density calculation
        if self.is_spike is True:  # turn on spike
            # Extract parameters
            self.rhos,self.rs,self.n,self.mBH,self.tBH,self.rh,self.alpha,self.sigv = params.merge('halo','spike',**kwargs).values()
            # Setup profile
            self._nxsp = HaloSpike(self.mBH,self.tBH,self.alpha)
            self._nx = lambda r,mx: self._nxsp(r,mx,self.sigv,self.rhos,self.rs,self.n)
        elif self.is_spike is False:  # turn off spike
            # Extract parameters
            self.rhos,self.rs,self.n = params.merge('halo',**kwargs).values()
            # Setup profile
            self._nx = lambda r,mx: rhox(r,self.rhos,self.rs,self.n)/mx
        else:
            raise FlagError('Argument \'is_spike\' must be a bool.')
        
    def nx(self,r,mx) -> float:
        """
        DM number density
        """
        return self._nx(r,mx)
    
    def sigma_xe(self,Tx,mx) -> float:
        """Obtain total sigma_xe for a given (Tx,mx), cm^2"""
        me = self.me
        Ex = Tx + mx                             # Total BDM energy, kinetic + mass
        s = 2*me*Ex + mx**2 + me**2
        p_squared = (s - (me + mx)**2) * (s - (me - mx)**2)/4/s
        t_min,t_max = get_tBound(Tx,mx,me) # Allowed t range for dsigma/dt integration
        sigma = quad(lambda t: self.amp2_xe(s,t,2 * (mx**2 + me**2) - s - t,mx),t_min,t_max)[0]
        sigma = sigma/32/s/p_squared             # integration over azimuthal angle 2pi is incoporated
        sigma = sigma*self.perMeV2cm**2          # convert 1/MeV^2 to cm^2
        return sigma

    def dsigma_xv(self,Tx,mx,psi) -> float:
        """Obtain diff sigma_xv for a given (Tx,mx,psi), cm^2"""
        varMandelstam = Mandelstam(Tx,0,mx,psi)  # Get the Mandelstam variables and dLips at psi for (Tx,mx,psi)
        s,t,u,dLips = varMandelstam.s,varMandelstam.t,varMandelstam.u,varMandelstam.dLips
        sigma = self.amp2_xv(s,t,u,mx)*dLips     # Evaluate the differential cross section at psi
        sigma = sigma*self.perMeV2cm**2          # convert 1/MeV^2 to cm^2
        return sigma
            
    def _emissivity_jx(self,Tx,Ev,dEv,mx,d,r,psi,d_cut=3.24e-15) -> float:
        """
        Emissivity jx of supernova-neutrino-boost dark matter at boost point.
    
        Parameters
        ----------
        Tx : float
            BDM kinetic energy, MeV.
        Ev : float
            The supernova neutrino energy, MeV.
        dEv : float
            The Jacobian (dEv/dTx)*(vx/c) that converts per netrino energy width, dEv,
            to per BDM kinetic energy width, dTx.
        mx : float
            Dark matter mass, MeV.
        d : float
            Distance from supernova to boost point, kpc.
        r : float
            Distance from galactic center to boost point, kpc.
        psi : float
            The scattering angle in lab frame at boost point, rad. 
        d_cut: scalar
            Terminating point for d. Below the value will return 0.
            Default is 3.24e-15 kpc, approximating 100 km, the size of neutrino sphere.

        Returns
        -------
        out : scalar
            BDM emissivity at boost point along the direction psi, 1/MeV/cm^3/s/sr
    
        See Eq. (13) in BDM Physics for detail.
        """
        # Mutiplied by 10 is becasue we assume total energy released in 10 seconds
        # thus L_tot = E_tot/10. Given we integrate all the contribution within neutrino
        # shell, we should multply it back. See BDM Physics for discussion.
        dfv = 10*sn_nu_spectrum(Ev,d,d_cut,is_density=False)   # SNv flux
        dsigma = self.dsigma_xv(Tx,mx,psi)   # DM-v diff. cross section, cm^2/sr
        ndx = self.nx(r,mx) # DM number density
        # Evaluate BDM emissivity
        return ndx * dfv * dsigma * dEv

    def _differential_flux(self,t,Tx,mx,theta,phi,d_cut=3.24e-15,r_cut=1e-8) -> float:
        """
        The differential supernova-neutrin-boosted dark matter flux at Earth at specific time t 
        and angular direction (theta,phi).
        
        Parameters
        ----------
        t : float
            Time t, relative to the SN neutrino's arrival.
        Tx : float
            BDM kinetic energy, MeV.
        mx : float
            DM mass, MeV.
        theta : float
            The zenith angle theta, rad.
        phi : float
            The azimuthal angle that centers SN, rad.
        d_cut : scalar
            Terminating point for d. Below the value will return 0.
            Default is 3.24e-15 kpc, approximating 100 km, the size of neutrino sphere.
        r_cut : float
            Terminating nx when r' < r_cut, kpc. If one needs to incorporate dark matter spike
            in the central region, r_cut cannot be too large. Otherwise, the spike effect will
            be chopped off before it has any noticeble consequence. Default is 1e-8 kpc.
        **kwargs
            Keyword arguments for characteristic parameters of NFW profile and
            spike halo. If 'is_spike = False', the parameters for configuring
            spiky halo will be deactivated. Default values assume Milky Way.
        
        Returns
        -------
        out : scalar
            The differential BDM flux at Earth, 1/MeV/cm^2/s/sr
    
        See the integrand of Eq. (18) in BDM Physics.
        """
        Rs,Rg,beta = self.Rs,self.Rg,self.beta
        # Dimensionless BDM velocity, vx/c
        vx = get_vx(Tx,mx)
        
        # Initializing Propagation class to account time-dependency
        # in propagation geometry
        bdmProgagation = Propagation(t,vx,theta,phi,Rs,Rg,beta)
        l = bdmProgagation.l  # l.o.s, to evaluate Jacobian J
        d = bdmProgagation.d  # Distance from SN to B, to evaluate neutrino flux
        r = bdmProgagation.rprime  # Distance from GC to B, to evaluate nx
        psi = arccos(bdmProgagation.cos_psi)  # Scattering angle that points Earth direction
        
        # Initializing Neutrino class to get required SNv properties
        snNu = Neutrino(Tx,mx,psi)
        Ev = snNu.Ev  # required SNv energy
        dEv = snNu.dEv * vx  # the Jacobian (dEv/dTx)*(vx/c)
        sanity = snNu.sanity  # is this scattring process allowed?

        # Evaluate differential flux, integrand of Eq. (18) in BDM Physics
        if sanity and r > r_cut: # if energy conservation holds
            # Evaluate the BDM emissivity
            jx = self._emissivity_jx(Tx,Ev,dEv,mx,d,r,psi,d_cut)
            # # Jacobian, it should not diverge as we already require d > d_trunct
            J = d * vx / (vx * (l - Rs * cos(theta)) + d) * constant.c
            # Differential flux 
            diff_flux = J * jx * sin(theta)
        else:
            diff_flux = 0
        return diff_flux
        
    def flux(self,t,Tx,mx,**kwargs) -> float:
        """
        The supernova-neutrino-boosted dark matter flux at time t on Earth after integrated over
        a field-of-view dOmega. Note that zenith angle theta is integrated up to thetaMax
        and azimuthal angle phi from 0 to 2pi.
    
        Parameters
        ----------
        t : float
            The BDM ToF, relative to the first SN neutrino's arrival
        Tx : float
            BDM kinetic energy, MeV
        mx : float
            DM mass, MeV
        **kwargs
            Optional parameters for mimimum distances and vegas
        
        Returns
        -------
        out : scalar
            The time-depenent boosted dark matter flux at Earth, 1/MeV/cm^2/s
    
        See Eq. (18) in BDM Physics 
        """     
        d_cut,r_cut,nitn,neval = params.merge('min_distance','vegas',**kwargs).values()
        Rs = self.Rs
        def diff_flux(x):
            """
            The integrand in (18) with only theta and phi are left as inputs.
            Note that theta = x[0] and phi = x[1]. This matches the vegas
            inputs.
            """
            theta, phi = x
            df = self._differential_flux(t=t,Tx=Tx,mx=mx,theta=theta,phi=phi,d_cut=d_cut,r_cut=r_cut)
            return df
        
        # Get vanishing time
        _,t_van = _get_tof(Tx,mx,Rs)
        # Evaluate flux
        if t <= t_van:
            theta_max = get_thetaMax(t,Tx,mx,Rs)  # get the thetaMax
            integ = vegas.Integrator([[0,theta_max],[0,2*pi]])  # (theta,phi)
            flux = integ(diff_flux,nitn=nitn,neval=neval).mean
            return flux
        else:
           # t > t_van will yield zero BDM
           return 0.0
    
    def event(self,mx,Tx_range=[5,30],t_range=[10,35*constant.year2Seconds],**kwargs) -> float:
        """
        The supernova-neutrino-boosted dark matter evnet per electron. To retrieve the correct
        event number, one should mutiply the total electron number Ne.

        For instance, if the BDM event rate obtained from this function is Nx0, then the total
        BDM event in a detector with electron number is

            Nx = Ne * Nx0
    
        Parameters
        ----------
        mx : float
            DM mass, MeV
        Tx_range : list
            Integration range for BDM kinetic energy [Tx_min,Tx_max], MeV
        t_range : list
            Integration range for exposure time [t_min,t_max], seconds
        **kwargs
            Optional parameters for mimimum distances and vegas
        
        Returns
        -------
        out : scalar
            Event number of supernova-neutrino-boosted dark matter per electron.
        """
        d_cut,r_cut,nitn,neval = params.merge('min_distance','vegas',**kwargs).values()
        Rs = self.Rs
        def diff_event(x):
            """
            The integrand in (18) with t, Tx, theta and phi are left as inputs.
            Note that t = x[0], Tx = [1], theta = x[2] and phi = x[3]. This
            matches the vegas inputs.
            """
            t, Tx, theta, phi = x
            df = self._differential_flux(t=t,Tx=Tx,mx=mx,theta=theta,phi=phi,d_cut=d_cut,r_cut=r_cut)
            return df * self.sigma_xe(Tx,mx)

        # Integration range for Tx
        Tx_min,Tx_max = Tx_range
        # Integration range for t
        t_min,t_max = t_range
        
        _,t_van = _get_tof(Tx_min,mx,Rs)   # get vanishing time
        theta_max = get_thetaMax(t_min,Tx_min,mx,Rs)   # get the thetaMax
        
        if t_van <= t_max: # check if user-input maximum exposure time t_max is smaller than the vanishing time
            t_max = t_van  # if so, reset t_max as t_van
        
        if t_min < t_max:  # sometimes the user-input beginning time could be larger than the vanishing time if DM mass is very light 
            integ = vegas.Integrator([[t_min,t_max],[Tx_min,Tx_max],[0,theta_max],[0,2*pi]])  #(t,Tx,theta,phi)
            event = integ(diff_event,nitn=nitn,neval=neval).mean
            return event
        else:
            return 0.0


def galactic_to_beta(l,b,GC_coord=[0,0]):
    """
    Transform galactic coordinate to off-center angle beta

    Parameters
    ----------
    l : array_like
        Galactic longitude, rad
    b : array_like
        Galactic latitude, rad
    GC_coord : list
        Galactic coordinate for arbitrary galactic center [lg,bg].
        Default is Milky Way center [lg,bg] = [0,0]

    Returns
    -------
    out : scalar/ndarray
        Off-center angle beta. The result is scalar if all
        inputs are scalars.
    """
    lg,bg = GC_coord
    l,b = asarray((l,b))
    cos_beta = cos(b) * cos(bg) * cos(l - lg) + sin(b) * sin(bg)
    cos_beta = clip(cos_beta,-1,1)
    return arccos(cos_beta)


def equatorial_to_beta(ra,dec,GC_coord=None):
    """
    Transform equatorial coordinate to off-center angle and
    galactic coordinate (beta,l,b)

    Parameters
    ----------
    ra : array_like
        Right ascension, hms in string type. Eg. '5h6.7m4.4s'.
    dec : array_like
        Declination, dms in string type.  Eg. '6d10.7m9.4s'.
    GC_coord : None/list
        The equatorial coordinate for arbitrary galactic center.
        Default is None, which automatically implements our
        Milky Way center. For a specific GC coordinate, it should
        have GC_coord = [RA,DEC] where RA and DEC are, similar to 
        'ra' and 'dec', in hms and dms units respectively. Additionally,
        they should be subject to ICRS J2000.0.

    Returns
    -------
    out : tuple
        Tuple of (beta,l,b). Each component is scalar if all
        inputs are scalars
    """
    # Galactic coordinate for GC
    if GC_coord is not None:
        ra_g, dec_g = GC_coord
        gc_coord = SkyCoord(ra = ra_g,dec = dec_g)
        lg = gc_coord.galactic.l.radian
        bg = gc_coord.galactic.b.radian
    else:
        lg,bg = 0,0
    # Galactic coordinate for SN
    eq_coords = SkyCoord(ra = ra,dec = dec)
    # Get galactic coordinate
    l,b = eq_coords.galactic.l.radian,eq_coords.galactic.b.radian
    beta = galactic_to_beta(l,b,GC_coord = [lg,bg])
    return beta,l,b