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

__all__ = ['sn_nu_spectrum',
           'dsigma_xv',
           'emissivity_jx',
           'differential_flux',
           'flux',
           'event',]


#---------- Import required utilities ----------#

from numpy import pi,exp,arccos,cos,sin,broadcast_arrays,atleast_1d,clip,where
import vegas
from halo import nx,nxSpike
from kinematics import Neutrino,get_gx,get_vx,get_thetaMax,_get_tof
from geometry import Propagation
from constants import constant
from sysmsg import FlagError


##########################################################################
#                                                                        #
#   General Classes and Functions for Numerics                           #
#                                                                        #
##########################################################################


"""
This script contains following classes and functions

Classes
------
None

Functions
------
1. sn_nu_spectrum
2. dsigma_xv
3. emissivity_jx
4. differential_flux
5. flux
6. event

The docstrings should be sufficient for their self-explanations
"""


def sn_nu_spectrum(Ev,d,d_cut=3.24e-15,is_density=False):
    """
    Supernova neutrino spectrum at distance d to supernova.
    
    Parameters
    ----------
    Ev : array_like
        Supernova neutrino energy, MeV
    d : array_like
        Distance from supernova to the boost point, kpc
    d_cut: float
        Terminating point for d. Below the value will return 0.
        Default is 3.24e-15 kpc, approximating 100 km, the size of
        neutrino sphere.
    is_density : bool
        Should convert the output to the unit of number density. Default
        is False and output has the unit of flux.
    
    Returns
    -------
    out : scalar/ndarray
        Supernova neutrino flux, #/MeV/cm^2/s, or number density, #/MeV/cm^3,
        depending on is_density. The output is scalar if all inputs are
        scalars.

    See Eq. (9-12) in BDM Physics for detail
    """
    # --- Preparing for vectorization --- #
    Ev = atleast_1d(Ev) # at least 1d Ev
    Ev,d = broadcast_arrays(Ev,d)
    # Truncate d at d_cut to avoid divergence because 1/d^2
    #d = clip(d,d_cut,None)
    
    # --- Evaluating flux --- #
    Lv = constant.Lv * constant.erg2MeV # convert erg/s to MeV/s
    
    def _fv(Ev,Tv):
        """
        Fermi-Dirac distribution
        """
        # clip exponent to avoid overflow, by clip we mean exponent must < 300
        exponent = clip(Ev/Tv - 3,None,300)
        fv = lambda exponent: (1/18.9686) * Tv**(-3) * (Ev**2 / (exp(exponent) + 1))
        return where(exponent < 300, fv(exponent), 0)
    
    # PDF for all neutrino species
    nue_dist = _fv(Ev,2.76) / 11  # Nu_e
    nueb_dist = _fv(Ev,4.01) / 16 # Nu_e bar
    nux_dist = _fv(Ev,6.26) / 25  # Rest of the 4 species
    # Evaluate the total flux
    lum = lambda d: Lv / (4 * pi * d**2)
    luminosity = where(d > d_cut, lum(d * constant.kpc2cm), 0)
    flux = luminosity * (nue_dist + nueb_dist + 4 * nux_dist)

    # Should the output be flux or number density    
    if is_density is True:
        # Output should be number density
        flux /= constant.c
    return flux if flux.size > 1 else flux.item()


def dsigma_xv(Ev,mx,psi,sigxv0=1e-45):
    """
    Differential DM-nu scattering cross section at angle psi in lab frame.
    
    Parameters
    ----------
    Ev : array_like
        Neutrino energy, MeV
    mx : array_like
        Dark matter mass, MeV
    psi : array_like
        Lab frame scattering angle in [0,pi/2]
    sigxv0 : array_like
        Energy-independent DM-nu cross section, cm^2
        Default is 1e-45 cm^2.
    
    Returns
    -------
    out : scalar/ndarray
        Differential DM-nu cross section at angle psi, cm^2/sr.
        Out is scalar if all inputs are scalars.

    See Eqs. (2) and (3) in BDM Physics.
    """ 
    dsigma_dOmega = sigxv0 * get_gx(Ev,mx,psi)
    return dsigma_dOmega


def emissivity_jx(Ev,dEv,mx,d,r,psi,
                  sigxv0=1e-45,profile='MW',d_cut=3.24e-15,
                  is_spike=False,sigv=None,tBH=1e10,alpha='3/2') -> float:
    """
    Emissivity jx of supernova-neutrino-boost dark matter at boost point.
    
    Parameters
    ----------
    Ev : array_like
        The supernova neutrino energy, MeV.
    dEv : array_like
        The Jacobian (dEv/dTx)*(vx/c) that converts per netrino energy width, dEv,
        to per BDM kinetic energy width, dTx.
    mx : array_like
        Dark matter mass, MeV.
    d : array_like
        Distance from supernova to boost point, kpc.
    r : array_like
        Distance from galactic center to boost point, kpc.
    psi : array_like
        The scattering angle in lab frame at boost point, rad. 
    sigxv0 : float
        Total DM-nu cross section, cm^2. It will be multiplied by snorer.get_gx
        to account for the angular distribution and makes it cm^2/sr.
    profile : str
        'MW' or 'LMC' stands for Milky Way or Large Magellanic Cloud profile in use.
    d_cut: scalar
        Terminating point for d. Below the value will return 0.
        Default is 3.24e-15 kpc, approximating 100 km, the size of neutrino sphere.
    is_spike : bool
        Is halo spike included? Default is False.
    sigv : float
        Dark matter annihilation cross section, in the unit of 1e-26 cm^3/s.
        None indicates no annihilation. It is disregarded if 'is_spike' = False.
    tBH : float
        Age of supermassive black hole in the galactic center, years.
        It is disregarded if 'is_spike' = False.
    alpha : str
        Slope of the spike, '3/2' or '7/3'. It is disregarded if 'is_spike' = False.
    
    Returns
    -------
    out : scalar/ndarray 
        BDM emissivity at boost point along the direction psi, 1/MeV/cm^3/s/sr

    See Eq. (13) in BDM Physics for detail.
    """
    dfv = sn_nu_spectrum(Ev,d,d_cut,is_density=False)   # SNv flux
    dsigma = dsigma_xv(Ev,mx,psi,sigxv0)   # DM-v diff. cross section, cm^2/sr
    
    # Incorporating spike feature?
    if is_spike is True:
        ndx = nxSpike(r,mx,profile,sigv,tBH,alpha)
    elif is_spike is False:
        ndx = nx(r,mx,profile)
    else:
        raise FlagError('\'is_spike\' must be either True or False.')
    
    # Evaluate BDM emissivity
    return ndx * dfv * dsigma * dEv


def differential_flux(t,Tx,mx,theta,phi,Rs,beta,
                      sigxv0=1e-45,profile='MW',Re=8.5,tau=10,d_cut=3.24e-15,r_cut=1e-8,
                      is_spike=False,sigv=None,tBH=1e10,alpha='3/2') -> float:
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
    Rs : float
        Distance from supernova to Earth, kpc.
    beta : float
        The off-center angle, characterizes how SN deviates from GC-Earth axis angularly, rad.
    sigxv0 : float
        Total DM-neutrino cross section, cm^2. Default is 1e-45 cm^2.
    profile : str
        'MW' or 'LMC' stands for Milky Way or Large Magellanic Cloud profile is used.
    Re : float
        The distance from GC to Earth, kpc. Default is 8.5 kpc.
    tau : float
        The duration of SN explosion, seconds. Default is 10 s.
    d_cut : scalar
        Terminating point for d. Below the value will return 0.
        Default is 3.24e-15 kpc, approximating 100 km, the size of neutrino sphere.
    r_cut : float
        Terminating nx when r' < r_cut, kpc. If one needs to incorporate dark matter spike
        in the central region, r_cut cannot be too large. Otherwise, the spike effect will
        be chopped off before it has any noticeble consequence. Default is 1e-8 kpc.
    is_spike : bool
        Whether spike feature is included in nx. Default is False.
    sigv : float 
        DM annihilation cross section, in the unit of 1e-26 cm^3/s.
        None stands for no annihilation. Default is None. It is disregarded if 'is_spike' = False.
    tBH : float
        Age of supermassive black hole, years. Default is 1e+10 years. It is disregarded if
        'is_spike' = False.
    alpha : str
        Slope of the DM spike. Default is '3/2'. It is disregarded if 'is_spike' = False.
    
    Returns
    -------
    out : scalar
        The differential BDM flux at Earth, 1/MeV/cm^2/s/sr

    See the integrand of Eq. (18) in BDM Physics.
    """
    # Dimensionless BDM velocity, vx/c
    vx = get_vx(Tx,mx)
    
    # Initializing Propagation class to account time-dependency
    # in propagation geometry
    bdmProgagation = Propagation(t,vx,theta,phi,Rs,Re,beta)
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
        jx = emissivity_jx(Ev,dEv,mx,d,r,psi,sigxv0,profile,d_cut,is_spike,sigv,tBH,alpha)
        # Jacobian, it should not diverge as we already require d > d_trunct
        J = d * vx / (vx * (l - Rs * cos(theta)) + d) * constant.c
        # Differential flux 
        diff_flux = tau * J * jx * sin(theta)
    else:
        diff_flux = 0
    return diff_flux


def flux(t,Tx,mx,Rs,beta,
         sigxv0=1e-45,profile='MW',Re=8.5,tau=10,d_cut=3.24e-15,r_cut=1e-5,
         is_spike=False,sigv=None,tBH=1e10,alpha='3/2',
         nitn=10,neval=30000) -> float:
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
    Rs : float
        Distance from supernova to Earth, kpc
    beta : float
        The off-center angle, characterizes how SN deviates from GC-Earth axis angularly, rad
    sigxv0 : float
        Total DM-neutrino cross section, cm^2. Default is 1e-45 cm^2
    profile : str
        'MW' or 'LMC' stands for Milky Way or Large Magellanic Cloud profile is used.
    Re : float
        The distance from GC to Earth, kpc. Default is 8.5 kpc.
    tau : float
        The duration of SN explosion, seconds. Default is 10 s.
    d_cut: scalar
        Terminating point for d. Below the value will return 0.
        Default is 3.24e-15 kpc, approximating 100 km, the size of neutrino sphere.
    r_cut : float
        Terminating nx when r' < r_cut, kpc. If one needs to incorporate dark matter spike
        in the central region, r_cut cannot be too large. Otherwise, the spike effect will
        be chopped off before it has any noticeble consequence. Default is 1e-8 kpc.
    is_spike : bool
        Whether spike feature is included in nx. Default is False.
    sigv : float 
        DM annihilation cross section, in the unit of 1e-26 cm^3/s.
        None stands for no annihilation. Default is None. It is disregarded if 'is_spike' = False.
    tBH : float
        Age of supermassive black hole, years. Default is 1e+10 years. It is disregarded if
        'is_spike' = False.
    alpha : str
        Slope of the DM spike. Default is '3/2'. It is disregarded if 'is_spike' = False.
    nitn : int
        Number of chains for vegas to evaluate the integral. Default is 10.
    neval : int
        Number of evaluation number in each chain in vegas. Default is 30000.
    
    Returns
    -------
    out : scalar
        The time-depenent boosted dark matter flux at Earth, 1/MeV/cm^2/s

    See Eq. (18) in BDM Physics 
    """
    
    def diff_flux(x):
        """
        The integrand in (18) with only theta and phi are left as inputs.
        Note that theta = x[0] and phi = x[1]. This matches the vegas
        inputs.
        """
        theta, phi = x[0], x[1]
        df = differential_flux(t=t, Tx=Tx, mx=mx, theta=theta, phi=phi, Rs=Rs, beta=beta,
                               sigxv0=sigxv0, profile=profile, Re=Re, tau=tau, d_trunct=d_cut, r_trunct=r_cut,
                               is_spike=is_spike, sigv=sigv, tBH=tBH, alpha=alpha)
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


def event(mx,Rs,beta,
          Tx_range=[5,30],t_range=[10,35*constant.year2Seconds],
          sigxv0=1e-45,profile='MW',Re=8.5,tau=10,d_cut=3.24e-15,r_cut=1e-5,
          is_spike=False,sigv=None,tBH=1e10,alpha='3/2',
          nitn=10,neval=30000) -> float:
    """
    The supernova-neutrino-boosted dark matter evnet per electron with DM-e cross section
    normalized to 1 cm^2 at Earth. The field-of-view (dOmega) is integrated over entirely
    and the kinetic energy Tx, exposure time t can be integrated within user-defined
    ranges.

    The result is normalized to one electron with DM-e cross section that is normalized to
    1 cm^2. To retrieve the correct event rate 'Nx' inside a detector with total electron
    number 'Ne' and DM-e cross section 'sigma_xe', one should multiply the correction factor

            Nx = Ne * (sigma_xe/1 cm^2).

    For example, SK contains roughly 7.34e+33 electrons and let 'sigma_xe' = 1e-45 cm^2, the correct
    Nx is

            Nx = 7.34e+33 * (1e-45/1)
        
    and note that by default we have used 'sigma_xv' = 1e-45 cm^2.

    Parameters
    ----------
    mx : float
        DM mass, MeV
    Rs : float
        Distance from supernova to Earth, kpc
    beta : float
        The off-center angle, characterizes how SN deviates from GC-Earth axis angularly, rad
    Tx_range : list
        Integration range for BDM kinetic energy [Tx_min,Tx_max], MeV
    t_range : list
        Integration range for exposure time [t_min,t_max], seconds
    sigxv0 : float
        Total DM-neutrino cross section, cm^2. Default is 1e-45 cm^2
    profile : str
        'MW' or 'LMC' stands for Milky Way or Large Magellanic Cloud profile is used.
    Re : float
        The distance from GC to Earth, kpc. Default is 8.5 kpc.
    tau : float
        The duration of SN explosion, seconds. Default is 10 s.
    d_cut: scalar
        Terminating point for d. Below the value will return 0.
        Default is 3.24e-15 kpc, approximating 100 km, the size of
        neutrino sphere.
    r_cut : float
        Terminating nx when r' < r_cut, kpc. If one needs to incorporate dark matter spike
        in the central region, r_cut cannot be too large. Otherwise, the spike effect will
        be chopped off before it has any noticeble consequence. Default is 1e-8 kpc.
    is_spike : bool
        Whether spike feature is included in nx. Default is False.
    sigv : float 
        DM annihilation cross section, in the unit of 1e-26 cm^3/s.
        None stands for no annihilation. Default is None. It is disregarded if 'is_spike' = False.
    tBH : float
        Age of supermassive black hole, years. Default is 1e+10 years. It is disregarded if
        'is_spike' = False.
    alpha : str
        Slope of the DM spike. Default is '3/2'. It is disregarded if 'is_spike' = False.
    
    Returns
    -------
    out : scalar
        Event number of supernova-neutrino-boosted dark matter per electron.
    """

    def diff_flux(x):
        """
        The integrand in (18) with t, Tx, theta and phi are left as inputs.
        Note that t = x[0], Tx = [1], theta = x[2] and phi = x[3]. This
        matches the vegas inputs.
        """
        t, Tx, theta, phi = x
        df = differential_flux(t=t, Tx=Tx, mx=mx, theta=theta, phi=phi, Rs=Rs, beta=beta,
                               sigxv0=sigxv0, profile=profile, Re=Re, tau=tau, d_trunct=d_cut, r_trunct=r_cut,
                               is_spike=is_spike, sigv=sigv, tBH=tBH, alpha=alpha)
        return df

    # Integration range for Tx
    Tx_min,Tx_max = Tx_range
    # Integration range for t
    t_min,t_max = t_range
    
    _,t_van = _get_tof(Tx_min,mx,Rs)   # get vanishing ti,e
    theta_max = get_thetaMax(t_min,Tx_min,mx,Rs)   # get the thetaMax
    
    if t_van <= t_max: # check if user-input maximum exposure time t_max is smaller than the vanishing time
        t_max = t_van  # if so, reset t_max as t_van
    
    if t_min < t_max:  # sometimes the user-input beginning time could be larger than the vanishing time if DM mass is very light 
        integ = vegas.Integrator([[t_min,t_max],[Tx_min,Tx_max],[0,theta_max],[0,2*pi]])  #(t,Tx,theta,phi)
        event = integ(diff_flux,nitn=nitn,neval=neval).mean
        return event
    else:
        return 0.0

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors    
if __name__ == "__main()__":
    Ev_vals = np.logspace(-3,2,100) # Ev values
    d_vals = np.logspace(-6,2,100) # d values

    # Setup meshgrid for (mx,Tx) plane
    Ev,D = np.meshgrid(Ev_vals,d_vals,indexing='ij')
    # Evaluating tvan and convert it to years
    DNvDEv = sn_nu_spectrum(Ev,D)

    # Plot
    fig, ax = plt.subplots()
    # log-scaler color
    norm = mcolors.LogNorm(vmin=DNvDEv.min(), vmax=DNvDEv.max())
    # Contour plot
    contour = ax.contourf(Ev, D, DNvDEv, levels=20, cmap="viridis", norm=norm)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$E_\nu$ [MeV]')
    ax.set_ylabel(r'$d$ [kpc]')
    # Color bar
    cbar = fig.colorbar(contour, ax=ax)
    cbar.set_label(r"$dN_\nu/dE_\nu$ [MeV$^{-1}$ cm$^{-2}$ s$^{-1}$]")
    plt.savefig('dNvdEv.svg',bbox_inches='tight')