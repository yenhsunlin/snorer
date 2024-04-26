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


from numpy import pi,exp,arccos,cos,sin,isclose
import vegas
from .halo import dmNumberDensity
from .kinematics import Kinematics,get_thetaRange,get_tof
from .geometry import Geometry
from .constant import Constants
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
None

Functions
------
1. snNuSpectrum
2. dsigma_xv
3. emissivity
4. diff_flux
5. flux
6. event

The docstrings should be sufficient for their self-explanations
"""


constant = Constants()


def snNuSpectrum(Ev,D,D_trunct=3.24e-15,is_density=False):
    """
    SN neutrino spectrum at D, flux or number density
    
    Input
    ------
    Ev: Neutrino energy in MeV
    D: Distance from the boosted point to the SN, kpc
    D_trunct: Truncation point for D to prevent divergence when D -> 0,
        default is 3.24e-15 kpc ~ 100 km
    is_density: output neutrino number denisty within shell
    
    Output
    ------
    dfv/dEv: SNv flux, #/MeV/cm^2/s, if is_density = False
    dnv/dEv: SNv number density, #/MeV/cm^3, if is_density = True
    """
    if D >= D_trunct:
        Lv = constant.Lv*constant.erg2MeV
        D = D*constant.kpc2cm
        #D = 5*kpc2cm
    
        # Fermi dirac distribution
        def _fv(Ev,Tv):
            exponent = Ev/Tv - 3
            if exponent <= 700:
                return (1/18.9686)*Tv**(-3)*(Ev**2/(exp(exponent) + 1))
            else:
                return 0
        nue_dist = _fv(Ev,2.76)/11
        nueb_dist = _fv(Ev,4.01)/16
        # total 4 species for x
        nux_dist = _fv(Ev,6.26)/25
    
        luminosity = Lv/(4*pi*D**2) 
        flux = luminosity*(nue_dist+nueb_dist+4*nux_dist)
        
        if is_density is False:
            return flux
        elif is_density is True:
            return flux/constant.c
        else:
            raise FlagError('Keyword argument \'is_density\' must be a boolean.')
    else:
        return 0


def dsigma_xv(Ev,mx,psi,sigxv0=1e-45):
    """
    Differential DM-nu scattering angle in lab frame
    
    Input
    ------
    Ev: Neutrino energy, MeV
    mx: DM mass, MeV
    psi: lab frame scattering angle in [0,Pi/2]
    sigxv0: Energy-independent DM-v cross section, cm^2
    
    Output
    ------
    dsigma: cm^2 1/sr
    """ 
    dndOmega = Kinematics.dmdOmega_lab(Ev,mx,psi) # angular distribution in lab frame
    return dndOmega*sigxv0


def emissivity(Ev,dEv,mx,psi,r,D,sigxv0=1e-45,is_spike=True,rh=None,sigv=None,tBH=1e10,profile='MW',alpha='3/2',gamma=1):
    """
    SNv BDM emissivity at upscattered point, note the returned result is divided by sigxv and dimensionless DM velocity
    
    Input
    ------
    Ev: The SN neutrino energy
    dEv: The Jacobian dEv/dTx
    mx: DM mass, MeV
    psi: the scattering angle in lab frame
    r: the distance from the scattering point to GC, in kpc
    D: Distance from the boosted point to the SN, kpc
    sigxv0: DM-nu cross section, cm^2
    is_spike: Turn on/off spike feature, bool
    rh: SMBH influence radius, kpc
        None indicates automatically calculated using the given mBH
    sigv: DM annihilation cross section, in the unit of 1e-26 cm^3/s
        None indicates no annihilation
    tBH: SMBH age, years
    profile: str, 'MW' or 'LMC'
    alpha: Slope of the spike, str type, '3/2' or '7/3'
    gamma: Slope of the initial profile
    
    Output
    ------
    BDM emissivity: # per cm^3 per second


    /*---------------------------*/
    
    Additional Description
    ------

    By default the snNuSpectrum in this function yields flux unit,
    cf. Eq. (5) in Rev. D 108, 083013 (2023). To convert a flux into
    number density we divide it flux by its velocity. In terms of
    neutrinos, it is the speed of light.

    This explains the factor c in front of Eq. (4) as the dnv/dEv is
    expected to be SNv number density within the SNv shell, cf. Fig. 4.
    It is simply (number density) * (cross section) * velocity.

    In numerics, the c in the denominator of SNv flux can be canceled
    by 
    
      c*dnv/dEv*sigma_xv = c*(df/dEv/c)*sigma_xv = df/dEv*sigma_xv.
    
    Thus we remove c in the computation to avoid the round-off error
    during multiplication and division. This does not undermine the
    physics picture behind it!
    """
    dnv = snNuSpectrum(Ev,D)
    dsigma = dsigma_xv(Ev,mx,psi,sigxv0)
    nx = dmNumberDensity(r,mx,is_spike,rh,sigv,tBH,profile,alpha,gamma)
    jx = dnv*dsigma*dEv*nx
    return jx


def diff_flux(t,Tx,mx,theta,phi,Rstar,beta,
              sigxv0=1e-45,Re=8.5,r_cut=1e-8,tau=10,
              is_spike=True,rh=None,sigv=None,tBH=1e10,profile='MW',alpha='3/2',gamma=1):
    """
    The differential SNv BDM flux at Earth
    
    Input
    ------
    t: The BDM ToF, relative to the first SN neutrino's arrival
    Tx: BDM kinetic energy, MeV
    mx: DM mass, MeV
    theta: The open angle theta
    phi: The azimuthal angle along the Earth-SN axis, rad
    Rstar: Distance from Earth to SN, kpc
    beta: The deviation angle, characterizing how SN deviates the GC, rad
    sigxv: DM-neutrino cross section, default 1e-45 cm^2
    Re: The distance from Earth to GC, default 8.5 kpc
    r_cut: Ignore the BDM contribution when r' < r_cut, default 1e-5 kpc
    tau: The duration of SN explosion, default 10 s
    
    Output
    ------
    scalar: The diff. BDM flux at Earth, # per MeV per cm^2 per second per sr

    See the integrand in Eq. (6) in Phys. Rev. D 108, 083013 (2023), cf. Eq. (14)
    """
    vx = Kinematics.get_vx(Tx,mx)        # Get BDM velocity
    
    # Initializing BDM Geometry class
    bdmGeometry = Geometry(t,theta,vx,Rstar,Re,beta,phi)
    d = bdmGeometry.d
    D = bdmGeometry.D 
    rprime = bdmGeometry.rprime
    psi = arccos(bdmGeometry.cosPsi) # Get the required scattering angle
    
    # Initializing BDM Kinematics class
    bdmKinematics = Kinematics(Tx,mx,psi)
    Ev = bdmKinematics.Ev                # Get the required Ev
    dEv = bdmKinematics.dEv              # Get Jacobian dEv/dTx
    is_sanity = bdmKinematics.is_sanity  # check if the combination (Tx,mx,psi) satisfies energy conservation
                                         # True means energy conservation satisfied
                                         # False means energy conservation violated
    
    if rprime >= r_cut and is_sanity:
        # BDM scattering angle
        
        # Evaluate the BDM emissivity
        jx = emissivity(Ev,dEv,mx,psi,rprime,D,sigxv0,is_spike,rh,sigv,tBH,profile,alpha,gamma)
        # Jacobian
        if ~isclose(0,D,atol=1e-100):
            J = constant.c/((d - Rstar*cos(theta))/D + 1/vx)
        else:
            J = constant.c*vx
        # BDM flux
        return J*jx*vx*sin(theta)*tau
    else:
        return 0


def flux(t,Tx,mx,Rstar,beta,
         sigxv0=1e-45,Re=8.5,r_cut=1e-8,tau=10,
         is_spike=True,rh=None,sigv=None,tBH=1e10,profile='MW',alpha='3/2',gamma=1,
         nitn=10,neval=30000) -> float:
    """
    The SNv BDM flux at Earth after integrated over zenith angle theta and
    azimuthal angle phi

    Input
    ------
    t: The BDM ToF, relative to the first SN neutrino's arrival
    Tx: BDM kinetic energy, MeV
    mx: DM mass, MeV
    Rstar: Distance from Earth to SN, kpc
    beta: The deviation angle, characterizing how SN deviates the GC, rad
    sigxv: DM-neutrino cross section, default 1e-45 cm^2
    Re: The distance from Earth to GC, default 8.5 kpc
    r_cut: Ignore the BDM contribution when r' < r_cut, default 1e-5 kpc
    tau: The duration of SN explosion, default 10 s
    
    Output
    ------
    scalar: The diff. BDM flux at Earth, # per MeV per cm^2 per second

    See Eq. (6) in Phys. Rev. D 108, 083013 (2023) 
    """
    _,t_van = get_tof(Tx,mx,Rstar)                           # get the vanishing time t_van 
    if t <= t_van:
        theta_min,theta_max = get_thetaRange(t,Tx,mx,Rstar)  # get the zenith angle range that contains non-zero BDM flux
        integ = vegas.Integrator([[theta_min,theta_max],[0,2*pi]])  # (theta,phi)
        flux = integ(lambda x: diff_flux(t=t,Tx=Tx,mx=mx,theta=x[0],phi=x[1],Rstar=Rstar,beta=beta,
                                                sigxv0=sigxv0,Re=Re,r_cut=r_cut,tau=tau,
                                                is_spike=is_spike,rh=rh,sigv=sigv,tBH=tBH,profile=profile,alpha=alpha,gamma=gamma),
                    nitn=nitn,neval=neval).mean
        return flux
    else:                                                    # t > t_van will yield zero BDM
        return 0


def event(mx,Rstar,beta,TxRange=[5,30],tRange=[10,35*constant.year2Seconds],
          sigxv0=1e-45,Re=8.5,r_cut=1e-8,tau=10,
          is_spike=True,rh=None,sigv=None,tBH=1e10,profile='MW',alpha='3/2',gamma=1,
          nitn=10,neval=30000) -> float:
    """
    The SNv BDM evnet at Earth after integrated over exposure time t, BDM
    kinetic energy Tx, zenith angle theta and azimuthal angle phi

    Input
    ------
    mx: DM mass, MeV
    Rstar: Distance from Earth to SN, kpc
    beta: The deviation angle, characterizing how SN deviates the GC, rad
    TxRange: BDM kinetic energy range of interest, [Tx_min,Tx_max], MeV
    tRange: Detector exposure time, [t_min,t_max], seconds
    sigxv: DM-neutrino cross section, default 1e-45 cm^2
    Re: The distance from Earth to GC, default 8.5 kpc
    r_cut: Ignore the BDM contribution when r' < r_cut, default 1e-5 kpc
    tau: The duration of SN explosion, default 10 s
    
    Output
    ------
    scalar: Total SNv BDM event in the given period and energy range

    See Eq. (16) in Phys. Rev. D 108, 083013 (2023) and the discussion in
    the maintext.
    """
    Tx_min,Tx_max = TxRange
    t_min,t_max = tRange
    
    _,t_van = get_tof(Tx_min,mx,Rstar)                             # get ToF
    theta_min,theta_max = get_thetaRange(t_min,Tx_min,mx,Rstar)    # get the theta range with non-zero BDM flux

    
    if t_van <= t_max: # check if user-input maximum exposure time t_max is smaller than the vanishing time
        t_max = t_van  # if so, reset t_max as t_van
    
    if t_min < t_max:  # sometimes the user-input beginning time could be larger than the vanishing time if DM mass is very light 
        integ = vegas.Integrator([[t_min,t_max],[Tx_min,Tx_max],[theta_min,theta_max],[0,2*_np.pi]])  #(t,Tx,theta,phi)
        event = integ(lambda x: diff_flux(t=x[0],Tx=x[1],mx=mx,theta=x[2],phi=x[3],Rstar=Rstar,beta=beta,
                                          sigxv0=sigxv0,Re=Re,r_cut=r_cut,tau=tau,
                                          is_spike=is_spike,rh=rh,sigv=sigv,tBH=tBH,profile=profile,alpha=alpha,gamma=gamma),
                    nitn=nitn,neval=neval).mean
        return event.mean
    else:
        return 0    