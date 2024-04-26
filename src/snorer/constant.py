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


class Constants:
    
    def __init__(self):
        pass
    
    @property
    def md2MeVperCubicCM(self):
        return 3.796e-05        # Msun/kpc^3 to MeV/cm^3
    
    @property
    def md2MeVperQuadCM(self):
        return 1.171e+17        # Msun/kpc^3 to MeV/cm^2
    
    @property
    def year2Seconds(self):
        return 31556926         # year to seconds
    
    @property
    def erg2MeV(self):
        return 6.241e+05        # erg to MeV
    
    @property
    def kpc2cm(self):
        return 3.085e+21        # kpc to cm
    
    @property
    def me(self):
        return 5.110e-01        # electron mass, MeV
    
    @property
    def mn(self):
        return 9.395e+02        # neutron mass, MeV
    
    @property
    def mp(self):
        return 9.382e+02        # proton mass, MeV
    
    @property
    def Msun(self):
        return 1.115e+60        # Solar mass, MeV
    
    @property
    def Msun_kg(self):
        return 1.981e+30        # Solar mass, kg
    
    @property
    def Mmw(self):
        return 5.290e+10        # MW stellar mass, Msun
    
    @property
    def Mhalo(self):
        return 1.290e+12        # MW halo mass, Msun
    
    @property
    def Rhalo(self):
        return 2.300e+02        # MW halo radius, kpc

    @property
    def sigma0(self):
        return 1.000e-35        # sigma_0, cm^2

    @property
    def c(self):
        return 29979245800      # light speed, cm/s
    
    @property
    def H0(self):
        return 7.300e+01        # Hubble constant, km/s/Mpc
    
    @property
    def rho_c(self):
        return 1.500e-07        # critical density, Msun/pc^3
    
    @property
    def Lv(self):
        return 3.000e+52/6      # Supernova neutrino luminosity (single specie), erg/s
    
    @property
    def Omega_0m(self):
        return 3.150e-01        # Cosmological matter fraction
    
    @property
    def Omega_0L(self):
        return 6.850e-01        # Cosmological dark energy fraction
    
    @property
    def Omega_0r(self):
        return 2.300e-03        # Cosmological radiation fraction
    @property
    def Omega_0(self):
        return 1.000            # Cosmological total energy
    
    @property
    def D_H0(self):
        return 4.280e+03        # Mpc
    
    @property
    def G(self):
        return 4.301e-03        # Newton gravitational constant, pc Msun^-1 (km/s)^2

    @property
    def MagicalNumber(self):
        return 2.572e-64        # Magician does Magic!
    