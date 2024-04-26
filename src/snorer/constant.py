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
    def md2MeVperCubicCM(self): # Msun/kpc^3 to MeV/cm^3
        return 3.796e-05        
    
    @property
    def md2MeVperQuadCM(self):  # Msun/kpc^3 to MeV/cm^2
        return 1.171e+17        
    
    @property
    def year2Seconds(self):     # year to seconds
        return 31556926         
    
    @property
    def erg2MeV(self):          # erg to MeV
        return 6.241e+05        
    
    @property
    def kpc2cm(self):           # kpc to cm
        return 3.085e+21        
    
    @property
    def me(self):               # electron mass, MeV
        return 5.110e-01        
    
    @property
    def mn(self):               # neutron mass, MeV
        return 9.395e+02        
    
    @property
    def mp(self):               # proton mass, MeV
        return 9.382e+02        
    
    @property
    def Msun(self):             # Solar mass, MeV
        return 1.115e+60        
    
    @property
    def Msun_kg(self):          # Solar mass, kg
        return 1.981e+30        
    
    @property
    def Mmw(self):              # MW stellar mass, Msun
        return 5.290e+10        
    
    @property
    def M_SgrA(self):           # Sgr A* mass, Msun
        return 4.29e+06         
    
    @property
    def Mhalo(self):            # MW halo mass, Msun
        return 1.290e+12        
    
    @property
    def Rhalo(self):            # MW halo radius, kpc
        return 2.300e+02        

    @property
    def sigma0(self):           # sigma_0, cm^2
        return 1.000e-35        

    @property
    def c(self):                # light speed, cm/s
        return 29979245800      
    
    @property
    def H0(self):               # Hubble constant, km/s/Mpc
        return 7.300e+01        
    
    @property
    def rho_c(self):            # critical density, Msun/pc^3
        return 1.500e-07        
    
    @property
    def Lv(self):               # Supernova neutrino luminosity (single specie), erg/s
        return 3.000e+52/6      
    
    @property
    def Omega_0m(self):         # Cosmological matter fraction
        return 3.150e-01        
    
    @property
    def Omega_0L(self):         # Cosmological dark energy fraction
        return 6.850e-01        
    
    @property
    def Omega_0r(self):         # Cosmological radiation fraction
        return 2.300e-03        
    
    @property
    def Omega_0(self):          # Cosmological total energy
        return 1.000            
    
    @property
    def D_H0(self):             # Mpc
        return 4.280e+03        
    
    @property
    def G(self):                # Newton gravitational constant, pc Msun^-1 (km/s)^2
        return 4.301e-03        
    
    @property
    def SgrA_coord(self):       # Sgr A*'s RA, DEC and distance
        return ['17h45m40.0383s','-29d00m28.069s',8.13]
    
    @property
    def LMC_coord(self):       # LMC's RA, DEC and distance
        return ['05h23m34.5264s', '-69d45m22.053s',49.97]
    
    @property
    def SN1987a_coord(self):   # SN1987a's RA, DEC and distance
        return ['05h35m27.8733s', '-69d16m10.478s',51.7]

    @property
    def MagicalNumber(self):    # Magician does Magic!
        return 2.572e-64        
    