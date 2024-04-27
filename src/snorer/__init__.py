# Created by Yen-Hsun Lin (Academia Sinica) in 03/2024.
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


"""
This is the proxy script
"""


__name__         = 'snorer'
__version__      = '0.0.0'
__description__  = 'This package evaluates the time-of-flight signatures of boosted dark matter due to supernova neutrinos from our Milky Way, SN1987a in LMC or from SN in arbitrary distant galaxy'
__author__       = 'Yen-Hsun Lin'
__email__        = 'yenhsun@phys.ncku.edu.tw'
__url__          = 'https://github.com/yenhsunlin/snorer'
__license__      = 'GNU GPL-3.0'
__all__          = ['snNuSpectrum',
                    'dsigma_xv',
                    'emissivity',
                    'diff_flux',
                    'flux',
                    'event',
                    'haloSpike',
                    'rhox',
                    'M_sigma',
                    'radiusInfluence',
                    'radiusSchwarzschild',
                    'dmNumberDensity',
                    'dmNumberDensity_general',
                    'Geometry',
                    'Neutrino',
                    'get_vx',
                    'get_maxPsi',
                    'get_thetaRange',
                    'get_tof',
                    'fx_lab',
                    'Constants',
                    'constant',
                    'sanity_check',
                    'get_Ev',
                    'get_dEv',]


#---------- Useful utilities for user ----------#

from .snorerMain import snNuSpectrum,dsigma_xv,emissivity,diff_flux,flux,event
from .halo import haloSpike,rhox,M_sigma,radiusInfluence,radiusSchwarzschild,dmNumberDensity,dmNumberDensity_general
from .geometry import Geometry
from .kinematics import Neutrino,get_vx,get_maxPsi,get_thetaRange,get_tof,fx_lab
from .constants import Constants,constant


#---------- Useful instances for user ----------#

sanity_check = Neutrino.sanity_check
get_Ev = Neutrino.get_Ev
get_dEv = Neutrino.get_dEv