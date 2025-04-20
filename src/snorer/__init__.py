# Created by Yen-Hsun Lin (Academia Sinica) in 03/2025.
# Copyright (c) 2025 Yen-Hsun Lin.
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
__version__      = '2.1.0'
__description__  = 'This package evaluates the time-of-flight signatures of boosted dark matter due to supernova neutrinos from our Milky Way, SN1987a in LMC or from SN in arbitrary distant galaxy'
__author__       = 'Yen-Hsun Lin'
__email__        = 'yenhsun@phys.ncku.edu.tw'
__url__          = 'https://github.com/yenhsunlin/snorer'
__license__      = 'GNU GPL-3.0'
__all__          = ['sn_nu_spectrum','dsigma_xv','emissivity_jx','differential_flux','flux','event',
                    'HaloSpike','rhox','M_sigma','radiusInfluence','radiusSchwarzschild','nx',
                    'Geometry','Propagation',
                    'Kinematics','Mandelstam','Neutrino','get_vx','get_psiMax','get_thetaMax','_get_tof','get_gx','KallenLambda','get_tBound','get_tvan',
                    'Constants','constant',
                    'params',
                    'BoostedDarkMatter','galactic_to_beta','equatorial_to_beta',]


#---------- Useful utilities for user ----------#

from .snorerMain import sn_nu_spectrum,dsigma_xv,emissivity_jx,differential_flux,flux,event
from .halo import HaloSpike,rhox,M_sigma,radiusInfluence,radiusSchwarzschild,nx
from .geometry import Geometry,Propagation
from .kinematics import Kinematics,Mandelstam,Neutrino,get_vx,get_psiMax,get_thetaMax,_get_tof,get_gx,get_tvan,get_tBound,KallenLambda
from .constants import Constants,constant
from .params import params
from .utils import BoostedDarkMatter,galactic_to_beta,equatorial_to_beta