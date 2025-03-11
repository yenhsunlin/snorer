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
This script stores the default settings for parameters
"""

__all__ = ['params',]

from typing import NamedTuple, Optional, Dict, Any


class MinimumDistanceConfig(NamedTuple):
    """
    This data class stores parameters for truncated distances
    for both dark matter halo and supernova neutrino spectrum.
    """
    # Trucation distance for supernova neutrino spectrum
    # Below this value the spectrum will output 0
    d_cut: float = 3.24e-15 # kpc, ~100 km
    
    # Trucation distance for dark matter halo number density
    # Below this value the density will output 0 when
    # calculates flux and event
    r_cut: float = 1e-8 # kpc


class HaloConfig(NamedTuple):
    """
    This data class stores parameters for configuring halo denisty
    based on NFW profile. Default values are for Milky Way.
    """
    # Characteristic density
    rhos: float = 184 # MeV/cm^3
    # Characteristic distance
    rs: float = 24.42 # kpc
    # profile index
    n: float = 2


class SpikeConfig(NamedTuple):
    """
    This data class stores parameters for configuring halo with
    spike. Only valid when 'is_spike = True'.
    """
    # Supermassive black hole mass
    mBH: float = 4.29e6 # M_sun
    # Supermassive black hole age
    tBH: float = 1e9 # yr
    # Influence radius
    rh: float = 0.002 # kpc
    # Spike index
    alpha: str = '3/2' # str
    # Dark matter annihilation cross section
    # Eg. if type 3 the it indicates 3e-26 cm^3/s
    # None indicates no annihilation
    sigv: Optional[float] = None # 1e-26 cm^3/s


class vegasConfig(NamedTuple):
    """
    This data class stores parameters for configuring vegas.
    """
    # Number of chains in each evaluation
    nitn: int = 10
    # Number of events in each chain
    neval: int = 10000


class ParamsConfig:
    def __init__(self):
        self.min_distance = MinimumDistanceConfig() # mimimum distances
        self.halo = HaloConfig() # halo profile
        self.spike = SpikeConfig() # halo with spike
        self.vegas = vegasConfig() # vegas setup

    def merge(self, *categories: str, **overrides) -> Dict[str, Any]:
        """Merge multiple parameters"""
        merged_config = {}
        for category in categories:
            if hasattr(self, category):
                merged_config.update(getattr(self, category)._asdict())
            else:
                raise ValueError(f"Configuration parameters: '{category}' does not exist!")

        # Ensure all user inputs are legal and exist
        if overrides:
            invalid_keys = set(overrides.keys()) - set(merged_config.keys())
            if invalid_keys:
                raise ValueError(f"Invalid keyword argument: {invalid_keys}, please check again.")
            merged_config.update(overrides)
        
        # Override with legal inputs, if exist
        merged_config.update(overrides)
        return merged_config


# Initialize parameter configuration
params = ParamsConfig()
