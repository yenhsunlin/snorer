<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>



# snorer.Constants


### *`class`* snorer.Constants() <!-- omit in toc -->

A data class stores many physical constants and coversion factors for convenient use.

<!-- **<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**
> None -->


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Attributes:</div>**
> `perMeV2cm` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Convert MeV to MeV<sup>−1</sup> to cm

> `md2MeVperCubicCM` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Convert $M_\odot$ kpc<sup>−3</sup> to MeV cm<sup>−3</sup>

> `year2Seconds` : *int* <br>&nbsp;&nbsp;&nbsp;&nbsp;Convert one year to seconds

> `erg2MeV` : *int* <br>&nbsp;&nbsp;&nbsp;&nbsp;Convert erg to MeV

> `kpc2cm` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Convert kiloparsec to cm

> `me` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Electron mass, MeV

> `mn` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Neutron mass, MeV

> `mp` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Proton mass, MeV

> `Msun` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Solar mass, MeV

> `Msun_kg` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Solar mass, kg

> `Mmw` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Stellar mass in Milky Way, Msun (halo not included)

> `Mhalo` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter halo mass in Milky Way, Msun

> `Rhalo` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Halo radius, kpc

> `sigma0` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Constant cross section, cm<sup>2</sup>

> `c` : *int* <br>&nbsp;&nbsp;&nbsp;&nbsp;Speed of light, cm s<sup>−1</sup>

> `sigma0` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Constant cross section, cm<sup>2</sup>

> `H0` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Hubble constant, km Mpc<sup>−1</sup> s<sup>−1</sup>

> `rho_c` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Critical densityt $\rho_c$ in cosmology, $M_\odot$ pc<sup>−3</sup>

> `Lv` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Neutrino luminosity, divided by 10 seconds, for single species, total six (3$\nu$ and 3$\bar{\nu}$), erg s<sup>−1</sup>

> `Omega_0m` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Cosmological matter fraction

> `Omega_0L` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Cosmological dark energy fraction

> `Omega_0r` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Cosmological radiation fraction

> `Omega_0` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Cosmological total energy, 1

> `D_H0` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Cosmological distance, Mpc

> `G` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Gravitational constant, pc $M_\odot^{-1}$ kms<sup>2</sup> s<sup>−1</sup>

> `SgrA_coord` : *list* <br>&nbsp;&nbsp;&nbsp;&nbsp;Sagittarius A\* coordinate and distance, [*str*,*str*,*float*] = [RAC,DEC,kpc]

> `LMC_coord` : *list* <br>&nbsp;&nbsp;&nbsp;&nbsp;Coordinate for Large Magellanic Cloud center and distance, [*str*,*str*,*float*] = [RAC,DEC,kpc]

> `SN1987a_coord` : *list* <br>&nbsp;&nbsp;&nbsp;&nbsp;SN1987a coordinate and distance, [*str*,*str*,*float*] = [RAC,DEC,kpc]


**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

Import `snorer` and do
```python
>>> import snorer as sn
>>> sn.Constants().perMeV2cm # Convert 1/MeV to cm
1.973e-11
>>> sn.Constants().SgrA_coord # Sgr A J2000 coordinate and its distance in kpc
['17h45m40.0383s', '-29d00m28.069s', 8.13]
```

For simplicity, one can also call the constants and conversion factors by
```python
>>> sn.constant.perMeV2cm
1.973e-11
```
Both `Constants` and `constant` are equivalent and is related by
`constant = Constants()`.
