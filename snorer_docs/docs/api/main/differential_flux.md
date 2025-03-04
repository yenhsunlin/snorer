<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>

<style>
.mono {
    font-family: monospace;
}
</style>


# snorer.differential_flux


###  <span class="mono">snorer.differential_flux(*t*,*Tx*,*mx*,*theta*,*phi*,*Rs*,*beta*,*sigxv0=1e-45*,*profile='MW'*,*d_trunct=3.24e-15*,*r_trunct=1e-5*,*is_spike=False*,*sigv=None*,*tBH=1e10*,*alpha='3/2'*)</span>

The differential supernova-neutrin-boosted dark matter flux at Earth at specific time $t$ and angular direction $(\theta,\varphi)$

$$
\left.\sin\theta\mathcal{J}j_\chi(d,r,T_\chi,\psi)\right|_{t=\frac{d}{c}+\frac{\ell}{v_\chi}-t_\nu}.
$$

This is the integrand of Eq. (18) in [BDM Physics](../../manual/overview.md#from-line-of-sight-to-time-dependency){:target="_blank"}, cf. [Fig. 1](../../manual/overview.md#snv_scheme){:target="_blank"} too.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `t` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Time $t$, relative to the SN$\nu$'s arrival

> `Tx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy, MeV.

> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

> `theta` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The zenith angle $\theta$, rad.

> `phi` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The azimuthal angle $\varphi$ that centers SN, rad.

> `Rs` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from supernova to Earth, kpc.

> `beta` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;The off-center angle, characterizes how SN deviates from GC-Earth axis angularly, rad.

> `sigxv0` : float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Total DM-$\nu$ cross section, cm<sup>2</sup>. It will be multiplied by `snorer.get_gx` to account for the angular distribution and makes it cm<sup>2</sup> sr<sup>−1</sup>.

> `profile` : *str* <br>&nbsp;&nbsp;&nbsp;&nbsp;**'MW'** or **'LMC'** stands for Milky Way or Large Magellanic Cloud profile in use.

> `d_trunct` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;Truncation point for $d$ to prevent supernova neutrino flux diverges at $d\to 0$. Default is $3.24\times10^{-15}$ kpc, approximating 100 km.

> `r_trunct` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;Truncating $n_\chi$ when $r < r_{\rm trunct}$, kpc. Default is $10^{-5}$ kpc.

> `is_spike` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Is halo spike included? Default is `False`.

> `sigv` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter annihilation cross section, in the unit of $10^{-26}$ cm<sup>3</sup> s<sup>−1</sup>. `None` indicates no annihilation. It is disregarded if `is_spike = False`.

> `tBH` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Age of supermassive black hole in the galactic center, years. It is disregarded if `is_spike = False`.

> `alpha` : *str* <br>&nbsp;&nbsp;&nbsp;&nbsp;Slope of the spike, `'3/2'` or `'7/3'`. It is disregarded if `is_spike = False`.


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;The differential BDM flux at Earth, MeV<sup>−1</sup> cm<sup>−2</sup> s<sup>−1</sup> sr<sup>−1</sup>.

