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


# snorer.emissivity_jx


###   <span class="mono">snorer.emissivity_jx(*Ev*,*dEv*,*mx*,*d*,*r*,*psi*,*sigxv0=1e-45*,*profile='MW'*,*d_trunct=3.24e-15*,*is_spike=False*,*sigv=None*,*tBH=1e10*,*alpha='3/2'*)</span>

Emissivity $j_\chi$ of supernova-neutrino-boost dark matter at boost point.
See Eq. (13) in [BDM Physics](../../manual/overview.md#emissivity-on-the-shell){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `Ev` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;The supernova neutrino energy, MeV.

> `dEv` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;The Jacobian $(dE_\nu/dT_\chi)(v_\chi/c)$ that converts per netrino energy width, $dE_\nu$, to per BDM kinetic energy width, $dT_\chi$.

> `mx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

> `d` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from supernova to boost point, kpc.

> `r` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from galactic center to boost point, kpc.

> `psi` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;The scattering angle in lab frame at boost point, rad. 

> `sigxv0` : float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Total DM-$\nu$ cross section, cm<sup>2</sup>. It will be multiplied by **snorer.get_gx** to account for the angular distribution and makes it cm<sup>2</sup> sr<sup>−1</sup>.

> `profile` : *str* <br>&nbsp;&nbsp;&nbsp;&nbsp;`'MW'` or `'LMC'` stands for Milky Way or Large Magellanic Cloud profile in use.

> `d_trunct` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;Truncation point for $d$ to prevent supernova neutrino flux diverges at $d\to 0$. Default is $3.24\times10^{-15}$ kpc, approximating 100 km.

> `is_spike` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Is halo spike included? Default is `False`.

> `sigv` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter annihilation cross section, in the unit of $10^{-26}$ cm<sup>3</sup> s<sup>−1</sup>. `None` indicates no annihilation. It is disregarded if `is_spike = False`.

> `tBH` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Age of supermassive black hole in the galactic center, years. It is disregarded if `is_spike = False`.

> `alpha` : *str* <br>&nbsp;&nbsp;&nbsp;&nbsp;Slope of the spike, `'3/2'` or `'7/3'`. It is disregarded if `is_spike = False`.


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM emissivity at boost point along the direction $\psi$, MeV<sup>−1</sup> cm<sup>−3</sup> s<sup>−1</sup> sr<sup>−1</sup>.

