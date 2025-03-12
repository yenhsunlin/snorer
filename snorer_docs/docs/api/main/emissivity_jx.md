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


###   <span class="mono">snorer.emissivity_jx(*Ev*,*dEv*,*mx*,*d*,*r*,*psi*,*sigxv0=1e-45*,*d_cut=3.24e-15*,*is_spike=False*,*\*\*kwargs*)</span>

Emissivity $j_\chi$ of supernova-neutrino-boost dark matter at boost point.
See Eq. (13) in [BDM Physics](../../manual/overview.md#emissivity-on-the-shell){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `Ev` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;The supernova neutrino energy, MeV.

> `dEv` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;The Jacobian $(dE_\nu/dT_\chi)(v_\chi/c)$ that converts per netrino energy width, $dE_\nu$, to per BDM kinetic energy width, $dT_\chi$.

> `mx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

> `d` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from supernova to boost point, kpc.

> `r` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from galactic center to boost point, kpc.

> `psi` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;The scattering angle in lab frame at boost point, rad. 

> `sigxv0` : float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Total DM-$\nu$ cross section, cm<sup>2</sup>. It will be multiplied by `snorer.get_gx` to account for the angular distribution and makes it cm<sup>2</sup> sr<sup>−1</sup>.

> `d_cut` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Terminating point for $d$. Below the value will return 0. Default is $3.24\times 10^{-15}$ kpc, approximating 100 km, the size of neutrino sphere.

> `is_spike` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Whether spike feature is included in $n_\chi$. Default is False.


> ***`**kwargs`***  <br>&nbsp;&nbsp;&nbsp;&nbsp; Keyword arguments for characteristic parameters of NFW profile and spike halo, . If `is_spike = False`, the parameters for configuring spiky halo will be deactivated. Default values assume Milky Way. See default arguments in [`snorer.params.halo`](../params/params.md#__attr__-snorerparamshalo){:target="_blank"} and [`snorer.params.spike`](../params/params.md#__attr__-snorerparamsspike){:target="_blank"}.


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM emissivity at boost point along the direction $\psi$, MeV<sup>−1</sup> cm<sup>−3</sup> s<sup>−1</sup> sr<sup>−1</sup>.

