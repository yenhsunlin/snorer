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


###  <span class="mono">snorer.differential_flux(*t*,*Tx*,*mx*,*theta*,*phi*,*Rs*,*beta*,*Re=8.5*,*sigxv0=1e-45*,*is_spike=False*,*\*\*kwargs*)</span>

The differential supernova-neutrin-boosted dark matter flux at Earth at specific time $t$ and angular direction $(\theta,\varphi)$

$$
\left.\sin\theta\mathcal{J}j_\chi(d,r,T_\chi,\psi)\right|_{t=\frac{d}{c}+\frac{\ell}{v_\chi}-t_\nu}.
$$

This is the integrand of Eq. (18) in [BDM Physics](../../manual/overview.md#from-line-of-sight-to-time-dependency){:target="_blank"}, cf. [Fig. 1](../../manual/overview.md#snv_bdm_scheme){:target="_blank"} too.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `t` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Time $t$, relative to the SN$\nu$'s arrival

> `Tx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy, MeV.

> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

> `theta` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The zenith angle $\theta$, rad.

> `phi` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The azimuthal angle $\varphi$ that centers SN, rad.

> `Rs` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from supernova to Earth, kpc.

> `beta` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;The off-center angle, characterizes how SN deviates from GC-Earth axis angularly, rad.

> `Re` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;The distance from GC to Earth, kpc. Default is 8.5 kpc.

> `sigxv0` : float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Total DM-$\nu$ cross section, cm<sup>2</sup>. It will be multiplied by `snorer.get_gx` to account for the angular distribution and makes it cm<sup>2</sup> sr<sup>−1</sup>.


> `is_spike` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Is halo spike included? Default is `False`.

> ***`**kwargs`*** <br>&nbsp;&nbsp;&nbsp;&nbsp; Keyword arguments for characteristic parameters of NFW profile and spike halo, . If `is_spike = False`, the parameters for configuring spiky halo will be deactivated. Default values assume Milky Way. See default arguments in [`snorer.params.min_distance`](../params/params.md#snorerparamsmin_distance){:target="_blank"}, [`snorer.params.halo`](../params/params.md#snorerparamshalo){:target="_blank"} and [`snorer.params.spike`](../params/params.md#snorerparamsspike){:target="_blank"}.



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;The differential BDM flux at Earth, MeV<sup>−1</sup> cm<sup>−2</sup> s<sup>−1</sup> sr<sup>−1</sup>.

