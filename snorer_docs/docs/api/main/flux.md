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


# snorer.flux


###  <span class="mono">snorer.flux(*t*,*Tx*,*mx*,*Rs*,*beta*,*Re=8.5*,*sigxv0=1e-45*,*is_spike=False*)</span>

The supernova-neutrino-boosted dark matter flux at time $t$ on Earth after integrated over
a field-of-view $d\Omega$. Note that zenith angle $\theta$ is integrated up to $\theta^*_M$ and azimuthal angle $\varphi$ from $0$ to $2\pi$.
See Eqs. (18) and (24) in [BDM Physics](../../manual/overview.md#snnu-bdm-flux){:target="_blank"}, cf. [Fig. 1](../../manual/overview.md#snv_scheme){:target="_blank"} too.
**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `t` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Time $t$, relative to the SN$\nu$'s arrival

> `Tx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy, MeV.

> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

> `Rs` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from supernova to Earth, kpc.

> `beta` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;The off-center angle, characterizes how SN deviates from GC-Earth axis angularly, rad.

> `Re` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;The distance from GC to Earth, kpc. Default is 8.5 kpc.

> `sigxv0` : float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Total DM-$\nu$ cross section, cm<sup>2</sup>. It will be multiplied by `snorer.get_gx` to account for the angular distribution and makes it cm<sup>2</sup> sr<sup>−1</sup>.

> `is_spike` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Is halo spike included? Default is `False`.

> ***`**kwargs`*** <br>&nbsp;&nbsp;&nbsp;&nbsp; Keyword arguments for characteristic parameters of NFW profile and spike halo, . If `is_spike = False`, the parameters for configuring spiky halo will be deactivated. Default values assume Milky Way. See default arguments in [`snorer.params.min_distance`](../params/params.md#__attr__-snorerparamsmin_distance){:target="_blank"}, [`snorer.params.halo`](../params/params.md#__attr__-snorerparamshalo){:target="_blank"}, [`snorer.params.spike`](../params/params.md#__attr__-snorerparamsspike){:target="_blank"} and [`snorer.params.vegas`](../params/params.md#__attr__-snorerparamsvegas){:target="_blank"}.



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;The differential BDM flux at Earth, MeV<sup>−1</sup> cm<sup>−2</sup> s<sup>−1</sup> sr<sup>−1</sup>.

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

We use [**vegas**](https://github.com/gplepage/vegas){:target="_blank"} to evaluate the integral Eq. (18). This explains why we incorporate `nitn` and `neval` in keyword arguments. Increasing these values will improve the accuracy but the computation time enhances too.
One may need to find a balance between acceptable accuracy and evaluation time.