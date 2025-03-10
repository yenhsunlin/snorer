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


###  <span class="mono">snorer.flux(*t*,*Tx*,*mx*,*Rs*,*beta*,*sigxv0=1e-45*,*profile='MW'*,*d_cut=3.24e-15*,*r_cut=1e-5*,*is_spike=False*,*sigv=None*,*tBH=1e10*,*alpha='3/2'*,*nitn=10*,*neval=30000*)</span>

The supernova-neutrino-boosted dark matter flux at time $t$ on Earth after integrated over
a field-of-view $d\Omega$. Note that zenith angle $\theta$ is integrated up to $\theta^*_M$ and azimuthal angle $\varphi$ from $0$ to $2\pi$.
See Eqs. (18) and (24) in [BDM Physics](../../manual/overview.md#snnu-bdm-flux){:target="_blank"}, cf. [Fig. 1](../../manual/overview.md#snv_scheme){:target="_blank"} too.
**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `t` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Time $t$, relative to the SN$\nu$'s arrival

> `Tx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy, MeV.

> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

> `Rs` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from supernova to Earth, kpc.

> `beta` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;The off-center angle, characterizes how SN deviates from GC-Earth axis angularly, rad.

> `sigxv0` : float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Total DM-$\nu$ cross section, cm<sup>2</sup>. It will be multiplied by `snorer.get_gx` to account for the angular distribution and makes it cm<sup>2</sup> sr<sup>−1</sup>.

> `profile` : *str* <br>&nbsp;&nbsp;&nbsp;&nbsp;`'MW'` or `'LMC'` stands for Milky Way or Large Magellanic Cloud profile in use.

> `d_cut` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Terminating point for $d$. Below the value will return 0. Default is $3.24\times 10^{-15}$ kpc, approximating 100 km, the size of neutrino sphere.

> `r_cut` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;Terminating $n_\chi$ when $r^\prime <$ `r_cut`, kpc. If one needs to incorporate dark matter spike in the central region, `r_cut` cannot be too large. Otherwise, the spike effect will be chopped off before it has any noticeble consequence. Default is $10^{-8}$ kpc.

> `is_spike` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Is halo spike included? Default is `False`.

> `sigv` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter annihilation cross section, in the unit of $10^{-26}$ cm<sup>3</sup> s<sup>−1</sup>. `None` indicates no annihilation. It is disregarded if `is_spike = False`.

> `tBH` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Age of supermassive black hole in the galactic center, years. It is disregarded if `is_spike = False`.

> `alpha` : *str* <br>&nbsp;&nbsp;&nbsp;&nbsp;Slope of the spike, `'3/2'` or `'7/3'`. It is disregarded if `is_spike = False`.

> `nitn` : *int* <br>&nbsp;&nbsp;&nbsp;&nbsp;Number of chains in [**vegas**](https://github.com/gplepage/vegas){:target="_blank"}  to evaluate the integral. Default is 10.

> `neval` : *int* <br>&nbsp;&nbsp;&nbsp;&nbsp;Number of evaluation number in each chain in [**vegas**](https://github.com/gplepage/vegas){:target="_blank"} . Default is 30000.


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;The differential BDM flux at Earth, MeV<sup>−1</sup> cm<sup>−2</sup> s<sup>−1</sup> sr<sup>−1</sup>.

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

We use [**vegas**](https://github.com/gplepage/vegas){:target="_blank"} to evaluate the integral Eq. (18). This explains why we incorporate `nitn` and `neval` as the parameters. Increasing these values will improve the accuracy but the computation time enhances too.
One may need to find a balance between acceptable accuracy and evaluation time.