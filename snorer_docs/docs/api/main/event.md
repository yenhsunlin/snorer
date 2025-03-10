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


# snorer.event


###   <span class="mono">snorer.event(*mx*,*Rs*,*beta*,*Tx_range=[5,30]*,<br>*t_range=[10,35\*snorer.constant.year2Seconds]*,*sigxv0=1e-45*,*profile='MW'*,*d_cut=3.24e-15*,*r_cut=1e-5*,*is_spike=False*,*sigv=None*,*tBH=1e10*,*alpha='3/2'*,*nitn=10*,*neval=30000*)</span>

The supernova-neutrino-boosted dark matter evnet per electron with DM-$e$ cross section $\sigma_{\chi e}$
normalized to 1 cm<sup>2</sup> at Earth. The field-of-view $d\Omega$ is integrated over entirely
and the kinetic energy $T_\chi$, exposure time $t$ can be integrated within user-defined ranges. Precisely speaking, the event $N_\chi$ is, using Eqs. (18) in [BDM Physics](../../manual/overview.md#snnu-bdm-flux){:target="_blank"},

$$
N_\chi = N_e \sigma_{\chi e} \int_{t_{\rm min}}^{t_{\rm max}} dt \int_{T_{\chi,{\rm min}}}^{T_{\chi,{\rm max}}} dT_\chi \frac{d\Phi_\chi}{dT_\chi}
$$

and `snorer.event` presumes $N_e=1$ and $\sigma_{\chi e}=1$ cm<sup>2</sup>.
One can restore the correct $N_\chi^{\rm correct}$ for any detector by multiplying the true $N_e^{\rm true}$ for that detector and $\sigma_{\chi e}^{\rm true}$,

$$
N_\chi^{\rm correct} = N_\chi \times \frac{N_e}{1}\times \frac{\sigma_{\chi e}^{\rm true}}{1 \,{\rm cm^2}} \times \frac{\sigma_{\chi\nu}^{\rm true}}{10^{-45}\,{\rm cm^2}}
$$

where we have set $\sigma_{\chi \nu}=10^{-45}$ cm<sup>2</sup> by default in the function.
**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

> `Rs` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from supernova to Earth, kpc.

> `beta` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;The off-center angle, characterizes how SN deviates from GC-Earth axis angularly, rad.

> `Tx_range` : *list* <br>&nbsp;&nbsp;&nbsp;&nbsp;Integration range for BDM kinetic energy `[Tx_min,Tx_max]`, MeV

> `t_range` : *list* <br>&nbsp;&nbsp;&nbsp;&nbsp;Integration range for exposure time `[t_min,t_max]`, seconds

> `sigxv0` : float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Total DM-$\nu$ cross section, cm<sup>2</sup>. It will be multiplied by <span class="mono">snorer.get_gx</span> to account for the angular distribution and makes it cm<sup>2</sup> sr<sup>−1</sup>.

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

> `out` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;Event number of supernova-neutrino-boosted dark matter per electron.

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

We also use [**vegas**](https://github.com/gplepage/vegas){:target="_blank"} to evaluate $N_\chi$. See **Notes** in [`snorer.flux`](flux.md){:target="_blank"}.