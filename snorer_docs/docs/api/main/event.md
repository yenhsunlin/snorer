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


###   <span class="mono">snorer.event(*mx*,*Rs*,*beta*,*Re=8.5*,*Tx_range=[5,30]*,<br>*t_range=[10,1.1045e+09]*,*sigxv0=1e-45*,*is_spike=False*,*\*\*kwargs*)</span>

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

> `Re` : *floate* <br>&nbsp;&nbsp;&nbsp;&nbsp;The distance from GC to Earth, kpc. Default is 8.5 kpc.

> `Tx_range` : *list* <br>&nbsp;&nbsp;&nbsp;&nbsp;Integration range for BDM kinetic energy `[Tx_min,Tx_max]`, MeV

> `t_range` : *list* <br>&nbsp;&nbsp;&nbsp;&nbsp;Integration range for exposure time `[t_min,t_max]`, seconds

> `sigxv0` : float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Total DM-$\nu$ cross section, cm<sup>2</sup>. It will be multiplied by <span class="mono">snorer.get_gx</span> to account for the angular distribution and makes it cm<sup>2</sup> sr<sup>−1</sup>.

> `is_spike` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Is halo spike included? Default is `False`.

> ***`**kwargs`*** <br>&nbsp;&nbsp;&nbsp;&nbsp; Keyword arguments for characteristic parameters of NFW profile and spike halo, . If `is_spike = False`, the parameters for configuring spiky halo will be deactivated. Default values assume Milky Way. See default arguments in [`snorer.params.min_distance`](../params/params.md#snorerparamsmin_distance){:target="_blank"}, [`snorer.params.halo`](../params/params.md#snorerparamshalo){:target="_blank"}, [`snorer.params.spike`](../params/params.md#snorerparamsspike){:target="_blank"} and [`snorer.params.vegas`](../params/params.md#snorerparamsvegas){:target="_blank"}.


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;Event number of supernova-neutrino-boosted dark matter per electron.

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

We also use [**vegas**](https://github.com/gplepage/vegas){:target="_blank"} to evaluate $N_\chi$. See **Notes** in [`snorer.flux`](flux.md){:target="_blank"}.