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



# snorer.nx


###  <span class="mono">snorer.nx(*r*,*mx*,*is_spike=False*,*\*\*kwargs**)</span>

Dark matter number density of Milky Way of Large Magellanic Cloud at distance $r$ to the galactic center.  Spike feature is not included.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `r` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance to galactic center $r$, kpc

> `mx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV

> `is_spike` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Is halo spike included? Default is `False`.

> ***`**kwargs`*** <br>&nbsp;&nbsp;&nbsp;&nbsp; Keyword arguments for characteristic parameters of NFW profile and spike halo, . If `is_spike = False`, the parameters for configuring spiky halo will be deactivated. Default values assume Milky Way. See default arguments in [`snorer.params.halo`](../params/params.md#snorerparamshalo){:target="_blank"} and [`snorer.params.spike`](../params/params.md#snorerparamsspike){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter number density at $r$, cm<sup>−3</sup>

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

Let's plot $n_\chi(r)$ for different profiles.

```python
import numpy as np
import matplotlib.pyplot as plt
import snorer as sn

mx = 0.01
# radius, kpc
r_vals = np.logspace(-3,2,100)
# profiles
profiles = [sn.constant.MW_profile,sn.constant.LMC_profile]
labels = ['MW','LMC']

# Make plot
for i in range(2):
    rhos,rs,n,_,_ = profiles[i].values()
    nx_vals = sn.nx(r_vals,mx,rhos=rhos,rs=rs,n=n)
    plt.plot(r_vals,nx_vals,label=labels[i])
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$r$ [kpc]')
plt.ylabel(r'$n_\chi(r)$ [cm$^{-3}$]')
plt.title(fr'$m_\chi = {mx:.2f}$ MeV')
plt.legend()
plt.show()
```
<figure>
<center><img src="../../../figs/nx2.svg" alt="nx2" style="width: 40%;">
</figure>

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

Number density is just density divided by mass,
$$
n_\chi(r)=\frac{\rho_\chi(r)}{m_\chi}.
$$
See also [**snorer.rhox**](rhox.md){:target="_blank"}.
