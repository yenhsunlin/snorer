<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>



# snorer.nx


###  snorer.nx(*r*, *mx*, *profile = 'MW'*)

Dark matter number density of Milky Way of Large Magellanic Cloud at distance $r$ to the galactic center.  Spike feature is not included.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `r` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance to galactic center $r$, kpc

> `mx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV

> `profile` : *str* <br>&nbsp;&nbsp;&nbsp;&nbsp;`'MW'` or `'LMC'`, stands for MW halo or LMC halo

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter number density at $r$, 1/cm^3

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

Let's plot $n_\chi(r)$ for different profiles.

```python
import numpy as np
import matplotlib.pyplot as plt
import snorer as sn

# DM mass, keV
mx = 0.01
# radius, kpc
r_vals = np.logspace(-3,2,100)
# profiles
profiles = ['MW','LMC']

# Make plot
for profile in profiles:
    nx_vals = nx(r_vals,mx,profile=profile)
    plt.plot(r_vals,nx_vals,label=profile)
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
