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

# snorer.nxSpike


###  <span class="mono">snorer.nxSpike(*r*,*mx*,*profile='MW'*,*sigv=None*,*tBH=1e+10*,*alpha='3/2'*)</span>

Dark matter number density of Milky Way of Large Magellanic Cloud at
    distance $r$ to the galactic center. Spike feature is included.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `r` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance to galactic center, kpc

> `mx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV

> `profile` : *str* <br>&nbsp;&nbsp;&nbsp;&nbsp;`'MW'` or `'LMC'`, stands for MW halo or LMC halo

> `sigv` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM annihilation cross section, in the unit of $10^{-26}$ cm<sup>3</sup> s<sup>−1</sup>. `None` indicates no annihilation

> `tBH` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Supermassive black hole age, years

> `alpha` : *str* <br>&nbsp;&nbsp;&nbsp;&nbsp;Slope of the spike, `'3/2'` or `'7/3'`

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter number density at r with spike in the center, cm<sup>−3</sup>

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

Let's plot $n_\chi$ for different $\langle \sigma v\rangle$.

```python
import numpy as np
import matplotlib.pyplot as plt
import snorer as sn

# DM mass, keV
mx = 0.01
# radius, kpc
r_vals = np.logspace(-5,2,100)
# profiles

sigv_vals = [None,0.01,0.1,3]

for sigv in sigv_vals:
    # calculate nx
    nx_vals = nxSpike(r_vals,mx,sigv=sigv,profile='LMC')
    if sigv is None: sigv = 0 # legend label
    plt.plot(r_vals,nx_vals,label=r'$\langle\sigma v\rangle=$' + str(sigv))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$r$ [kpc]')
plt.ylabel(r'$n_\chi(r)$ [cm$^{-3}$]')
plt.title(fr'LMC with spike and $m_\chi = {mx:.2f}$ MeV')
plt.legend()
plt.show()
```
<figure>
<center><img src="../../../figs/nx_spike.svg" alt="nx2" style="width: 40%;">
</figure>

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

To realize $n_\chi$ with spike feature we initialized a `snorer.HaloSpike` instance inside the function `snorer.nxSpike` and utilize the callable feature. However, such callable function does not support vectorization. To mimic vectorized inputs/outputs, we employ `numpy.nditer`. It could become clumsy if the points to be calculated are massive.