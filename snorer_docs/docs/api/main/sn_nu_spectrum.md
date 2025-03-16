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


# snorer.sn_nu_spectrum


###  <span class="mono">snorer.sn_nu_spectrum(*Ev*,*d*,*d_cut=3.24e-15*,*is_density=False*)</span>

Supernova neutrino spectrum at distance $d$ to supernova,
$$
\frac{dN_\nu}{dE_\nu}=\sum_i\frac{L_{\nu_i}}{4\pi d^2\langle E_{\nu_i}\rangle} f_{\nu_i}(E_\nu).
$$
See Eqs. (9-12) in [BDM Physics](../../manual/overview.md#snnu-spectrum){:target="_blank"} for detail.
**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `Ev` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Supernova neutrino energy, MeV.

> `d` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from supernova to the boost point, kpc.

> `d_cut` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Terminating point for $d$. Below the value will return 0. Default is $3.24\times 10^{-15}$ kpc, approximating 100 km, the size of neutrino sphere.

> `is_density` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Should convert the output to the unit of number density. Default is `False` and output has the unit of flux.


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Outputs flux [MeV<sup>−1</sup> cm<sup>−2</sup> s<sup>−1</sup>] when `is_density = False`, or number density [MeV<sup>−1</sup> cm<sup>−3</sup>] when `is_density = True`. The output is scalar if all inputs are scalars.

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

In this example, we show $dN_\nu/dE_\nu$ over $(E_\nu,d)$ plane. One can clearly see
that $d<$ `d_cut` the flux is 0.

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import snorer as sn

Ev_vals = np.logspace(-3,2,100) # Ev values
d_vals = np.logspace(-16,2,200) # d values

# Setup meshgrid for (Ev,d) plane
Ev,D = np.meshgrid(Ev_vals,d_vals,indexing='ij')
# Evaluate SNv flux
DNvDEv = sn.sn_nu_spectrum(Ev,D)

# Plot
fig, ax = plt.subplots()
# log-scaler color
norm = mcolors.LogNorm(vmin=DNvDEv.min() + 1, vmax=DNvDEv.max())
# Contour plot
contour = ax.contourf(Ev, D, DNvDEv, levels=20, cmap="viridis", norm=norm)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$E_\nu$ [MeV]')
ax.set_ylabel(r'$d$ [kpc]')
# Color bar
cbar = fig.colorbar(contour, ax=ax)
cbar.set_label(r"$dN_\nu/dE_\nu$ [MeV$^{-1}$ cm$^{-2}$ s$^{-1}$]")
plt.show()
```
<figure id="gx">
<center><img src="../../../figs/dNvdEv.svg" alt="scheme" style="width: 50%;">
</figure>
