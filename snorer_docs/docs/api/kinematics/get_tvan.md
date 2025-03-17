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

# snorer.get_tvan


###  <span class="mono">snorer.get_tvan(*Tx*,*mx*,*Rs*)</span>

Get the BDM vanishing time. The time-zero is set as the arrival of SN$\nu$ at Earth.
See Eqs. (22) and (23) in [BDM Physics](../../manual/overview.md#time-dependent-feature){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `Tx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy $T_\chi$, MeV

> `mx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM mass $m_\chi$, MeV

> `Rs` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance to supernova, $R_s$, kpc


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM vanishing time $t_{\rm van}$, seconds

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

In this example, we show $t_{\rm van}$ on $(m_\chi,T_\chi)$ plane.

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import snorer as sn

Rs = 9.6 # SN distance, kpc
Tx_vals = np.logspace(-6,2,100) # Tx values
mx_vals = np.logspace(-6,3,100) # mx values

# Setup meshgrid for (mx,Tx) plane
MX,TX = np.meshgrid(mx_vals,Tx_vals,indexing='ij')
# Evaluating tvan and convert it to years
TVAN = sn.get_tvan(TX,MX,Rs)/sn.constant.year2Seconds

# Plot
fig, ax = plt.subplots()
# log-scaler color
norm = mcolors.LogNorm(vmin=TVAN.min(), vmax=TVAN.max())
# Contour plot
contour = ax.contourf(MX, TX, TVAN, levels=20, cmap="viridis", norm=norm)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$m_\chi$ [MeV]')
ax.set_ylabel(r'$T_\chi$ [MeV]')
# Color bar
cbar = fig.colorbar(contour, ax=ax)
cbar.set_label(r"$t_{\rm van}$ [yrs]")
plt.show()
```
<figure id="gx">
<center><img src="../../../figs/tvan.svg" alt="scheme" style="width: 50%;">
</figure>

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

Practically speaking, the underlying algorithm of `snorer.get_tvan` is not vectorized. It relies on `numpy.nditer` to support vectorized inputs/outputs.
The kernel of `snorer.get_tvan` is the internal function `snorer._get_tof` which has the same **Parameters** and **Returns** as `snorer.get_tvan` but only accepts scalar inputs/outputs. It could become clumsy when the points to be calcuated are massive.
