<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>



# snorer.get_gx


###  snorer.get_gx(*Ev*, *mx*, *psi*)

Calculate the probability density for cross section at scattering
angle psi and averaged over azimuthal angle in lab frame. This is
for energy-independent cross section.
See Eq. (3) in [User Manual/Physics Overview <i class="fa-regular fa-bookmark"></i>](../../manual/overview.md#particle-kinematics){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `Ev` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;The incoming neutrino energy $E_\nu$, MeV

> `mx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM mass $m_\chi$, MeV

> `psi` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Lab frame scattering angle $\psi$, rad


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Probability density for cross section at $\psi$ and averaged over azimuthal angle $2\pi$. The result is a scalar if the three inputs are all scalars. The unit is sr<sup>−1</sup>

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

In this example, we show $2\pi g_\chi\sin\psi$ vs. $\psi$ for various $m_\chi$.

```python
import numpy as np
import matplotlib.pyplot as plt
import snorer as sn

# Neutrino energy, mx values and psi range
Ev = 10
mx_vals = np.logspace(-3,0,4)
psi_vals = np.linspace(0,np.pi/2,500)

# Draw gx plots for various mx
for mx in mx_vals:
    dOmega = 2*np.pi*np.sin(psi_vals)
    gx_vals = get_gx(Ev,mx,psi_vals)*dOmega
    plt.plot(psi_vals,gx_vals,label=fr'$m_\chi={1000*mx:.0f}$ keV')
plt.yscale('log')
plt.ylim(9.5e-3,)
plt.xlabel(r'$\psi$ [rad]')
plt.ylabel(r'$2\pi g_\chi\sin\psi$')
plt.legend()
plt.show()
```
<figure id="gx">
<center><img src="../../../figs/gx.svg" alt="scheme" style="width: 50%;">
</figure>


