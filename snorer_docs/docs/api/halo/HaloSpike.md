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


# snorer.HaloSpike


### *`class`*  <span class="mono">snorer.HaloSpike(*mBH*,*tBH*,*alpha*)</span>

Superclass: `snorer.Constants`

Class for constructing dark matter halo with spike due to supermassive black hole
(SMBH) in the galactic center.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `mBH` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;SMBH mass, $M_\odot$


> `tBH` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;SMBH age, years


> `alpha` : *str* <br>&nbsp;&nbsp;&nbsp;&nbsp;Slope of the spike profile, only `'3/2'` or `'7/3'` is acceptable




**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Attributes:</div>**

> `mBH` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;SMBH mass, user's input


> `tBH` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;SMBH age, user's input


> `alpha` : *obj* <br>&nbsp;&nbsp;&nbsp;&nbsp;Slope of the spike profile, user's input and is a `Fraction` object

> `rh` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp; SMBH influence radius $r_h$, kpc

> `Rsp` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp; Spike radius $R_{\rm sp}$, kpc

### *`__call__`*  <span class="mono">(*r*,*mx*,*sigv*,*rhos*,*rs*,*n*)</span>

After initializing `snorer.HaloSpike` instance, it is callable like normal function with the following required inputs. 

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `r` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance to GC, kpc


> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM mass, MeV


> `sigv` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp; DM annihilation cross section $\langle\sigma v\rangle$ in the unit of cm<sup>3</sup> s<sup>−1</sup>. `None` means no annihilation

> `rhos` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Charasteristic density $\rho_s$, MeV s<sup>−3</sup>

> `rs` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Characteristic radius $r_s$, kpc

> `n` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Slope of the DM profile outside $R_{\rm sp}$


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM number density at $r$, cm<sup>−3</sup>


**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

Initializing instance and check its attributes.

```python
>>> import snorer as sn
>>> nx = sn.HaloSpike(mBH=1e7,tBH=1e10,alpha='3/2') # initializing instance
>>> nx # print instance information
       SMBH mass: 1.000e+07 M_sun
     Spike slope: 3/2
   Initial slope: 1.000e+00
    Spike radius: 9.504e-10 kpc
Influence radius: 3.412e-03 kpc
>>> nx.alpha # check alpha
Fraction(3,2)
```
The influence radius $r_h$ is auto generated but can be replaced by user defined number. It can be reset to the default value by giving `None`.
```python
>>> nx.rh # default value of rh
0.003411783398329804
>>> nx.rh = 2.54e-3 # replace rh with user-defined value
>>> nx.rh
2.54e-03
>>> nx.rh = None # reset it to the default value
>>> nx.rh
0.003411783398329804
```
Make it a callable function that can calculate DM number density at different $r$.

```python
import numpy as np
import matplotlib.pyplot as plt
import snorer as sn

# Get MW rhos, rs, n, mBH and rh
rhos,rs,n,mBH,rh = constant.MW_profile.values()
# Assuming BH age is 1 Gyr
tBH = 1e10
# DM mass, MeV
mx = 0.1
# Annihilation cross section
sigv_list = [3,0.03,None]
sigv_label = ['3','0.03','0']
# Initializing instances with two different alphas
nx = HaloSpike(mBH=mBH,tBH=tBH,alpha='3/2')  # alpha = 3/2

# radius, kpc
r_vals = np.logspace(-5,2,100)

for i in range(3):
    sigv = sigv_list[i]
    nx_vals = [nx(r,mx,sigv,rhos,rs,n) for r in r_vals]
    plt.plot(r_vals,nx_vals,label=sigv_label[i] + r'$\times10^{-26}\,{\rm cm^3~s^{-1}}$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$r$ [kpc]')
plt.ylabel(r'$n_\chi(r)$ [cm$^{-3}$]')
plt.title(fr'$m_\chi = {mx:.1f}$ MeV')
plt.legend()
plt.show()
```
<figure id="22scat">
<center><img src="../../../figs/nx.svg" alt="nx" style="width: 50%;">
</figure>
### References
1. P. Gondolo and J. Silk, *Phys. Rev. Lett.* **83**, 1719 (1999)
2. J. Cline and M. Puel, *JCAP* **06**, 004 (2023)