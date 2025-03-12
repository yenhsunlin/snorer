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

# snorer.BoostedDarkMatter


### *`class`* <span class="mono">snorer.BoostedDarkMatter(*Rs*,*Rg*,*beta*,*amp2_xv*,*amp2_xe*,<br>*is_spike=False*,*\*\*kwargs*)</span>

Superclass: `Constants`
    
Class with medoths that evaluate supernova-neutrino-boosted dark matter (SN$\nu$ BDM) coming from SN in arbitrary distant
galaxy with DM-$\nu$ and DM-$e$ interaction cross sections descrbied by a specific particle model. This class integrates functions like `snorer.flux`, `snorer.event`
as methods for user to calculate the SN$\nu$ BDM flux and event associated to any models.

To construct the scattering amplitude from the specific model, we take the model discussed in Ref. [[1](#bib_Lin2023PRD)] for instance. Both
DM-$\nu$ and DM-$e$ have the same amplitude square show by Eq. (3) in Ref. [[1](#bib_Lin2023PRD)] in terms of Mandelstam variables $s$, $t$, and $u$,

$$
|\mathcal{M}|^2 = 2 \left(\frac{\mathcal{Q}}{t} - m_V^2\right)^2 (s^2 + u^2 + 4t (m_1^2 + m_2^2) - 2 (m_1^2 + m_2^2)^2)
$$

where $\mathcal{Q}$ is the multiplication of coupling constants, $m_1$ and $m_2$ are the masses of
incident and target particles respectively and $m_V$ the mediator mass.

Thus we can construct DM-v amplitude square by letting, $m_1 = 0$ and assume $m_V = m_\chi/3$, $g_V = 10^{-6}$ and $g_\chi = 10^{-2}$,

```python
def amp2_xv(s,t,u,mx) -> float:
    mV = mx/3
    gV,gx = 1e-06,1e-02
    Q = gV*gx
    return (s**2 + u**2 + 4*t*(mx**2) - 2*(mx**2)**2)*(Q/(t - mV**2))**2
```

Similarily, for DM-$e$ scattering, $m_1 = m_\chi$, $m_2 = m_e$ and kinetic mixing $\varepsilon = 10^{-6}$,

```python
def amp2_xe(s,t,u,mx) -> float:
    mV = mx/3
    me = constant.me
    gx,eps = 1e-02,1e-06
    Q = gx*eps
    return 2*(s**2 + u**2 + 4*t*(me**2 + mx**2) - 2*(me**2 + mx**2)**2)*(Q/(t - mV**2))**2
```

These are the desired amplitudes and serve as the inputs in the class.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `Rs` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from Earth to SN, kpc.


> `Rg` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from Earth to the center of a distant galaxy, kpc.


> `beta` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Off-center angle, rad.


> `amp2_xv` : *func* <br>&nbsp;&nbsp;&nbsp;&nbsp;Amplitude squared for DM-$\nu$ interaction with 4 positioning arguments. `amp2_xv = func(s,t,u,mx)`: the first 3 are Mandelstam variables and the last one is the DM mass.

> `amp2_xe` : *func* <br>&nbsp;&nbsp;&nbsp;&nbsp;Identical to `amp2_xv`, but is for DM-$e$ interaction.

> `is_spike` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Whether spike feature is included in $n_\chi$. Default is False.

> ***`**kwargs`***  <br>&nbsp;&nbsp;&nbsp;&nbsp; Keyword arguments for characteristic parameters of NFW profile and spike halo, . If `is_spike = False`, the parameters for configuring spiky halo will be deactivated. Default values assume Milky Way. See default arguments in [`snorer.params.halo`](../params/params.md#snorerparamshalo){:target="_blank"} and [`snorer.params.spike`](../params/params.md#snorerparamsspike){:target="_blank"}.

####  <span class="mono">nx(*r*,*mx*)</span>

Method.

Yields DM number density at place distant $r$ to GC. Whether spike feature is on depending the initial setting of `is_spike` when the class instance is initialized.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**
> `r` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from Earth to GC, kpc.

> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter number density at $r$, cm<sup>−3</sup>.

####  <span class="mono">dsigma_xv(*Tx*,*mx*,*psi*)</span>

Method.

Yields differential DM-$\nu$ cross section for a given $(T_\chi,m_\chi,\psi)$ associated with `amp2_xv`.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**
> `Tx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy, MeV.

> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

> `psi` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Scattering angle in lab frame, rad.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Differential DM-$\nu$ cross section, cm<sup>2</sup> sr<sup>−1</sup>.


####  <span class="mono">sigma_xe(*Tx*,*mx*)</span>

Method.

Yields total DM-$e$ cross section for a given $(T_\chi,m_\chi)$ associated with `amp2_xe`.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**
> `Tx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy, MeV.

> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Total DM-$e$ cross section, cm<sup>2</sup>.


####  <span class="mono">flux(*t*,*Tx*,*mx*,*\*\*kwargs*)</span>

Method.

The SN$\nu$ BDM flux at time $t$ on Earth after integrated over
a field-of-view $d\Omega_{\rm lab}$. Note that zenith angle $\theta$ is integrated up to $\theta_{\rm M}^*$
and azimuthal angle $\varphi$ from $0$ to $2\pi$.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `t` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The BDM ToF, relative to the first SN neutrino's arrival.

> `Tx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy, MeV.

> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

> ***`**kwargs`*** <br>&nbsp;&nbsp;&nbsp;&nbsp;Keyword arguments for min distances and vegas. See default arguments in [`snorer.params.min_distance`](../params/params.md#snorerparamsmin_distance){:target="_blank"} and [`snorer.params.vegas`](../params/params.md#snorerparamsvegas){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The time-depenent boosted dark matter flux at Earth, MeV<sup>−1</sup> cm<sup>2</sup> s<sup>−1</sup>.

####  <span class="mono">event(*mx*,*Tx_range=[5,30]*,*t_range=[10,1.1045e+09]*,*\*\*kwargs*)</span>

Method.

The SN$\nu$ BDM evnet *per electron*, $N_{\chi,0}$. To retrieve the correct
event number, one should mutiply the total electron number $N_e$.
For instance, if the BDM event rate obtained from this function is $N_\chi$, then the total BDM event in a detector with electron number  $N_e$ is

$$
N_\chi^{\rm correct} = N_e\times N_{\chi,0}.
$$


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

> `Tx_range` : *list* <br>&nbsp;&nbsp;&nbsp;&nbsp;Integration range for BDM kinetic energy `[Tx_min,Tx_max]`, MeV. Default is `[5,30]`.

> `t_range` : *list* <br>&nbsp;&nbsp;&nbsp;&nbsp;Integration range for exposure time `[t_min,t_max]`, seconds. Default is `[10,1.1045e+09]` and implies `t_max` is around 35 years.

> ***`**kwargs`*** <br>&nbsp;&nbsp;&nbsp;&nbsp;Keyword arguments for min distances and vegas. See default arguments in [`snorer.params.min_distance`](../params/params.md#snorerparamsmin_distance){:target="_blank"} and [`snorer.params.vegas`](../params/params.md#snorerparamsvegas){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Event number of SN$\nu$ BDM per electron.

#### References
1. <p id="bib_Lin2023PRD">Y.-H. Lin *et al.*, *Phys. Rev. D* **108**, 083013 (2023)</p>