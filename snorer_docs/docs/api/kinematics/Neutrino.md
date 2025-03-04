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


# snorer.Neutrino


### *`class`* <span class="mono">snorer.Neutrino(*Tx*,*mx*,*psi*)</span>

Superclass: `snorer.Kinematics`

This class constructs the required neturino energy to have BDM with
$(T_\chi,m_\chi,\psi)$. See [Fig. 2](../../manual/overview.md#lab_scatt){:target="_blank"} in [BDM Physics <i class="fa-regular fa-bookmark"></i>](../../manual/overview.md).
We have assumed neutrino mass $m_\nu=0$.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `Tx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy $T_\chi$, MeV

> `mx` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM mass $m_\chi$, MeV

> `psi` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Lab frame scattering angle $\psi$, rad



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Attributes:</div>**
> `Ev` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The required neutrino energy $E_\nu$ to boost DM with $m_\chi$ to $T_\chi$, MeV

> `dEv` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The Jacobian $dE_\nu/dT_\chi$, dimensionless

> `x` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;$x:=\cos\psi \in [1,-1]$

> `sanity` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Is the reaction physically plausible? `True` for plausible and `False` for physically impossible.

> `dLips` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Value for differential Lorentz invariant phase space



**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

Import `snorer` and do
```python
>>> import snorer as sn
>>> Tx,mx,psi = 15,1e-3,0.05 # BDM kinetic energy, mx, scattering angle
>>> snv = sn.Neutrino(Tx,mx,psi)
>>> snv.Ev # required Ev
-0.8451953159962898
>>> snv.dEv # Jacobian
0.0031707324661873464
>>> snv.sanity # is this physically possible?
False
```
This example is identical to the example conducted in [**snorer.Kinematics**](Kinematics.md){:target="_blank"} as `snorer.Kinematics` is the superclass of `snorer.Neutrino`. One understands that $T_1=E_\nu$, $T_2=T_\chi$, $m_1=m_\nu=0$ and $m_2=m_\chi$.
