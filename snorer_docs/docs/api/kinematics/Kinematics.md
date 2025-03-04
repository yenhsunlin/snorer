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

# snorer.Kinematics


### *`class`* <span class="mono">snorer.Kinematics(*T2*,*m1*,*m2*,*psi*)</span>

This class constructs the required kinetic energy $T_1$ of incoming particle with
mass $m_1$ to boost the target with mass $m_2$ to kinetic energy $T_2$ along the direction
$\psi$. See [Fig. 1](22scat.md/#22scat){:target="_blank"} in  [2-2 elastic scattering](22scat.md){:target="_blank"}.
**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `T2` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Kinetic energy $T_2$ received by the particle 2, MeV


> `m1` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Mass of particle 1 (incident) $m_1$, MeV


> `m2` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Mass of particle 2 (target) $m_2$, MeV


> `psi` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Lab frame scattering angle $\psi$, rad



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Attributes:</div>**
> `T1` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The required kinetic energy $T_1$ of particle 1, MeV

> `dT1` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The Jacobian $dT_1/dT_2$, dimensionless

> `x` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;$x:=\cos\psi \in [1,-1]$

> `sanity` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Are the parameters physically plausible? `True` for plausible and `False` for physically impossible.

> `dLips` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Value for differential Lorentz invariant phase space



**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

Import `snorer` and do
```python
>>> import snorer as sn
>>> T2,m1,m2,psi = 15,0,1e-3,0.05 # kinetic energy, m1, m2, scattering angle
>>> snv = sn.Kinematics(T2,m1,m2,psi)
>>> snv.T1 # required kinetic energy T1 for particle 1
-0.8451953159962898
>>> snv.dT1 # Jacobian
0.0031707324661873464
>>> snv.sanity # is this physically possible?
False
```
It is clear that massless particle 1 is no way to upscatter particle 2 with $m_2=10^{-3}$ MeV to $T_2=15$ MeV at angular direction $\psi$. Becasue the required $T_1$ (`snv.T1`) is negative. 
The attribute `snv.sanity` is `False`, which implies this reaction is physically impossible.
