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

# snorer.Mandelstam


### *`class`* <span class="mono">snorer.Mandelstam(*T2*,*m1*,*m2*,*psi*)</span>

Superclass: `snorer.Kinematics`

This class constructs the associated Mandelstam variables $s$, $t$ and $u$ associated with the
scattering process depicted in [Fig. 1](22scat.md/#22scat){:target="_blank"} in [2-2 elastic scattering <i class="fa-regular fa-bookmark"></i>](22scat.md){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `T2` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Kinetic energy $T_2$ received by the particle 2, MeV


> `m1` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Mass of particle 1 (incident) $m_1$, MeV


> `m2` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Mass of particle 2 (target) $m_2$, MeV


> `psi` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Lab frame scattering angle $\psi$, rad



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Attributes:</div>**

> `s` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The $s$-channel in this scattering process, MeV<sup>2</sup>

> `t` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The $t$-channel in this scattering process, MeV<sup>2</sup>

> `u` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The $u$-channel in this scattering process, MeV<sup>2</sup>

> `T1` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The required kinetic energy $T_1$ of particle 1, MeV

> `dT1` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The Jacobian $dT_1/dT_2$, dimensionless

> `x` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;$x:=\cos\psi \in [1,-1]$

> `sanity` : *bool* <br>&nbsp;&nbsp;&nbsp;&nbsp;Are the parameters physically plausible? `True` for plausible and `False` for physically impossible.

> `dLips` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Value for differential Lorentz invariant phase space


**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

Given $T_1$ is obtained from its superclass `snorer.Kinematics`, we can evaluate all Mandelstam variables easily. Thus,

$$
\begin{align*}
s &= (p_1+p_2)^2 = m_1^2+m_2^2 + 2 E_1 m_1 \\
&= m_1^2+m_2^2 + 2(T_1+m_1)m_2
\end{align*}
$$

for $s$-channel, and

$$
\begin{align*}
t &= (p_2^\prime - p_2)^2 = 2m_2^2 - 2E_2 E_2^\prime \\
&= 2m_2^2 - 2(T_2+m_2)m_2.
\end{align*}
$$

For $u$-channel, we use the identity
$$
s+t+u = \sum_i m_i^2 = 2(m_1^2+m_2^2)
$$
where $i$ indicates all particle masses before and after the reaction.

### References
1. M. Peskin and D. Schroeder, *An Introduction To Quantum Field Theory*, Westview (1995)