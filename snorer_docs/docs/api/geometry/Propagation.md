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

# snorer.Propagation


### *`class`* <span class="mono">snorer.Propagation(*t*,*vx*,*theta*,*phi*,*Rs*,*Re*,*beta*)</span>

Superclass: `snorer.Geometry`

The class constructs the dynamical geomatrical relations for $d$, $r^\prime$ and $\cos\psi$
when $(l,\theta,\phi)$ and $(R_s,R_e,\beta)$ are specified.

Unlike its superclass Geometry, the class parameter `l` is now replaced by a specific
time `t` and dimensionless BDM velocity `vx`. This allows it to incorporate time-dependent feature when evaluating the geometrical quantities during propagation.  

This class is also not exclusively for SN in MW or LMC, it can be generalized to SN
in arbitrary distant galaxy as long as the aforementioned inputs are determined.
The BDM emissivity along the line-of-sight then can be determined when calculate
the BDM flux and event at Earth associated to that particular SN.

See [Positioning](Positioning.md){:target="_blank"} for more detail.
**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `t` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The BDM at specific time $t$, seconds. Time-zero is set to be the arrival of SN$\nu$ at Earth

> `vx` : *float*  <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM dimesionless velocity $v_\chi/c$

> `theta` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The zenith angle $\theta$ at Earth, centers SN, rad


> `phi` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Azimuthal angle $\varphi$ at Earth, centers SN, rad


> `Rs` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from Earth to SN $R_s$, kpc

> `Re` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from Earth to GC $R_e$, kpc

> `beta`: *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Off-center angle $\beta$, rad



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Attributes:</div>**

> `l` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;The line-of-sight distance $\ell$, kpc

> `d` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from SN to boost point $d$, kpc

> `rprime` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from GC to boost point $r^\prime$, kpc

> `cos_psi` : *float* <br>&nbsp;&nbsp;&nbsp;&nbsp;$\cos\psi$ at boost point where $\psi$ is the direction for BDM at B pointing Earth




**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

Import `snorer` and do
```python
>>> bdmProp = Propagation(t=59,vx=0.9,theta=1e-4,phi=0,Rs=11.6,Re=8.5,beta=0.71)
>>> print(bdmProp.l)  # The corresponding l.o.s. distance
5.160120751743069e-09
>>> print(bdmProp.d)  # The distace from SN to boost point
11.59999999483988
>>> print(bdmProp.rprime)  # The distance from GC to boost point
8.49999999608676
>>> print(bdmProp.cos_psi)  # Scattering angle that points Earth at B
1.0
```