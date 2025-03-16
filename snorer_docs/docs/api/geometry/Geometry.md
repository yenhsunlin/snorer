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



# snorer.Geometry


### *`class`* <span class="mono">snorer.Geometry(*l*,*theta*,*phi*,*Rs*,*Re*,*beta*)</span>

The class constructs the static geomatrical relations for $d$, $r^\prime$ and $\cos\psi$
when $(l,\theta,\phi)$ and $(R_s,R_e,\beta)$ are specified. 
See [Positioning](Positioning.md){:target="_blank"} for more detail.
**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `l` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;The line-of-sight distance, kpc


> `theta` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;The zenith angle $\theta$ at Earth, centers SN, rad


> `phi` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Azimuthal angle $\varphi$ at Earth, centers SN, rad


> `Rs` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from Earth to SN, kpc

> `Re` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from Earth to GC, kpc

> `beta`: *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Off-center angle $\beta$, rad



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Attributes:</div>**
> `d` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from SN to boost point $d$, kpc

> `rprime` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance from GC to boost point $r^\prime$, kpc

> `cos_psi` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;$\cos\psi$ at boost point where $\psi$ is the direction for BDM at B pointing Earth




**<div style="border-bottom: 1px solid lightgray; width: 100%;">Examples</div>**

Import `snorer` and do
```python
>>> l,theta,phi,Rs,Re,beta = 5.160e-9,1e-4,0,11.6,8.5,0.71  # specify quantities
>>> bdmGeo = Geometry(l,theta,phi,Rs,Re,beta)
>>> print(bdmGeo.d)  # SN to boost point
11.59999999483988
>>> print(bdmGeo.rprime) # GC to boost point
8.49999999608676
>>> print(bdmGeo.cos_psi) # cos(psi)
0.9999999721078604
```