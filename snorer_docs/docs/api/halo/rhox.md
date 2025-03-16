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

# snorer.rhox


###  <span class="mono">snorer.rhox(*r*,*rhos*,*rs*,*n*)</span>

Dark matter density at $r$.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `r` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance to GC, kpc

> `rhos` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Characteristic density $\rho_s$, MeV cm<sup>−3</sup>

> `rs` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Characteristic radius $r_s$, kpc

> `n` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Slope of the DM profile

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM density at $r$, MeV cm<sup>−3</sup>. Out is scalar if all inputs are scalars


**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

This function evaluates the DM density profile
$$
\rho_\chi(r)=\frac{\rho_s}{\frac{r}{r_s}(1+\frac{r}{r_s})^n}
$$
where rhos and rs are characteristic density and radius respectively.
When $(\rho_s,r_s,n) =$ (184 MeV cm<sup>−3</sup>, 24.42 kpc, 2), it is the famous
NFW profile. If divided by $m_\chi$, it becomes DM number density.

### References
1. G. Bertone *et al.*, *Phys. Rept.* **405**, 279 (2005)