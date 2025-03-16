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

# snorer.radiusInfluence


###  <span class="mono">snorer.radiusInfluence(*mBH*)</span>

Influence radius of a supermassive black hole


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `mBH` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Supermassive black hole mass, $M_\odot$



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Influence radisu, kpc. Out is scalar if the input is scalar.


**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

This function evaluates the influence radius of a supermassive black hole
$$
r_h = \frac{GM_{\rm BH}}{\sigma_s^2}
$$
where sigma_s is the stellar dispersion near SMBH. See also [`snorer.M_sigma`](M_sigma.md){:target="_blank"}.
