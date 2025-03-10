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


# snorer.galactic_to_beta


###  <span class="mono">snorer.galactic_to_beta(*l*,*b*,*GC_coord=[0,0]*)</span>

Transform galactic coordinate $(\ell,b)$ to off-center angle $\beta$.
See Eqs. (2) in [Coordinate Transformations](coord_transf.md){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `l` : *array_like*  <br>&nbsp;&nbsp;&nbsp;&nbsp;Galactic longitude, rad.

> `b` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Galactic latitude, rad.

> `GC_coord` : *list* <br>&nbsp;&nbsp;&nbsp;&nbsp;Galactic coordinate for arbitrary galactic center $(\ell_g,b_g)$. Default is Milky Way center `GC_coord = [0,0]`.


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Off-center angle $\beta$, rad. The result is scalar if all inputs are scalars.
