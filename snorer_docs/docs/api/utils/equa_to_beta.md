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


# snorer.equatorial_to_beta


###  <span class="mono">snorer.equatorial_to_beta(*ra*,*dec*,*GC_coord=None*)</span>

Transform equatorial coordinate to off-center angle and galactic coordinate $(\beta,\ell,b)$.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `ra` : *array_like*  <br>&nbsp;&nbsp;&nbsp;&nbsp;Right ascension, hms in string type. Eg. `'5h6.7m4.4s'`.

> `dec` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Declination, dms in string type.  Eg. `'6d10.7m9.4s'`.

> `GC_coord` : *None/list*  <br>&nbsp;&nbsp;&nbsp;&nbsp;The equatorial coordinate for arbitrary galactic center. Default is `None`, which automatically implements our
Milky Way center. For a specific GC coordinate, it should
have `GC_coord = [RA,DEC]` where RA and DEC are, similar to `ra` and `dec`, in hms and dms units respectively. Additionally, they should be subject to ICRS J2000.0.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *tuple* <br>&nbsp;&nbsp;&nbsp;&nbsp;Tuple of $(\beta,\ell,b)$ in rad. Each component is scalar if all inputs are scalars.

**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

We rely on [`astropy.coordinates.SkyCoord`](https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html){:target="_blank"} in [**Astropy**](https://www.astropy.org/){:target="_blank"} to resolve $(\ell,b)$ from $(\alpha,\delta)$ and obtain $\beta$ by
[`snorer.galactic_to_beta`](gal_to_beta.md){:target="_blank"}.
