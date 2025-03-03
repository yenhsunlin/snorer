<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>



# snorer.M_sigma


###  snorer.M_sigma(*mBH*)

Stellar dispersion relation under the influence of black hole. Also known
as $M-\sigma$ relation.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `mBH` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Supermassive black hole mass, $M_\odot$



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Stellar velocity dispersion, km s<sup>−1</sup>. Out is scalar if the input is scalar too.



**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

This function evaluates the stellar dispersion $\sigma_s$ near supermassive black hole
$$
\log_{10}\left(\frac{M_{\rm BH}}{M_\odot}\right)
= 8.29 + 5.12\log_{10}\left(\frac{\sigma_s}{200\,{\rm km\,s^{-1}}}\right).
$$

### References

1.  N. McConnell *et al.*, *Nature* **480**, 215 (2011)
