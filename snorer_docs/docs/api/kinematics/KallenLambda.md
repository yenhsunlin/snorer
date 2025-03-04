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

# snorer.KallenLambda


###  <span class="mono">snorer.KallenLambda(*x*,*y*,*z*)</span>

K&auml;llen lambda function
$$
\lambda(x,y,z)=x^2+y^2+z^2-2(xy+yz+zx),
$$
a useful function for evaluating kinetical quantities in particle physics.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `x` : *array_like*

> `y` : *array_like*

> `z` : *array_like*


**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;`out` is scalar if all inputs are scalars.



### References
1. <p id="bib_ConceptQFT">V. Ilisie, *Concepts in quantum field theory: A practitioner's toolkit*, Springer (2016)</p> 


