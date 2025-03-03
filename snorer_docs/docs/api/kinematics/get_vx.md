<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>



# snorer.get_vx


###  snorer.get_vx(*Tx*, *mx*)

Get dimensionless BDM velocity $v_\chi/c$.
See [User Manual/Physics Overview <i class="fa-regular fa-bookmark"></i>](../../manual/overview.md#from-line-of-sight-to-time-dependency){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `Tx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy $T_\chi$, MeV

> `mx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM mass $m_\chi$, MeV



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dimensionless BDM velocity. `out` is scalar if all inputs are scalars.




