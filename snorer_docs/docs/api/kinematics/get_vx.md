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


# snorer.get_vx


###  <span class="mono">ssnorer.get_vx(*Tx*,*mx*)</span>

Get dimensionless BDM velocity $v_\chi/c$.
See [BDM Physics](../../manual/overview.md#from-line-of-sight-to-time-dependency){:target="_blank"} for detail.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `Tx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy $T_\chi$, MeV

> `mx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM mass $m_\chi$, MeV



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dimensionless BDM velocity. `out` is scalar if all inputs are scalars.




