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


# snorer.get_psiMax


###  <span class="mono">ssnorer.get_psiMax(*Tx*,*mx*)</span>

Get the maximumly allowed scattering angle $\psi_{\rm max}$.
Se Eq. (6) [BDM Physics](../../manual/overview.md#constraint-by-positive-definite-e_nu){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `Tx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy $T_\chi$, MeV

> `mx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM mass $m_\chi$, MeV



**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Maximum allowed scattering angle $\psi_{\rm max}$ [rad]. `out` is scalar if all inputs are scalars.




