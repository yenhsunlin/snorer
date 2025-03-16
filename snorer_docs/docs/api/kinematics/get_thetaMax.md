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


# snorer.get_thetaMax


###  <span class="mono">snorer.get_thetaMax(*t*,*Tx*,*mx*,*Rs*)</span>

Find the maximum BDM field-of-view, $\theta_{M}^*$, that centers SN at particular time $t^*$.
See Eq. (24) in [BDM Physics](../../manual/overview.md#field-of-view-across-the-sky){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `t` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;The BDM at particular time $t^*$, seconds. If $t^* > t_{\rm van}$, the result is unphysical.

> `Tx` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;BDM kinetic energy $T_\chi$, MeV


> `mx` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;DM mass $m_\chi$, MeV

> `Rs` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;Distance to supernova, $R_s$, kpc.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar* <br>&nbsp;&nbsp;&nbsp;&nbsp;Maximum field-of-view centers supernova, $\theta^*_M$ [rad].




