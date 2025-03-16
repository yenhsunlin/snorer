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


# snorer.dsigma_xv


###  <span class="mono">snorer.dsigma_xv(*Ev*,*mx*,*psi*,*sigxv0='1e-45'*)</span>

Differential DM-$\nu$ scattering cross section at angle $\psi$ in lab frame,
$$
\frac{d\sigma_{\chi\nu}}{d\Omega_{\rm lab}}=\sigma_0 \times g_\chi(\psi).
$$
See Eqs. (2) and (3) in [BDM Physics](../../manual/overview.md#particle-kinematics){:target="_blank"}.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Parameters:</div>**

> `Ev` : *array_like*  <br>&nbsp;&nbsp;&nbsp;&nbsp;Neutrino energy, MeV.</pre>

> `mx` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Dark matter mass, MeV.

> `psi` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Lab frame scattering angle $\psi \in [0,\pi/2]$.

> `sigxv0` : *array_like* <br>&nbsp;&nbsp;&nbsp;&nbsp;Energy-independent DM-$\nu$ cross section $\sigma_0$, cm<sup>2</sup>. Default is $10^{-45}$ cm<sup>2</sup>.

**<div style="background-color: lightgrey; padding: 5px; width: 100%;">Returns:</div>**

> `out` : *scalar/ndarray* <br>&nbsp;&nbsp;&nbsp;&nbsp;Differential DM-$\nu$ cross section, cm<sup>2</sup> sr<sup>−1</sup>. Out is scalar if all inputs are scalars.


**<div style="border-bottom: 1px solid lightgray; width: 100%;">Notes</div>**

The result is simply `sigxv0 * snorer.get_gx`. See also [`snorer.get_gx`](../kinematics/get_gx.md){:target="_blank"}.
