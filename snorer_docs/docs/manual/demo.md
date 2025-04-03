<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>

# Demonstration

In this page, we show three demonstrations for how BDM evolving in time
on the celestial sphere. Practically, we present the $T_\chi$-integrated
BDM flux by
\begin{equation}
\frac{d\Phi_\chi(t)}{d\Omega}=\int_{T_{\chi,{\rm min}}}^{T_{\chi,{\rm max}}}\frac{d\Phi_\chi(t)}{dT_\chi d\Omega}
\end{equation}
where $(T_{\chi,{\rm min}},T_{\chi,{\rm max}})=(5,50)$, $(300,1000)$ and
$(1,1000)$ keV.

We assume the Milky Way (MW) core-collapse supernova (CCSN) rate approximates 1.63 per century
and the simulation time is 100,000 years.
In average, the expected CCSN number is around 1,630 in this period.
Suppos the CCSN location follows the MW baryonic density distribution
and the birth dates are unifomly distributed in 100,000 years.

The animations are shown below and the time progresses in logarithmic scale.
The dataset contains MW CCSNe locations and birth dates in the following animations
can be [downloaded here](../../dat/MW_sne.txt).

## Animations

### $(T_{\chi,{\rm min}},T_{\chi,{\rm max}})=(5,50)~{\rm keV}$

<video controls width="100%">
  <source src="../../anis/5_50keV.mp4" type="video/mp4">
  Your browser does not suppor MP4 playback.
</video>

[Download video](../../anis/5_50keV.mp4)

### $(T_{\chi,{\rm min}},T_{\chi,{\rm max}})=(300,1000)~{\rm keV}$

<video controls width="100%">
  <source src="../../anis/300_1000keV.mp4" type="video/mp4">
  Your browser does not suppor MP4 playback.
</video>

[Download video](../../anis/300_1000keV.mp4)