<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>

# Demonstration


In this page, we present three demonstrations showing how BDM evolves over time on the celestial sphere. 
Specifically, we display the **\( T_\chi \)-integrated BDM flux** for \( m_\chi = 5 \, \text{keV} \), defined by 
\begin{equation}
\frac{d\Phi_\chi(t)}{d\Omega} = \int_{T_{\chi,\mathrm{min}}}^{T_{\chi,\mathrm{max}}} \frac{d\Phi_\chi(t)}{dT_\chi\, d\Omega}
\end{equation} 
with integration ranges \( (T_{\chi,\mathrm{min}}, T_{\chi,\mathrm{max}}) = (5, 50) \), \( (300, 1000) \), and \( (1, 1000) \) keV.

We assume the Milky Way (MW) core-collapse supernova (CCSN) rate is approximately 1.63 per century, 
and the simulation spans a total duration of 100,000 years. 
On average, this corresponds to about 1,630 CCSNe occurring during the simulation period. 
The CCSN locations are sampled from the MW baryonic density distribution, 
and their birth times are assumed to be uniformly distributed over the 100,000-year period.

The animations below illustrate the time evolution in **logarithmic time scale**. 
The dataset containing the simulated CCSN positions and birth dates used in these animations 
can be [downloaded here](../../dat/MW_sne.txt){:target="_blank"}.


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