<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>

# Positioning

To position the BDM signals, cf. [Fig. 1](../../manual/overview.md#snv_bdm_scheme){:target="_blank"} in [User Manual/Physics Overview <i class="fa-regular fa-bookmark"></i>](../../manual/overview.md#general-picture){:target="_blank"}, we have to analyze the geometry of propagation in order to get how large the emissivity is at B.
Once the geometrical relations are understood, we can proceed to construct a class that solve this matter on-the-fly. See also Ref. [[1](#yhl_prd)] and note that we slightly modify the notations in this document for transparency.

## Geometry

We replot [Fig. 1](../../manual/overview.md#snv_bdm_scheme){:target="_blank"} in [User Manual/Physics Overview <i class="fa-regular fa-bookmark"></i>](../../manual/overview.md#general-picture){:target="_blank"} and leaves only relevant information in the following [Fig. 1](#geometry1). Three coplanar points S, G and E remained the same. Repeated quantities are the same while we also mark addditional auxiliary terms which their purposes will be clear soon.

<figure id="geometry1">
<center><img src="../../../figs/geometry1.svg" alt="geometry1" style="width: 60%;">
<figcaption>Figure 1. The geometry of BDM propagation.
</figure>

By examining the blue and green triangles, we have

$$
\begin{align}
h &= \ell \sin\theta \\
b &= \ell \cos\theta.
\end{align}
$$

Moreover, blue and brown triangles are identical, as brown one simply rotates around axis $\overline{\mathsf{SE}}$ by $\varphi$. See  [Fig. 2](#geometry2).

<figure id="geometry2">
<center><img src="../../../figs/geometry2.svg" alt="geometry1" style="width: 60%;">
<figcaption>Figure 2. Blue and brown triangles are identical.
</figure>

Suppose B is distant $r$ away from G, then B′ is $r^\prime$. It says that $r$ is a special case of $r^\prime$ when $\varphi=0$.
This is curcial that when SN is not at the GC, DM number density $n_\chi$ is not spherical symmetric relative to SN.
To resolve the correct $j_\chi$ at boost point, one needs to know $n_\chi(r^\prime(\varphi))$.



<figure id="geometry3">
<center><img src="../../../figs/geometry3.svg" alt="geometry1" style="width: 90%;">
<figcaption>Figure 3. Another set of auxiliary triangles.
</figure>

We draw another set of auxiliary triangles in [Fig. 3](#geometry3), also cf. [Fig. 1](#geometry1), and immediately see that, left figure,

$$
\begin{equation}\label{eq:rprime}
r^{\prime2} = a^2 + h^2\cos^2\varphi
\end{equation}
$$

while 

$$
\begin{equation}
a^2 = (\rho\sin\delta - h\sin\varphi)^2 + \rho^2\cos^2\delta
\end{equation}
$$

and, from law of cosine, the right figure gives,


$$
\begin{equation}
\rho = \sqrt{b^2+R_e^2 - 2b R_e \cos\beta}.
\end{equation}
$$

To determine $\delta$, we apply law of cosine (cause we need to know whether $\delta> \pi/2$ or not) on right figure,

$$
R_e^2 =\rho^2 + b^2  - 2\rho b \cos(\pi-\delta)
$$

which yields


$$
\begin{equation}
 \cos\delta =\frac{R_e^2-\rho^2-b^2}{2\rho b}
\end{equation}
$$

One can check that we already determine $r^\prime$, Eq. $\eqref{eq:rprime}$, by known quantities $(d,\ell,\theta,\varphi)$ and $(R_s,R_e,\beta)$.
The first set will be specified during evaluating BDM signatures and the second set of parameters gives the SN location.

Last thing to be calculated is the scattering angle $\psi$,
which can be easily obtained by law of cosine again,


$$
R_s^2 = d^2 + \ell^2 -2d \ell \cos(\pi-\psi) 
$$

such that

$$
\begin{equation}
\cos\psi = \frac{R_s^2-d^2-\ell^2}{2d\ell}
\end{equation}
$$

where
$$
\begin{equation}\label{eq:d}
d = \sqrt{\ell^2 + R_s^2 - 2\ell R_s\cos\theta}.
\end{equation}
$$

This indicates that $d$ is not an independent quantity but subject to $\ell$ and $\theta$.

## Static to time-dependent

From  Eq. (15) in [User Manual/Physics Overview <i class="fa-regular fa-bookmark"></i>](../../manual/overview.md#from-line-of-sight-to-time-dependency){:target="_blank"} and offseting it by $t_\nu=R_s/c$, we have

$$
t = \frac{d}{c} + \frac{\ell}{v_\chi}- t_\nu
$$

hence

$$
\begin{equation}\label{eq:t_dependent}
d + \frac{\ell}{\beta_\chi}  = R_s + ct.
\end{equation}
$$

where $\beta_\chi = v_\chi/c$. For convenience, we define

$$
\begin{equation}
\zeta(t) := R_s + ct
\end{equation}
$$

and plug Eq. $\eqref{eq:d}$ into Eq. $\eqref{eq:t_dependent}$ and solve for $\ell$,
$$
\begin{equation}\label{eq:ell_t}
\ell(t) = -\frac{\beta_\chi}{1-\beta^2_\chi}\left(\alpha
+ \gamma-\zeta \right).
\end{equation}
$$

where

$$
\begin{align*}
\alpha &= \sqrt{(R_s^2-\zeta^2)(1-\beta_\chi^2)+(R_s\beta_\chi\cos\theta-\zeta)^2},\\
\gamma &= R_s\beta_\chi\cos\theta.
\end{align*}
$$

We now can determine $\ell$ at any time $t$ from Eq. $\eqref{eq:ell_t}$ and $d$ is also subject to the change of $\ell$ when $\theta$ is specified..
The geometry for BDM becomes time-dependent from a static profile.

### References
1. <p id="yhl_prd">Y.-H. Lin *et al.*, *Phys. Rev. D.* **108**, 083013 (2023)</p> 