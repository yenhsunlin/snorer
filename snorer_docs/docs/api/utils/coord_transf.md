<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>

# Coordinate Transformations

We briefly introduce how to extract the off-center angle $\beta$, cf. [Fig. 1](../../manual/overview.md#snv_bdm_scheme){:target="_blank"} in [BDM Physics](../../manual/overview.md#general-picture){:target="_blank"}, from two famous astronomical coordinate systems, galactic and equatorial coordinates.
This aids us to input arbitrary supernova (SN) loactions from astrophysical database, eg. [SRcat](http://snrcat.physics.umanitoba.ca/){:target="_blank"},  into **snorer** and evaluates its BDM signature.

## Galactic coordinate

The galactic coordinate is shown in [Fig. 1](#galcoord) with its longitute $\ell$ and latitude $b$.
Once $(\ell,b)$ for a SN is specified, its position on the celestial sphere is marked.

<figure id="galcoord">
<center><img src="../../../figs/galactic_coord.svg" alt="galcoord" style="width: 80%;">
<figcaption>Figure 1. The galactic coordinate, including longitutde \(\ell\) and latitude \(b\).
</figure>

Given that the distance from Earth to SN, $R_s$, will be provided by the database, the last task is to retrieve $\beta$.
Before proceed, there are two conventions in this coordinate. In terms of Cartesian representation, the Earth is at the origin and the galactic center (GC) at $(x,y,z)=(0,R_e,0)$.
See the lateral view to this geometrical system in [Fig. 2](#gal_to_bR). 

<figure id="gal_to_bR">
<center><img src="../../../figs/gal_to_betaRs.svg" alt="gal_to_bR" style="width: 60%;">
<figcaption>Figure 2. The Cartesian representation of galactic coordinate with distance included.
</figure>

With two unit vectors that one points from Earth to SN $\hat{\mathbf{s}} = (x_s,y_s,z_s)/R_s$ and the other points from Earth to GC $\hat{\mathbf{g}} = (0,1,0)$, we immediately obtain

$$
\begin{equation}
y_s = R_s \cos b \cos\ell
\end{equation}
$$

The reason why we don't need $x_s$ and $z_s$ is clear that when we evaluate $\beta$,

$$
\begin{equation}\label{eq:cos_beta_MW}
\cos\beta = \hat{\mathbf{g}}\cdot \hat{\mathbf{s}} = \cos b\cos\ell
\end{equation}
$$

there is no need for these two.

Note that when we retrieve $\beta$ from galactic coordinate, there is impossible to transform it back.
It is obvious that $\beta$ to $(\ell,b)$ is a one to two transformation. Without one additional constraint,
the reverse is ill-posed. Reconstructing the SN location relative to galactic plane
is implausible simply by $\beta$.

For completeness, to have the exact location, the missing information one needs is the azimuthal angle, say, around the 
Earth-GC axis. This requires both $x_s$ and $z_s$. We omit the detail here due to its irrelevance.

### Supernova in arbitrary distant galaxy

We now consider a more general case that SN lies in a arbitrary distant galaxy and we want to know the BDM signature from that SN.
Again, no matter how this scene chages, the underlying concept is to obtain three things: $R_g$, $R_s$ and $\beta$. The last two
are quite familiar but the first one indicates the distance between Earth and the center of the distant galaxy.
On the other hand, we actually replace $R_e$ by $R_g$. See [Fig. 3](#arb_gal).

<figure id="arb_gal">
<center><img src="../../../figs/arbitrary_galx.svg" alt="arb_gal" style="width: 60%;">
<figcaption>Figure 3. SN in arbitrary distant galaxy.
</figure>

Suppose the galactic coordinates for the two stellar objects are $(\ell_s,b_s)$ and $(\ell_g,b_g)$ respectively, we can follow
the previous derivation to obtain their Cartesian representations $(x_s,y_s,z_s)$ and $(x_g,y_g,z_g)$.
Thus

$$
\begin{align*}
x_s &=  R_s \cos b_s \sin (2\pi- \ell_s), \\
y_s & = R_s \cos b_s \cos \ell_s, \\
z_s & = R_s \sin b_s
\end{align*}
$$

for SN and

$$
\begin{align*}
x_g &=  R_g \cos b_g \sin (2\pi- \ell_g), \\
y_g & = R_g \cos b_g \cos \ell_g, \\
z_g & = R_g \sin b_g
\end{align*}
$$

for the distant galaxy.

We define two unit vectors for the two stellar objects
$$
\hat{\mathbf{s}}=(-\cos b_s \sin \ell_s, \cos b_s \cos \ell_s,\sin b_s)
$$
and
$$
\hat{\mathbf{g}}=(-\cos b_g \sin \ell_g, \cos b_g \cos \ell_g,\sin b_g)
$$
respectively. The off-center angle $\beta$ can be retrieved by the same formula

$$
\begin{equation}
\cos \beta = \hat{\mathbf{g}}\cdot\hat{\mathbf{s}} = \cos b_s \cos b_g \cos(\ell_s - \ell_g) + \sin b_s \sin b_g.
\end{equation}
$$

In terms of our Milky Way, $\ell_g = b_g=0$, we recover Eq. $\eqref{eq:cos_beta_MW}$.

## Equatorial coordinate

Another commonly used coordinate system is the equatorial coordinate and is specified by right ascension $\alpha$ (RA) and declination $\delta$ (DEC).
Due to its complexity, we do not provide any detail mathematical conversion from $(\alpha,\delta)$ to $\beta$ here.

To tackle this task, we employ [**Astropy**](https://www.astropy.org/){:target="_blank"} to resolve $(\ell,b)$ from $(\alpha,\delta)$.
Then proceed the discussion in the last section to recover $\beta$. Note that we will assume the input $(\alpha,\delta)$ in **snorer** is expressed in terms of ICRS J2000.0.