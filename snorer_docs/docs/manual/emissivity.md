<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>

# Emissivity

We systematically derive the BDM emissivity $j_\chi$, given by Eq. (13) in [BDM Physics](overview.md#emissivity-on-the-shell){:target="_blank"},

\begin{equation}
j_\chi = c n_\chi \frac{dn_\nu}{dE_\nu} \frac{d\sigma_{\chi\nu}}{d\Omega} \left( \frac{dE_\nu}{dT_\chi} \frac{v_\chi}{c} \right).
\end{equation} 
Below, we provide two approaches: one that is more familiar to particle physicists, and another that is more heuristic.
All terms can be expressed in vector notations generally, but we omitt this in the following discussion for simplicity.

## Particle approach

### Specific number intensity

First, recall that the particle number $N$ is defined as

\begin{equation}
N\equiv \int d^3 p d^3x\, f(x,p)
\end{equation}

where \( f(x, p) \) is the phase space distribution and \( p  \) the particle momentum.
Hence the differential particle flux is given by

\begin{equation}\label{eq:diff_flux}
d\Phi = \int d^3p f(x, p)v(p) = \int dp d\Omega\, p^2 f(x, p)\, v(p)
\end{equation}
with \( v \) the particle velocity. 
This distribution is defined in momentum space because \( d^3p/E \) is Lorentz invariant, with the dispersion relation \( E = \sqrt{p^2 + m^2} \).
Transforming to other phase space variables is straightforward via an appropriate Jacobian factor. From now on, we omit \(x \) for simplicity.

We refer the integrand in Eq. \eqref{eq:diff_flux} as the *specific number intensity* \( I(p) \),
\begin{equation}\label{eq:I}
I(p) \equiv p^2 f(p)v(p) = \frac{d\\Phi}{dp d\Omega}.
\end{equation} 
This quantity has the dimension of per unit momentum (energy) per unit area per steradian.

### Emissivity on the SN$\nu$ shell

Since the flux equals number density times velocity, we have  

\begin{equation}\label{eq:flux_I}
\frac{d \Phi}{dp d\Omega} = \frac{dn}{dp d\Omega}v(p).
\end{equation}

Alternatively, the flux is also defined as  

\begin{equation}\label{eq:flux_j}
\frac{d\Phi}{dp d\Omega} = \int ds \, j_\chi
\end{equation}
where $j_\chi$ is the emissivity
and we only look for the component such that $v$ is parallel to $s$. Combining the two expressions above, we obtain  

\begin{equation}\label{eq:jx_def}
\frac{dn}{dp d\Omega} v = \int ds\,  j_\chi = v \int dt \, j_\chi
\end{equation} 
where \( ds = v dt \) is used.

Now, we turn to the BDM case and rewrite the relation above as  

\begin{equation}
j_\chi = \frac{dn_\chi}{dp_\chi d\Omega dt} .
\end{equation}

Here, $d\Omega$ is exclusively at the boost point B.
Thus, we obtain

\begin{equation}\label{eq:jx}
 j_\chi = \frac{dn_\chi}{dp_\chi d\Omega dt} = c n_\chi(r) \frac{dn_\nu}{dE_\nu} \frac{d\sigma_{\chi \nu}}{d\Omega} \frac{dE_\nu}{dp_\chi}
\end{equation}
 
and one can further convert \( dE_\nu/dp_\chi \) to \( dE_\nu/dT_\chi \) by 

\[
\frac{dT_\chi}{dp_\chi} = \frac{v_\chi}{c}.
\]

Note that in the earlier step, we considered only the $v$-parallel-to-$\ell$ component, which has been accounted for in the definition of \( d\sigma_{\chi\nu}/d\Omega \).




## Heuristic approach

In this section, we re-derive the same expression for \( j_\chi \) using a heuristic method inspired by astrophysical domain. 

Suppose the number of BDM particles \( dN_\chi \) passing through a surface element \( dA \), in a direction specified by an angle \( \psi \), is proportional to \( dt \), \( dp_\chi \), \( d\Omega \), and the projected surface area element \( \cos\psi\, dA \equiv dA_\perp \). 
See [Fig. 1](#emissivity) below. 
The *proportional factor* is referred to as the specific intensity \( I_\chi \).

<figure id="emissivity">
<center><img src="../../figs/emissivity.svg" alt="scheme" style="width: 45%;">
<figcaption>Figure 1. Overall BDM emitted at boost point B.
</figure>

We can therefore write

\begin{equation}
dN_\chi = I_\chi \, dA_\perp \, d\Omega \, dp_\chi \, dt
\end{equation}

which implies

\begin{equation}
I_\chi = \frac{dN_\chi}{dp_\chi\, d\Omega\, dA_\perp\, dt} = \frac{d\Phi_\chi}{dp_\chi\, d\Omega} \equiv \frac{dn_\chi}{dp_\chi\, d\Omega} v_\chi
\end{equation}

where one can verify that \( dN_\chi / dA_\perp dt \) indeed has the units of flux, and \( n_\chi \) denotes the BDM number density.
The radiative transfer equation implies:

\begin{equation}
\frac{dI_\chi}{ds} = -\alpha_\chi I_\chi + j_\chi
\end{equation}

where \( \alpha_\chi \) characterizes the absorption rate of BDM particles during propagation. However, we assume this effect is highly suppressed and thus negligible in our analysis. 
Therefore, we obtain
\begin{equation}
I_\chi = \int ds\, j_\chi.
\end{equation}
Repeating the same reasoning applied in Eq. \eqref{eq:jx_def} then leads us directly to the desired result in Eq. \eqref{eq:jx}.
