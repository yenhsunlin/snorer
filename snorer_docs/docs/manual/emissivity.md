<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>

# Emissivity

We systematically derive how the BDM emissivity  
\begin{equation}
j_\chi = c n_\chi \frac{dn_\nu}{dE_\nu} \frac{d\sigma_{\chi\nu}}{d\Omega} \left( \frac{dE_\nu}{dT_\chi} \frac{v_\chi}{c} \right)
\end{equation} 
was obtained in the previous note. Below, we provide two approaches: one that is more familiar to particle physicists, and another that is more heuristic.

## Specific number intensity

First, recall that the particle flux at position \( \mathbf{x} \) is defined as  
\begin{equation}
\boldsymbol{\Phi}(\mathbf{x}) = \int d^3p\, f(\mathbf{x}, \mathbf{p}) \mathbf{v}(\mathbf{p}) = \int dp d\Omega\, p^2 f(\mathbf{x}, \mathbf{p})\, \mathbf{v}(\mathbf{p})
\end{equation}
where \( f(\mathbf{x}, \mathbf{p}) \) is the phase space distribution, \( p = |\mathbf{p}| \), and \( \mathbf{v} \) is the particle velocity. 
This distribution is defined in momentum space because \( d^3p/E \) is Lorentz invariant, with \( E = \sqrt{p^2 + m^2} \). Transforming to other phase space variables is straightforward via an appropriate Jacobian factor. From now on, we omit \( \mathbf{x} \) for simplicity.

We then define the specific number intensity \( \mathbf{I}(\mathbf{p}) \) as  
\begin{equation}\label{eq:I}
\mathbf{I}(\mathbf{p}) \equiv p^2 f(\mathbf{p})\, \mathbf{v}(\mathbf{p}) = \frac{d\boldsymbol{\Phi}(\mathbf{p})}{dp d\Omega}
\end{equation} 
in the phase space \( dp d\Omega \). This quantity has the dimension: per unit energy, per steradian, per volume.

## Emissivity on the SN$\nu$ shell

Since the flux equals number density times velocity, we have  

\begin{equation}\label{eq:flux_I}
\frac{d\boldsymbol{\Phi}}{dp d\Omega} = \frac{dn}{dp d\Omega}\, \mathbf{v}(\mathbf{p}).
\end{equation}

Alternatively, the flux can be expressed in terms of the emissivity \( j_\chi \) as  

\begin{equation}\label{eq:flux_j}
\frac{d\boldsymbol{\Phi}}{dp d\Omega} = \int d\boldsymbol{\ell} \, j_\chi.
\end{equation}

Combining the two expressions above and projecting only the component along the direction \( \hat{\boldsymbol{\ell}} = \hat{\mathbf{v}} \), we obtain  

\begin{equation}\label{eq:jx_def}
\frac{dn}{dp d\Omega}\, v = \int d\ell \, j_\chi = v \int dt \, j_\chi
\end{equation} 
where we used  \( d\ell = v dt \), with \( v = |\mathbf{v}| \).

Now, recall that the BDM emissivity is restricted to the SN$\nu$ shell, bounded by two Heaviside theta functions. We rewrite the relation above as  

\begin{equation}
\frac{dn_\chi}{dp_\chi d\Omega} \approx \tau_s \int dt\, j_\chi\, \underbrace{ \delta\left( t - \frac{d}{c} - \frac{\ell}{v_\chi} \right) }_{\text{on the shell}} = \tau_s \left. j_\chi \right|_{\text{on the shell}}.
\end{equation}

Here, the subscript *on the shell* means the condition \( t = d/c + \ell/v_\chi \) must be satisfied everywhere. 
Thus, we obtain 

\begin{equation}\label{eq:jx}
\left. j_\chi \right|_{\text{on the shell}} = \frac{dn_\chi}{dp_\chi d\Omega dt} = c n_\chi(r) \frac{dn_\nu}{dE_\nu} \frac{d\sigma_{\chi \nu}}{d\Omega} \frac{dE_\nu}{dp_\chi}.
\end{equation}

Here, the factor \( dt \) is effectively replaced by \( \tau_s \). 
One can further convert \( dE_\nu/dp_\chi \) to \( dE_\nu/dT_\chi \) by the Jacobian  

\[
\frac{dT_\chi}{dp_\chi} = \frac{v_\chi}{c}.
\]

Note that in the earlier step, we considered only the component along \( \hat{\boldsymbol{\ell}} = \hat{\mathbf{v}} \), which has been accounted for in the definition of \( d\sigma_{\chi\nu}/d\Omega \).




## Heuristic approach

In this section, we re-derive the same expression for \( j_\chi \) using a heuristic method inspired by astrophysics. 

Suppose the number of BDM particles \( dN_\chi \) passing through a surface element \( dA \), in a direction specified by an angle \( \psi \), is proportional to \( dt \), \( dp_\chi \), \( d\Omega \), and the projected surface area element \( \cos\psi\, dA \equiv dA_\perp \). 
See [Fig. 1](#emissivity) below. 
The proportionality constant is referred to as the specific intensity \( I_\chi \).

<figure id="emissivity">
<center><img src="../../figs/emissivity.svg" alt="scheme" style="width: 45%;">
<figcaption>Figure 1. Overall BDM emitted at boost point B.
</figure>

We can therefore write:

\begin{equation}
dN_\chi = I_\chi \, dA_\perp \, d\Omega \, dp_\chi \, dt
\end{equation}

which implies:

\begin{equation}
I_\chi = \frac{dN_\chi}{dp_\chi\, d\Omega\, dA_\perp\, dt} = \frac{d\Phi_\chi}{dp_\chi\, d\Omega} \equiv \frac{dn_\chi}{dp_\chi\, d\Omega} v_\chi
\end{equation}

where one can verify that \( dN_\chi / (dA_\perp dt) \) indeed has the units of flux, and \( n_\chi \) denotes the BDM number density.
The radiative transfer equation implies:

\begin{equation}
\frac{dI_\chi}{d\ell} = -\alpha_\chi I_\chi + j_\chi
\end{equation}

where \( \alpha_\chi \) characterizes the absorption rate of BDM particles during propagation. However, we assume this effect is highly suppressed and thus negligible in our analysis. 
Therefore, we obtain:
\begin{equation}
I_\chi = \int d\ell\, j_\chi.
\end{equation}
Repeating the same reasoning applied in Eq. \eqref{eq:jx_def} then leads us directly to the desired result in Eq. \eqref{eq:jx}.
