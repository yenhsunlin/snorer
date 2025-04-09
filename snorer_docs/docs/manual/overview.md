<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>

# BDM Physics

We review all the details provided in Refs. [[1](#bib_Lin2023PRL), [2](#bib_Lin2023PRD)]. 
While these papers present mathematical expressions suitable for general readers, the content in this section is tailored for programming purposes. 
The goal of this document is not to replace Refs. [[1](#bib_Lin2023PRL), [2](#bib_Lin2023PRD)] but to bridge the gap between traditional journal writing and the technical details required for programming. 
Typos and misstatements are corrected here to ensure clarity and to facilitate understanding of the subsequent Python code.



## General picture

Given a non-zero cross section between DM ($\chi$) and neutrino ($\nu$), it is expected that SN$\nu$ can scatter with halo $\chi$ particles as they propagate outward from the explosion site. 
We illustrate the general setup in [Fig. 1](#snv_bdm_scheme). In this figure, the supernova (SN), galactic center (GC), and Earth are labeled as S, G, and E respectively.

<figure id="snv_bdm_scheme">
<center><img src="../../figs/scheme.svg" alt="scheme" style="width: 60%;">
<figcaption>Figure 1. The 3D schematic of SN\(\nu\)-induced BDM.</figcaption>
</figure>

In [Fig. 1](#snv_bdm_scheme), S, G, and E lie on the same plane, and a $\chi$ particle is boosted at point B on the SN$\nu$ shell. 
Here, $d$ denotes the distance $\overline{\mathsf{SB}}$, $\ell$ the distance $\overline{\mathsf{BE}}$, $R_s$ the distance between the SN and Earth, and $R_e$ the distance between the GC and Earth. 
The DM number density at point B is given by $n_\chi(r)$, where $r$ is the radial distance from the GC to B. 
We adopt the NFW profile [[3](#bib_NFW)] for $n_\chi$:

\begin{equation}
n_\chi(r) = \frac{\rho_s}{m_\chi} \frac{1}{\frac{r}{r_s} ( 1 + \frac{r}{r_s} )^2}
\end{equation}

where $\rho_s = 184$ MeV cm<sup>−3</sup> and $r_s = 24.4$ kpc.
Due to the spherically symmetric nature of the DM halo, there is no need to specify the orientation of the galactic plane. 
However, unless S lies precisely along the axis $\overline{\mathsf{GE}}$, the resulting BDM distribution will not be azimuthally symmetric in $\varphi$.

We define $\beta$ as the off-center angle, which quantifies the angular displacement of the SN relative to the GC. 
This quantity is further explained in the API documentation:  
[Positioning](../api/geometry/Positioning.md){:target="_blank"} and 
[Coordinate Transformations](../api/utils/coord_transf.md){:target="_blank"}.




### Particle kinematics 

When a $\nu$ carries energy $E_\nu$ and scatters with a static $\chi$ at point B, the DM particle receives kinetic energy $T_\chi$ given by
$$
T_\chi = \frac{E_\nu^2}{E_\nu + m_\chi / 2} \left( \frac{1 + \cos\theta_c}{2} \right)
$$
where $\theta_c$ is the scattering angle in the center-of-mass (cm) frame (not to be confused with the $\theta$ in [Fig. 1](#snv_bdm_scheme)). 
Assuming $\theta_c$ is uniformly distributed in $[0, \pi]$, it can be related to the lab frame scattering angle $\psi$ by
$$
\theta_c = 2 \tan^{-1}(\gamma \tan\psi)
$$
where $\psi \in [0, \pi/2]$ and the boost factor is defined as
$$
\gamma = \frac{E_\nu + m_\chi}{\sqrt{m_\chi(2E_\nu + m_\chi)}}.
$$
Thus, the differential cross section in the lab frame becomes 
\begin{equation}
\frac{d\sigma_{\chi\nu}}{d\Omega_{\rm lab}} = \frac{1}{2\pi} \frac{d\sigma_{\chi\nu}}{d\cos\psi} = \sigma_0 \times g_\chi(\psi)
\end{equation} 
where 
\begin{equation}\label{eq:gx}
g_\chi(\psi) = \frac{\gamma^2 \sec^3\psi}{\pi (1 + \gamma^2 \tan^2\psi)^2}
\end{equation} 
is the angular distribution of the $\chi\nu$ scattering cross section in the lab frame. 
One can verify that Eq. \eqref{eq:gx} satisfies the normalization condition
$$
2\pi \int_0^{\pi/2} g_\chi(\psi)\, \sin\psi\, d\psi = 1
$$
and that the result is independent of $\gamma$.
When $d$, $\ell$, and $R_s$ are specified, the scattering angle $\psi$ can be determined using the law of cosines:
$$
R_s^2 = d^2 + \ell^2 - 2d\ell \cos(\pi - \psi)
$$
which gives 
\begin{equation}\label{eq:geo_psi}
\psi = \cos^{-1} \left( \frac{R_s^2 - d^2 - \ell^2}{2d\ell} \right).
\end{equation}




### Constraint by positive-definite $E_\nu$
Note that although the valid range of $\psi$ from Eq. \eqref{eq:geo_psi} lies in $[0, \pi]$, Eq. \eqref{eq:gx} imposes a stronger constraint. 
While Eq. \eqref{eq:geo_psi} is purely a geometric relation, Eq. \eqref{eq:gx} arises from the kinematics of the scattering process.
In fact, further analysis reveals that the true valid range for $\psi$ is even narrower. 
Let us consider the $\chi\nu$ scattering in the lab frame, illustrated in [Fig. 2](#lab_scatt).

<figure id="lab_scatt">
<center><img src="../../figs/lab_scattering.svg" alt="scheme" style="width: 40%;">
<figcaption>Figure 2. \(\chi\nu\) scattering in lab frame.</figcaption>
</center>
</figure>

We write the four-momenta for each particle as

$$
\begin{align*}
p_{\nu} &= (E_{\nu}, 0, 0, E_{\nu}), \\
p_{\chi} &= (m_{\chi}, 0, 0, 0), \\
p_{\nu}^{\prime} &= (E_{\nu}^{\prime}, x, y, z), \\
p_{\chi}^{\prime} &= \left(E_{\chi}, -|\mathbf{p}_{\chi}|\sin\psi, 0, |\mathbf{p}_{\chi}|\cos\psi \right),
\end{align*}
$$

where $(x, y, z)$ are irrelevant to the calculation. 
Assuming $m_\nu = 0$ and using the metric signature $g^{\mu\nu}=g_{\mu\nu} = \mathrm{diag}(1, -1, -1, -1)$, relativistic kinematics gives 

$$
\begin{align*}
(p_\nu - p_\chi^\prime)^2 &= (p_\chi - p_\nu^\prime)^2, \\
m_\chi^2 - 2E_\nu (E_\chi - |\mathbf{p}_{\chi}|\cos\psi) &= m_\chi^2 - 2m_\chi E_\nu^\prime.
\end{align*}
$$

Using the relations $E_\chi = T_\chi + m_\chi$ and $T_\chi = E_\nu - E_\nu^\prime$, the equation becomes

\begin{equation}\label{eq:Ev}
E_\nu = \frac{T_\chi m_\chi}{|\mathbf{p}_{\chi}|\cos\psi - T_\chi}.
\end{equation}

Given the BDM kinetic energy $T_\chi$ and scattering angle $\psi$ at point B, Eq. \eqref{eq:Ev} determines the required $E_\nu$ to produce such a configuration. 
Since $E_\nu$ must be positive and finite, this imposes a stronger constraint:

\begin{equation}\label{eq:psi_max}
|\mathbf{p}_\chi| \cos\psi - T_\chi > 0 \quad \Rightarrow \quad \psi < \psi_{\rm max} := \cos^{-1} \left( \frac{T_\chi}{|\mathbf{p}_{\chi}|} \right).
\end{equation}

It is also straightforward to derive

$$
|\mathbf{p}_\chi| = \sqrt{T_\chi (T_\chi + 2m_\chi)}.
$$

Hence, the physically allowed range for $\psi$ is

\begin{equation}
\psi \in [0, \psi_{\rm max})
\end{equation}

and $E_\nu$ diverges as $\psi \to \psi_{\rm max}$.
Moreover, by differentiating Eq. \eqref{eq:Ev} with respect to $T_\chi$, we obtain

\begin{equation}\label{eq:dEv/dTx} 
\frac{dE_\nu}{dT_\chi} = \left( \frac{m_\chi}{|\mathbf{p}_{\chi}|\cos\psi - T_\chi} \right)^2 \frac{T_\chi}{|\mathbf{p}_{\chi}|} \cos\psi.
\end{equation}

This expression will be useful in later derivations.


## Dark emissivity

The key factor is determining how many $\chi$ particles are boosted at point B. 
This is generally characterized by the **emissivity**, which has units of cm<sup>−3</sup> s<sup>−1</sup>, or MeV<sup>−1</sup> cm<sup>−3</sup> s<sup>−1</sup> when expressed in terms of the energy spectrum.


### SN$\nu$ spectrum

We begin by writing down the SN$\nu$ energy spectrum on the shell at distance $d$:
\begin{equation}\label{eq:snv_spectrum}
\frac{dN_\nu}{dE_\nu} = \sum_i \frac{L_{\nu_i}}{4\pi d^2 \langle E_{\nu_i} \rangle} E_\nu^2 f_{\nu_i}(E_\nu)
\end{equation}
where $L_{\nu_i} = L_{\rm tot} / 6$ is the luminosity for each neutrino species ($\nu_{e,\mu,\tau}$ and their antiparticles), 
and $f_{\nu_i}(E_\nu)$ is the neutrino phase space distributin that obeys Fermi-Dirac statistics [[4](#bib_Duan2006PRD)]:
\begin{equation}
E_\nu^2 f_{\nu_i}(E_\nu) = \frac{1}{F_2(\eta_\nu)} \frac{1}{T_\nu^3} \frac{E_\nu^2}{e^{E_\nu/T_\nu - \eta_\nu} + 1}
\end{equation}
where
\begin{equation}
F_k(\eta) = \int_0^\infty dx\, \frac{x^k}{e^{x - \eta} + 1}.
\end{equation}
We list all the numerical values for the parameters mentioned above in [Tab. 1](#numerical_values).

<div style="display: flex; justify-content: center;">
    <table id ="numerical_values" style="border-collapse: collapse; text-align: center;">
        <tr>
            <th>Parameters</th>
            <th>Values</th>
            <th>Parameters</th>
            <th>Values</th>
        </tr>
        <tr>
            <td>\(\langle E_{\nu_e}\rangle\)</td>
            <td>11 MeV</td>
            <td>\(T_{\nu_e}\)</td>
            <td>2.76 MeV</td>
        </tr>
        <tr>
            <td>\(\langle E_{\bar{\nu}_e}\rangle\)</td>
            <td>16 MeV</td>
            <td>\(T_{\bar{\nu}_e}\)</td>
            <td>4.01 MeV</td>
        </tr>
        <tr>
            <td>\(\langle E_{\nu_x,\bar{\nu}_x}\rangle\)</td>
            <td>25 MeV</td>
            <td>\(T_{\nu_x,\bar{\nu}_x}\)</td>
            <td>6.26 MeV</td>
        </tr>
        <tr>
            <td>\(L_{\rm tot}\)</td>
            <td>\(3\times10^{52}\) erg s<sup>−1</sup></td>
            <td>\(\eta_i\)</td>
            <td>3</td>
        </tr>
    <caption style="font-style: normal;">Table 1. Numerical values for parameters.</caption>
    </table>
</div>

Suppose the duration of the SN explosion is approximately $\tau_s = 10$ s. 
The total energy released in the form of neutrinos from a single explosion is roughly $10^{53}$ erg.

### Number density on the shell

After converting erg to MeV, one can verify that Eq. \eqref{eq:snv_spectrum}, when multiplied by $\tau_s$, has units of MeV<sup>−1</sup> cm<sup>−2</sup>. 
Since neutrinos travel at the speed of light, the thickness of the SN$\nu$ shell can be estimated as $h = c \tau_s$.
The number density of neutrinos on the shell is then given by
\begin{equation}\label{eq:snv_nd}
\frac{dn_\nu}{dE_\nu} = \frac{dN_\nu}{dE_\nu} \frac{\tau_s}{h} = \sum_i \frac{L_{\nu_i}}{4\pi d^2 \langle E_{\nu_i} \rangle c} E_\nu^2 f_{\nu_i}(E_\nu).
\end{equation}
One can confirm that $dn_\nu/dE_\nu$ has units of MeV<sup>−1</sup> cm<sup>−3</sup>, 
which has the correct dimension.
<!-- By multiplying by $E_\nu$, the MeV<sup>−1</sup> unit is removed, yielding
$$
n_\nu = \frac{E_\nu}{c} \frac{dN_\nu}{dE_\nu},
$$
which gives the SN$\nu$ number density. -->


### Emissivity on the shell

To obtain the emissivity $j_\chi$ at point B on the SN$\nu$ shell, one simply computes
$$
j_\chi = c\, n_\nu\, n_\chi\, \frac{d\sigma_{\chi\nu}}{d\Omega_{\rm lab}},
$$
which has units of cm<sup>−3</sup> s<sup>−1</sup> sr<sup>−1</sup>, where sr<sup>−1</sup> indicates per steradian.
It is often more convenient to restore the energy spectrum form and express it in terms of $T_\chi$, yielding
\begin{equation}\label{eq:jx}
j_\chi(d, r, T_\chi, \psi) = c n_\chi \frac{dn_\nu}{dE_\nu} \left( \frac{1}{2\pi} \frac{d\sigma_{\chi\nu}}{d\cos\psi} \right) \left( \frac{dE_\nu}{dT_\chi} \frac{v_\chi}{c} \right)
\end{equation}
where $dE_\nu / dT_\chi$ is given in Eq. \eqref{eq:dEv/dTx}.
The above equation describes the BDM emissivity at any point on the shell and has units of MeV<sup>−1</sup> cm<sup>−3</sup> s<sup>−1</sup> sr<sup>−1</sup>.
A more detailed derivation of Eq. \eqref{eq:jx} can be found in 
[Emissivity](emissivity.md){:target="_blank"}.



## SN$\nu$ BDM flux

When the emissivity at any point is known, one can calculate the BDM flux at Earth by integrating $j_\chi$ along the line of sight (l.o.s.) $\ell$. 
Thus, cf. [Figs. 1](#snv_bdm_scheme) and [3](#2d_shell),
\begin{equation}\label{eq:total_Phi}
\frac{d\Phi_\chi}{dT_\chi\, d\Omega} = \int d\ell\, j_\chi\, \Theta(r_\nu - d)\, \Theta(d + h - r_\nu)
\end{equation}
where $d\Omega$ denotes the field-of-view (f.o.v.) centered on the SN. 
The two $\Theta$ functions ensure that $j_\chi$ is only non-zero within the shell.
Without loss of generality, we assume $h \ll d$, so that
$$
\Theta(r_\nu - d)\, \Theta(d + h - r_\nu) \approx h\, \delta(r_\nu - d) = c\tau_s\, \delta(r_\nu - d),
$$
where $\delta(x)$ is the Dirac delta function. 
Thus, Eq. \eqref{eq:total_Phi} describes the total BDM flux observable on Earth.
Note that $\delta(r_\nu - d)$ has units of [L<sup>−1</sup>], which cancels the length dimension from $c\tau_s$. 
This ensures the approximation is dimensionless and self-consistent.

<figure id="2d_shell">
<center><img src="../../figs/2D_scheme.svg" alt="scheme" style="width: 45%;">
<figcaption>Figure 3. BDM production on the shell with thickness \(h\) at a distance \(d\) from the SN.  
When \(\ell\) and \(\theta\) are specified, \(j_\chi\) is non-zero at \(r_\nu\) only if \(d \leq r_\nu \leq d + h\).</figcaption>
</center>
</figure>


### From line-of-sight to time-dependency

However, Eq. \eqref{eq:total_Phi} is inadequate because it does not take the DM velocity into account.
Given that DM is massive, BDM cannot propagate at the same velocity as SN$\nu$. 
Depending on where the $\chi$ particle is upscattered, it will arrive at Earth at different times. 
This results in a broadened flux-vs-time profile compared to that of the SN$\nu$ burst.
BDM will arrive earlier if it was upscattered closer to Earth, and significantly later if it was upscattered farther away.

To incorporate the time-dependent behavior, let the time of the SN explosion define time-zero. 
Then, the time required for SN$\nu$ to propagate from S to E is
$$
t_\nu = \frac{R_s}{c}.
$$
For BDM produced at point B, the arrival time at Earth is
\begin{equation}\label{eq:tprime}
t^\prime = \frac{d}{c} + \frac{\ell}{v_\chi},
\end{equation}
where the first term accounts for the time taken by SN$\nu$ to travel from S to B, and the second term accounts for the BDM's time of flight to Earth. 
The BDM velocity is given by
$$
v_\chi = \frac{\sqrt{T_\chi (T_\chi + 2m_\chi)}}{T_\chi + m_\chi} c.
$$
By applying the law of cosines, we express $d$ as
\begin{equation}\label{eq:d}
d = \sqrt{\ell^2 + R_s^2 - 2\ell R_s \cos\theta}.
\end{equation}
Taking the total differential of Eq. \eqref{eq:tprime}, we find
\begin{equation}\label{eq:dtprime}
dt^\prime = \left( \frac{\ell - R_s \cos\theta}{c d} + \frac{1}{v_\chi} \right) d\ell = \mathcal{J}^{-1} d\ell.
\end{equation}
Note that $\mathcal{J}$ has the same dimension as velocity, i.e., [L T<sup>−1</sup>].
Before recasting Eq. \eqref{eq:total_Phi} into a time-dependent form, we first manipulate
$$
c\tau_s \delta(r_\nu - d) = \tau_s \delta\left( \frac{r_\nu}{c} - \frac{d}{c} \right) = \tau_s \delta\left( \frac{r_\nu}{c} - t^\prime + \frac{\ell}{v_\chi} \right) = \tau_s \delta\left( t^\prime - \frac{r_\nu}{c} - \frac{\ell}{v_\chi} \right),
$$
using the identity $\delta(ax) = \delta(x)/|a|$.
Under the thin-shell approximation ($h \ll d$), we write
$$
d \leq r_\nu \leq d + h \approx 1 \leq \frac{r_\nu}{d} \leq 1 + \frac{h}{d} \quad \Rightarrow \quad r_\nu \approx d.
$$
Thus, the delta function becomes
$$
c\tau_s \delta(r_\nu - d) \approx \tau_s \delta\left( t^\prime - \frac{d}{c} - \frac{\ell}{v_\chi} \right).
$$
Hence, we can recast Eq. \eqref{eq:total_Phi} from a line-of-sight integration into a time-domain integration:

$$
\frac{d\Phi_\chi}{dT_\chi d\Omega} = \tau_s \int dt^\prime\, \mathcal{J} j_\chi \delta\left( t^\prime - \frac{d}{c} - \frac{\ell}{v_\chi} \right) = \left. \tau_s \mathcal{J} j_\chi \right|_{t^\prime = \frac{d}{c} + \frac{\ell}{v_\chi}}.
$$

In the last step, the $\delta$ function constrains the integrand to be non-zero only at $t^\prime = d/c + \ell/v_\chi$.
Finally, we integrate over the f.o.v. as seen from Earth:

\begin{equation}\label{eq:BDM_flux}
\frac{d\Phi_\chi(t^\prime)}{dT_\chi} = \left. \tau_s \int_0^{2\pi} d\varphi \int_0^{\pi/2} d\theta\, \sin\theta \mathcal{J} j_\chi(d, r, T_\chi, \psi) \right|_{t^\prime = \frac{d}{c} + \frac{\ell}{v_\chi}}.
\end{equation}

This is the **time-dependent SN$\boldsymbol{\nu}$-induced BDM flux** observed at Earth.
To verify dimensional consistency:
$$
\left[ \tau_s \cdot d\varphi\, d\theta\, \sin\theta \cdot \mathcal{J} \cdot j_\chi \right] = \mathrm{s} \cdot \mathrm{sr} \cdot \frac{\mathrm{cm}}{\mathrm{s}} \cdot \frac{1}{\mathrm{MeV\;cm^3\;s\;sr}} = \mathrm{MeV^{-1}\;cm^{-2}\;s^{-1}},
$$
showing that Eq. \eqref{eq:BDM_flux} has the correct units for differential flux per energy.
Instead of setting $t^\prime = 0$ at the moment of SN explosion, it is more convenient to shift the reference point to the arrival of SN$\nu$ at Earth:
\begin{equation}\label{eq:t}
t = t^\prime - t_\nu = \frac{d}{c} + \frac{\ell}{v_\chi} - t_\nu,
\end{equation}
which leads to the final expression for the observed BDM flux:

\begin{equation}\label{eq:BDM_flux_earth}
\frac{d\Phi_\chi(t)}{dT_\chi} = \left. \tau_s \int_0^{2\pi} d\varphi \int_0^{\pi/2} d\theta\, \sin\theta \mathcal{J} j_\chi(d, r, T_\chi, \psi) \right|_{t = \frac{d}{c} + \frac{\ell}{v_\chi} - t_\nu}.
\end{equation}

Note that the prefactor $\tau_s$ is actually redundant. 
We have defined
$$
L_{\rm tot} = \frac{\mathcal{E}_{\rm tot}}{\tau_s},
$$
where $\mathcal{E}_{\rm tot} \approx 3 \times 10^{53}$ erg is the total energy released in the form of neutrinos (see [Tab. 1](#numerical_values)). 
Hence,

$$
\tau_s \times L_{\rm tot} = \mathcal{E}_{\rm tot}.
$$

In numerical calculations, we can therefore drop the explicit factor of $\tau_s$ and replace $L_{\rm tot}$ by $\mathcal{E}_{\rm tot}$ in Eq. \eqref{eq:snv_spectrum}.


### Time-dependent feature

Judging from Eq. \eqref{eq:BDM_flux_earth}, it is clear that the BDM flux from an individual SN is not everlasting. 
The final portion of the signal arrives from the BDM particles with the longest propagation time. 
Hence, we define the **vanishing time** as
$$
t_{\rm van} = \max \left[ \frac{d(\theta)}{c} + \frac{\ell(\theta)}{v_\chi} - t_\nu \right].
$$
Using the law of sines (cf. [Fig. 1](#snv_bdm_scheme)),
$$
\frac{\ell}{\sin(\psi - \theta)} = \frac{d}{\sin\theta} = \frac{R_s}{\sin\psi},
$$
we obtain the time-delay formula
\begin{equation}\label{eq:t_psi}
t = \frac{R_s}{c} \frac{\sin\theta}{\sin\psi} + \frac{R_s}{v_\chi} \frac{\sin(\psi - \theta)}{\sin\psi} - t_\nu.
\end{equation}
It is evident that, for fixed $\theta$, $t$ increases monotonically with $\psi \in [0, \pi/2]$. 
Therefore, to globally maximize $t$, we set $\psi = \psi_{\rm max}$ (cf. Eq. \eqref{eq:psi_max}), which depends only on $T_\chi$ and $m_\chi$.
To find the corresponding $\theta$ that yields the maximum $t$, we solve
\begin{equation}\label{eq:theta_maximum_t}
\frac{\cos\theta}{c} = \frac{\cos(\psi_{\rm max} - \theta)}{v_\chi}.
\end{equation}
This equation can be solved numerically for $\theta$. 
Suppose the solution is $\theta_{\rm MAX}$, then the vanishing time is given by
\begin{equation}\label{eq:tvan}
t_{\rm van} = \frac{d(\theta_{\rm MAX})}{c} + \frac{\ell(\theta_{\rm MAX})}{v_\chi} - t_\nu.
\end{equation}
We emphasize that Eq. \eqref{eq:tvan} provides the exact solution for $t_{\rm van}$, 
and it remains valid in both the relativistic and non-relativistic regimes. 
In contrast, the approximation given by Eq. (13) in Ref. [[2](#bib_Lin2023PRD)] is only suitable for the relativistic case.


### Field-of-view across the sky

From the previous subsection, we learned that $\theta_{\rm MAX}$ determines the maximum field of view f.o.v. angle at the vanishing time $t_{\rm van}$. 
Now, we consider the inverse problem: given a particular time $t^* < t_{\rm van}$, we ask what is the corresponding angular extent of the observable BDM flux.
Let $\theta^*_M$ satisfy Eq. \eqref{eq:t_psi}:

\begin{equation}
t^* = \frac{R_s}{c} \frac{\sin\theta^*_M}{\sin\psi_{\rm max}} + \frac{R_s}{v_\chi} \frac{\sin(\psi_{\rm max} - \theta^*_M)}{\sin\psi_{\rm max}} - t_\nu,
\end{equation}

then $\theta^*_M$ determines the maximal f.o.v. angle at observation time $t = t^*$.
Thus, at any given time, the region in the sky containing non-zero BDM flux corresponds to the angular range
\begin{equation}
\theta \in [0, \theta^*_M).
\end{equation}


#### References

1. <p id="bib_Lin2023PRL">Y.-H. Lin *et al.*, *Phys. Rev. Lett.* **130**, 111002 (2023)</p>
2. <p id="bib_Lin2023PRD">Y.-H. Lin *et al.*, *Phys. Rev. D.* **108**, 083013 (2023)</p>
3. <p id="bib_NFW">J. F. Navarro *et al.*, *Astrophys. J.* **462**, 563 (1996)</p>
4. <p id="bib_Duan2006PRD">H. Duan *et al.*, *Phys. Rev. D* **74**, 105014 (2006)</p>