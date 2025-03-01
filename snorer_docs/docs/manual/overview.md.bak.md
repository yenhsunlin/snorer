<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>

# Physics Overview

We review all the details provided in Refs. [[1](#bib_Lin2023PRL), [2](#bib_Lin2023PRD)].
While these papers present mathematical expressions suitable for general readers, the content in this section is tailored for programming purposes.
The goal of this document is not to replace Refs. [[1](#bib_Lin2023PRL), [2](#bib_Lin2023PRD)] but to bridge the gap between traditional journal writing and the technical details required for programming.
Typos and misinformation are corrected here to ensure clarity and facilitate understanding of the subsequent Python code.

## General picture

Given non-zero cross section between DM ($\chi$) and neutrino ($\nu$), it is expected that when SN$\nu$ could scatter with halo $\chi$ when propagates outward
from the explosion site. We display the general scheme in [Fig. 1](#snv_bdm_scheme). In this figure, SN, galactic center (GC) and Earth are denoted as S, G and E respectively.

<figure id="snv_bdm_scheme">
<center><img src="/figs/scheme.svg" alt="scheme" style="width: 60%;">
<figcaption>Figure 1. The 3D scheme for SN\(\nu\) BDM.
</figure>

In [Fig. 1](#snv_bdm_scheme), S, G and E are coplanar and $\chi$ is boosted at B on the SN$\nu$ shell. $d$ is the distance $\overline{\mathsf{SB}}$, $\ell$ the distance $\overline{\mathsf{BE}}$, $R_s$ the distance between SN and Earth and $R_e$ the distance between GC and Earth. The DM number density at B is determined by $n_\chi(r)$ where $r$ is the radial distance between GC and B. We adopt NFW profile  [[3](#bib_NFW)] for $n_\chi$,
$$
\begin{equation}
n_\chi(r)=\frac{\rho_s}{m_\chi}\frac{1}{\frac{r}{r_s}(1+\frac{r}{r_s})^2}
\end{equation}
$$
where $\rho_s=184$ MeV cm<sup>−3</sup> and $r_s=24.4$ kpc.
Due to its spherical symmetric nature, one does not need to know where the galactic plane lies. Unless S is on the axis $\overline{\mathsf{GE}}$, BDM is not azimuthally symmetric around $\varphi$.

### Particle kinematics

When a $\nu$ carries energy $E_\nu$ and scatters with a static $\chi$ at B, DM will receive kinetic energy $T_\chi$ by
$$
T_\chi = \frac{E_\nu^2}{E_\nu+m_\chi/2}\left(\frac{1+\cos\theta_c}{2}\right)
$$
where $\theta_c$ is the scattering angle at center-of-mass (cm) frame (do not confuse with the $\theta$ in [Fig. 1](#snv_bdm_scheme)). Assuming $\theta_c$ uniformly distributes in $[0,\pi]$, we can relate it with lab frame scattering $\psi$ by
$$
\theta_c =2\tan^{-1}(\gamma \tan\psi)
$$
where $\psi\in [0,\pi/2]$ and $\gamma=(E_\nu+m_\chi)/\sqrt{m_\chi(2E_\nu+m_\chi)}$.
Thus we have
$$
\begin{equation}
\frac{d\sigma_{\chi\nu}}{d\Omega_{\rm lab}}=\frac{1}{2\pi}\frac{d\sigma_{\chi\nu}}{d\cos\psi}=\sigma_0\times g_\chi(\psi)
\end{equation}
$$
where
$$
\begin{equation}\label{eq:gx}
g_\chi(\psi) = \frac{\gamma^2\sec^3\psi}{\pi(1+\gamma^2\tan^2\psi)^2} 
\end{equation}
$$
is the angular distribution for $\chi\nu$ scattering cross section in lab frame. One can easily verify that Eq. $\eqref{eq:gx}$ satisfies
$$
2\pi \int_0^{\pi/2} g_\chi(\psi) \sin\psi d\psi =1
$$
and is independent of $\gamma$.



When $d$, $\ell$ and $R_s$ are specified, one can determined $\psi$ by law of cosine
$$
R_s^2 = d^2+\ell^2-2d\ell \cos(\pi-\psi)
$$
where
$$
\begin{equation}\label{eq:geo_psi}
 \psi =\cos^{-1}\left( \frac{R_s^2-d^2-\ell^2}{2d\ell}\right).
\end{equation}
$$

### Constraint by positive-definite $E_\nu$

Note that though the valid range of $\psi$ from Eq. $\eqref{eq:geo_psi}$ is in $[0,\pi]$ but Eq. $\eqref{eq:gx}$ puts stronger constraint. As Eq. $\eqref{eq:geo_psi}$ is merely a geometrical relation but Eq. $\eqref{eq:gx}$ is required by the kinematics.
However, $[0,\pi/2]$ is simply the full range allowed by kinematics, if we examine further, we will find that the true valid range for $\psi$ is narrower.
Let's consider the $\chi\nu$ scattering in lab frame in [Fig. 2](#lab_scatt)

<figure id="lab_scatt">
<center><img src="/figs/lab_scattering.svg" alt="scheme" style="width: 40%;">
<figcaption>Figure 2. \(\chi\nu\) scattering in lab frame.
</figure>

We write down the four-momenta for every particles,

$$
\begin{align*}
p_{\nu} & =(E_{\nu},0,0,E_{\nu}),\\
p_{\chi} & =(m_{\chi},0,0,0),\\
p_{\nu}^{\prime} & =(E_{\nu}^{\prime},x,y,z),\\
p_{\chi}^{\prime} & =(E_{\chi},-|\mathbf{p}_{\chi}|\sin\psi,0,|\mathbf{p}_{\chi}|\cos\psi),
\end{align*}
$$

where $(x,y,z)$ indicates their values are irrelevant. Relativistic mechanics implies that, assuming $m_\nu=0$ and the metric tensor $g_{\mu\nu}={\rm diag}(1,-1,-1,-1)$,

$$
\begin{align*}
(p_\nu-p_\chi^\prime)^2 &= (p_\chi-p_\nu^\prime)^2, \\
m_\chi^2 - 2E_\nu (E_\chi -|\mathbf{p}_{\chi}|\cos\psi) &= m_\chi^2-2m_\chi E_\nu^\prime.
\end{align*}
$$

We can rewrite the last line by $E_\chi = T_\chi + m_\chi$ and $T_\chi = E_\nu - E_\nu^\prime$ which yields

$$
\begin{equation}\label{eq:Ev}
E_\nu = \frac{m_\chi T_\chi}{ |\mathbf{p}_{\chi}|\cos\psi - T_\chi}.
\end{equation}
$$

Once BDM kinetic energy $T_\chi$ and scattering angle $\psi$ at B are given, Eq. $\eqref{eq:Ev}$ specifies the required $E_\nu$ to have such properties.
Without ambiguity, $E_\nu$ is positive-definite and non-divergent, hence an even stronger constraint is put by

$$
\begin{equation}\label{eq:psi_max}
|\mathbf{p}_\chi| \cos\psi - T_\chi > 0 \to\psi < \psi_{\rm max}:= \cos^{-1}\left(\frac{T_\chi}{|\mathbf{p}_{\chi}|}\right).
\end{equation}
$$

It is also not hard to deduce that $|\mathbf{p}_{\chi}|=\sqrt{T_\chi (2m_\chi + T_\chi)}$.
Hence the realistic range for $\psi$ is
$$
\begin{equation}
\psi \in [0,\psi_{\rm max})
\end{equation}
$$
while $E_\nu$ diverges at $\psi=\psi_{\rm max}$.

## Dark emissivity

The most crucial thing is to know how many $\chi$s are boosted at B. This is generally characterized by the emissivity which has a unit of cm<sup>−3</sup> s<sup>−1</sup> and of MeV<sup>−1</sup> cm<sup>−3</sup> s<sup>−1</sup> in terms of energy spectrum.

### SN$\nu$ spectrum

We begin by writing down the SN$\nu$ energy spectrum on the shell at $d$
$$
\begin{equation}\label{eq:snv_spectrum}
\frac{dN_\nu}{dE_\nu} = \sum_i \frac{L_{\nu_i}}{4\pi d^2 \langle E_{\nu_i}\rangle} f_{\nu_i}(E_\nu)
\end{equation}
$$
where $L_{\nu_i}=L_{\rm tot}/6$ is the luminosity for each neutrino ($\nu_{e,\mu,\tau}$ and their anti-particles) and $f_{\nu_i}(E_\nu)$ is the Fermi-Dirac distribution  [[4](#bib_Duan2006PRD)],
$$
\begin{equation}
f_{\nu_i}(E_\nu) = \frac{1}{F_2(\eta_\nu)}\frac{1}{T_\nu^3}\frac{E_\nu^2}{e^{E_\nu/T_\nu-\eta_\nu}+1}
\end{equation}
$$
where
$$
\begin{equation}
F_k(\eta) = \int_0^\infty dx \frac{x^k}{e^{x-\eta}+1}.
\end{equation}
$$
We list all the numerical values for the parameters mentioned above in [Tab. 1](#numerical_values)

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

Suppose the duration of the SN explosion is around $\tau_s=10$ s, the total energy released in the form of neutrinos approximates $10^{53}$ erg for individual explosion.

### Number density on the shell

After converting erg to MeV, one can examine that Eq. $\eqref{eq:snv_spectrum}$, after mutiplying $\tau_s$, has the unit MeV<sup>−1</sup> cm<sup>−2</sup>.
Becasue $\nu$ travels at light speed, we can estimate the shell thickness by $h=c\tau_s$.
The SN$\nu$ shell would have a number density
$$
\begin{equation}\label{eq:snv_nd}
\frac{dn_\nu}{dE_\nu} = \frac{dN_\nu}{dE_\nu}\frac{\tau_s}{h}= \sum_i \frac{L_{\nu_i}}{4\pi d^2 \langle E_{\nu_i}\rangle c} f_{\nu_i}(E_\nu).
\end{equation}
$$
One examines that $dn_\nu/dE_\nu$ has MeV<sup>−1</sup> cm<sup>−3</sup>, which is exactly the energy spectrum of SN$\nu$ number density on the shell.

By mutiplying $E_\nu$, one can eliminate MeV<sup>−1</sup> and receive
$$
n_\nu =\frac{E_\nu}{c} \frac{dN_\nu}{dE_\nu}
$$
the SN$\nu$ number density.

### Emissivity on the shell

To obtain emissivity at B on the SN$\nu$ shell, one simply do
$$
j_\chi = c n_\nu n_\chi \frac{d\sigma_{\chi\nu}}{d\Omega_{\rm lab}}
$$
which has a unit cm<sup>−3</sup> s<sup>−1</sup> sr<sup>−1</sup> and sr<sup>−1</sup> indicates per steradian.
It would be more convenient to restore the energy spectrum form and translate it into $T_\chi$ expression,
$$
\begin{equation}\label{eq:jx}
j_\chi(d,r,T_\chi,\psi) = cn_\chi \frac{dn_\nu}{dE_\nu}  \left(\frac{1}{2\pi}\frac{d\sigma_{\chi\nu}}{d\cos\psi}\right)\frac{dE_\nu}{dT_\chi}
\end{equation}
$$
where $dE_\nu/dT_\chi$ can be obtained by differentiating Eq. $\eqref{eq:Ev}$ w.r.t $T_\chi$ directly.
Eq. $\eqref{eq:jx}$ now describes the BDM emissivity at any point on the shell
and has a unit MeV<sup>−1</sup> cm<sup>−3</sup> s<sup>−1</sup> sr<sup>−1</sup>.



## SN$\nu$ BDM flux

When the emissivity at any point is known, one can calculate the BDM flux at Earth by integrating $j_\chi$ over the ling-of-sight (l.o.s) $\ell$.
Thus, cf. [Figs. 1](#snv_bdm_scheme) and [3](#2d_shell),
$$
\begin{equation}\label{eq:total_Phi}
\frac{d\Phi_\chi}{dT_\chi d\Omega}=\int d\ell ~  j_\chi \Theta(r_\nu - d)\Theta(d+h-r_\nu)
\end{equation}
$$
where $d\Omega$ is the field-of-view centering SN.
The two $\Theta$ functions restrict $j_\chi$ is only non-zero within the shell.
Without loss of generality, $h\ll d$ always holds such that 
$$
\Theta(r_\nu-d)\Theta(d+h -r_\nu)\approx h \delta(r_\nu-d) =c\tau_s \delta(r_\nu-d)
$$
where $\delta(x)$ is the Dirac-delta function. Thus, Eq. $\eqref{eq:total_Phi}$ describes the total BDM flux that could be observed on the Earth.
Note that $\delta(r_\nu-d)$ actually carries \[L<sup>−1</sup>\] and cancels the dimensiton from $c\tau_s$.
This approximation remains dimensionless and is self-consistent.

<figure id="2d_shell">
<center><img src="/figs/2D_scheme.svg" alt="scheme" style="width: 45%;">
<figcaption>Figure 3. BDM production on the shell with thickness \(h\) at \(d\) distant to SN.
When \(\ell\) and \(\theta\) are specified, \(j_\chi\) is non-zero at \(r_\nu\) only when \(d\leq r_\nu\leq d+h\).
</figure>

### From line-of-sight to time-dependency

However, Eq. $\eqref{eq:total_Phi}$ is inadequate because it does not take DM velocity into account.
Given massive DM, BDM cannot have the same velocity as SN$\nu$. Depending on where $\chi$ is upscattered, it will arrive on Earth at different times. This results in a prolonged flux vs. time compared to the SN$\nu$ flux. BDM arrives earlier if it was upscattered closer and much later if it was upscattered farther away.

To incorporate time-dependent feature, we first set the time of SN explosion as the time-zero. Then SN$\nu$ needs
$$
t_\nu = \frac{R_s}{c}
$$
to propagate from S to E. For BDM produced at B, it will arrive at
$$
\begin{equation}\label{eq:tprime}
t^\prime = \frac{d}{c} + \frac{\ell}{v_\chi}
\end{equation}
$$
where the first term accounts for SN$\nu$ traveling from S to B and the BDM velocity
$$
v_\chi =\frac{\sqrt{T_\chi(T_\chi+2m_\chi )}}{T_\chi+m_\chi}c.
$$
Applying law of cosine, we have
$$
\begin{equation}\label{eq:d}
d = \sqrt{\ell^2 + R_s^2 - 2\ell R_s \cos\theta}
\end{equation}
$$
and take total differentiation on both sides of Eq. $\eqref{eq:tprime}$,
$$
\begin{equation}\label{eq:dtprime}
dt^\prime =\left( \frac{\ell - R_s \cos\theta}{c d} + \frac{1}{v_\chi} \right) d\ell := \mathcal{J}^{-1}d\ell.
\end{equation}
$$
Note that $\mathcal{J}$ has the same dimension as velocity \[L T<sup>−1</sup>\].
Before recasting Eq. $\eqref{eq:total_Phi}$ into time-dependent form, we firstly manage
$$
c\tau_s \delta(r_\nu-d) = \tau_s \delta\left(\frac{r_\nu}{c}- \frac{d}{c}\right) = \tau_s\delta\left(\frac{r_\nu}{c}- t^\prime +\frac{\ell}{v_\chi}\right)=\tau_s\delta\left(t^\prime - \frac{r_\nu}{c}-\frac{\ell}{v_\chi}\right).
$$
Using thin-shell approximation, $h\ll d$, we approximate
$$
d \leq r_\nu \leq d+h \approx  1 \leq \frac{r_\nu}{d} \leq 1+\frac{h}{d} \to r_\nu\approx d.
$$
The $\delta$-function becomes
$$
c\tau_s\delta(r_\nu-d) \approx \tau_s\delta \left(t^\prime - \frac{d}{c}-\frac{\ell}{v_\chi} \right).
$$

Hence we can recast Eq. $\eqref{eq:total_Phi}$ into

$$
\frac{d\Phi_\chi}{dT_\chi d\Omega} = \tau_s \int dt^\prime  \mathcal{J}  j_\chi \delta \left(t^\prime - \frac{d}{c}-\frac{\ell}{v_\chi} \right) = \left. \tau_s\mathcal{J}  j_\chi \right|_{t^\prime = \frac{d}{c}+\frac{\ell}{v_\chi}},
$$

and finally,

$$
\begin{equation}\label{eq:BDM_flux}
\frac{d\Phi_\chi(t^\prime)}{dT_\chi} =  \left. \tau_s\int_0^{2\pi} d\varphi \int_0^{\pi/2} d\theta~\sin\theta\mathcal{J}  j_\chi (d,r,T_\chi,\psi) \right|_{t^\prime = \frac{d}{c}+\frac{\ell}{v_\chi}}
\end{equation}
$$

which is the **time-dependent SN$\boldsymbol{\nu}$ BDM flux** at Earth.
One last thing is to check the dimension,
$$
\left[\tau_s \times  d\varphi d\theta \sin\theta \times  \mathcal{J} \times j_\chi \right] = {\rm s}\cdot{\rm sr}\cdot\frac{\rm cm}{\rm s} \cdot \frac{1}{\rm MeV~cm^3~s~sr}= {\rm MeV^{-1}~cm^{-2}~s^{-1}}
$$
that shows Eq. $\eqref{eq:BDM_flux}$ has the same unit as flux per energy width.

The time-zero for $t^\prime$ is set to the SN explosion. For our convenience, we can shift it by $t_\nu$,
$$
\begin{equation}\label{eq:t}
t := t^\prime - t_\nu = \frac{d}{c} + \frac{\ell}{v_\chi} - t_\nu
\end{equation}
$$
and this offsets the time-zero to the arrival of SN$\nu$ at Earth,

$$
\begin{equation}\label{eq:BDM_flux_earth}
\frac{d\Phi_\chi(t)}{dT_\chi} =  \left. \tau_s\int_0^{2\pi} d\varphi \int_0^{\pi/2} d\theta~\sin\theta\mathcal{J}  j_\chi (d,r,T_\chi,\psi) \right|_{t = \frac{d}{c}+\frac{\ell}{v_\chi}-t_\nu}.
\end{equation}
$$

### Time-dependent feature

Judging by Eq. $\eqref{eq:BDM_flux_earth}$, it is easy to realize that the BDM flux from an individual SN is not everlasting.
The last bit of signal comes from the BDM with maximum propagation time.
Hence we define a vanishing time,
$$
t_{\rm van} = \max \left[ \frac{d(\theta)}{c} + \frac{\ell(\theta)}{v_\chi} - t_\nu\right].
$$
Using law of sine, cf. [Fig. 1](#snv_bdm_scheme),
$$
\frac{\ell}{\sin(\psi-\theta)}=\frac{d}{\sin\theta}=\frac{R_s}{\sin\psi}
$$
then
$$
\begin{equation}\label{eq:t_psi}
t =\frac{R_s}{c}\frac{\sin\theta}{\sin\psi}+\frac{R_s}{v_\chi}\frac{\sin(\psi-\theta)} {\sin\psi}.
\end{equation}
$$
We omit $t_\nu$ since it is a constant.
Obviously, when we fix $\theta$ in Eq. $\eqref{eq:t_psi}$, $t$ increases monotonically with $\psi \in [0,\pi/2]$.
Thus to find the corresponding $\theta$ that maximizes $t$, we restrict ourselves at $\psi=\psi_{\rm max}$, cf. Eq. $\eqref{eq:psi_max}$.
Hence we can do $dt/d\theta=0$ and solve
$$
\begin{equation}\label{eq:theta_maximum_t}
\frac{\cos\theta}{c}=\frac{\cos(\psi_{\rm max}-\theta)}{v_\chi}
\end{equation}
$$
for $\theta$. Suppose $\theta^*$ satisfies Eq. $\eqref{eq:theta_maximum_t}$, then we have

$$
\begin{equation}\label{eq:tvan}
t_{\rm van} = \frac{d(\theta^*)}{c} + \frac{\ell(\theta^*)}{v_\chi} - t_\nu.
\end{equation}
$$

We note that Eq. $\eqref{eq:tvan}$ is the exact solution for $t_{\rm van}$ and valid for both relativistic and non-relativistic cases. However, an approximation given by Eq. (13) in Ref. [[2](#bib_Lin2023PRD)] is only suitable for relativistic case.


1. <p id="bib_Lin2023PRL">Y.-H. Lin *et al.*, *Phys. Rev. Lett.* **130**, 111002 (2023)</p>
2. <p id="bib_Lin2023PRD">Y.-H. Lin *et al.*, *Phys. Rev. Lett.* **108**, 083013 (2023)</p>
3. <p id="bib_NFW">J. F. Navarro *et al.*, *Astrophys. J.* **462**, 563 (1996)</p>
4. <p id="bib_Duan2006PRD">H. Duan *et al.*, *Phys. Rev. D* **74**, 105014 (2006)</p>