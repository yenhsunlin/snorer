<script>
window.MathJax = {
  tex: {
    tags: "ams"  // Auto-numbering, AMS based
  }
};
</script>

# 2-2 elastic scattering

In [User Manual/Physics Overview <i class="fa-regular fa-bookmark"></i>](../../manual/overview.md#particle-kinematics){:target="_blank"} we have reviewed the 2-2 elastic scattering with one is massless neutrino.
In order to facilitate a broader case, we relax the massless assumption and do not distinguish them as any particle species.
This aids us to build an independent class that is suitable for any 2-2 scattering with non-zero masses.

## General expressions

The scheme for scattering of two particles, labeled 1 and 2, is shown in [Fig. 1](#22scat). Each particle carries non-zero mass $m_1$ and $m_2$ respectively.
Four-momenta for evary particles are shown on the sides. The scattering angle measures the degrees of deflection relative to the incoming direction of particle 1.
After scattering, particles 1 and 2 have scattering angles $\vartheta$ and $\psi$ respectively. 

<figure id="22scat">
<center><img src="../../../figs/22scat.svg" alt="22scat" style="width: 40%;">
<figcaption>Figure 1. The 2-2 particle scattering in lab frame.
</figure>

<!-- We can write down the 4-momenta in c.m. frame

$$
\begin{align*}
p_1^c &= (E_1^c,\mathbf{p}),\\
p_2^c &= (E_2^c,-\mathbf{p}),\\
p_1^{c\prime} &= (E_1^c,\mathbf{p}^\prime),\\
p_2^{c\prime} &= (E_2^c,-\mathbf{p}^\prime).\\
\end{align*}
$$

The energies do not change due to elastic scattering only deflects momenta directions.
We also have the following relations in c.m. frame [[1](#bib_ConceptQFT)],

$$
\begin{align}
E_1^c &=\frac{1}{2\sqrt{s}}(s+m_1^2-m_2^2)=E_1^{c\prime},\\
E_2^c &=\frac{1}{2\sqrt{s}}(s+m_2^2-m_1^2)=E_2^{c\prime},\\
|\mathbf{p}| & =\frac{1}{2\sqrt{s}} \lambda^{1/2}(s,m_1^2,m_2^2)=|\mathbf{p}^\prime|,
\end{align}
$$

where
$$
\begin{equation}
\lambda(x,y,z)= x^2+y^2+z^2-2(xy+yz+xz)
\end{equation}
$$
is the K&auml;llen lambda function.
 -->

We can write down the 4-momenta in lab frame,

$$
\begin{align*}
p_1 = (E_1,\mathbf{p}_1),&\quad p_2 = (m_2,\mathbf{0}),\\
p_1^{\prime} = (E_1^\prime,\mathbf{p}_1^\prime),&\quad p_2^{\prime} = (E_2^\prime,\mathbf{p}_2^\prime),
\end{align*}
$$

<!-- and the following identities [[1](#bib_ConceptQFT)],

$$
\begin{align}
E_1 &=\frac{1}{2m_2}(s-m_1^2-m_2^2),\\
E_1^\prime &=\frac{1}{2m_2}(m_1^2+m_2^2-u),\\
E_2^\prime &=\frac{1}{2m_2}(2m_2^2-t),
\end{align}
$$

where $s$, $t$ and $u$ are the Mandelstam variables and statisfy
$$
\begin{equation}\label{eq:Man_relation}
s+t+u = 2(m_1^2+m_2^2).
\end{equation}
$$ -->

and the corresponding $u$-channel 

$$
\begin{gather*}
(p_1-p_2^\prime)^2 = m_1^2+m_2^2-2E_1E_2^\prime +2 |\mathbf{p}_1| |\mathbf{p}_2^\prime| \cos\psi,\\
(p_2-p_1^\prime)^2 = m_1^2+m_2^2 - 2 E_1^\prime m_2.
\end{gather*}
$$

The two identities are equivalent due to Lorentz invariant nature.
Similar to the BDM case, suppose we know $E_2^\prime$ and its kinetic energy such that
$$
T_2 = E^\prime_2 - m_2 = E_1 - E_1^\prime.
$$
We thus have 
$$
\begin{equation}
E_1 (T_2+m_2) -|\mathbf{p}_1||\mathbf{p}_2^\prime| x = (E_1-T_2)m_2,
\end{equation}
$$
where $x:= \cos\psi$. By letting $|\mathbf{p}_1|=\sqrt{E_1^2-m_1^2}$ and $|\mathbf{p}_2^\prime| =\sqrt{T_2(T_2+2m_2)}$, the only unknown in the above equation is $E_1$, which can be solved analytically,

$$
\begin{equation}\label{eq:E1}
E_1=\frac{T_2^2 m_2 + |\mathbf{p}_2^\prime|x\sqrt{m_1^2 |\mathbf{p}_2^\prime|^2 x^2 + T_2^2 (m_2^2-m_1^2)}}{|\mathbf{p}_2^\prime|^2 x^2 - T_2^2}.
\end{equation}
$$

Note the $E_1$ gives the total energy of particle 1, thus $T_1=E_1 - m_1$. It is true that $E_1 = T_1$ only when $m_1=0$.


Moreover, one can differentiate $E_1$ w.r.t. $T_2$ and it yields
$$
\begin{equation}\label{eq:dE1/dT2}
\frac{dE_1}{dT_2}=m_2x^2\times \frac{\alpha +\beta + \gamma}{\eta}
\end{equation}
$$
where

$$
\begin{align*}
\alpha &= m_1^2 \delta, \\
\beta & = m_2^2 (2T_2+\delta),\\
\gamma &= 2m_2 x \kappa,\\
\eta &= \delta^2  x \kappa,
\end{align*}
$$

with $\delta=-T_2 + (T_2 + 2m_2)x^2$ and $\kappa=\sqrt{(T_2+2m_2)(\alpha+ T_2 m_2^2 )}$.

In most of the situations, elastic scattering does not change particle mass, thus
$$
\begin{equation}
\frac{dE_1}{dT_2} = \frac{d}{dT_2}(T_1 + m_1) = \frac{dT_1}{dT_2}
\end{equation}
$$
and this allows to express every quantities in terms of kinetic energies $T_i$.
When we construct the corresponding `class`, the namespace will always refers to $T_i$ instead of $E_i$.


Though the angle $\vartheta$ for particle 1 is irrelevant to our study, we can still determine it by 3-momentum conservation,
$$
\begin{equation}
\vartheta =\sin^{-1}\left( \frac{|\mathbf{p}_2^\prime|}{|\mathbf{p}_1^\prime|}\sin\psi\right)
\end{equation}
$$
where $|\mathbf{p_1^\prime}|=\sqrt{E_1^{\prime 2}-m_1^2}$ and $E_1^\prime = E_1-T2$.

## Validation

Now recall the $\nu\chi$ scattering 
with $E_1=E_\nu$, $T_2=T_\chi$, $m_1=m_\nu=0$ and $m_2=m_\chi$.
With $|\mathbf{p}_2^\prime|=|\mathbf{p}_\chi|=\sqrt{T_\chi(T_\chi+m_\chi)}$, Eq. $\eqref{eq:E1}$ becomes

$$
\begin{equation}\label{eq:Ev}
E_\nu = T_\chi m_\chi \frac{T_\chi + |\mathbf{p}_\chi|x }{( |\mathbf{p}_\chi|x-T_\chi)(|\mathbf{p}_\chi|x+T_\chi)}=\frac{T_\chi m_\chi}{|\mathbf{p}_\chi|x-T_\chi}.
\end{equation}
$$

Additionally,

$$
\begin{align*}
\alpha &= 0, \\
\beta & = m_\chi^2 (T_\chi+(T_\chi +2m_\chi)x^2),\\
\gamma &= 2m_\chi^2 |\mathbf{p}_\chi|x,\\
\eta &= \delta^2 m_\chi |\mathbf{p}_\chi|x ,
\end{align*}
$$

and with some tedious math one obtains

$$
\begin{align*}
m_\chi x^2\frac{\alpha+\beta+\gamma}{\eta}&= \frac{x}{\delta^2}\frac{m_\chi^2}{|\mathbf{p}_\chi|}(T_\chi+(T_\chi+2m_\chi)x^2+2|\mathbf{p}_\chi|x) \\
&=\frac{x}{\delta^2}\frac{m_\chi^2}{|\mathbf{p}_\chi|}\left(T_\chi +\frac{|\mathbf{p}_\chi|^2}{T_\chi}x^2+2|\mathbf{p}_\chi|x\right)\\
& = \frac{x}{\delta^2}\frac{m_\chi^2}{|\mathbf{p}_\chi|} \frac{(T_\chi + |\mathbf{p}_\chi|x)^2}{T_\chi }
\end{align*}
$$

we further use

$$
\begin{align*}
\delta^2 & =\left(-T_\chi + \frac{|\mathbf{p}_\chi|^2}{T_\chi}x^2\right)^2=\frac{(|\mathbf{p}_\chi|^2x^2 - T_\chi^2)^2}{T_\chi^2}\\
&= \frac{(|\mathbf{p}_\chi|x - T_\chi)^2(|\mathbf{p}_\chi|x + T_\chi)^2}{T_\chi^2}
\end{align*}
$$

then finally reaches

$$
\begin{equation}\label{eq:dEv/dTx}
m_\chi x^2\frac{\alpha+\beta+\gamma}{\eta}=\left(\frac{m_\chi}{|\mathbf{p}_\chi|x - T_\chi}\right)^2\frac{T_\chi}{|\mathbf{p}_\chi|}x = \frac{dE_\nu}{dT_\chi}.
\end{equation}
$$


We discover that Eqs. $\eqref{eq:Ev}$ and $\eqref{eq:dEv/dTx}$ are exactly Eqs. (5) and (8) in [User Manual/Physics Overview <i class="fa-regular fa-bookmark"></i>](../../manual/overview.md#constraint-by-positive-definite-e_nu){:target="_blank"} respectively.
Thus conclude that
Eqs. $\eqref{eq:E1}$ and $\eqref{eq:dE1/dT2}$ are the general expressions for massive 2-2 elastic scattering.



<!-- 1. <p id="bib_ConceptQFT">V. Ilisie, *Concepts in quantum field theory: A practitioner's toolkit*, Springer (2016)</p>  -->