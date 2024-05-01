<a href = "https://python.org" target = "_blank">![Python](https://img.shields.io/badge/python-3.8-blue.svg)</a>
<a href = "https://choosealicense.com/licenses/gpl-3.0/"  target = "_blank">![License](https://img.shields.io/badge/License-GPL_3.0-blue.svg)</a>
<a href = "https://arxiv.org/abs/2206.06864"  target = "_blank">![ArXiv](https://img.shields.io/badge/arXiv-2206.06864-yellowgreen.svg)</a>
<a href = "https://arxiv.org/abs/2307.03522"  target = "_blank">![ArXiv](https://img.shields.io/badge/arXiv-2307.03522-yellowgreen.svg)</a>

# snorer: *S*upernova-*N*eutrino-b*O*osted da*R*k matt*ER*


`snorer` is a package for evaluating time-of-flight signatures of supernova-neutrino-boosted dark matter (SN*ν* BDM) from our Milky Way (MW) and SN1987a in Large Magellanic Cloud (LMC) based on
<a href = "https://doi.org/10.1103/PhysRevLett.130.111002" target = "_blank">*Phys. Rev. Lett.* **130**, 111002 (2023)</a> [<a href = "https://arxiv.org/abs/2206.06864" target = "_blank">arXiv:2206.06864</a>]
and
<a href = "https://doi.org/10.1103/PhysRevD.108.083013" target = "_blank">*Phys. Rev. D* **108**, 083013 (2023)</a>
[<a href = "https://arxiv.org/abs/2307.03522" target = "_blank">arXiv:2307.03522</a>].

## Installation

To install, excute the following command on the prompt

    $ pip install snorer

and everything should be processed on-the-fly.

### Dependency

`snorer` requires python >= 3.8 and the following packages these external packages

- `numpy` >= 1.20.0
- `scipy` >= 1.10.0
- `vegas` >= 6.0.1
- `astropy` >= 6.0.0

where `vegas` is a the backend engine for evaluating multidimensional integrals based on adaptive Monte Carlo vegas algorithm, see its homepage: <a heref = "https://pypi.org/project/vegas/" target = "_blank">https://pypi.org/project/vegas/</a>.

Other packages maybe required by these dependencies during the installation, see `requirements.txt` for details.
The versions of these dependencies are not strict, but are recommended to update to the latest ones to avoid incompatibility.

## Usage

We briefly summarize the usage in this section and a comprehensive tutorial can be found in the jupyter notebook in `examples/tutorial.ipynb`.

To import, do

    >>> import snorer

in python terminal and is similar in the jupyter notebook.
All functions and classes can be accessed by typing `snorer.foo()` where `foo` is the name.




### Useful Constants

We document various useful physical constants and conversion factors...etc as the *attributes* of an instance `snorer.constant`.
For example, we can retrieve

    >>> snorer.constant.me      # electron mas, MeV
    0.511
    >>> snorer.constant.kpc2cm  # kpc to cm
    3.085e+21

We did this instead of naming them as constant variables in some module to prevent them from being edited.

### BDM Velocity

A boosted dark matter (BDM) with mass $m_\chi$ and kinetic energy $T_\chi$ has the velocity $v_\chi$,

$$
\frac{v_\chi}{c} = \frac{\sqrt{T_\chi(2m_\chi+T_\chi)}}{m_\chi+T_\chi}.
$$

Let $(T_\chi,m_\chi)=(15,0.075)$ MeV, we can use
`snorer.get_vx()` to evaluate

    >>> Tx,mx = 15,0.075
    >>> snorer.get_vx(Tx,mx)
    0.9999876239921284


### SN*ν* BDM Flux 

The BDM flux, or *afterglow* to the ${\rm SN}\nu$, due to SN that is $R_\star$ distant away from us can be evaluated by 

$$
\frac{d\Phi_{\chi}(T_\chi, t^\prime)}{dT_{\chi}dt} =
\tau\int_0^{2\pi} d\phi\int_{0}^{\pi/2}\sin\theta d\theta~ \mathcal{J} j_{\chi}(r(\phi),D,T_{\chi},\psi)
$$

where $t$ is the BDM ToF with time-zero at the discovery of SN*ν* on Earth and $t^\prime$ is the total time. We focus on $t$ instead of $t^\prime$.
Zenith angle $\theta$ and azimuthal angle $\phi$ are relative to the SN-Earth line-of-sight. The default DM-*ν* cross section is $\sigma_{\chi\nu}=10^{-45}$ cm<sup>2</sup>.

The function to evaluate this flux is `snorer.flux()` with $(t,T_\chi,m_\chi,R_\star,\beta)$ are the necessary inputs. 
Suppose SN's location is at GC, we have $R_\star=8.5$ kpc and $\beta=0$, and examine the flux with turning on DM spike feature

    >>> t,Tx,mx,Rstar,beta = 100,15,1e-2,8.5,0
    >>> snorer.flux(t,Tx,mx,Rstar,beta,neval=15000)
    4.572295175982701e-16

Users can turn off the spike feature by inserting `is_spike=False` and will find both numerical results are similar. It implies the contribution to the BDM flux is due to the place outside the spike's influencial region.

### SN*ν* BDM Event

The BDM event number in a detector after exposing to the flux for a period of time, $(t_{\rm min},t_{\rm max})$, can be evaluated by

$$
N_{\rm BDM} = \int_{t_{\rm min}}^{\rm t_{\rm max}} dt \int_{T_{\chi,{\rm min}}}^{T_{\chi,{\rm max}}}
\frac{d\Phi_{\chi}(T_\chi, t^\prime)}{dT_{\chi}dt} \times N_e \sigma_{\chi e}
$$

where $N_e$ is the total electron number in the detector and $\sigma_{\chi e}$ the DM-*e* cross section.
This can be accomplished by `snorer.event()` with $(m_\chi,R_\star,\beta)$ the necessary inputs.
The default $(t_{\rm min},t_{\rm max})=(10~{\rm s},35~{\rm yrs})$ and $(T_{\chi,{\rm min}},T_{\chi,{\rm max}})=(5,30)$ MeV.

Note that this function is normalized to $N_e=1$ and $\sigma_{\chi e}=1$ cm<sup>2</sup>.
To have the correct $N_{\rm BDM}$ in a specific detector, users have to mutiply the corresponding $N_e$ and $\sigma_{\chi e}$.

Now let $m_\chi=0.015$ MeV,

    >>> mx,Rstar,beta = 0.015,8,0
    >>> N_BDM = snorer.event(mx,Rstar,beta,is_spike=False,neval=50000)
    >>> N_BDM
    1.662174035857532e-06

Suppose it happened in Super-Kamiokande with $N_e\approx 7\times 10^{33}$ and assume $\sigma_{\chi e}=10^{-35}$ cm<sup>2</sup>. The correct $N_{\rm BDM}$ would be

    >>> Ne,sigma_xe = 7e33,1e-35
    >>> N_BDM*Ne*sigma_xe
    1.1635218251002724e-07

### *Experimental* :: Implementation of Particle Physics Model and SN in Arbitrary Distant Galaxy

The aforementioned functions for evaluating BDM signatures are based on model-agnostic picture. It means the cross sections between dark and visble sectors are generally independent of any physical quantities, eg. energy, mass and coupling constants.

The most important feature of `snorer` is that it offers a general interface for users to implement their favorite particle models.
Furthermore, SN is not necessary residing in MW or LMC. As long as users can provide these celetial objects' coordinates expressed in *ICRS J2000.0* system, `snorer` can do the calculation.
`snorer` also allows users to customize the halo shape, by manipulating $\rho_s$, $r_s$ and $n$...etc, and including or excluding spike feature, of such distant galaxy.

This will be done by introducing a *class* `snorer.GeneralInterface`.
All these user-specified features will compose an instance of `snorer.GeneralInterface`.
The BDM signatures can be evaluated by calling the associated *methods* within it.
We have an example in `examples/tutorial.ipynb`, also see the in-class docstring for more information.


### Other Useful Functions Classes

We also provide many useful functions and classes at users' disposal. See `examples/tutorial.ipynb` for details.

## Known Issue

To evaluate BDM event, `snorer` uses `vegas` to handle the multidimensional integration.
The sampling method of `vegas` cannot manipulate`snorer.event()` as well as the method in the instance of `snorer.GeneralInterface` properly, when SN is exactly at GC with spike feature turning on and no DM self-annihilation.

Since the spike is a highly singular behavior, the sampling method may miss the substantial DM contribution from the inner galactic region and causes underestimate of $N_{\rm BDM}$ plus unstable results. 
To avoid this, users may try to displace the SN from GC when evaluating $N_{\rm BDM}$ whith DM sipke turning on and no DM annihilation.
For BDM flux evaluation, there is no such issue.

To be fair, the probability of a very cuspy DM spike surving the gravitational disturbance without annihilating away and SN happening exactly at the GC might be very rare.

This issue is scheduled to fix in the next version of `snorer`.

## Bugs and troubleshooting

Please report to the author, Yen-Hsun Lin, via [yenhsun@phys.ncku.edu.tw](mailto:yenhsun@phys.ncku.edu.tw).

