<a href = "https://python.org" target = "_blank">![Python](https://img.shields.io/badge/python-3.8-blue.svg)</a>
<a href = "https://choosealicense.com/licenses/gpl-3.0/"  target = "_blank">![License](https://img.shields.io/badge/License-GPL_3.0-blue.svg)</a>
<a href = "https://arxiv.org/abs/2206.06864"  target = "_blank">![ArXiv](https://img.shields.io/badge/arXiv-2206.06864-yellowgreen.svg)</a>
<a href = "https://arxiv.org/abs/2307.03522"  target = "_blank">![ArXiv](https://img.shields.io/badge/arXiv-2307.03522-yellowgreen.svg)</a>

# snorer: *S*upernova-*N*eutrino-b*O*osted da*R*k matt*ER*


`snorer` is a package for evaluating time-of-flight signatures of supernova-neutrino-boosted dark matter (SN$\nu$ BDM)from our Milky Way (MW) and SN1987a in Large Magellanic Cloud (LMC) based on
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




### Useful constants

We document various useful physical constants and conversion factors...etc as the *attributes* of an instance `snorer.constant`.
For example, we can retrieve

    >>> snorer.constant.me      # electron mas, MeV
    0.511
    >>> snorer.constant.kpc2cm  # kpc to cm
    3.085e+21

We did this instead of naming them as constant variables in some module to prevent them from being edited.

### BDM velocity

A boosted dark matter (BDM) with mass $m_\chi$ and kinetic energy $T_\chi$ has the velocity $v_\chi$,

$$
\frac{v_\chi}{c} = \frac{\sqrt{T_\chi(2m_\chi+T_\chi)}}{m_\chi+T_\chi}.
$$

Let $(T_\chi,m_\chi)=(15,0.075)$ MeV, we can use
`snorer.get_vx()` to evaluate

    >>> Tx,mx = 15,0.075
    >>> snorer.get_vx(Tx,mx)
    0.9999876239921284


### SN 

$$
\frac{d\Phi_{\chi}(T_\chi, t^\prime)}{dT_{\chi}dt} =
\left.\tau\int_0^{2\pi} d\phi\int_{0}^{\pi/2}\sin\theta d\theta~ \mathcal{J} j_{\chi}(r(\phi),D,T_{\chi},\psi)\right|_{t^{\prime}=\frac{D}{c}+\frac{d}{v_{\chi}}}
$$

## Bugs and troubleshooting

Please report to the author, Yen-Hsun Lin, via [yenhsun@phys.ncku.edu.tw](mailto:yenhsun@phys.ncku.edu.tw).

