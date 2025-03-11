<a href = "https://python.org" target = "_blank">![Python](https://img.shields.io/badge/python-3.8-blue.svg)</a>
<a href = "https://choosealicense.com/licenses/gpl-3.0/"  target = "_blank">![License](https://img.shields.io/badge/License-GPL_3.0-blue.svg)</a>
<a href = "https://arxiv.org/abs/2206.06864"  target = "_blank">![ArXiv](https://img.shields.io/badge/arXiv-2206.06864-yellowgreen.svg)</a>
<a href = "https://arxiv.org/abs/2307.03522"  target = "_blank">![ArXiv](https://img.shields.io/badge/arXiv-2307.03522-yellowgreen.svg)</a>

# snorer: *S*upernova-*N*eutrino-b*O*osted da*R*k matt*ER*


`snorer` is a package for evaluating time-of-flight signatures of supernova-neutrino-boosted dark matter (SN*ν* BDM) from our Milky Way (MW), SN1987a in Large Magellanic Cloud (LMC) and SN in arbitrary distant galaxy based on
<a href = "https://doi.org/10.1103/PhysRevLett.130.111002" target = "_blank">*Phys. Rev. Lett.* **130**, 111002 (2023)</a> [<a href = "https://arxiv.org/abs/2206.06864" target = "_blank">arXiv:2206.06864</a>]
and
<a href = "https://doi.org/10.1103/PhysRevD.108.083013" target = "_blank">*Phys. Rev. D* **108**, 083013 (2023)</a>
[<a href = "https://arxiv.org/abs/2307.03522" target = "_blank">arXiv:2307.03522</a>].

### Citation

If you use this package or part of the code in your research, please cite the followings:

1. Y.-H. Lin *et al.*, *Phys. Rev. Lett.* **130**, 111002 (2023), arXiv:2206.06864
2. Y.-H. Lin *et al.*, *Phys. Rev. D* **108**, 083013 (2023), arXiv:2307.03522
3. `snorer`: https://github.com/yenhsunlin/snorer


## Installation

To install, excute the following command on the prompt

    $ pip install snorer

and everything should be processed on-the-fly.

### Dependency

`snorer` requires python >= 3.8 and the following external packages

- `numpy` >= 1.20.0
- `scipy` >= 1.10.0
- `vegas` >= 6.0.1
- `astropy` >= 6.0.0

where `vegas` is a the backend engine for evaluating multidimensional integrals based on adaptive Monte Carlo vegas algorithm, see its homepage: <a heref = "https://pypi.org/project/vegas/" target = "_blank">https://pypi.org/project/vegas/</a>.

Other packages maybe required by these dependencies during the installation, see `requirements.txt` for details.
The versions of these dependencies are not strict, but are recommended to update to the latest ones to avoid incompatibility.

## Documentation

See the ***snorer*** documentation: https://yenhsunlin.github.io/snorer for more details.

## Known Issue

To evaluate BDM event, `snorer` uses `vegas` to handle the multidimensional integration.
The sampling method of `vegas` cannot manipulate event calculation, e.g. `snorer.event` and the method in the instance of `snorer.BoostedDarkMatter`, properly, when SN is exactly at GC with spike and no DM self-annihilation.

Since the spike is a highly singular behavior, the sampling method may miss the substantial DM contribution from the inner galactic region and causes underestimate of $N_{\rm BDM}$ plus unstable results. 
To avoid this, users may try to displace the SN from GC a little bit when evaluating $N_{\rm BDM}$ with DM sipke and no DM annihilation.
For BDM flux evaluation, there is no such issue.

To be fair, the probability of a very cuspy DM spike surving the gravitational disturbance without annihilating away and SN happening exactly at the GC might be very rare.

This issue is scheduled to fix in the future update.

## Bugs and troubleshooting

Please report to the author, Yen-Hsun Lin, via [yenhsun@phys.ncku.edu.tw](mailto:yenhsun@phys.ncku.edu.tw).

