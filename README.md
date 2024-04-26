[![Python](https://img.shields.io/badge/python-3.8-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/License-GPL_3.0-blue.svg)](https://choosealicense.com/licenses/gpl-3.0/)
[![ArXiv](https://img.shields.io/badge/arXiv-2206.06864-yellowgreen.svg)](https://arxiv.org/abs/2206.06864) 
[![ArXiv](https://img.shields.io/badge/arXiv-2307.03522-yellowgreen.svg)](https://arxiv.org/abs/2307.03522) 

# snorer: *S*upernova-*N*eutrino-b*O*osted da*R*k matt*ER*


`snorer` is a package for evaluating time-of-flight signatures of supernova-neutrino-boosted dark matter from our Milky Way and SN1987a in Large Magellanic Cloud based on [*Phys. Rev. Lett.* **130**, 111002 (2023)](https://doi.org/10.1103/PhysRevLett.130.111002) [[arXiv:2206.06864](https://arxiv.org/abs/2206.06864)] and [*Phys. Rev. D* **108**, 083013 (2023)](https://doi.org/10.1103/PhysRevD.108.083013) [[arXiv:2307.03522](https://arxiv.org/abs/2307.03522)].

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

where `vegas` is a the backend engine for evaluating multidimensional integrals based on adaptive Monte Carlo vegas algorithm, see its homepage: [https://pypi.org/project/vegas/](https://pypi.org/project/vegas/).

Other packages maybe required by these dependencies during the installation, see `requirements.txt` for details.
The versions of these dependencies are not strict, but are recommended to update to the latest ones to avoid incompatibility.




## Bugs and troubleshooting

Please report to the author, Yen-Hsun Lin, via [yenhsun@phys.ncku.edu.tw](mailto:yenhsun@phys.ncku.edu.tw).