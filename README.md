# snorer


`snorer` is a package for evaluating the signatures of **S**upernova-**N**eutrino-b**O**osted da**R**k matt**ER** from our Milky Way based on *Phys. Rev. Lett.* **130**, 111002 (2023) and *Phys. Rev. D* **108**, 083013 (2023).

## Installation

To install, excute the following command on the prompt

    $ pip install snorer

and everything should be processed on-the-fly.

### Dependency

`dukes` requires these external packages

- `numpy` >= 1.20.0
- `scipy` >= 1.10.0
- `vegas` >= 6.0.1

where `vegas` is a the backend engine for evaluating multidimensional integrals based on adaptive Monte Carlo vegas algorithm, see its homepage: [https://pypi.org/project/vegas/](https://pypi.org/project/vegas/).

Other packages, e.g. `gvar`, maybe required by these dependencies during the installation.
The versions of these dependencies are not strict, but are recommended to update to the latest ones to avoid incompatibility. 



## Misc

Bug report and troubleshooting please contact the author Yen-Hsun Lin via [yenhsun@phys.ncku.edu.tw](mailto:yenhsun@phys.ncku.edu.tw).