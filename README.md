
<p align="center">
<img src="images/logo.png" />
</p>

---

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8251777.svg)](https://doi.org/10.5281/zenodo.8251777) [![CI](https://github.com/brian-sinquin/FiberNlse.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/brian-sinquin/FiberNlse.jl/actions/workflows/CI.yml) [![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/N/NamedDims.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html) [![codecov](https://codecov.io/gh/curio-sitas/FiberNlse.jl/branch/main/graph/badge.svg?token=O9L9P064J1)](https://codecov.io/gh/curio-sitas/FiberNlse.jl) [![code style blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle) [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

Documentation : [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://curio-sitas.github.io/FiberNlse.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://curio-sitas.github.io/FiberNlse.jl/dev)

---

 A Generalised Non-Linear Schrödinger Equation (GNLSE) solver for julia aimed towards fiber optics or other waveguides. 

## Features

The *FiberNlse.jl* package simulates the propagation of an optical field envelope signal in a waveguide showing *Kerr* nonlinearity, group velocity dispersion, Raman scattering and linear losses.

The algorithm used is ERK4(5)-IP with adaptative step-size. The **API** permits to include the Raman scattering or/and self-steepening in the integration process.

---

## Installation
To install you can simply type :
```
] add FiberNlse
```
in your julia terminal or clone this repository and include the `src/FiberNlse.jl` file in your project.

## ⚙️ Current development state

- [x] Implement ERK4(5)-IP integration
- [x] Register DOI
- [ ] Document code
- [ ] New Logo / Visuals
- [ ] Benchmark / Comparison with other language packages (Python/Matlab/C)
- [ ] Documentation + sources
- [ ] Examples (NL Compression, MI, Super-Continuum, Flaticons, Similaritons, Solitons, Dispersion spreading)
- [x] Setup continuous integration
- [ ] Add progress bar option
- [x] higher order dispersion
- [ ] Add variable dispersion, loss and nonlinearity
- [x] Add Self-steepening
- [x] Add Raman scattering
- [ ] Utils and helper functions for visualization and field manipulation
- [ ] Phase sign convention

     
## 💡 Future ideas
- [ ] Many solvers (as in [pychi](https://github.com/pychi-code/pychi/tree/main))
- [ ] Make it compatible with **AD** for optimisation or inverse problems      

## 🤝 Contributions

Any help to optimize, refactor or enrich the package is welcome, whatever it is on the core side or practical examples.
Don't hesitate to open an issue if something goes wrong with the package.

## 📝 Citation

Please cite this repository if you use it to publish data in a research paper.


```
@software{sinquin_fibernlse,
  author       = {Sinquin Brian},
  title        = {FiberNlse.jl},
  month        = mar,
  year         = 2022,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.8251777},
  url          = {https://doi.org/10.5281/zenodo.8251777}
}
```
