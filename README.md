# SmoQyDEAC.jl Paper Scripts
Repository of code used to generate plots for SmoQyDEAC.jl's paper [arXiv:TBD](https://arxiv.org/abs/FIXME). This is a reimplementation of the DEAC algorithm developed by Nathan S. Nichols, Paul Sokol, and Adrian Del Maestro [arXiv:2201.04155](https://arxiv.org/abs/2201.04155).

Other software needed for these scripts:
- [`ACFlow`](https://github.com/huangli712/ACFlow): Li Huang's Analytic Continuation package
- [`SAC`](https://github.com/JefferyWangSH/sac): Jeffrey Wang's implementation of Stochastic Analytic Continuation
- [`SynthAC`]( https://github.com/sandimas/SynthAC.jl): Package to convert spectral functions to Matsubara space 

For questions and comments regarding these scripts, please email James Neuhaus at [jneuhau1@utk.edu](mailto:jneuhau1@utk.edu).
# Funding
This work was supported by the U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences, under Award Number DE-SC0022311. N.S.N. was supported by the Argonne Leadership Computing Facility, which is a U.S. Department of Energy Office of Science User Facility operated under contract DE-AC02-06CH11357. 

# Comparison Plots
Our paper included a few representative cases for usage of SmoQyDEAC.jl. All 10 tested distributions are accounted for here, for completeness. See "1_MakeDistributions.jl" in comparisons subfolders for true distributions tested.
## Additional Fermion Plots
![Single Gaussian](img/f1.png)
![Single Lorentzian](img/f2.png)
![Double Gaussian](img/f3.png)
![Double Lorentzian](img/f4.png)
![Lorentzian and Gaussian](img/f5.png)
![Double Plateau](img/f6.png)
![Triple Peak 1](img/f7.png)
![Triple Peak 2](img/f8.png)
![Chaos](img/f9.png)
![Batman](img/f10.png)

## Additional Boson (symmetric) Plots
![Gaussian at 1.0](img/b1.png)
![Gaussian at 0.1](img/b2.png)
![Gaussian at 0.01](img/b3.png)
![Lorentzian at 1.0](img/b4.png)
![Double Gaussian 1](img/b5.png)
![Double Gaussian 2](img/b6.png)
![Gaussian and Singularity 1](img/b7.png)
![Gaussian and Singularity 2](img/b8.png)
![Chaos](img/b9.png)
![Chaos with Singularity](img/b10.png)
