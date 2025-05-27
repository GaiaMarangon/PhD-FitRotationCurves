# PhD-FitRotationCurves

This repository contains the MATLAB code to fit the Schroedinger-Poisson eigenstates to the experimental rotation curves from the [SPARC database](http://astroweb.cwru.edu/SPARC/) ([Lelli et al., 2016](https://iopscience.iop.org/article/10.3847/0004-6256/152/6/157)).  

The repository includes three main files:
- the `main_estimatedSampling` computes the predicted curve using a priori estimates to select the set of parameters to be tested;
- the `main_elaborateResults` analyzes the computed predicted curves to plot the optimal predictions, the predictions as the parameters are varied one-by-one, the $\chi^2$ metric as the parameters are varied one-by-one, and the reference features in the rotation curves used to constrain the free parameters.
- the `main_uniformSampling` computes the predicted curve for a uniform grid of free parameters, centered around a reference value. Usually, the optimal set of parameters obtained from the `main_estimatedSampling` is used as a reference, and the `main_uniformSampling` computations are used as an additional analysis, to construct better plots (in particular that of the $\chi^2$ metric as the parameters are varied one-by-one)


Computing the predicted curves require some input data:
- the rotation curves provided by the [SPARC database](http://astroweb.cwru.edu/SPARC/)
- the density parameters $a_D$, $r_{0,D}$, $a_B$, $r_{0,B}$, $a_G$, $r_{0,G}$ defining the exponential density profiles $\rho_D(r)=a_D e^{r/r_{0,D}}$, $\rho_B(r)=a_B e^{r/r_{0,B}}$, $\rho_G(r)=a_G e^{r/r_{0,G}}$ for disk, bulge and gas. See the [SPARC_analysis.git](https://github.com/GaiaMarangon/SPARC_analysis.git) repository for instructions on how to derive these parameters from SPARC data.

Computed data, both for the `main_estimatedSampling` and for the `main_uniformSampling`, are organized by galaxy and by model type. For each galaxy and each model type, the results include:
- a table `paramTable` in the `paramTable.mat` file, collecting all tested combination of parameters $(Q,n,\xi,m)$. When a new combination needs to be computed, the table is checked to see if the combination has already be tested. Only if it's not we do perform the computations. Once a combination is tested, all the results are stored in a specific file and the `paramTable` in the `paramTable.mat` file is updated.
- a `.mat` file for each tested combination of parameters, including all the results. This file is named after the galaxy, model type and tested combination of free parameters, and it is loaded in the `main_elaborateResults` when the results for that specific combination need to be plotted or elaborated.

Available galaxies include all those for which the input data are available.

Available model types include the following options:
- the Schroedinger-Poisson model defining the eigenstates may or may not include a $4\pi$ factor in the right-hand side of the Poisson equation, when expressed in physical units;
- the density profiles may be of exponential type, of hard-sphere type or of truncated-Plummer type (see the [SPARC_analysis.git](https://github.com/GaiaMarangon/SPARC_analysis.git) repository for detailed definitions).

The `main_elaborateResults` file accesses the results of the `main_estimatedSampling` and of the `main_uniformSampling` with no distinction. It saves the plots in dedicated folders and subfolders.