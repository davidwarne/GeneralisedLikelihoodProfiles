# Generalised Likelihood Profiles  

This repository contains reference MATLAB functions and scripts to demonstrate the generalised likelihood profile approach for frequentist style parameter estimation when the likelihood is intractable. 

## Developers

David Warne$^{1,2}$ (david.warne@qut.edu.au), https://scholar.google.com.au/citations?user=t8l-kuoAAAAJ&hl=en
Christopher Drovandi$^{1,2}$
Elliot Carr$^{1}$
Matthew Simpson$^{1}$

1. School of Mathematical Sciences, Faculty of Science, Queensland Univeristy of Technology, Australia
2. Centre for Data Science, Queensland University of Technology, Australia

## Citation Information

This code is provided as supplementary information to the paper,

David J Warne, Oliver J. Maclaren, Elliot J. Carr, Matthew J. Simpson, and Christopher Drovandi. Generalised likelihood profiles for models with intractable likelihoods. ArXiv preprint (TBA) 

## Licensing
This source code is licensed under the GNU General Public License Version 3.
Copyright (C) 2023 David J. Warne

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Contents

```bash
The directory structure is as follows
|-- init.m                                  Adds all functions to the MATLAB Path
|-- LossFuctions/
    |-- discrete_fisher_divergence.m        
    |-- moments_distance.m
|-- GeneralisedLikelihood/
    |-- generalisedLikelihoodProfiles.m    Main algorithm to calibrate profiles
    |-- compute_coverage.m
|-- Examples/
    |-- ConwayMaxwellPoisson/
        |-- glp_cmp_model.m                Run Conway-Maxwell-Poisson example
        |-- com_rnd.m
        |-- com_pdf_unnorm.m
        |-- com_pdf.m
        |-- com_logpdf.m
        |-- com_cdf.m
    |-- BiasedStochasticDiffusion/
        |-- glp_biased_diffusion.m        Run biased stochastic transport example
        |-- bootstrapMoments.m
        |-- Stochastic_Model.m
        |-- Exact_Moments.m
        |-- Exact_Moments_x0.m
        |-- Numerical_Moments.m
        |-- Numerical_Moments_x0.m
```
## Usage

1. Start MATLAB
2. In MATLAB browse to the repository folder GeneralisedLikelihoodProfiles/
3. In the MATLAB command prompt, run 
   `>> init` 
   to set up the paths of all the example implementations.
4. Run an example either  
   `>> glp_cmp_model` to generate profiles like Fig 1 for the Conway-Maxwell-Poisson model.
   or
   `>> glp_biased_diffusion` to generate Figs 3 and 4 for th biased stochastic transport model.

5. To obtain coverage plots set `validateTF = 1` at line  19 in `glp_biased_diffusion.m`

