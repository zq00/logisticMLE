This repository contains raw code that reproduces results in the paper "The Asymptotic Distribution of the MLE in High-dimensional Logistic Models: Arbitrary Covariance".

Each directory contains code and output in the corresponding section. Folders in the directories are organized as follows:
- **code** : R code to generate parameters, simulate data, estimate parameters and fit logistic regressions
- **param**: raw parameters 
- **output**: raw output. Sometimes the outputs are splitted to conform to size restriction.
- **results**: R code to analyze results and generate table/figure.

A r package that streamlines procedures described in the paper to adjust high-dimensional logistic MLE is still under construction. In the meantime, if you are interested in the phenomenon, you can look at the example data analysis in *real* folder. We will provide a quick demonstration with simulated data here soon. 
