This repository contains raw code that reproduces results in the paper "The Asymptotic Distribution of the MLE in High-dimensional Logistic Models: Arbitrary Covariance".

Each directory contains code and output in the corresponding section. Folders in the directories are organized as follows:
- **code** : R code to generate parameters, simulate data, estimate parameters and fit logistic regressions
- **param**: raw parameters 
- **output**: raw output. Sometimes the outputs are splitted to conform to size restriction.
- **results**: R code to analyze results and generate table/figure.

If you are interested in the phenomenon in high-dimensional logistic MLE, you can take a look at the quick demonstration *demo.Rmd* here, or check out the example data analysis in *real* folder. You can also install the R package implementing methodology in this paper with

```R
install.packages("devtools")
devtools::install_github("zq00/glmhd")
```

