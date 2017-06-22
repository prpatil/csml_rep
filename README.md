# csml_rep
## Code to replicate analyses from "Cross-Study Machine Learning"

### File contents

- csml_runner.R: Main function to (1) generate dataset list from curatedOvarianData (2) run main simulation to produce Figure 3 of paper + supplemental results
- curovsim.R: Code to produce Figure 2 of main text
- debulksim.R: Code to produce Table 1 of main text
- multistudysim.R: Simulation functions to accompany csml_runner.R
- mutlistudysim_utils.R: Load libraries, learner functions, helper functions
- signatures.Rda: List of published signatures used for Figure 2 of main text

### R libraries needed

RcppArmadillo, e1071, ranger, glmnet, rpart, nnet, nnls, foreach, parallel

### Bioconductor libraries needed:

curatedOvarianData, genefilter

### Comments

Main simulations are run parallelized to save some time. If you'd rather run them in serial, let me know and I can send you modified code.

Start with csml_runner.R to generate the list of datasets, and save that into a .Rda object.

Let me know if you have any issues! I intend to produce a .Rmd with all results in it, but generating this will take some time.
