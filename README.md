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

### Instructions

First, source multistudysim_utils.R and multistudysim.R. Begin with csml_runner.R, which contains at the top code to generate a list of datasets from curatedOvarianData.
Then you may proceed with the remainder of that file which runs the simulations for Figure 3 in the paper and figures and tables in the supplement. 
curovsim.R and debulksim.R produce Figure 2 and Table 1 of the main text, respectively, and rely on the list of datasets created in csml_runner.R.

### Comments

Main simulations are run in parallel to save some time. If you'd rather run them serially, let me know and I can send you modified code.

Please let me know if you have any issues. I intend to produce a .Rmd with all results in it, but generating this will take some time.