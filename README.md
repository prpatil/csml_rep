# csml_rep
## Code to replicate analyses from "Training Replicable Predictors in Multiple Studies"

### File s

- csml_runner.R: Main function to (1) generate dataset list from curatedOvarianData (2) run main simulation to produce Figure 2 of main text + supplemental results
- curovsim.R: Code to produce Figure 1 of the supplement
- debulksim.R: Code to produce Table 1 of the supplement
- multistudysim.R: Simulation functions to accompany csml_runner.R
- mutlistudysim_utils.R: Load libraries, learner functions, helper functions
- ovariansurvival.R: Code to produce Figure 3 of the main text
- signatures.Rda: List of published signatures used for Figure 2 of main text
- training-replicable-supplement.pdf: PDF file of supplement to "Training Replicable Predictors in Multiple Studies"

### R libraries needed

RcppArmadillo, e1071, ranger, glmnet, rpart, nnet, nnls, foreach, parallel, doParallel, metafor

### Bioconductor libraries needed:

curatedOvarianData, genefilter

### Other packages used

survHD

### Instructions

First, source multistudysim_utils.R and multistudysim.R. Begin with csml_runner.R, which contains at the top code to generate a list of datasets from curatedOvarianData.
Then you may proceed with the remainder of that file which runs the simulations for Figure 2 in the paper and figures and tables in the supplement.
Next, run ovariansurvival.R to produce Figure 3 of the main text.
curovsim.R and debulksim.R produce Supplemental Figure and Table 1 and rely on the list of datasets created in csml_runner.R.

### Comments

Main simulations are run in parallel to save some time. If you'd rather run them serially, please contact us for modified code.
Please contact us if you have any issues.