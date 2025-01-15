## de Barros et al. (in review) "Quantifying diets using Bayesian inference: comparing juvenile fish diets among restored seascapes"

This is a repository that contains data and code to reproduce all analyses and figures in the manuscript. Each file is as follows:

### Folders
- `data/`: contains all data used in the analyses
- `code/`: contains all code used to run the analyses

### Stan Files
- `stan/logitReg.stan`: Stan code for the Bayesian logistic regression model used to estimate PoE
- `stan/betareg.stan`: Stan code for the Bayesian mixed beta regression model used to estimate gut fullness

### R functions
- `plotPostPrior.R`: function to plot priors vs. posteriors for the mixed beta regression model
- `plotPostPriorLogit.R`: function to plot priors vs. posteriors for the logistic regression model
- `plotOntoShifts.R`: function to plot ontogenetic shifts in PoE

### R files
- `analysesLogitReg.R`: contains the framework to gather data, run, and plot logistic model results
- `analysesBetaReg.R`: contains the framework to gather data, run, and plot mixed beta regression model results



