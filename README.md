## de Barros et al. (in review) "Quantifying diets using Bayesian inference: comparing juvenile fish diets among restored seascapes"

This is a repository that contains data and code to reproduce all analyses and figures in the manuscript. Each file is as follows:

### Stan Files
- `stan/logitReg.stan`: Stan code for the Bayesian logistic regression model used to estimate PoE
- `stan/betareg.stan`: Stan code for the Bayesian mixed beta regression model used to estimate gut fullness

### R functions
- `functions/plotPostPrior.R`: function to plot priors vs. posteriors for the mixed beta regression model
- `functions/plotPostPriorLogit.R`: function to plot priors vs. posteriors for the logistic regression model
- `functions/plotOntoShifts.R`: function to plot ontogenetic shifts in PoE

### R files
- `analyses/analysesLogitReg.R`: contains the framework to gather data, run, and plot logistic model results
- `analyses/analysesBetaReg.R`: contains the framework to gather data, run, and plot mixed beta regression model results

### Model specification

#### Mixed effects beta regression

Gut fullness as a % body weight is assumed to follow a beta distribution with parameters A and B.
The A and B parameters are defined by a precision term ($\phi$) that is shared across groups and a linear
predictor ($\mu_{[i]}$). The linear predictor contains a GLM whose parameters describe the overall mean
($\beta_{0}$), the effect of site (site-specific random effects) ($\beta_{1[s]}$), and the effect of the total length ($\beta_{TL}$).

$$ GF_{[i]} \sim beta(A,B)  $$

$$ A = \phi*\mu_{i} $$

$$ B = \phi*(1-\mu_{[i]})  $$

$$ \mu_{[i]} = inverse-logit(\beta_{0} + \beta_{1[s]} + \beta_{TL}TL_{[i]})  $$

#### Logistic regression model

Prey Probabilities of Encounter (PoEs) are modelled using a random effects logistic regression model.
We assume $Y_{[i]}$ is a Bernoulli random variable indicating presence (1) or absence (0) of a given prey 
with probability $PoE$. The probability of encounter is modelled as a logistic function of a linear predictor
containing an intercept $\alpha$, site-specific random effects $\beta_{[s]}$, and the effect of total length $\beta_{TL}$.

$$ Y_{[i]} \sim Bernoulli(PoE) $$

$$ log(\frac{PoE}{1-PoE}) = \alpha + \beta_{[s]} + \beta_{TL_{[i]}}  $$













