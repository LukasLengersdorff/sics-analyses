
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##+++++++++++ SINGLE SUBJECT ANALYSES +++++++++++++++++++++++++++++++++++++++++
##+++++++++++  conditional EP Method  +++++++++++++++++++++++++++++++++++++++++
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### !!! Make sure that you are in the correct working directory! ####
###                     analysis/01-single-subject               ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Packages:
library(tidyverse) # for convenience
library(osfr)
library(sics)

###############
#### Data ##### --------------------------------------------------------------------------------------------------------------------------------------------------
###############

source("../get_data.R") #-> loads "data"

#++++++++++++++++++++++++++
#+++ Single participant +++ -------------------------------------------------------------------------------------------------------------
#++++++++++++++++++++++++++

# Settings ------------------------------------------------------------------------------------------------------------------------------

subject = 1;
max_k = 8;
# Maximum dimension of learning rate parameter
# max_k = 8 means that we have hypotheses
# H1: a1
# H2: a1 > a2
# H3: a1 > a2 > a3
# ...
# H8: a1 > a2 > a3 > ... > a8

R = 100 # Number of repetitions of each analysis, to investigate precision and efficiency
N = 30000 # Number of posterior samples to be taken in total.
# Exact number of posterior samples per estimated quantity varies with the method and contrast, to increase comparability

set.seed(2)

res.dir = "results/cem"
dir.create(res.dir, showWarnings = FALSE, recursive = TRUE)
res = list()
data.used = filter(data, id == subject)

### Unconditional encompassing method ------------------------------------------------------------------------------------------------

## Sub-settings for CEM
# Prior spec
exp_lambda = 1
beta_shape1 = 1
beta_shape2 = 1

# MCMC
chains = 4
cores = chains

res$cem = data.frame(bf = numeric(R),
                     pC = numeric(R),
                     pCD = numeric(R),
                     pD = numeric(R),
                     t.user = numeric(R),
                     t.system = numeric(R),
                     t.elapsed = numeric(R))

for (j in 2:max_k) { # Level of constraints/hypotheses
  write_csv(res$cem[NULL,], file.path(res.dir,sprintf("con%i.csv", j)), append=FALSE)
  samples.per.chain = ceiling(N/(2*(j-1)*chains))
  frame = with(data.used,
               SICSFrame_RL(
                 choice,
                 reward,
                 reset,
                 k,
                 max_k = j,
                 exp_lambda = exp_lambda,
                 beta_shape1 = beta_shape1,
                 beta_shape2 = beta_shape2
               ))
  for (r in seq_len(R)) { # Level of repetitions
    t0 = proc.time()
    pCi = numeric(j-1)
    pCDi = numeric(j-1)

    for (m in 1:(j-1)) { # Level of partial bayes factors
      fit = model(frame,
                  constraint = if (m > 1) {EmptyConstraint(j) + OrdinalConstraint(1:m)} else {EmptyConstraint(j)},
                  chains_posterior = chains,
                  chains_prior = chains,
                  iter_posterior = samples.per.chain*2, # Warmup takes half of the samples by rstan default
                  iter_prior = samples.per.chain*2, # Warmup takes half of the samples by rstan default
                  force_prior_sampling = TRUE, # For cEM, we actually want the prior of an empty constraint! :-)
                  cores = cores)
      # Prior
      a = rstan::extract(fit$fit0, "alpha_t")[[1]]
      pCi[m] = mean(a[,m] > a[,m+1])

      # Posterior
      a = rstan::extract(fit$fit1, "alpha_t")[[1]]
      pCDi[m] = mean(a[,m] > a[,m+1])

      # Compute p(D) if m = 1 (encompassing model)
      if (m == 1) {
        pD = logml(fit)[1]
      }
    }
    # P(C)
    pC = prod(pCi)
    # P(C|D)
    pCD = prod(pCDi)
    # BF
    BF = pCD/pC
    # Time
    t = summary(proc.time()-t0)

    res$cem[r,] = c(BF, pC, pCD, pD, t)
    write_csv(res$cem[r,,drop=FALSE], file.path(res.dir,sprintf("con%i.csv", j)), append=TRUE)
  }
}
