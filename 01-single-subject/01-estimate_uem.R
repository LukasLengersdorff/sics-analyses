
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##+++++++++++ SINGLE SUBJECT ANALYSES +++++++++++++++++++++++++++++++++++++++++
##+++++++++++ unconditional EP Method +++++++++++++++++++++++++++++++++++++++++
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
N = 30000 # Number of samples drawn per estimated term

set.seed(1)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
res = list()
data.used = filter(data, id == subject)

### Unconditional encompassing method ------------------------------------------------------------------------------------------------

## Sub-settings for UEM
# Prior spec
exp_lambda = 1
beta_shape1 = 1
beta_shape2 = 1

# MCMC
chains = 4
cores = chains

samples.per.chain = ceiling(N/(2*chains))

res.dir = "results/uem"
dir.create(res.dir, showWarnings = FALSE, recursive = TRUE)
res$uem = data.frame(bf = numeric(R),
                        pC = numeric(R),
                        pCD = numeric(R),
                        pD = numeric(R),
                        t.user = numeric(R),
                        t.system = numeric(R),
                        t.elapsed = numeric(R))

### For convenience, we also use SICS functions for sampling from the priors and posteriors, with only "empty" constraints

for (j in 2:max_k) { # Level of constraints/hypotheses
  write_csv(res$uem[NULL,], file.path(res.dir,sprintf("con%i.csv", j)), append=FALSE)
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
  con = EmptyConstraint(j)
  for (r in seq_len(R)) { # Level of repetitions
    t0 = proc.time()
    fit = model(frame,
                constraint = con,
                chains_posterior = chains,
                chains_prior = chains,
                iter_posterior = samples.per.chain*2, # Warmup takes half of the samples by rstan default
                iter_prior = samples.per.chain*2, # Warmup takes half of the samples by rstan default
                force_prior_sampling = TRUE, # For UEM, we actually want the prior of an empty constraint! :-)
                cores = cores)
    # P(C)
    a = rstan::extract(fit$fit0, "alpha_t")[[1]]
    boolmat = sapply(1:(j-1), function(i) {a[,i] > a[,i+1]})
    # This is a [samples] x [constraints] matrix
    pC = apply(boolmat, 1, all) %>% mean

    # P(C|D)
    a = rstan::extract(fit$fit1, "alpha_t")[[1]]
    boolmat = sapply(1:(j-1), function(i) {a[,i] > a[,i+1]})
    # This is a [samples] x [constraints] matrix
    pCD = apply(boolmat, 1, all) %>% mean

    # BF
    BF = pCD/pC

    # P(D)
    pD = logml(fit)[1]

    # Time
    t = summary(proc.time()-t0)

    res$uem[r,] = c(BF, pC, pCD, pD, t)
    write_csv(res$uem[r,,drop=FALSE], file.path(res.dir,sprintf("con%i.csv", j)), append=TRUE)
  }
}
