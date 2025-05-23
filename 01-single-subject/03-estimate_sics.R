
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##+++++++++++ SINGLE SUBJECT ANALYSES +++++++++++++++++++++++++++++++++++++++++
##+++++++++++         SICS            +++++++++++++++++++++++++++++++++++++++++
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

set.seed(3)

res.dir = "results/sics"
dir.create(res.dir, showWarnings = FALSE, recursive = TRUE)

data.used = filter(data, id == subject)

### SICS ------------------------------------------------------------------------------------------------

## Sub-settings
# Prior spec
exp_lambda = 1
beta_shape1 = 1
beta_shape2 = 1

# MCMC
chains = 4
cores = chains

samples.per.chain = ceiling(N/(2*chains))

res = data.frame(pDC = numeric(R),
                 t.user = numeric(R),
                 t.system = numeric(R),
                 t.elapsed = numeric(R))

for (j in 1:max_k) { # Level of constraints/hypotheses
  write_csv(res[NULL,], file.path(res.dir,sprintf("con%i.csv", j)), append=FALSE)
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

    # P(D|C)
    fit = model(frame,
                constraint = if (j > 1) OrdinalConstraint(1:j) else EmptyConstraint(1),
                chains_posterior = chains,
                chains_prior = chains,
                iter_posterior = samples.per.chain*2, # Warmup takes half of the samples by rstan default
                iter_prior = samples.per.chain*2, # Warmup takes half of the samples by rstan default
                #force_prior_sampling = TRUE
                cores = cores)
    pDC = logml(fit)[1]

    # Time
    t = summary(proc.time()-t0)

    res[r,] = c(pDC, t)
    write_csv(res[r,,drop=FALSE], file.path(res.dir,sprintf("con%i.csv", j)), append=TRUE)
  }
}
