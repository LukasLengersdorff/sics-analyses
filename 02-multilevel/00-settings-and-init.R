
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##+++++++++++ SETTINGS FOR MULTILEVEL ANALYSIS ++++++++++++++++++++++++++++++++
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### !!! Make sure that you are in the correct working directory! ####
###                     analysis/02-multilevel                   ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Settings --------------------------------------------------------------------------------------------------------------------------------------------------

ver = "1"     # Version string - if you'd like to try different things, like MCMC sample sizes, or a different set of participants

subs = 1:20
# Subject IDs included in analysis
# This is an ARBITRARY choice! Again, just demonstrating things here...

max_k = 8
# Maximum dimension of learning rate parameter
# max_k = 8 means that we have hypotheses
# H1: a1
# H2: a1 > a2
# H3: a1 > a2 > a3
# ...
# H8: a1 > a2 > a3 > ... > a8
#
# Note the following non-trivial decision: In the analyses, we do not
# actually define a max_k - dimensional parameter, and then use equality constraints;
# Rather, we define a 1-dimensional parameter for H1; a 2-dimensional parameter with
# (a1 > a2) for H2; and so on.
# This does not make a great difference for single level estimations. But using the first
# approach has unwanted consequences in multilevel estimation:
#   - In approach ML1 (constraint on population parameter), every participant gets max_k parameters
#     that are not constrained other than being sampled from the constrained population parameters.
#     But the (in our eyes) sensible thing to do is to also constrain ai = aj when mu_ai = mu_aj.
#   - In approach ML2 (constraint on every single participant parameter), population parameters for
#     max_k-dimensional participant parameters are sampled for every hypothesis, although effectively,
#     the participant parameters are only i-dimensional for hypothesis Hi. Conceptually, this is not a problem,
#     but it plays VERY badly with sampling (especially the correlation matrix parameter doesn't like this).

## STAN settings
cores = 4

## Hyperpriors
hyperprior_lkj_eta = 1
hyperprior_mu_alpha = c(0,0.5)
hyperprior_sigma_alpha = c(0,0.5)
hyperprior_mu_beta = c(0,0.5)
hyperprior_sigma_beta = c(0,0.5)

# Number of posterior samples
n_posterior = 2000

## --------------------------------------------------------------------------------------------------------------------------------------------------
# Packages:
library(dplyr) # for convenience
#library(osfr)
library(sics)

# ++++++++++++++++++++
# ++++++ Data: +++++++
# ++++++++++++++++++++
#
# Schaaf et al. (2023), https://osf.io/preprints/psyarxiv/chq5a
#
# Data at https://osf.io/k9quw
#
# The dataset contains more information and conditions than we actually use in our analyses.
# We run these analyses primarily to demonstrate our method (SICS). We do not claim that this
# is the "one way" this dataset should be analysed.



## Download data from OSF to working directory (if necessary)
if (!file.exists("data_RL.xls")) {
  osfr::osf_retrieve_file("k9quw") %>% osfr::osf_download()
}

# File is saved as "data_RL.xls" in working directory
# It's not in valid XLS, however, just tab delimited
data0 = readr::read_delim("data_RL.xls", delim = "\t", col_names = TRUE)

data = data0[order(data0$pp, data0$pair1),]
# Reorders data:
# In this RL task, in each block, two symbol pairs (identified by pair1) were presented
# in an interspersed manner. The line above reorders the data in such a way that all the
# trials for one data pair are grouped together. The temporal order is preserved, of course!
# We do this for simplicity, because the analyses conducted here do not take these aspects
# of the data into account (no time effects across blocks, no effects of pairs changing within a block,
# or similar). For analyses that take these aspects into account, the data would need to be
# organized in a different way.

data = dplyr::filter(data, !is.na(hit1)) # Remove missing trials

data = dplyr::mutate(data, reward = dplyr::if_else(condition1 == 1, outcome1, outcome1 + 1))
# We pool data from the two conditions "reward learning" and "punishment learning":
#   condition1 == 1: outcome is no reward (0) vs. reward (1);
#   condition1 == 2: outcome is punishment (-1) vs. no punishment (0);
# So `reward` is 0 for the "worse" outcome and 1 for the "better" outcome.
# Note that the used models are agnostic to this coding.

data = dplyr::select(data, id = pp, trial = newtrial1, pair = pair1, choice = hit1, reward = reward)
# Simplify data set. Note that we completely disregard the data from time point 2 of this data set

data$choice[is.na(data$choice)] = -9 # code for missing values
data$reward[is.na(data$reward)] = -9 # -"-

data$reset = c(1, diff(data$pair)); data$reset = as.numeric(data$reset != 0)
# 1 at first trial of new pair, 0 otherwise. Resets the subjective values in the model code.

data$k0 =  tapply(data$choice == 0, list(data$id,data$pair), cumsum) %>% t %>% c %>% {do.call(c,.)}
# k0: How often was symbol 0 chosen up until and including this trial?
data$k1 =  tapply(data$choice == 1, list(data$id,data$pair), cumsum) %>% t %>% c %>% {do.call(c,.)}
# k1: How often was symbol 1 chosen up until and including this trial?
data = dplyr::mutate(data, k = dplyr::if_else(choice == 0, k0, k1))
# k : How often was the chosen symbol chosen up until and including this trial?


## Constraints ------------------------------------------------------------------------------------------------------------------------------------------------

# Constraints:
# The constraints are the same for all approaches, so we can create them once at the beginning

cons = vector(mode = "list", length = max_k)
for (k in seq_len(max_k)) {
  if (k == 1) {
    con = EmptyConstraint(1)
  } else {
    con = OrdinalConstraint(1:k, n = k)
  }
  cons[[k]] = con
}

## Data wrangling ---------------------------------------------------------------------------------------------------------------------------------------------

.data = dplyr::filter(data, id %in% subs)
max_T = with(.data, tapply(id, id, length)) %>% max

.datalist = dplyr::group_by(.data, id) %>% group_split

pad = function(v, max_T, pad_value = -9) {
  out = rep(pad_value, max_T)
  out[1:length(v)] = v
  out
}

choice = do.call(cbind,
                 lapply(1:length(subs), function(i) pad(.datalist[[i]]$choice, max_T, 0)))
reward = do.call(cbind,
                 lapply(1:length(subs), function(i) pad(.datalist[[i]]$reward, max_T)))
reset = do.call(cbind,
                lapply(1:length(subs), function(i) pad(.datalist[[i]]$reset, max_T, 0)))
k = do.call(cbind,
            lapply(1:length(subs), function(i) pad(.datalist[[i]]$k, max_T, 1)))

## Create directory

dir.create("results", showWarnings = FALSE)
