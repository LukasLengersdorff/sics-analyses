#############################################################
### Load and prepare example data ###########################
#############################################################


## --------------------------------------------------------------------------------------------------------------------------------------------------
# Packages:
library(tidyverse) # for convenience
library(osfr) # to load data from OSF

###############
#### Data ##### --------------------------------------------------------------------------------------------------------------------------------------------------
###############
#
# Schaaf et al. (2023), https://osf.io/preprints/psyarxiv/chq5a
#
# Data at https://osf.io/k9quw
#
# The dataset contains more information and conditions than we actually use in our analyses.
# We run these analyses primarily to demonstrate our method (SICS). We do not claim that this
# is the "one way" this dataset should be analysed.


## Download data from OSF to working directory (if necessary) ---------------------------------------------------------------------------------------
if (!file.exists("data_RL.xls")) {
  osf_retrieve_file("k9quw") %>% osf_download()
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

data = filter(data, !is.na(hit1)) # Remove missing trials

data = mutate(data, reward = if_else(condition1 == 1, outcome1, outcome1 + 1))
# We pool data from the two conditions "reward learning" and "punishment learning":
#   condition1 == 1: outcome is no reward (0) vs. reward (1);
#   condition1 == 2: outcome is punishment (-1) vs. no punishment (0);
# So `reward` is 0 for the "worse" outcome and 1 for the "better" outcome.
# Note that the used models are agnostic to this coding.

data = select(data, id = pp, trial = newtrial1, pair = pair1, choice = hit1, reward = reward)
# Simplify data set. Note that we completely disregard the data from time point 2 of this data set

data$choice[is.na(data$choice)] = -9 # code for missing values
data$reward[is.na(data$reward)] = -9 # -"-

data$reset = c(1, diff(data$pair)); data$reset = as.numeric(data$reset != 0)
# 1 at first trial of new pair, 0 otherwise. Resets the subjective values in the STAN model code.

data$k0 =  tapply(data$choice == 0, list(data$id,data$pair), cumsum) %>% t %>% c %>% {do.call(c,.)}
# k0: How often was symbol 0 chosen up until and including this trial?
data$k1 =  tapply(data$choice == 1, list(data$id,data$pair), cumsum) %>% t %>% c %>% {do.call(c,.)}
# k1: How often was symbol 1 chosen up until and including this trial?
data = mutate(data, k = if_else(choice == 0, k0, k1))
# k : How often was the chosen symbol chosen up until and including this trial?
