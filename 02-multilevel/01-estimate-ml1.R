
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### !!! Make sure that you are in the correct working directory! ####
###                     analysis/02-multilevel                   ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## ML1 --------------------------------------------------------------------------------------------------------------------------------------------------
# Approach "ML1": Multilevel model with constraint on population parameter

save_fits = TRUE  #Set TRUE if you'd like to inspect single fits. But will create quite a lot of data

source("00-settings-and-init.R")

res.dir = sprintf("results/ml1/v_%s", ver)
dir.create(res.dir, showWarnings = FALSE, recursive = TRUE)
ML = matrix(nrow = max_k, ncol = 3)

for (j in seq_len(max_k)) {
  frame = SICSFrame_RL_ml1(choice, reward, reset, k, max_k = j,
                           hyperprior_lkj_eta = hyperprior_lkj_eta,
                           hyperprior_mu_alpha = hyperprior_mu_alpha,
                           hyperprior_sigma_alpha = hyperprior_sigma_alpha,
                           hyperprior_mu_beta = hyperprior_mu_beta,
                           hyperprior_sigma_beta = hyperprior_sigma_beta)

  fit = model(frame, constraint = cons[[j]],
              cores = cores,
              iter_posterior=2*n_posterior,
              refresh = 1)
  ML[j,] = logml(fit)
  if (save_fits) saveRDS(fit, file = file.path(res.dir, sprintf("con%i", j)))
  readr::write_csv(data.frame(ML[j,,drop=FALSE]), file.path(res.dir,"logml.csv"), append=TRUE)
}
