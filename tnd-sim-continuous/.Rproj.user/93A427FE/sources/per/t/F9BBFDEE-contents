library(numDeriv, lib.loc = c(.libPaths(), "/home/qijunli/R_libraries"))
library(sandwich, lib.loc = c(.libPaths(), "/home/qijunli/R_libraries"))
library(parallel, lib.loc = c(.libPaths(), "/home/qijunli/R_libraries"))
library(dplyr, lib.loc = c(.libPaths(), "/home/qijunli/R_libraries"))

source("param.R")
source("tnd-nc-continuous.R")
source("tnd-iptw.R")
source("logit_reg.R")

NSIM <- 1000

set.seed(2022)



ii <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

beta0 <- param_grid$beta[ii]
mu_0y <- param_grid$mu_0y[ii]



RES <- mclapply(1:NSIM, function(jj) {
  print(paste(ii, jj))
  
  source("data-gen.R", local = TRUE)  
  
  ## oracle method using the true treatment bridge function
  tnd_oracle <- try(tnd_nc_continuous_oracle(study_data, 
                                             mu_0a, mu_ua, mu_0z, mu_uz, 
                                             sigma_z, mu_xa, mu_xz), silent = TRUE)
  
  
  tnd_nc_gmm <- try(tnd_nc_gmm_continuous(study_data), silent = TRUE)
  
  logit_reg_res <- try(logit_reg(study_data), silent = TRUE)
  iptw_res <- try(tnd_iptw(study_data), silent = TRUE)
  
  
  if (class(tnd_oracle) == "try-error" | 
      class(tnd_nc_gmm) == "try-error" |
      class(logit_reg_res) == "try-error" |
      class(iptw_res) == "try-error") {
    return(c(oracle_bias = NA, 
             oracle_coverage = NA, 
             oracle_se = NA,
             rr_gmm_bias = NA,
             rr_gmm_coverage = NA,
             rr_gmm_se = NA,
             logit_reg_bias = NA,
             logit_reg_coverage = NA,
             logit_reg_se = NA,
             iptw_bias = NA,
             iptw_coverage = NA,
             iptw_se = NA,
             ave_n = NA))
  } else {
    return(c(oracle_bias = tnd_oracle$beta - beta0, 
             oracle_coverage = (beta0 > tnd_oracle$CI[1]) & (beta0 < tnd_oracle$CI[2]), 
             oracle_se = tnd_oracle$SE,
             gmm_bias =  tnd_nc_gmm$beta - beta0,
             gmm_coverage = (beta0 > tnd_nc_gmm$CI[1]) & (beta0 < tnd_nc_gmm$CI[2]),
             gmm_se = tnd_nc_gmm$SE,
             logit_reg_bias = logit_reg_res$beta - beta0,
             logit_reg_coverage = (beta0 > logit_reg_res$CI[1]) & (beta0 < logit_reg_res$CI[2]),
             logit_reg_se = logit_reg_res$SE,
             iptw_bias = iptw_res$beta - beta0,
             iptw_coverage = (beta0 > iptw_res$CI[1] ) & (beta0 < iptw_res$CI[2]),
             iptw_se = iptw_res$SE,
             ave_n = nrow(study_data)))
  } 
}, mc.cores = 5)

RES <- as.data.frame(bind_rows(RES))
save(RES, file = paste0("./results/RES_setting_", ii, ".RData"))



RES <- mclapply(1:3, function(jj) {
  print(paste(ii, jj))
  
  source("data-gen.R", local = TRUE)  
  
  ## oracle method using the true treatment bridge function
  tnd_oracle <- try(tnd_nc_continuous_oracle(study_data, 
                                             mu_0a, mu_ua, mu_0z, mu_uz, 
                                             sigma_z, mu_xa, mu_xz), silent = TRUE)
  
  
  tnd_nc_gmm <- try(tnd_nc_gmm_continuous(study_data), silent = TRUE)
  
  logit_reg_res <- try(logit_reg(study_data), silent = TRUE)
  iptw_res <- try(tnd_iptw(study_data), silent = TRUE)
  
  
  if (class(tnd_oracle) == "try-error" | 
      class(tnd_nc_gmm) == "try-error" |
      class(logit_reg_res) == "try-error" |
      class(iptw_res) == "try-error") {
    return(c(oracle_bias = NA, 
             oracle_coverage = NA, 
             oracle_se = NA,
             rr_gmm_bias = NA,
             rr_gmm_coverage = NA,
             rr_gmm_se = NA,
             logit_reg_bias = NA,
             logit_reg_coverage = NA,
             logit_reg_se = NA,
             iptw_bias = NA,
             iptw_coverage = NA,
             iptw_se = NA,
             ave_n = NA))
  } else {
    return(c(oracle_bias = tnd_oracle$beta - beta0, 
             oracle_coverage = (beta0 > tnd_oracle$CI[1]) & (beta0 < tnd_oracle$CI[2]), 
             oracle_se = tnd_oracle$SE,
             gmm_bias =  tnd_nc_gmm$beta - beta0,
             gmm_coverage = (beta0 > tnd_nc_gmm$CI[1]) & (beta0 < tnd_nc_gmm$CI[2]),
             gmm_se = tnd_nc_gmm$SE,
             logit_reg_bias = logit_reg_res$beta - beta0,
             logit_reg_coverage = (beta0 > logit_reg_res$CI[1]) & (beta0 < logit_reg_res$CI[2]),
             logit_reg_se = logit_reg_res$SE,
             iptw_bias = iptw_res$beta - beta0,
             iptw_coverage = (beta0 > iptw_res$CI[1] ) & (beta0 < iptw_res$CI[2]),
             iptw_se = iptw_res$SE,
             ave_n = nrow(study_data)))
  } 
})
