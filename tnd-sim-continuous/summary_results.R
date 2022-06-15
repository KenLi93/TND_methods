source("param.R")
library(dplyr)
library(xtable)
library(ggplot2)
library(gridExtra)
library(ggpubr)

NSIM <- 500

results_panel <- data.frame(ave_n = NA,
                            oracle_bias = NA, oracle_SE = NA, oracle_mean_SE = NA, oracle_coverage = NA, oracle_power = NA,
                            gmm_bias = NA, gmm_SE = NA, gmm_mean_SE = NA, gmm_coverage = NA, gmm_power = NA,
                            logit_reg_bias = NA, logit_reg_SE = NA, logit_reg_mean_SE = NA, logit_reg_coverage = NA, logit_reg_power = NA,
                            iptw_bias = NA, iptw_SE = NA, iptw_mean_SE = NA, iptw_coverage = NA, iptw_power = NA,
                            na = NA)

results_grid <- cbind(param_grid, results_panel)

all_RES <- data.frame(beta0 = NULL, mu_0y = NULL, bias = NULL, se = NULL, method = NULL)
## tabulate the results
for (ii in 1:nrow(results_grid)) {
  load(paste0("results/RES_setting_", ii, ".RData"))
  RES <- as.data.frame(t(RES))
  names(RES)[7] <- "logit_reg_bias"
  beta0 <- results_grid$beta[ii]
  
  results_grid$oracle_power[ii] <- mean(pchisq((beta0 + RES$oracle_bias) ^ 2 / RES$oracle_se ^ 2, df = 1, lower.tail = F) < 0.05)
  results_grid$gmm_power[ii] <- mean(pchisq((beta0 + RES$gmm_bias) ^ 2 / RES$gmm_se ^ 2, df = 1, lower.tail = F) < 0.05)
  results_grid$logit_reg_power[ii] <- mean(pchisq((beta0 + RES$logit_reg_bias) ^ 2 / RES$logit_reg_se ^ 2, df = 1, lower.tail = F) < 0.05)
  results_grid$iptw_power[ii] <- mean(pchisq((beta0 + RES$iptw_bias) ^ 2 / RES$iptw_se ^ 2, df = 1, lower.tail = F) < 0.05)
  
  
  results_grid$na[ii] <- NSIM - nrow(RES)
  
  results_grid$ave_n[ii] <- mean(RES$ave_n, na.rm = TRUE)
  results_grid$oracle_bias[ii] <- mean(RES$oracle_bias, na.rm = TRUE)
  results_grid$oracle_SE[ii] <- sd(RES$oracle_bias, na.rm = TRUE)
  results_grid$oracle_mean_SE[ii] <- median(RES$oracle_se, na.rm = TRUE)
  results_grid$oracle_coverage[ii] <- mean(RES$oracle_coverage, na.rm = TRUE)
  
  
  results_grid$gmm_bias[ii] <- mean(RES$gmm_bias, na.rm = TRUE)
  results_grid$gmm_SE[ii] <- sd(RES$gmm_bias, na.rm = TRUE)
  results_grid$gmm_mean_SE[ii] <- mean(RES$gmm_se, na.rm = TRUE)
  results_grid$gmm_coverage[ii] <- mean(RES$gmm_coverage, na.rm = TRUE)
  
  results_grid$logit_reg_bias[ii] <- mean(RES$logit_reg_bias, na.rm = TRUE)
  results_grid$logit_reg_SE[ii] <- sd(RES$logit_reg_bias, na.rm = TRUE)
  results_grid$logit_reg_mean_SE[ii] <- mean(RES$logit_reg_se, na.rm = TRUE)
  results_grid$logit_reg_coverage[ii] <- mean(RES$logit_reg_coverage, na.rm = TRUE)
  
  results_grid$iptw_bias[ii] <- mean(RES$iptw_bias, na.rm = TRUE)
  results_grid$iptw_SE[ii] <- sd(RES$iptw_bias, na.rm = TRUE)
  results_grid$iptw_mean_SE[ii] <- mean(RES$iptw_se, na.rm = TRUE)
  results_grid$iptw_coverage[ii] <- mean(RES$iptw_coverage, na.rm = TRUE)
  
  
  add_df <- data.frame(beta0 = param_grid$beta[ii],
                       mu_0y = param_grid$mu_0y[ii],
                       bias = with(RES, c(gmm_bias, oracle_bias, logit_reg_bias, iptw_bias)),
                       se = with(RES, c(gmm_se, oracle_se, logit_reg_se, iptw_se)),
                       method = rep(c("NC", "NC-Oracle", "Logit Reg.", "IPTW"), each = nrow(RES)))
  
  all_RES <- rbind(all_RES, add_df)
}

results_tab <- results_grid[, c("beta", "mu_0y", "ave_n",
                                "gmm_bias", "gmm_SE", "gmm_mean_SE", "gmm_coverage", "gmm_power",
                                "oracle_bias", "oracle_SE", "oracle_mean_SE", "oracle_coverage", "oracle_power",
                                "logit_reg_bias", "logit_reg_SE", "logit_reg_mean_SE", "logit_reg_coverage", "logit_reg_power",
                                "iptw_bias", "iptw_SE", "iptw_mean_SE", "iptw_coverage", "iptw_power")] %>%
  arrange(beta, mu_0y)

results_xtab <- xtable(results_tab, digits = c(0, 3, 3,  0, rep(3, 20)))
print.xtable(results_xtab, include.rownames = F, file = "results_tab.txt")

results_tab_inf <- data.frame(beta0 = rep(round(results_tab$beta, 3), 4),
                              mu_0y = rep(round(results_tab$mu_0y, 3), 4),
                              bias = with(results_tab, c(gmm_bias, oracle_bias, logit_reg_bias, iptw_bias)),
                              true_se = with(results_tab, c(gmm_SE, oracle_SE, logit_reg_SE, iptw_SE)),
                              coverage = with(results_tab, c(gmm_coverage, oracle_coverage, logit_reg_coverage, iptw_coverage)),
                              power = with(results_tab, c(gmm_power, oracle_power, logit_reg_power, iptw_power)),
                              method = rep(c("NC", "NC-Oracle", "Logit Reg.", "IPTW"), each = nrow(results_tab)))
results_tab_inf$method <- factor(results_tab_inf$method, levels = c("NC", "NC-Oracle", "Logit Reg.", "IPTW"),
                                 labels = c("NC", "NC-Oracle", "Logit Reg.", "IPTW"))

## plot the bias
for (ee in c(log(0.01), log(0.1), log(0.2))) {
  plot_dat <- all_RES %>% filter(mu_0y == ee) 
  
  plot_dat$beta0 <- round(plot_dat$beta0, 3)
  
  plot_dat$method <- factor(plot_dat$method, levels = c("NC", "NC-Oracle", "Logit Reg.", "IPTW"),
                            labels = c("NC", "NC-Oracle", "Logit Reg.", "IPTW"))
  
  gg_bias <- ggplot(data = plot_dat, aes(x = method, y = bias, color = method)) +
    geom_boxplot(lwd = 1, width = 0.7) +
    geom_hline(yintercept = 0, color = "black", lty = 2) +
    facet_wrap(~beta0, nrow = 1,
               labeller = label_bquote(beta[0] == .(beta0))) +
    theme_pubr() +
    ylab("Bias") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_line(),
          strip.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = "bottom") 
  ggsave(filename = paste0("plots/tnd_continuous_bias_mu_0y_", round(ee * 100),"_.png"), plot = gg_bias,
         width = 6, height = 4)
  
  ## plot SE
  plot_dat_inf <- results_tab_inf %>% filter(mu_0y == round(ee, 3)) 
  
  gg_se <- ggplot(data = plot_dat, aes(x = method, y = se, color = method)) +
    geom_boxplot(lwd = 1, width = 0.7) +
    geom_point(data = plot_dat_inf, aes(x = method, y = true_se), pch = 17, color = "black", size = 3) +
    facet_wrap(~beta0, nrow = 1,
               labeller = label_bquote(beta[0] == .(beta0))) +
    theme_pubr() +
    ylab("SE") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_line(),
          strip.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = "bottom")
  ggsave(filename = paste0("plots/tnd_continuous_se_mu_0y_", round(ee * 100),"_.png"), plot = gg_se,
         width = 6, height = 4)
  
  ## plot bar chart of coverage probability
  gg_coverage <- ggplot(data = plot_dat_inf, aes(x = method, y = coverage, fill = method)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.95, color = "black", lty = 2) +
    facet_wrap(~beta0, nrow = 1,
               labeller = label_bquote(beta[0] == .(beta0))) +
    theme_pubr() +
    ylab("Coverage rate") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_line(),
          strip.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = "bottom")
  ggsave(filename = paste0("plots/tnd_continuous_coverage_mu_0y_", round(ee * 100),"_.png"), plot = gg_coverage,
         width = 6, height = 4)
  
  
  
  ## plot bar chart of coverage probability
  gg_power <- ggplot(data = plot_dat_inf, aes(x = method, y = power, fill = method)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.05, color = "black", lty = 2) +
    facet_wrap(~beta0, nrow = 1,
               labeller = label_bquote(beta[0] == .(beta0))) +
    theme_pubr() +
    ylab("Power") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_line(),
          strip.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = "bottom")
  ggsave(filename = paste0("plots/tnd_continuous_power_mu_0y_", round(ee * 100),"_.png"), plot = gg_power,
         width = 6, height = 4)
}



for (bb in round(c(log(0.2), 0), 3)) {
  plot_dat <- results_tab_inf %>% filter(beta0 == bb) %>% mutate(delta = exp(mu_0y)) %>%
    filter(delta < 0.55)
    
  
  gg_bias_tr <- ggplot(data = plot_dat, aes(x = delta, y = bias, color = method)) +
    geom_point(shape = 4, size = 2, stroke = 1.5) +
    geom_line(lwd = 1.2) +
    #ylim(c(-0.6, 0.2)) +
    xlab(bquote(delta ~ "=exp(" * eta[0] * ")")) + 
    ylab("Bias") +
    theme_pubr() +
    theme(axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_line(),
          strip.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = "bottom")
  ggsave(filename = paste0("plots/tnd_continuous_bias_tr_beta0_", round(bb * 100),"_.png"), plot = gg_bias_tr,
         width = 6, height = 5)
  
  
  gg_coverage_tr <- ggplot(data = plot_dat, aes(x = delta, y = coverage, color = method)) +
    geom_point(shape = 4, size = 2, stroke = 1.5) +
    geom_line(lwd = 1.2) +
    ylim(c(-0.05, 1)) +
    geom_hline(yintercept = 0.95, color = "red", lty = 2, lwd = 1.2) +
    xlab(bquote(delta ~ "= exp(" * eta[0] * ")")) + 
    ylab("Coverage rates") +
    theme_pubr() +
    theme(axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_line(),
          strip.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = "bottom")
  ggsave(filename = paste0("plots/tnd_continuous_coverage_tr_beta0_", round(bb * 100),"_.png"), 
         plot = gg_coverage_tr,
         width = 6, height = 5)
}

