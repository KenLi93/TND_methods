trt_bridge_oracle <- function(p_0a, p_ua, p_0z, p_uz) {
  q_z0 <- solve(matrix(c((1 - p_0a) * (1 - p_0z),
                         (1 - p_0a) * p_0z,
                         (1 - (p_0a + p_ua)) * (1 - (p_0z + p_uz)),
                         (1 - (p_0a + p_ua)) * (p_0z + p_uz)),
                       nrow = 2, ncol = 2, byrow = TRUE),
                c(1, 1)) 
  q_z1 <- solve(matrix(c(p_0a * (1 - p_0z),
                         p_0a * p_0z,
                         (p_0a + p_ua) * (1 - (p_0z + p_uz)),
                         (p_0a + p_ua) * (p_0z + p_uz)),
                       nrow = 2, ncol = 2, byrow = TRUE),
                c(1, 1)) 
  
  q_value <- data.frame(Z = c(0, 1, 0, 1),
                        A = c(0, 0, 1, 1),
                        q = c(q_z0, q_z1))
  
  return(q_value)
}

## estimate the bridge function in the control sample
trt_bridge_binary <- function(data) {
  n0 <- sum(data$Y == 0)
  
  ## parameters in the saturated model of the bridge function
  ## that solves the estimating equation
  tau <- with(data[data$Y == 0,], 
              solve(matrix(c(n0, sum(Z), sum(A), sum(Z * A),
                             sum(W), sum(W * Z), sum(W * A), sum(W * A * Z),
                             sum(A), sum(Z * A), sum(A), sum(Z * A),
                             sum(W * A), sum(W * A * Z), sum(W * A), sum(W * A * Z)),
                           nrow = 4, byrow = TRUE), 
                    c(2 * n0, 2 * sum(W), n0, sum(W))))
  
  q_value <- data.frame(Z = c(0, 1, 0 ,1),
                        A = c(0, 0, 1, 1))
  
  q_value$q <- with(q_value,
                    cbind(1, Z, A, Z * A) %*% tau)
  
  return(list(tau = tau, q_value = q_value))
}




tnd_nc_rr <- function(data, tau) {
  ## study sample size 
  nn <- nrow(data)
  
  ## estimate the treatment bridge function  
  tau_est <- tau
  
  ## given the parameters of the saturated model, calculate the values of the treatment bridge function
  q_val <- with(data, cbind(1, Z, A, Z * A) %*% tau_est)
  
  
  beta_hat <- log(sum(q_val[data$Y == 1 & data$A == 1]) / 
                    sum(q_val[data$Y == 1 & data$A == 0]))
  
  ## sandwich variance matrix for (beta, tau)
  Omega11 <- with(data, mean((1 - 2 * A) * q_val * Y * A * exp(-beta_hat * A)))
  
  
  Omega12 <- with(data, matrix(c(mean((2 * A - 1) * Y * exp(-beta_hat * A)),
                                 mean((2 * A - 1) * Y * Z * exp(-beta_hat * A)),
                                 mean((2 * A - 1) * Y * A * exp(-beta_hat * A)),
                                 mean((2 * A - 1) * Y * Z * A * exp(-beta_hat * A))),
                               nrow = 1))
  
  Omega21 <- matrix(rep(0, 4), ncol = 1)
  
  Omega22a <- array(NA, dim = c(4, 4, nn))
  
  for (ii in 1:nn) {
    Omega22a[, , ii] <- with(data, (1 - Y[ii]) * c(1, W[ii], A[ii], W[ii] * A[ii]) %o% 
                               c(1, Z[ii], A[ii], Z[ii] * A[ii]))
  }
  Omega22 <- apply(Omega22a, c(1, 2), mean)
  
  Omega <- rbind(cbind(Omega11, Omega12),
                 cbind(Omega21, Omega22))
  
  
  
  Mi2 <- array(NA, dim = c(5, 5, nn))
  
  for (ii in 1:nn) {
    Mi <- matrix(nrow = 5, ncol = 1)
    Mi[1,  1] <- with(data, (2 * A[ii] - 1) * q_val[ii] * Y[ii] * exp(-beta_hat * A[ii]))
    Mi[2:5, 1] <- with(data, (1 - Y[ii]) * (c(1, W[ii], A[ii], W[ii] * A[ii]) * q_val[ii] -
                                              c(2, 2 * W[ii], 1, W[ii])))
    Mi2[, , ii] <- Mi %*% t(Mi)
  }
  
  VAR_M <- apply(Mi2, c(1, 2), mean)
  
  VAR <- solve(Omega) %*% VAR_M %*% t(solve(Omega)) / nn
  
  VAR_beta <- as.numeric(VAR[1, 1])
  
  return(list(beta = beta_hat, VAR_beta = VAR_beta, SE = sqrt(VAR_beta),
              tau = tau_est, varcov = VAR, CI = beta_hat + qnorm(c(0.025, 0.975)) * sqrt(VAR_beta)))
}


## The risk ratio estimate using the oracle estimate

tnd_nc_rr_oracle <- function(data, p_0a, p_ua, p_0z, p_uz) {
  ## study sample size 
  nn <- nrow(data)
  
  ## obtain the true value of the q function
  q_val_mat <- trt_bridge_oracle(p_0a, p_ua, p_0z, p_uz)
  
  q_val <- left_join(data, q_val_mat, by = c("Z", "A"))$q
  
  beta_hat <- log(sum(q_val[data$Y == 1 & data$A == 1]) / 
                    sum(q_val[data$Y == 1 & data$A == 0]))
  
  ## sandwich variance matrix for (beta, tau)
  Omega11 <- with(data, mean((1 - 2 * A) * q_val * Y * A * exp(-beta_hat * A)))
  VAR_M <- mean(with(data, q_val ^ 2 * Y * exp(-2 * beta_hat * A)))
  
  VAR_beta <- VAR_M / (Omega11 ^ 2 * nn)
  
  
  return(list(beta = beta_hat, VAR_beta = VAR_beta, SE = sqrt(VAR_beta),
              CI = beta_hat + qnorm(c(0.025, 0.975)) * sqrt(VAR_beta)))
}




logit_reg <- function(data) {
  reg_mod <- glm(Y ~ A, family = "binomial", data = data)
  beta_hat <- coef(reg_mod)["A"]
  VAR_beta <- sandwich::vcovHC(reg_mod, "HC0")["A", "A"]
  return(list(beta = beta_hat, VAR_beta = VAR_beta, SE = sqrt(VAR_beta),
              CI = beta_hat + qnorm(c(0.025, 0.975)) * sqrt(VAR_beta)))
}
