
qfunc <- function(tau, data){
  q_est <- with(data, c(1 + exp((-1) ^ A * cbind(1, A, Z, X) %*% tau)))
  return(q_est)
}



U_q <- function(par, data){ ## dim = 4
  tau <- par
  q_est <- qfunc(tau, data[data$Y == 0, ])
  U <- with(data[data$Y == 0, ], 
            q_est * cbind(1, A, W, X) - cbind(2, 1, 2 * W, 2 * X))
  return(U)
}


GMMF <- function(mrf, param, data, ...){
  g0 <- mrf(par = param, data = data, ...)
  g <- apply(g0, 2, mean)
  gmmf <- sum(g ^ 2)
  
  return(gmmf)
}

G <- function(bfun, para, data){
  G <- numDeriv::jacobian(func = G1, bfun = bfun, x = para, data = data)
  return(G)
}

G1 <- function(bfun, param, data){
  G1 <- apply(bfun(param, data), 2, mean, na.rm=T)
  return(G1)
}

## sandwich variance
var_gmmf <- function(param, data){
  bG <- solve(G(bfun, param, data))
  bg <- bfun(param, data)
  spsz <- dim(bg)[1]
  Omega <- t(bg) %*% bg / spsz
  Sigma <- bG %*% Omega %*% t(bG) / spsz
  return(Sigma)
}





tnd_nc_gmm_continuous <- function(data) {


  
  nn <- nrow(data)
  
  ## first estimate q
  m_q1 <- nlm(p = rep(0, 4), f = GMMF,
                mrf = U_q, data = data)
  

  tau_est <- m_q1$estimate
  q_est <- qfunc(tau_est, data)
  
  ## use the estimated q to estimate the risk ratio

  beta_hat <- log(sum(q_est[data$Y == 1 & data$A == 1]) / 
                    sum(q_est[data$Y == 1 & data$A == 0]))
  
  ## sandwich variance matrix for (beta, tau)
  Omega11 <- with(data, mean((1 - 2 * A) * q_est * Y * A * exp(-beta_hat * A)))

  
  Omega12 <- with(data, matrix(c(mean(- Y * exp(-beta_hat * A) * (q_est - 1)),
                                 mean(- Y * exp(-beta_hat * A) * A * (q_est - 1)),
                                 mean(- Y * exp(-beta_hat * A) * Z * (q_est - 1)),
                                 mean(- Y * exp(-beta_hat * A) * X * (q_est - 1))),
                               nrow = 1))
  
  Omega21 <- matrix(rep(0, 4), ncol = 1)
  
  
  Omega22a <- array(NA, dim = c(4, 4, nn))
  
  for (ii in 1:nn) {
    Omega22a[, , ii] <- with(data, (1 - Y[ii]) * (1 - 2 * A[ii]) * (q_est[ii] - 1) *
                                c(1, A[ii], W[ii], X[ii]) %o% c(1, A[ii], Z[ii], X[ii]))
  }
  Omega22 <- apply(Omega22a, c(1, 2), mean)
  
  Omega <- rbind(cbind(Omega11, Omega12),
                 cbind(Omega21, Omega22))
  
  
  
  Mi2 <- array(NA, dim = c(5, 5, nn))
  
  for (ii in 1:nn) {
    Mi <- matrix(nrow = 5, ncol = 1)
    Mi[1,  1] <- with(data, (2 * A[ii] - 1) * q_est[ii] * Y[ii] * exp(-beta_hat * A[ii]))
    Mi[2:5, 1] <- with(data, (1 - Y[ii]) * (c(1, A[ii], W[ii], X[ii]) * q_est[ii] -
                                              c(2, 1, 2 * W[ii], 2 * X[ii])))
    Mi2[, , ii] <- Mi %*% c(Mi)
  }
  
  VAR_M <- apply(Mi2, c(1, 2), mean)
  
  VAR <- solve(Omega) %*% VAR_M %*% t(solve(Omega)) / nn
  
  VAR_beta <- as.numeric(VAR[1, 1])
  
  return(list(beta = beta_hat, VAR_beta = VAR_beta, SE = sqrt(VAR_beta),
              tau = tau_est, varcov = VAR, CI = beta_hat + qnorm(c(0.025, 0.975)) * sqrt(VAR_beta)))
} 



tnd_nc_continuous_oracle <- function(data, mu_0a, mu_ua, mu_0z, mu_uz, sigma_z, mu_xa, mu_xz) {
  nn <- nrow(data)
  
  ## first estimate q
  
  tau_true <- c(mu_0a - mu_ua * mu_0z / mu_uz - sigma_z ^ 2 * mu_ua ^ 2 / (2 * mu_uz ^ 2),
                sigma_z ^ 2 * mu_ua ^ 2 / mu_uz ^ 2 - mu_ua * mu_az / mu_uz,
                mu_ua / mu_uz,
                mu_xa - mu_xz * mu_ua / mu_uz)
  
  q_val <- qfunc(tau_true, data)
  
  ## use the estimated q to estimate the risk ratio
  
  beta_hat <- log(sum(q_val[data$Y == 1 & data$A == 1]) / 
                    sum(q_val[data$Y == 1 & data$A == 0]))

  
  ## sandwich variance matrix for beta
  Omega11 <- with(data, mean((1 - 2 * A) * q_val * Y * A * exp(-beta_hat * A)))
  VAR_M <- mean(with(data, q_val ^ 2 * Y * exp(-2 * beta_hat * A)))
  
  VAR_beta <- VAR_M / (Omega11 ^ 2 * nn)
  
  return(list(beta = beta_hat, VAR_beta = VAR_beta, SE = sqrt(VAR_beta),
              CI = beta_hat + qnorm(c(0.025, 0.975)) * sqrt(VAR_beta)))
} 

