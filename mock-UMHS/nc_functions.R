#============================================
# Functions for implementing the NC estimator
# and the IPTW estimator in Schnitzer 2022
# Author: Kendrick Li (qijunli@umich.edu)
# Date: 11/10/2022
#============================================
expit <- function(x) exp(x) / (1 + exp(x))

U_lr <- function(par, data) {
  tau <- par
  A <- data$A
  X <- data$X
  XX <- as.matrix(cbind(1, X))
  U <- c(A - exp(XX %*% tau)) * XX
}

U_rr_i <- function(par, data) {
  beta <- par[1]
  tau <- par[-1]
  A <- data$A
  X <- data$X
  Y <- data$Y
  XX <- as.matrix(cbind(1, X))
  pscore <- c(A / expit(XX %*% tau) + (1 - A) * (1 + exp(XX %*% tau)))
  U <- matrix((2 * A - 1) * Y * exp(- beta * A) * pscore, ncol = 1)
  return(U)
}

U_iptw <- function(par, data) {
  U <- cbind(U_lr(par[-1], data), U_rr_i(par, data))
  return(U)
}

GMMF <- function(mrf, param, data, ...) {
  g0 <- mrf(par = param, data = data, ...)
  g <- apply(g0, 2, mean)
  gmmf <- sum(g ^ 2)
  
  return(gmmf)
}

G <- function(bfun, para, data) {
  G <- numDeriv::jacobian(func = G1, bfun = bfun, x = para, data = data)
  return(G)
}

G1 <- function(bfun, param, data) {
  G1 <- apply(bfun(param, data), 2, mean, na.rm = T)
  return(G1)
}

var.gmmf <- function(bfun, param, data) {
  bG <- MASS::ginv(G(bfun, param, data))
  bg <- bfun(param, data)
  spsz <- dim(bg)[1]
  Omega <- t(bg) %*% bg / spsz
  Sigma <- bG %*% Omega %*% t(bG) / spsz
}

tnd_iptw <- function(data) {
  nn <- nrow(data)
  A <- data$A
  Y <- data$Y
  X <- as.matrix(data$X)
  tau_init <- coef(glm(A ~ X, family = binomial))
  
  beta_est <- log(mean(A * Y / expit(cbind(1, X) %*% tau_init)) /
                    mean((1 - A) * Y * (1 + exp(cbind(1, X) %*% tau_init))))
  
  m_q1 <- nlm(p = c(beta_est, tau_init), f = GMMF,
              mrf = U_iptw, data = data)
  
  param_est <- m_q1$estimate
  beta_est <- param_est[1]
  
  var_all <- var.gmmf(bfun = U_iptw, param = param_est, data)
  var_beta <- var_all[1, 1]
  
  return(list(beta = beta_est, var_beta = var_beta, se = sqrt(var_beta),
              ci = beta_est + qnorm(c(0.025, 0.975)) * sqrt(var_beta),
              ve = 1 - exp(beta_est),
              ve_ci = 1 - exp(beta_est + qnorm(c(0.975, 0.025)) * sqrt(var_beta))))
}

nc_covid_ve <- function(data, trt = "vaccinated", outcome = "COVID",
                        nce = "IMM", nco = "nco",
                        x = c("age_geq_18_le_60", "age_geq_60", "Gender", "caucasian", "wscore_geq_3",
                              "april", "may", "june", "july", "aug", "sep", "oct")) {
  nn <- nrow(data)
  
  A <- as.matrix(data[, trt])
  Z <- as.matrix(data[, nce])
  W <- as.matrix(data[, nco])
  Y <- as.matrix(data[, outcome])
  X <- as.matrix(data[, x])
  
  pp <- ncol(cbind(1, A, Z, Z * A, X))
  
  qfunc <- function(tau, A, Z, X) {
    q_est <- cbind(1, A, Z, Z * A, X) %*% tau
  }
  
  MM <- cbind(1, A, W, W * c(A), X, X * c(A))
  
  U_q <- function(par, A, Z, X, Y) {
    tau <- par
    q_est <- qfunc(tau, A[Y == 0,], Z[Y == 0,], X[Y == 0,])
    
    U <- MM[Y == 0,] * c(q_est) - cbind(2, 1, 2 * W, W, 2 * X, X)[Y == 0,]
    return(U)
  }
  
  U_rr <- function(par, data, q_est) {
    beta <- par
    U <- matrix((2 * A - 1) * Y * exp(-beta * A) * q_est, ncol = 1)
    return(U)
  }
  
  GMMF <- function(mrf, param, ...) {
    g0 <- mrf(par = param, ...)
    g <- apply(g0, 2, mean)
    gmmf <- sum(g ^ 2)
    
    return(gmmf)
  }
  
  m_q1 <- nlm(p = rep(0, pp), f = GMMF, mrf = U_q, A = A, Z = Z, X = X, 
              Y = Y, iterlim = 1000)
  
  tau_est <- m_q1$estimate
  
  q_est <- qfunc(tau_est, A, Z, X)
  
  beta_hat <- log(sum(q_est[Y == 1 & A == 1]) / sum(q_est[Y == 1 & A == 0]))
  
  Omega11 <- mean((1 - 2 * A) * q_est * Y * A * exp(-beta_hat * A))
  
  Omega12 <- matrix(c(mean((2 * A - 1) * Y * exp(-beta_hat * A)),
                      mean((2 * A - 1) * Y * exp(-beta_hat * A) * A),
                      mean((2 * A - 1) * Y * exp(-beta_hat * A) * Z),
                      mean((2 * A - 1) * Y * exp(-beta_hat * A) * A * Z),
                      colMeans(matrix(rep((2 * A - 1) * Y * exp(-beta_hat * A), ncol(X)), ncol = ncol(X)) * X)),
                    nrow = 1)
  Omega21 <- matrix(rep(0, ncol(MM)), ncol = 1)
  Omega22a <- array(NA, dim = c(ncol(MM), pp, nn))
  
  
  for (ii in 1:nn) {
    Omega22a[, , ii] <- (1 - Y[ii]) * MM[ii,] %o% c(1, A[ii], Z[ii], A[ii] * Z[ii], X[ii,])
  }
  
  Omega22 <- apply(Omega22a, c(1, 2), mean)
  
  Omega <- rbind(cbind(Omega11, Omega12),
                 cbind(Omega21, Omega22))
  
  Mi2 <- array(NA, dim = c(1 + ncol(MM), 1 + ncol(MM), nn))
  
  for (ii in 1:nn) {
    Mi <- matrix(nrow = 1 + ncol(MM), ncol = 1)
    Mi[1, 1] <- (2 * A[ii] - 1) * q_est[ii] * Y[ii] * exp(-beta_hat * A[ii])
    Mi[2:(1 + ncol(MM)), 1] <- (1 - Y[ii]) * MM[ii, ] * q_est[ii] - c(2, 1, 2 * W[ii, ], W[ii, ], 2 * X[ii, ], X[ii, ])
    
    Mi2[, , ii] <- Mi %*% c(Mi)
  }
  
  VAR_M <- apply(Mi2, c(1, 2), mean)
  
  VAR <- ginv(Omega) %*% VAR_M %*% t(ginv(Omega)) / nn
  
  VAR_beta <- as.numeric(VAR[1, 1])
  VAR_tau <- diag(VAR)[-1]
  
  logit_reg <- glm(Y ~ A + X, family = "binomial")
  
  beta_logit <- coef(logit_reg)["A"]
  sd_logit_reg <- summary(logit_reg)$coefficients["A", "Std. Error"]
  
  return(list(beta = beta_hat,
              VAR_beta = VAR_beta,
              tau_hat = tau_est,
              VAR_tau = VAR_tau,
              CI = beta_hat + qnorm(c(0.025, 0.975)) * sqrt(VAR_beta),
              VE = 1 - exp(beta_hat),
              VE_CI = 1 - exp(beta_hat + qnorm(c(0.975, 0.025)) * sqrt(VAR_beta)),
              beta_logit = beta_logit,
              CI_logit = beta_logit + qnorm(c(0.025, 0.975)) * sd_logit_reg,
              VE_logit = 1 - exp(beta_logit),
              VE_CI_logit = 1 - exp(beta_logit + qnorm(c(0.975, 0.025)) * sd_logit_reg),
              logit_more = summary(logit_reg)$coefficients))
}
