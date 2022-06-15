

expit <- function(x) {exp(x) / (1 + exp(x))}


## propensity score for logistic regression
U_lr <- function(par, data){ ## dim = 4
  tau <- par
  U <- with(data, 
            cbind(A - expit(tau)))
  return(U)
}



U_rr_i <- function(par, data){ ## dim=1
  
  beta <- par[1]
  tau <- par[-1]
  pscore <- with(data, c(A / expit(tau) + (1 - A) * (1 + exp(tau))))
  U <- with(data,
            matrix((2 * A - 1) * Y * exp(-beta * A) * pscore, ncol = 1))
  return(U)
}

U_iptw <- function(par, data){ ## dim=1
  U <- cbind(U_lr(par[-1], data), U_rr_i(par, data))
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
var.gmmf <- function(bfun, param, data){
  bG <- solve(G(bfun, param, data))
  bg <- bfun(param, data)
  spsz <- dim(bg)[1]
  Omega <- t(bg) %*% bg / spsz
  Sigma <- bG %*% Omega %*% t(bG) / spsz
  return(Sigma)
}




tnd_iptw <- function(data) {
  nn <- nrow(data)
  
  ## first estimate q
  tau_est <- coef(glm(A ~ 1, data = data, family = binomial))
  beta_est <- with(data,
                   log(mean(A * Y / expit(tau_est)) / 
                     mean((1 - A) * Y * (1 + exp(tau_est)))))

  ## use the estimated q to estimate the risk ratio
  param_est <- c(beta_est, tau_est)
  
  var_all <- var.gmmf(bfun = U_iptw, para = param_est, data)
  var_beta <- var_all[1, 1]
  
  return(list(beta = beta_est, VAR_beta = var_beta, SE = sqrt(var_beta),
              CI = beta_est + qnorm(c(0.025, 0.975)) * sqrt(var_beta)))
} 


