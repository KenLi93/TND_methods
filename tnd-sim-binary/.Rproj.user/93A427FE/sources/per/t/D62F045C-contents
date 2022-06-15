qfunc <- function(tau, data){
  q_est <- with(data, c(cbind(1, Z, A, Z * A) %*% tau))
  return(q_est)
}



U_q <- function(par, data){ ## dim = 4
  tau <- par
  q_est <- qfunc(tau, data[data$Y == 0, ])
  U <- with(data[data$Y == 0, ], 
            q_est * cbind(1, W, A, W * A) - cbind(2, 2 * W , 1, W))
  return(U)
}



U_rr <- function(par, data, q_est){ ## dim=1
  
  beta <- par
  U <- with(data,
            matrix((2 * A - 1) * Y * exp(-beta * A) * q_est), ncol = 1)
  return(U)
}


U_all <- function(par, data){
  tau <- par[2:5]
  beta <- par[1]
  q_est <- qfunc(tau, data)
  U_q <- with(data,
              c((1 - Y) * q_est) * cbind(1, W, A, W * A) - (1 - Y) * cbind(2, 2 * W, 1, W))
  U_rr <- with(data,
               (2 * A - 1) * Y * exp(- beta * A) * q_est)
  U <- cbind(U_q, U_rr)
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



tnd_nc_gmm_binary <- function(data) {
  ## first estimate q
  m_q1 <- nlm(f = GMMF, p = rep(0,4),
              mrf = U_q, data=data)
  
  tau_est <- m_q1$estimate
  q_est <- qfunc(tau_est, data)
  
  ## use the estimated q to estimate the risk ratio
  m_rr <- nlm(f = GMMF, p = 0, mrf = U_rr,
              data = data, q_est = q_est)
  rr <- m_rr$estimate
  
  ## compute the variance using the joint estimating equation
  var_all <- var.gmmf(bfun = U_all, para = c(tau_est, rr), data)
  var_rr <- var_all[4,4]
} 

