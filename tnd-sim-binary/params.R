p_u <- 0.5

p_0z <- 0.2; p_uz <- 0.4

p_0a <- 0.2; p_ua <- 0.4

param_grid <- rbind(expand.grid(beta = log(c(0.2, 0.5, 0.7, 1)),
                                eta_0y = log(c(0.01))),
                    expand.grid(beta = log(c(0.2, 1)),
                                eta_0y = log(c(0.02, 0.05, seq(0.1, 0.7, 0.1)))))


eta_uy <- log(0.5)

p_0w <- 0.02; p_uw <- -0.01
p_0d <- 0.02; p_ud <- -0.015

p_ys <- 0.1; p_uys <- 0.4
