
logit <- function(x) log(x / (1 - x))

expit <- function(x) exp(x) / (1 + exp(x))

N <- 7000000

mu_0a <- -1; mu_ua <- 1; mu_xa <- 0.25

mu_0z <- 0; mu_az <- 0.25; mu_xz <- 0.25; mu_uz <- 4; sigma_z <- 0.25

mu_0w <- 0; mu_xw <- 0.25; mu_uw <- 2; sigma_w <- 0.25


param_grid <- rbind(expand.grid(beta = log(c(0.2, 0.5, 0.7, 1)),
                                mu_0y = log(c(0.01))),
                    expand.grid(beta = log(c(0.2, 1)),
                                mu_0y = log(c(0.02, 0.05, seq(0.1, 0.7, 0.1)))))

mu_uy <- -2

mu_xy <- -0.25; mu_uxy <- -0.75

mu_0d <- log(0.01); mu_xd <- 0.25; mu_ud <- -0.2

mu_0s <- -1.4; mu_xs <- 0.5; mu_us <- 2; mu_uxs <- 1

# c(mu_0a - mu_ua * mu_0z / mu_uz - sigma_z ^ 2 * mu_ua ^ 2 / (2 * mu_uz ^ 2),
#  sigma_z ^ 2 * mu_ua ^ 2 / mu_uz ^ 2 - mu_ua * mu_az / mu_uz,
#  mu_ua / mu_uz,
#  mu_xa - mu_xz * mu_ua / mu_uz)