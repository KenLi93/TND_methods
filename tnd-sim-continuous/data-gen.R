
U <- runif(N, 0, 1)

X <- runif(N, 0, 1)

A <- rbinom(N, 1, expit(mu_0a + mu_ua * U + mu_xa * X))

Z <- rnorm(N, mu_0z + mu_az * A + mu_xz * X + mu_uz * U, sigma_z)

W <- rnorm(N, mu_0w + mu_xw * X + mu_uw * U, sigma_w)

Y <- rbinom(N, 1, exp(mu_0y + beta0 * A + mu_uy * U + mu_xy * X + mu_uxy * U * X))

D <- rbinom(N, 1, exp(mu_0d + mu_xd * X + mu_ud * U))

S <- pmax(Y, D) * rbinom(N, 1, expit(mu_0s + mu_xs * X + mu_us * U + mu_uxs * U * X))

full_data <- data.frame(X = X, A = A, Z = Z, W = W, Y = Y, U = U)
study_data <- full_data[S == 1, ]
