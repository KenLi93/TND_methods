N <- 7000000


U <- rbinom(N, 1, p_u)  ## p_u = 0.5
Z <- rbinom(N, 1, p_0z + p_uz * U)  ## p_0z = 0.4, p_uz = 0.2
A <- rbinom(N, 1, p_0a + p_ua * U)  ## p_0a = 0.4, p_ua = 0.2

Y <- rbinom(N, 1, exp(eta_0y + beta0 * A + eta_uy * U))  

W <- rbinom(N, 1, p_0w + p_uw * U)  ## p_0w = 0.4, p_uw = 0.2

D <- rbinom(N, 1, p_0d + p_ud * U)  ## p_
S <- rbinom(N, 1, p_ys * pmax(Y, D, W) + p_uys * U * pmax(Y, D, W))


full_data <- data.frame(A = A, Y = Y, Z = Z, W = W, S = S, U = U)

study_data <- full_data[full_data$S==1, c("A", "Y", "Z", "W")]

                          