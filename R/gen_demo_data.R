#' genenate a list of data with missing observation
#' @description genenate a list of data with missing observation for demo tests
#' @param N: number of observations
#' @param R: repeating times
#' @param T: time periods
#' @param K: something you know
#' @param rho1,rho2,sigu2,zeta: parameters you give to test
#' @return a list:
#' @return list$b is a initial parameter vector;
#' @return list$data is a list: data.frame X, Y with missing data, K1, K2,T = T,D xi,R,N.
#' @return list$real_para is a vector of the real paramaters.
#' @export
#' @examples
#' input_data = gen_demo_data()
#' input_data = gen_demo_data(N=20000,R=50,T=5,K=5,rho1 = -0.3, rho2 = 0.4, sigu2 = 0.6,zeta = 0.4)
#' b = input_data$b
#' data_list = input_data$data_list

gen_demo_data = function(N = 20000, R = 50, T = 5,K = 5,
                         rho1 = -0.3, rho2 = 0.4, sigu2 = 0.6,
                         zeta = 0.4 ){
  suppressWarnings(suppressPackageStartupMessages({
    library(MASS)
    library(stats)
  }))
  set.seed(2468)

  sig_u2 = 0.6#runif(1,0.2,0.8)
  zeta = 0.4#runif(1,-0.8,0.8)
  r = sig_u2 / (1 - zeta ^ 2)
  rho1 = -0.3
  rho2 = 0.4
  c1 = rho1 * sqrt(r)
  c2 = rho2 * sqrt(r)

  real_para <- list(rho1 = rho1, rho2=rho2,sigu2=sig_u2,zeta=zeta)

  pnormc = function(x) {
    k_G = pnorm(x)
    #k_G = 1/(1+exp(-pi*x/sqrt(3)))
    k_G[k_G > 0.999] = 0.999
    k_G[k_G < 0.001] = 0.001
    return(k_G)
  }

  SIGMA = diag(1, 2 * T - 1)

  beta1 = c(0.2009768,  0.4320059,  0.1347381,-0.1420730,  0.8136340)
  betar = c(-0.3589200,  0.6910340,  0.9173814,-0.3036093,  1.0749171)
  beta3 = c(0.6943055, 0.4320059,  0.1347381,-0.3420730,  0.8136340)
  X = matrix(0, N * T, K)

  for (i in 1:(2 * T - 1)) {
    for (j in 1:(2 * T - 1)) {
      if (i %% 2 == 1) {
        if (abs(j - i) %% 2 == 1) {
          SIGMA[i, j] = 0
        }
        if (abs(j - i) %% 2 == 0) {
          SIGMA[i, j] = r * zeta ^ (abs(i - j) / 2) + 1 - r
        }
        if ((j - i) == 1) {
          SIGMA[i, j] = c1
        }
        if ((j - i) == -1) {
          SIGMA[i, j] = c2
        }
      }
      if (i %% 2 == 0) {
        SIGMA[i, i + 1] = c2
        SIGMA[i, i - 1] = c1
      }
    }
  }

  #chol(SIGMA)

  eps = mvrnorm(N, mu = rep(0, 2 * T - 1), Sigma = SIGMA)

  X1_sig = diag(runif(K, 0.5, 1.5))
  Y = matrix(0, N, T)
  Re = matrix(0, N, T - 1)

  Xp1 = mvrnorm(N, rnorm(5), X1_sig)

  X[seq(1, N, 1), ] = Xp1
  mu1 = Xp1 %*% beta1 + eps[, 1]
  zeta = c(-Inf, quantile(mu1, probs = seq(0, 2 / 3, 1 / 3))[2:3], Inf)
  for (j in 1:3) {
    Y[, 1][mu1 >= zeta[j] & mu1 < zeta[j + 1]] = j
  }

  for (t in 2:T) {
    Xp1 = mvrnorm(N, rnorm(5), X1_sig)
    Xp1[, 1] = mu1
    X[seq(1 + (t - 1) * N, t * N, 1), ] = Xp1

    mu1 = Xp1 %*% beta1 + eps[, (2 * t - 1)]
    zeta = c(-Inf, quantile(mu1, probs = seq(0, 2 / 3, 1 / 3))[2:3], Inf)

    for (j in 1:3) {
      Y[, t][mu1 >= zeta[j] & mu1 < zeta[j + 1]] = j
    }
  }


  for (t in 1:(T - 1)) {
    #
    R_star =  X[seq(1 + (t - 1) * N, t * N, 1), ] %*% betar + eps[, t * 2]
    Re[, t] = I(R_star > 0)
    #
  }

  Y[, 2][Re[, 1] == 0] = NA
  Y[, 3][Re[, 2] == 0] = NA
  Y[, 4][Re[, 3] == 0] = NA
  Y[, 5][Re[, 4] == 0] = NA
  #summary(Y[, 2])
  #summary(Y[, 3])
  #summary(Y[, 4])


  D = 1 - is.na(Y[, -1])

  ###################model###############

  xi = list()
  for (i in 1:(2 * (T) - 1)) {
    xi[[i]] = runif(R)
  }
  Oprobit1 = MASS::polr(as.factor(Y[, 1]) ~ X[1:N, ],
                  , method = "probit")
  #summary(Oprobit1)


  Retention = glm(as.factor(Re[, 1]) ~ X[1:N, ] - 1,
                  , family = binomial(link = "probit"))
  #summary(Retention)

  a = list()

  a[[1]] = coef(Oprobit1)
  a[[2]] = coef(Retention)
  a[[3]] = rep(acos(0.2), 2)
  a[[4]] = rep(acos(0.2), 2)
  a[[5]] = c(Oprobit1$zeta[1], Oprobit1$zeta[2] - Oprobit1$zeta[1])


  b = unlist(a)

  K1 = length(coef(Oprobit1))
  K2 = length(coef(Retention))

  init_list <- list( X = X, #mat
                     Y = Y, #mat
                     K1 = K1, #int
                     K2 = K2, #int
                     T = T, # int
                     D = D, #mat
                     xi = xi, #List
                     R = R, # int
                     N= N # int
  )
  return(list(b=b, data_list = init_list, real_para = real_para))

}
