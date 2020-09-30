#' the main function to estimate models
#' @param b: a vector of initial parameters
#' @param data_list: a list of given data
#' @param real_para: a vector of real parameters you want to test
#' @return result: list of coef of MLE, LIST OF PARAMETERS.
#' @export
#' @examples
#' input_data <- gen_demo_data(N=20000,R=50,T=5,K=5,rho1 = -0.3, rho2 = 0.4, sigu2 = 0.6,zeta = 0.4)
#' b <- input_data$b
#' data_list <- input_data$data
#' real_para <- input_data$real_para
#' #main function to produce the result
#' jenkins_main(b=b, data_list = data_list, real_para = real_para)

jenkins_main <- function(b, data_list, real_para = NULL){
  suppressWarnings(suppressPackageStartupMessages({
    library(stats)
  }))
  MLE =  maxLik::maxLik(
    logLik = dprobit_rcpp ,
    print.level = 2,
    start = b,
    data_list = data_list,
    method = "BHHH")

  print(summary(MLE))

  par = split(coef(MLE), rep(1:5, c(5, 5, 2, 2, 2)))

  beta = par[[1]]
  beta_r = par[[2]]
  sigu2 = cos(par[[3]][1])
  zeta = cos(par[[3]][2])
  rho1 = cos(par[[4]][1])
  rho2 = cos(par[[4]][2])

  ## print real parameters
  if (is.null(real_para) == FALSE) {
  print("===real parameters:===\n")
  print(paste("rho1 = ", real_para$rho1))
  print(paste("rho2 = ", real_para$rho2))
  print(paste("sigu2 = ", real_para$sigu2))
  print(paste("zeta = ", real_para$zeta))
  }
  print("\n===estimated parameters:===\n")

  print(paste("rho1: ", rho1))
  print(paste("rho2: ", rho2))
  print(paste("sigu2: ",sigu2))
  print(paste("zeta: ",zeta))
  print(paste("beta: ", beta))
  print(paste("beta_r: ", beta_r))

  return(list(MLEcoef=coef(MLE),rho1= rho2,rho2=rho2,zeta=zeta,sigu2=sigu2,beta=beta,beta_r= beta_r))

}

