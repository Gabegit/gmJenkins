\name{gmJenkins-package}
\alias{gmJenkins-package}

\docType{package}
\title{
\packageTitle{gmJenkins}
}

\arguments{
  \item{b}{a vector of initial parameters for maxLik in MLE}
  \item{data_list}{a list of observed or generated data }
  \item{read_para}{a vector of real parameters you want to test}
}
\value{
  \code{gen_demo_data()} genenate a list of data with missing observation to demo.

  \code{jenkins_main()} Estimate jenkins' models using Xun's MLE approach.

  \code{dprobit_rcpp()} a log-likelihood function used in \code{maxLik}

}


\description{
\packageDescription{gmJenkins} It is used to estimate the parameters from the missing data.
}
\details{

The DESCRIPTION file:

 First, using \code{gen_demo_data()} function to generate a list data including initial parameter \code{b} and real data \code{input_data}. Of course, you can use your own data as the argments of the function \code{jenkins_main(b, data_list)}.

 Then, using \code{jenkins_main(b,b, data_list)} to estimate the model with maxLik, which will take a couple of minutes, which depends on the computer performance you use, such as about 7 mins in my Macbook pro 2014.
}
\author{
\packageAuthor{gmJenkins}

Maintainer: \packageMaintainer{gmJenkins}
}
\references{
~~Jenkins et al.(2004), Xun et al.(2020) ~~
}

\examples{
input_data <- gen_demo_data(N=20000,R=50,T=5,K=5,rho1 = -0.3, rho2 = 0.4, sigu2 = 0.6,zeta = 0.4)
b <- input_data$b

data_list <- input_data$data

real_para <- input_data$real_para

sum(dprobit_rcpp(b, data_list))

jenkins_main(b=b, data_list = data_list)
jenkins_main(b=b, data_list = data_list, real_para = real_para)

}
