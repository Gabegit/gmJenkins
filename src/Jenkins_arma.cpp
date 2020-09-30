// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include "./test_h.h"

//using namespace Rcpp; //Rcout, Environment,Funcion
//using namespace arma; //mat,vec,randn,repmat,diagmat

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// seq for arma::uvec index
arma::uvec seq_uvec(int from, int to, int by =1){
  int len = (((to-from)/by)+1);
  arma::uvec res(len);
  for (int i = 0;i<len;i++){
    res[i] = from + i*by;
  }
  return res;
}

Rcpp::NumericVector seq_rcpp(double from, double to, double by =1.0){
  int len = (((to-from)/by)+1);
  Rcpp::NumericVector res(len);
  for (int i = 0;i<len;i++){
    res[i] = from + i*by;
  }
  return res;
}

Rcpp::NumericVector pnormc_rcpp(Rcpp::NumericVector x){
  Rcpp::NumericVector K_G = pnorm(x);
  Rcpp::LogicalVector l1 = (K_G > 0.999);
  Rcpp::LogicalVector l2= (K_G < 0.001);
  //K_G(l1) = 0.990;
  for(int i =0; i < l1.size();i++){
    if(l1[i]) K_G(i) = 0.999;
    if(l2[i]) K_G(i) = 0.001;
  }

  return K_G;
}

arma::vec qnorm_rcpp(arma::vec x){
  arma::vec res;
  res = qnorm(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(x)));
  return res;
}

// [[Rcpp::export]]
arma::vec pnorm_arma(arma::vec x){
  arma::vec K_G = normcdf(x);
  K_G.elem(find(K_G > 0.999)).fill(0.999); //
  K_G.elem(find(K_G < 0.001)).fill(0.001);
  return K_G;
}

arma::vec rowSums(arma::mat M){
  arma::vec rs;
  rs = sum(M, 1); //rowSums(M,1),colSums(M,0)
  return rs;
}

// log-likelihood function
// [[Rcpp::export]]
arma::vec dprobit_rcpp(arma::vec b,Rcpp::List data_list){
  int N = data_list["N"];
  int T = data_list["T"];
  int R = data_list["R"];
  int K1 = data_list["K1"];
  int K2 = data_list["K2"];
  // float rho1 = init_list["rho1"];
  // float rho2 = init_list["rho2"];
  Rcpp::List xi = data_list["xi"]; //Rcpp::List
  arma::mat X = data_list["X"];    //arma::mat
  arma::mat Y = data_list["Y"];
  arma::mat D = data_list["D"];
  //mat zeta = init_list["zeta"];

  arma::vec beta; beta = b(arma::span(0,K1-1)); //subvec of b:b(span(start,end)) or b.subvec(s,e)
   //beta.print("beta:");

  arma::vec beta_r; beta_r = b(arma::span(K1,K1+K2-1));
   //beta_r.print("beta_r");

  // sigu2 = cos(a[[3]][1])
  // zeta = cos(a[[3]][2])
  // rho1 = cos(a[[4]][1])
  // rho2 = cos(a[[4]][2])
  //
  float sigu2 = cos(b(K1+K2));
  float zeta = cos(b(K1+K2+1));
  float rho1 = cos(b(K1+K2+2));
  float rho2 = cos(b(K1+K2+3));
   //Rcout << "sigu2:" << sigu2 << "\n zeta:" << zeta << "\nrho1:" << rho1 << "\nrho2:" << rho2<< endl;

  // tau = rep(0, 4)
  //   tau[1] = -Inf
  //   tau[2] = a[[5]][1]
  // tau[3] = a[[5]][1] + abs(a[[5]][2])
  //   tau[4] = Inf
  //
  //   taur = c(-Inf, 0, Inf)
 arma::vec tau; tau = {-datum::inf,0,0,datum::inf};
 tau[1] = b(K1+K2+2+2+0);
 tau[2] = b(K1+K2+2+2+0) + fabs(b(K1+K2+2+2+1));
 arma::vec taur; taur = {-datum::inf,0,datum::inf};
   //tau.print("tau");
   //taur.print("taur");

 // SIGMA_ = diag(1, 2 * T - 1) # T
 //  r = sigu2 / (1 - zeta ^ 2) #0.2083333
 // c1 = rho1 * sqrt(r) #-0.1369306
 // c2 = rho2 * sqrt(r) #0.1825742
 arma::vec v = arma::ones<vec>(2*T-1);
 arma::mat SIGMA_; SIGMA_ = diagmat(v);
 float r = sigu2/(1 - pow(zeta,2));
 float c1 = rho1 * sqrt(r);
 float c2 = rho2 * sqrt(r);
 //SIGMA_.print("sigma_");
   //Rcout << "r:" << r << "\nc1:" << c1 << "\nc2:" << c2 << endl;

 // for (i in 1:(2 * T - 1)) {
 //   for (j in 1:(2 * T - 1)) {
 //     if (i %% 2 == 1) { # %% 取余 13 %% 5 =3
 //       if (abs(j - i) %% 2 == 1) {
 //         SIGMA_[i, j] = 0
 //       }
 //       if (abs(j - i) %% 2 == 0) {
 //         SIGMA_[i, j] = r * zeta ^ (abs(i - j) / 2) + 1 - r
 //       }
 //       if ((j - i) == 1) {
 //         SIGMA_[i, j] = c1
 //       }
 //       if ((j - i) == -1) {
 //         SIGMA_[i, j] = c2
 //       }
 //     }
 //     if (i %% 2 == 0) {
 //       SIGMA_[i, i + 1] = c2
 //       SIGMA_[i, i - 1] = c1
 //     }
 //   }
 // }

 for (int i =1; i < (2*T -1 + 1 ); i++){
   for (int j = 1; j < (2*T -1 + 1); j++){
     if (i % 2 ==1){
       if ( abs(j-i) %2 == 1){
         SIGMA_(i-1,j-1) = 0;
       }
       if (abs(j-i) % 2 == 0){
         SIGMA_(i-1,j-1) = r * pow(zeta,abs(i-j)/2) + 1 -r;
       }
       if ((j-i) == 1){
         SIGMA_(i-1,j-1) = c1;
       }
       if ((j-i) == -1){
         SIGMA_(i-1,j-1) = c2;
       }
     }
     if ( i % 2 == 0){
       SIGMA_(i-1,i+1-1) = c2;
       SIGMA_(i-1,i-1-1) = c1;
     }
   }
 }

   // SIGMA_.print("SIGMA_");

 // A = t(chol(SIGMA_))
 //   L_hat_r = matrix(0, N, R)
 //   L_hat = matrix(0, N, 1)
 //
 //   eta = matrix(0, N, 2 * T - 1)
 arma::mat A; A = chol(SIGMA_).t();
   //A.print("mat A");
 arma::mat L_hat_r(N,R);
 arma::mat L_hat(N,1);
 arma::mat eta (N, 2*T-1);

//  sum_A_eta1 = matrix(0, N, T - 1)
//  sum_A_eta2 = matrix(0, N, T - 1)
 arma::mat sum_A_eta1(N,T-1);
 arma::mat sum_A_eta2(N,T-1);

// ###initial period###
//     XX = X[1:N, ] # X
//     mu = XX %*% as.matrix(beta, 1, K1)
 arma::mat XX; XX = X.rows(0,N-1);
 arma::mat m_beta(K1,1); m_beta.col(0) = beta;
 arma::mat mu; mu = XX * m_beta;
   //mu.rows(0,4).print("mu[1,]");

   //Rcout << "mu.n_cols:" << mu.row(0) << "\nnu.n_cols:" << mu.n_rows << endl;

// ** armadillo tips
//  R::pnorm = arma::normcdf
//  R::dnorm = arma::normpdf

//     kernel1_upper = pnorm((tau[Y[, 1] + 1] - mu) / A[1, 1]) #tau: vec
//     kernel1_lower = pnorm((tau[Y[, 1]] - mu) / A[1, 1])
//     f1 = kernel1_upper - kernel1_lower
// RcppArmadillo can't deal NA like R,use Rcpp NumericMatrix
 arma::mat kernel1_upper, kernel1_lower,f1;

 arma::uvec ind_y1 = arma::conv_to<arma::uvec>::from(Y.col(0));

 kernel1_upper = normcdf((repmat(tau.elem(ind_y1), 1, 1)-mu)/A(0,0));
 kernel1_lower = normcdf((repmat(tau.elem(ind_y1-1), 1, 1)-mu)/A(0,0));

 f1 = kernel1_upper - kernel1_lower; //correct

 //kernel1_upper.rows(0,4).print("kernel1_upper:");
 //kernel1_lower.rows(0,4).print("kernel1_lower:");
  //f1.rows(0,4).print("f1");


  arma::mat each_A_eta1,each_A_eta2;
  arma::mat Xr,mur,temp;
  Rcpp::NumericVector xi1,xit,xitp1;
  arma::mat beta_r_mat(K2,1);
  arma::uvec xridx, xxidx, ind_D1, ind_Y;
  arma::mat kernel2_upper, kernel2_lower, f2;
  arma::mat kernel3_upper(N,1), kernel3_lower(N,1), f3,ff;

  for (int ri = 0; ri < R; ri++){
    ff = f1;
    //eta[, 1] = qnorm((1 - xi[[1]][r]) * kernel1_lower + xi[[1]][r] * kernel1_upper)
    xi1 = xi[0];
    eta.col(0) = qnorm_rcpp((1-xi1[ri])*kernel1_lower.col(0)+xi1[ri]*kernel1_upper.col(0));
    //eta.rows(0,4).print("eta(0,)"); //checked

//     for (t in seq(2, 2 * (T - 1), 2)) {
//       Xr =   X[seq(1 + (t / 2 - 1) * N, t / 2 * N, 1), ]
//       XX =   X[seq(1 + (t / 2) * N, (t / 2 + 1) * N, 1), ]
//       XX[, 1] = mu
//       mur = Xr %*% as.matrix(beta_r, 1, K2)
//       mu  = XX %*% as.matrix(beta, 1, K1)
//

   //mat each_A_eta1;
   //mat sum_A_eta1(N, 4);
   for (int t =2; t <=2*(T-1); t = t+2){
      // Rcout << t << "\n";
      xridx = arma::regspace<arma::uvec>(1+(t/2-1)*N -1, t/2 *N -1);
      //Rcout << xridx(span(0,4)) << "\n";
      xxidx = arma::regspace<arma::uvec>(1+(t/2)*N -1, (t/2 +1)*N -1);

      Xr = X.rows(xridx); //cheecked
      //Xr.rows(0,3).print("Xr[1:4,]");
      //Rcout << "ok here: mat Xr" << endl;

      XX = X.rows(xxidx);
      XX.col(0) = mu.col(0);
      //Rcout << "ok here: XX.col(0)" << endl;

      beta_r_mat.col(0) = beta_r;

      mur = Xr * beta_r_mat;
      mu = XX * m_beta;

      // ##########
      //       each_A_eta1 = matrix(0, N, t - 1)
      //         for (k in 1:(t - 1)) {
      //           each_A_eta1[, k] = A[t, k] * eta[, k]
      //         }
      //         sum_A_eta1[, t / 2] = rowSums(each_A_eta1)

      each_A_eta1.set_size(N,t-1);

      for (int k =1; k <= t-1;k++){
        each_A_eta1.col(k-1) = A(t-1,k-1) * eta.col(k-1);

      }
      // each_A_eta1.rows(0,3).print("each_A_eta1[0:3,]");

      //each_A_eta1.row(0).print("first row of each_A_eta1");
      sum_A_eta1.col(t/2-1) = arma::sum(each_A_eta1,1); //rowSums returns a col
      //sum_A_eta1.rows(0,3).print("sum_A_eta1[0:3,]");

     // ###########
     //  kernel2_upper = pnorm((taur[D[, t / 2] + 2] - mur - sum_A_eta1[, t /
     //    2]) / A[t, t])
     //    kernel2_lower = pnorm((taur[D[, t / 2] + 1] - mur - sum_A_eta1[, t /
     //      2]) / A[t, t])
     //    f2 = kernel2_upper - kernel2_lower
     //    f2[f2 == 0] = 0.001

     // ##########

     ind_D1 =arma:: conv_to<arma::uvec>::from(D.col(t/2-1));
     kernel2_upper = normcdf((repmat(taur.elem(ind_D1+1), 1, 1) - mur -sum_A_eta1.col(t/2-1))/A(t-1,t-1));
     kernel2_lower = normcdf((repmat(taur.elem(ind_D1), 1, 1) -mur -sum_A_eta1.col(t/2-1))/A(t-1,t-1));
     f2 = kernel2_upper - kernel2_lower;
     f2.elem(find(f2 == 0)).fill(0.01); //checked

     //f2.rows(0,5).print("f2[1:5,]:");

     //  temp = qnorm((1 - xi[[t]][r]) * kernel2_lower + xi[[t]][r] * kernel2_upper)
     // #temp = qnorm((1-xi[[2]][r])*kernel2_lower+xi[[2]][r]*kernel2_upper)
     // temp[is.na(temp)] = 0
     //  eta[, t] = temp
     //
     xit = xi[t-1]; //NumericVector
     temp = qnorm_rcpp((1-xit[ri-1])*kernel2_lower.col(0) + xit[ri-1]*kernel2_upper.col(0));
     temp.elem(find_nonfinite(temp)).zeros(); //checked
     //temp.rows(0,4).print("temp[0:4,]");

     eta.col(t-1) = temp.col(0);
     //Rcout << "ok here: eta[0:5,]:\n" << eta.rows(0,3) << endl;

       // ##########
       // each_A_eta2 = matrix(0, N, t)
       // for (k in 1:t) {
       //   each_A_eta2[, k] = A[t + 1, k] * eta[, k]
       // }
       // sum_A_eta2[, t / 2] = rowSums(each_A_eta2)
      each_A_eta2.set_size(N,t);

      for (int kk=1;kk <=t;kk++ ){
        each_A_eta2.col(kk-1) = A(t, kk-1)*eta.col(kk-1);
      }
      sum_A_eta2.col(t/2 -1) = sum(each_A_eta2, 1); //rowSums
        //Rcout << "ok here: sum_A_eta2[0:4,]\n" << sum_A_eta2.rows(0,4) << endl;

      // kernel3_upper = pnormc((tau[Y[, t / 2 + 1] + 1] -
      //   mu - sum_A_eta2[, t /2]) / A[t + 1, t + 1])
      //   kernel3_lower = pnormc((tau[Y[, t / 2 + 1]] - mu -
      //     sum_A_eta2[, t /	2]) / A[t + 1, t + 1])
      //   f3 = kernel3_upper - kernel3_lower
      //   f3[f3 == 0] = 0.001
      // f3[is.na(f3)] = 1
      // temp = qnorm((1 - xi[[t + 1]][r]) * kernel3_lower + xi[[t + 1]][r] *
      //   kernel3_upper)
      // #temp = qnorm((1-xi[[3]][r])*kernel3_lower+xi[[3]][r]*kernel3_upper)
      //   temp[is.na(temp)] = 0
      // eta[, t + 1] = temp
      //   ff = ff * f2 * f3

      //attention::arma doesn't deal Nan correctly!!!

      ind_Y =arma:: conv_to<arma::uvec>::from(Y.col(t/2 + 1 -1));

      for (int id=0; id < Y.n_rows; id++){
        if (ind_Y(id)==0){  // 0 means nan in uvec
          //Rcout << "is nan" << endl;
          kernel3_upper(id,0) = arma::datum::nan;
          kernel3_lower(id,0) = arma::datum::nan;
        } else {
          kernel3_upper(id,0) = arma::normcdf((tau(ind_Y(id)) - mu(id) - sum_A_eta2(id, t/2-1))/A(t,t));
          kernel3_lower(id,0) = arma::normcdf((tau(ind_Y(id)-1)- mu(id) - sum_A_eta2(id, t/2-1))/A(t,t));
         // Rcout << "not nan" << endl;
        }

      }

       kernel3_upper.elem(find(kernel3_upper > 0.999)).fill(0.999); //
       kernel3_upper.elem(find(kernel3_upper < 0.001)).fill(0.001);
       kernel3_lower.elem(find(kernel3_lower > 0.999)).fill(0.999); //
       kernel3_lower.elem(find(kernel3_lower < 0.001)).fill(0.001);

      //Y.rows(0,4).print("Y[1:4,1:4]");


      //ind_Y.elem(find_nonfinite(ind_Y)).fill(4); //change NA into 4

      //kernel3_upper = pnorm_arma((repmat(tau.elem(ind_Y+1-1),1,1) - mu- sum_A_eta2.col(t/2-1))/A(t,t));

      //Rcout << "kernel3_upper\n" << kernel3_upper.rows(0,15) << endl;

      //Rcout << "ind_Y:" << ind_Y(span(0,5)) << "ind_Y-1" <<ind_Y(span(0,5)) -1 << endl;

      //kernel3_lower = pnorm_arma((repmat(tau.elem(ind_Y-1),1,1) - mu- sum_A_eta2.col(t/2-1))/A(t,t));
     //Rcout << "ok here kernel3_lower" << endl;

      f3 = kernel3_upper - kernel3_lower;
      //f3.rows(0,4).print("f3[1,]");

      f3.elem(find(f3==0)).fill(0.001);
      f3.elem(find_nonfinite((f3))).ones();

      xitp1 = xi[t];
      temp = qnorm_rcpp((1- xitp1[ri-1])*kernel3_lower.col(0) + xitp1[ri-1]*kernel3_upper.col(0));
      temp.elem(find_nonfinite(temp)).zeros();
      eta.col(t) = temp.col(0);

      ff = ff % f2 % f3;

   }
    L_hat_r.col(ri) = ff;
   //sum_A_eta1.row(0).print("first row of sum_A_eta1");
  // A.print("A");
   //Rcout << "ok here: sum_A_eta1" << sum_A_eta1.row(0) <<  "\n";

  }
  L_hat = arma::sum(L_hat_r,1) / R;
	L_hat.elem(arma::find(L_hat == 0)).ones();
	arma::vec LL = log(L_hat);

	return(LL);


}
