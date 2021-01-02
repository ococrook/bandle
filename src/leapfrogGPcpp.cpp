# include <RcppArmadillo.h>
# include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//' @title Computes gradients of hyperaparamters with respect to GP
//' 
//' @param Y The data
//' @param drv matrix derivative
//' @param nk number of observations
//' @param D number of samples
//' @param Z internal matrix required for computations
//' @param A internal matrix required for computations
//' @param sigmak standard deviation parameter
arma::vec gradienthyp(arma::vec Y,
                      arma::mat drv,
                      double nk,
                      double D,
                      arma::mat Z,
                      arma::mat A,
                      double sigmak){
  
  arma::mat B;
  arma::mat J = arma::ones(nk, nk);
  arma::mat eye = arma::eye(D, D);
  arma::mat C;
  arma::vec gradienthyp;
  arma::mat invZdrv = Z * drv * Z;
  
  //B = (((Xk.t() * covInv) * arma::kron(J, drv)) * (covInv * Xk)) / 2 ;
  B = (((Y.t() * invZdrv) * Y)) / (2 * sigmak * sigmak) ;
  //C = arma::accu(arma::kron(J, drv) % covInv)/2; //use trace is sum of entrywise product
  C = nk * trace(drv)/(2 * sigmak) - nk * arma::accu((eye - Z) % drv)/(2 * sigmak);
  
  gradienthyp = B - C;
  
  return(gradienthyp); 
}
//' @title Computes gradients of  GP
//' @param Xk The data
//' @param tau indexing term
//' @param h vector of hyperparamters
//' @param nk number of observations
//' @param D number of samples
// [[Rcpp::export]]
arma::vec gradientGPcpp(arma::mat Xk,
                        arma::vec tau,
                        arma::vec h,
                        int nk,
                        int D
){
  
  //terms are made unconstrained
  double sigmak; sigmak = exp(2 * h(2));
  double a; a = exp(2 * h(1));
  double l; l = exp(h(0));
  double negloglike;
  double logdet;
  arma::mat Z(D,D);
  arma::mat S(D,D); S = arma::zeros(D,D);
  arma::mat A(D,D);
  arma::mat R(D,D);
  arma::mat covInv;
  arma::mat drvl;
  arma::mat drva;
  double drvsigma;
  arma::vec temp;
  arma::vec v;
  arma::vec grdl;
  arma::vec grda;
  arma::vec grdsigma;
  arma::vec grad(3);
  arma::vec Y(D);
  double Y2;
  List trenchres = List::create(Rcpp::Named("logdet"),
                                Rcpp::Named("z"),
                                Rcpp::Named("v"));
  
  S = S.each_col() + tau;
  A = a * exp(- pow(S - S.t(), 2)/l);
  R = arma::eye(D,D) + (nk * A)/sigmak;
  
  trenchres = trenchDetcpp(R.row(0).t());
  logdet = trenchres["logdet"];
  v = as<arma::vec>(trenchres["v"]);
  
  if(R_IsNA(logdet)){
    negloglike = R_NegInf;
    grad.row(0) = negloglike;
    return(grad);
  }
  
  Z = trenchInvcpp(v);
  //covInv = covInvcpp(Z, A, sigmak, nk, D);
  
  drvl = A % (pow(S - S.t(), 2)/l);
  drva = 2 * A ;
  drvsigma = 2 * sigmak ;
  
  Y2 = arma::accu(pow(Xk, 2));  
  Y = sum(Xk, 1);
  
  grdl = gradienthyp(Y, drvl, nk, D, Z, A, sigmak);
  grda = gradienthyp(Y, drva, nk, D, Z, A, sigmak);
  grdsigma = gradientNoise(Y, drvsigma, nk, sigmak, Z , D, Y2);
  
  //return gradient of negative log like
  grad(0) = -grdl(0);
  grad(1) = -grda(0);
  grad(2) = -grdsigma(0);
  
  return(grad);
}

//' @description Leapfrog routine
//' 
// [[Rcpp::export]]
List LeapfrogGPcpp(arma::mat Xk,
                        arma::vec tau,
                        arma::vec p,
                        arma::vec x,
                        arma::vec m,
                        int nk,
                        int D,
                        int L,
                        double delta){
  
  arma::vec gradprior(3);
  gradprior(0) = x(0); gradprior(1) = 2*x(1); gradprior(2) = 2*x(2);

for (int t=0; t < L; t++){
  //half step for momenta
  p = p - delta * (gradientGPcpp(Xk, tau, x , nk, D) + gradprior)/2 ;
  //full step for position
  x = x + delta * p / m  ;
  //another half step for momenta
  p = p - delta * (gradientGPcpp(Xk, tau, x , nk, D) + gradprior)/2 ;
} 


return List::create(Rcpp::Named("p") = p,
                    Rcpp::Named("x") = x
                    );
  
}  
  