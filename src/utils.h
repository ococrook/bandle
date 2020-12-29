#include <RcppArmadillo.h>
using namespace Rcpp;

List trenchDetcpp(arma::vec c);
arma::mat trenchInvcpp(arma::vec v);
arma::vec gradientNoise(arma::vec Y,
                        double drv,
                        double nk,
                        double sigmak,
                        arma::mat Z,
                        double D,
                        double Y2);