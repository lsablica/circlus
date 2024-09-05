#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "Rcpp.h"
using namespace Rcpp;


//' @title Random Sampling from PKBD Distributions using ACG Envelopes
//' @description \code{rPKBD_ACG} generates a random sample from PKBD distributions.
//' @param n number of random draws.
//' @param rho a numeric giving the concentration parameter.
//' @param mu a numeric vector giving the mean direction parameter.
//' @return  A vector the generated values.
//' @details The function generates samples from PKBD using ACG envelopes
//' @rdname rPKBD_ACG
//' @useDynLib circlus
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
arma::mat rPKBD_ACG(int n, double rho, arma::vec &mu){
  double lambda = 2*rho/(1+rho*rho);
  double norm = as_scalar(sum(pow(mu,2)));
  int p = mu.n_elem;
  arma::mat A(n, p);
  if(lambda == 0 || norm == 0){/*uniform*/
    return normalise(A.randn(),2,1);
  }
  mu = mu/sqrt(norm);
  int count = 0;
  int Nt = 0;
  double unif, mACG, PKBD, mutz, ratio, mutx;
  arma::vec candidate;
  
  double pp = (double)p;
  arma::vec coe = { -4*(pp-1) , 4*pp-lambda*lambda*(pp-2)*(pp-2), 2*pp*(pp-2)*lambda*lambda, -pp*pp*lambda*lambda};
  arma::vec RO = arma::sort(arma::real(arma::roots(coe)));
  double b = RO(1);
  
  double minuslogM = log((1+sqrt(1-lambda*lambda/b))/2);
  double b2 = -1 + sqrt(1/(1-b));
  double b1 = b/(1-b);  
  
  while(count<n){
    candidate = arma::randn<arma::vec>(p);
    mutz = arma::dot(mu, candidate) ;
    norm = sqrt(arma::dot(candidate,candidate) + b1*mutz*mutz);
    mutx = mutz*(1+b2)/norm ;  
    PKBD = -log(1-lambda*mutx);
    mACG =  log(1-b*mutx*mutx);
    unif = arma::randu<double>();
    ratio = 0.5*p*(PKBD + mACG + minuslogM);
    if(log(unif)<ratio){
      candidate = (candidate + b2*mutz*mu)/norm;
      A.row(count) = arma::trans(candidate);
      count += 1;
    }
    Nt += 1;
    if(Nt % 1000000 == 0) Rcpp::checkUserInterrupt();
  }
  return A;
}
