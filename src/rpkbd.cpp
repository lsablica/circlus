#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "Rcpp.h"
using namespace Rcpp;

#ifdef __cplusplus
extern "C" {
#endif
  
#include "rdist.h"
  
#ifdef __cplusplus
}
#endif


NumericVector Tinflexsampler_sampler_from_c(int n,
                                            double lambda,
                                            double d,
                                            double cT,
                                            double rho) {
  NumericVector sexp_params = {lambda,d};
  NumericVector sexp_c = {cT};
  NumericVector sexp_rho = {rho};
  IntegerVector sexp_n = {n};
  NumericVector sexp_ib = {0,1};
  NumericVector sexp_max_intervals = {1001};
  SEXP ab = Tinflexsampler_sampler2(sexp_n, sexp_params, sexp_ib, sexp_c,
                                   sexp_rho, sexp_max_intervals);
  Rcpp::NumericVector result( ab );
  return result;
}

//' @title Random Sampling from PKBD Distributions
//' @description \code{rmwat} generates a random sample from PKBD distributions.
//' @param lambda a numeric giving the concentration parameter.
//' @param mu a numeric vector giving the mu parameter.
//' @return  A vector the generated values.
//' @details The function generates samples from PKBD using Tinflex package
//' @rdname rpkbd_Tinflex
//' @export
// [[Rcpp::export]]
arma::mat rpkbd_Tinflex(int n, double lambda, arma::vec &mu, double cT, double rho){
  double norm = as_scalar(sum(pow(mu,2)));
  int p = mu.n_elem;
  arma::mat A(n, p, arma::fill::randn);
  if(lambda == 0 || norm == 0){/*uniform*/
    return normalise(A,2,1);
  }
  mu = mu/sqrt(norm);
  NumericVector v = Tinflexsampler_sampler_from_c(n, lambda, p, cT, rho);
  arma::vec w = as<arma::vec>(wrap(v));
  arma::vec choice = {-1, 1};
  arma::vec index = RcppArmadillo::sample(choice, n, true); 
  w = w % index;
  A = A - A*mu*mu.t() ;
  A = arma::normalise(A, 2, 1);
  A = A.each_col() % sqrt(1 - w%w);
  A = w*mu.t() + A;
  return A;
} 



//' @title Random Sampling from PKBD Distributions
//' @description \code{rmwat} generates a random sample from PKBD distributions.
//' @param lambda a numeric giving the concentration parameter.
//' @param mu a numeric vector giving the mu parameter.
//' @return  A vector the generated values.
//' @details The function generates samples from PKBD using ACG envelopes
//' @rdname rpkbd_ACG
//' @export
// [[Rcpp::export]]
arma::mat rpkbd_ACG(int n, double lambda, arma::vec &mu){
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
