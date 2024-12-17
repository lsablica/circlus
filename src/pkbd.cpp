#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;


double hybridnewton(double c1, double c2, double c3, double tol, int maxiter) {
  double x,a,b;
  a = 0;
  b = 1;
  x = 0.5;
  
  double g, dg, oldx = 1;
  int niter = 0;
  
  while(tol < std::abs(x-oldx) && niter < maxiter){
    oldx = x; 
    g = x*c1+x*c2/(1-x*x)-c3;
    dg = c1+c2*(1+x*x)/pow(1-x*x,2);
    x -= g/dg;
    if(x < a || x > b){
      if(g>0){
        b = oldx;
      }
      else{
        a = oldx;
      }
      x= (b+a)/2;
    }
    niter +=1;
  }
  return x;
}

// [[Rcpp::export]]
List M_step_PKBD(const arma::mat &data, arma::vec weights, arma::vec mu_vec, double rho,
                 int n, int d, double tol = 1e-6, int maxiter = 100){ 
  double alpha = arma::mean(weights);
  
  arma::vec crossmat = data*mu_vec;
  arma::vec wscale =  1 + pow(rho, 2) - 2*rho*crossmat;
  arma::vec scaled_weight = weights/wscale;
  arma::vec mu = data.t() * scaled_weight; 
  double mu_norm = arma::vecnorm(mu);
  mu_vec = mu/mu_norm;
  double sums_scaled_weight = sum(scaled_weight);
  rho = hybridnewton(d*sums_scaled_weight, 2*n*alpha, d*mu_norm, tol, maxiter); 
  return List::create(Named("mu") = mu_vec, Named("rho") = rho);
  //Rcout << "rho: " << rho << "\n";
} 


// [[Rcpp::export]]
arma::vec logLik_PKBD(const arma::mat &data, arma::vec mu_vec, double rho){ 
  
  double d = data.n_cols;
  return log(1-rho*rho) - d*arma::log(1 + rho*rho -2*rho*data*mu_vec)/2; 
} 



