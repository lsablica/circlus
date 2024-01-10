#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double hybridnewton(double c1, double c2, double c3, double tol = 1e-6, int maxiter = 100) {
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
void M_step(const arma::mat &data, const arma::mat &beta_matrix, 
            arma::mat &mu_matrix, arma::vec &rho_vector,
            int k, int n, int d, double tol, int maxiter){ 
  arma::rowvec alpha = sum(beta_matrix)/n;
  
      
  arma::mat crossmat = data*mu_matrix;
  arma::mat rho_mat(n, k);
  rho_mat.each_row() = rho_vector.t();
  arma::mat wscale_mat =  1 + pow(rho_mat, 2) - 2*rho_mat%crossmat;
  arma::mat scaled_weight_matrix = beta_matrix/wscale_mat;
  arma::mat mu = scaled_weight_matrix.t() * data; 
  arma::vec mu_norms = arma::vecnorm(mu, 2, 1);
  mu_matrix = mu.each_col()/mu_norms;
  Rcout << "mu_matrix : " << mu_matrix << "\n";
  arma::rowvec sums_scaled_weight_matrix = sum(scaled_weight_matrix, 0);
  //standardize each
  double c1, c2, c3;
  for(int i = 0; i < k; i++){
    c1 = d*sums_scaled_weight_matrix(i);
    c2 = 2*n*alpha(i);
    c3 = d*mu_norms(i);
    rho_vector(i) = hybridnewton(c1,c2,c3,tol,maxiter); 
  }
  Rcout << "rho_vector : " << rho_vector << "\n";
}

