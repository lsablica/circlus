#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

 


double hyper2F1(const double b, const double c, const double x){
  double sum = 1.0;
  double element = 1.0;
  double i = 0.0;
  
  while(fabs(element/sum) > 1e-7){
    
    element *= (b+i) * x / (c+i);  
    sum +=  element;
    //Rcout << "element : " << element << "\n";
    //Rcout << "sum : " << sum << "\n";
    i += 1.0;
  }
  
  return sum;
} 


double n1d(const double d, const double x){
  double z = 4*x/((1+x)*(1+x));
  return 1 + 2*(1-z)*(1-hyper2F1(d/2, d, z))/z;
} 


double n1d_deriv(const double d, const double x){
  double z = 4*x/((1+x)*(1+x));
  double F = hyper2F1(d/2, d, z);
  return (1/2-1/(2*x*x))*(1-F) - (1-z)*((d/2+1)*hyper2F1(d/2+2, d+1, z) + d*(1-F)/z )/z;
} 

double hybridnewton(double d, double target, double tol = 1e-6, int maxiter = 100) {
  double x,a,b;
  a = 0;
  b = 1;
  x = 0.5;
  
  double g, dg, oldx = 1;
  int niter = 0;
  
  while(tol < std::abs(x-oldx) && niter < maxiter){
    oldx = x; 
    g = n1d(d,x) - target;
    dg = n1d_deriv(d,x);
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


arma::mat Moebius_S(arma::mat X, arma::vec mu, double rho){
  
  arma::mat Y = (1-rho*rho)*(X.each_row() + rho*mu.t());
  Y = Y.each_col()/(1+2*rho*X*mu+rho*rho);
  Y = Y.each_row() + rho*mu.t();
  
  return Y;
}

//' @title Random Sampling from Spherical Cauchy Distributions
//' @description Generates a random sample from spherical Cauchy distributions.
//' @param n The number of random draws.
//' @param rho A numeric value giving the rho parameter.
//' @param mu A numeric vector giving the mu direction parameter.
//' @return  A matrix with the generated values.
//' @rdname rspcauchy
//' @examples
//' rspcauchy(10, 0.95, c(1, 0, 0))
//' @export
// [[Rcpp::export]]
arma::mat rspcauchy(int n, double rho, arma::vec &mu){
  double norm = arma::as_scalar(sum(pow(mu,2)));
  int p = mu.n_elem;
  arma::mat A(n, p);
  A = normalise(A.randn(),2,1);
  if(rho == 0 || norm == 0){/*uniform*/
    return A;
  }
  return Moebius_S(A, mu, rho);
} 


// [[Rcpp::export]]
List M_step_sCauchy(const arma::mat &data, arma::vec weights,
                    int n, int d, double tol = 1e-6, int maxiter = 100){ 
  
  d = d - 1;
  weights =  arma::normalise(weights, 1);
  arma::vec weighted_means = data.t() * weights;
  int niter = 0;
  double norm, rho0, results_rho;
  arma::vec mu0, psi, psiold;
  arma::mat weighted_trans_data(n, d);
  psiold = 2*weighted_means;
  norm = arma::norm(weighted_means);
  
  mu0 = weighted_means/norm;
  rho0 = hybridnewton(d, norm, tol, maxiter);
  psi = rho0*mu0;
  
  while(arma::norm(psi-psiold, 2) > tol && niter < maxiter){
    psiold = psi;
    weighted_trans_data = Moebius_S(data, - mu0, rho0).t() * weights;
    psi = psiold + ((d+1)*(1-rho0*rho0)/(2*d))*weighted_trans_data; 
    rho0 = arma::norm(psi, 2);
    mu0 = psi/rho0;
    niter += 1;
  } 
  arma::rowvec results_mu = mu0.t();
  results_rho = rho0;
  return List::create(Named("mu") = results_mu, Named("rho") = results_rho);
  //Rcout << "rho_vector : " << results_rho << "\n";
  //Rcout << "results_mu : " << results_mu << "\n";
}  

// [[Rcpp::export]]
arma::vec logLik_sCauchy(const arma::mat &data, arma::vec mu_vec, double rho){ 
  
  double d = data.n_cols;
  return (d-1)*log(1-rho*rho) - (d-1)*arma::log(1 + rho*rho -2*rho*data*mu_vec); 
} 


