#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

 

// [[Rcpp::export]]
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

// [[Rcpp::export]]
double n1d(const double d, const double x){
  double z = 4*x/((1+x)*(1+x));
  return 1 + 2*(1-z)*(1-hyper2F1(d/2, d, z))/z;
} 


// [[Rcpp::export]]
double n1d_deriv(const double d, const double x){
  double z = 4*x/((1+x)*(1+x));
  double F = hyper2F1(d/2, d, z);
  return (1/2-1/(2*x*x))*(1-F) - (1-z)*((d/2+1)*hyper2F1(d/2+2, d+1, z) + d*(1-F)/z )/z;
} 
// [[Rcpp::export]]
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

// [[Rcpp::export]]
arma::mat Moebius_S(arma::mat X, arma::vec mu, double rho){
  
  arma::mat Y = (1-rho*rho)*(X.each_row() + rho*mu.t());
  Y = Y.each_col()/(1+2*rho*X*mu+rho*rho);
  Y = Y.each_row() + rho*mu.t();
  
  return Y;
}

// [[Rcpp::export]]
arma::mat rsCauchy(int n, double rho, arma::vec &mu){
  double norm = as_scalar(sum(pow(mu,2)));
  int p = mu.n_elem;
  arma::mat A(n, p);
  A = normalise(A.randn(),2,1);
  if(rho == 0 || norm == 0){/*uniform*/
    return A;
  }
  return Moebius_S(A, mu, rho);
} 


// [[Rcpp::export]]
void M_step_sCauchy(const arma::mat &data, arma::vec weights,
                    int n, int d, double tol = 1e-6, int maxiter = 100){ 
  
  d = d - 1;
  weights =  arma::normalise(weights, 1);
  arma::vec weighted_means = data.t() * weights;
  int niter = 0;
  double norm, rho0, results_rho;
  arma::vec mu0, psi, psiold, results_mu;
  arma::mat weighted_trans_data(n, d);
  psiold = 2*weighted_means;
  norm = arma::norm(weighted_means);
  
  mu0 = weighted_means/norm;
  rho0 = hybridnewton(d, norm, tol = tol, maxiter = maxiter);
  psi = rho0*mu0;
  
  while(arma::norm(psi-psiold, 2) > tol && niter < maxiter){
    psiold = psi;
    weighted_trans_data = Moebius_S(data, - mu0, rho0).t() * weights;
    psi = psiold + ((d+1)*(1-rho0*rho0)/(2*d))*weighted_trans_data; 
    rho0 = arma::norm(psi, 2);
    mu0 = psi/rho0;
    niter += 1;
  } 
  results_mu = mu0;
  results_rho = rho0;
  Rcout << "rho_vector : " << results_rho << "\n";
  Rcout << "results_mu : " << results_mu << "\n";
}  

// [[Rcpp::export]]
arma::vec logLik_PKBD(const arma::mat &data, arma::vec mu_vec, double rho){ 
  
  double d = data.n_rows;
  return (d-1)*log(1-rho) - (d-1)*arma::log(1 + rho*rho -2*rho*data*mu_vec); 
} 


//old code
/*
 
 
 // [[Rcpp::export]]
 void M_step_sCauchy2(const arma::mat &data, const arma::mat &beta_matrix,
 int k, int n, int d, double tol = 1e-6, int maxiter = 100){ 
 arma::rowvec sums = sum(beta_matrix);
 arma::rowvec alpha = sums/n;
 arma::mat weights =  beta_matrix.each_row()/sums;
 arma::mat weighted_means = data.t() * weights;
 
 // we iterate over all clusters k, calculate method of moments estimate and use it for MLE estimate
 int niter;
 double norm, rho0;
 arma::vec mu0, psi, psiold, w, results_rho(k);
 arma::mat weighted_trans_data(n, d), results_mu(d, k);
 for(int i = 0; i < k; i++){
 niter = 0;
 mu0 = weighted_means.col(i);
 w = weights.col(i);
 psiold = 2*mu0;
 norm = arma::norm(mu0, 2);
 mu0 = mu0/norm;
 rho0 = hybridnewton(d, norm, tol = tol, maxiter = maxiter);
 psi = rho0*mu0;
 Rcout << "psi0 : " << psi << "\n";
 
 while(arma::norm(psi-psiold, 2) > tol && niter < maxiter){
 psiold = psi;
 weighted_trans_data = Moebius_S(data, - mu0, rho0).t() * w;
 psi = psiold + ((d+1)*(1-rho0*rho0)/(2*d))*weighted_trans_data; 
 Rcout << "psi : " << psi << "\n";
 rho0 = arma::norm(psi, 2);
 mu0 = psi/rho0;
 niter += 1;
 }
 results_mu.col(i) = mu0;
 results_rho(i) = rho0;
 }
 Rcout << "rho_vector : " << results_rho << "\n";
 Rcout << "results_mu : " << results_mu << "\n";
 } 
 
 
 
 double n1d(const double d, const double x){
 return 1 + (1-x)*(1-x)*(1-hyper2F1(d/2, d, 4*x/((1+x)*(1+x))))/(2*x);
 }  
 
 
 
 double hyperg_2F1_series(const double a, const double b, const double c, const double x){
 double sum = 1.0;
 double del = 1.0;
 double del_prev;
 double k = 0.0;
 
 while(fabs(del/sum) > 1e-7){
 
 del_prev = del;
 del *= (a+k)*(b+k) * x / ((c+k) * (k+1.0));  
sum +=  del;
Rcout << "del : " << del << "\n";
Rcout << "sum : " << sum << "\n";
if (fabs(del_prev/sum) < 1e-7 && fabs(del/sum) < 1e-7) break;
k += 1.0;
}

return sum;
} 
 

 
 
 double n1d_deriv(const double d, const double x){
 double z = 4*x/((1+x)*(1+x));
 double F = hyper2F1(d/2, d, z);
 return (1/2-1/(2*x*x))*(1-F) - (1-x)*(1-x)*((d/2+1)*hyper2F1(d/2+2, d+1, z) + d*(1-F)/z )/(4*x);
 } 
 
 
 
 
 */

