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




//old code
/*
 
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

