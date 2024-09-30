
#' @title Spherical Cauchy routine using neural networks for flexmix
#' @description  \code{SCauchyNN_clust} offers a flexmix routine for spherical Cauchy distribution using neural networks. 
#' @param formula formula
#' @param EPOCHS number of epochs in the M-step estimation (default: 1)
#' @param LR learning rate used in the M-steo estimation (default: 0.5)
#' @param max_iter maximum number of iterations of the LBFGS optimizer (default: 200)
#' @param line_search_fn method used for line search in LBFGS (default: "strong_wolfe")
#' @param free_iter number of initial iterations for which the model in M-step is fully reseted (default: 3)
#' @return Object of type FLXMC for flexmix estimation.
#' @rdname SCauchyNN_clust
#' @import flexmix
#' @import torch
#' @importFrom methods new
SCauchyNN_clust <- function(formula = .~. , EPOCHS = 1, LR = 0.5, max_iter = 200, line_search_fn = "strong_wolfe", free_iter = 5){
  retval <- new ("FLXMC" , weighted = TRUE , formula = formula , dist = " SCauchy " ,
                 name = " Spherical Cauchy - based clustering using neural networks")
  retval@defineComponent <- function (para, df) {
    NNmodel = para$NNmodel
    logLik <- function (x , y ) {
      X = torch_tensor(x)
      Y = torch_tensor(y)
      
      NNmodel$eval()
      with_no_grad({ 
        para_new <- NNmodel(X)
      })
      scauchy_log_likelihood(mu = para_new$mu , rho = para_new$rho, Y)
      #(mu = para$mu , rho = para$rho, Y)
    }
    predict <- function ( x ) {
      X = torch_tensor(x)
      NNmodel$eval()
      with_no_grad({ 
        para_new <- NNmodel(X)
      })
      para_new$mu <- torch::as_array(para_new$mu)
      para_new$rho <- torch::as_array(para_new$rho)
      para_new
    }
    new ("FLXcomponent" , parameters = list(mu = torch::as_array(para$mu) , rho = torch::as_array(para$rho), model = NNmodel),
         df = para$df , logLik = logLik , predict = predict )
  }
  retval@fit <- function (x , y , w , component) {
    
    iteration <- eval(quote(get("iter")),parent.frame(n=8))
    n <- nrow(y)
    d <- ncol(y)
    input_dim = ncol(x)
    output_dim = d
    EPOCHS = EPOCHS
    LR = LR
    Y = torch_tensor(y)
    X = torch_tensor(x)
    W = torch_tensor(matrix(w/sum(w), ncol = 1))
    
    if(iteration <= free_iter){
      component$model = Spherical(input_dim, output_dim)
    } 
    NNmodel = component$model
    NNmodel = Spherical(input_dim, output_dim)
    optimizer = optim_lbfgs(NNmodel$parameters, lr = LR, max_iter = max_iter, line_search_fn = line_search_fn)
    NNmodel$train()
    
    calc_loss <- function() {
      optimizer$zero_grad()
      res = NNmodel(X)
      loss = scauchy_weighted_neg_log_likelihood(res$mu, res$rho, Y, W)
      loss$backward()
      loss
    }
    for(epoch in seq_len(EPOCHS)){
      optimizer$step(calc_loss)
    }
    
    NNmodel$eval()
    with_no_grad({ 
      para <- NNmodel(X)
    })
    df <- (d+1)
    retval@defineComponent(c(para , df = df, NNmodel = NNmodel))
  }
  retval
}



#' @title PKBD routine using neural networks for flexmix
#' @description  \code{PKBDNN_clust} offers a flexmix routine for PKBD distribution using neural networks. 
#' @param formula formula
#' @param EPOCHS number of epochs in the M-step estimation (default: 1)
#' @param LR learning rate used in the M-steo estimation (default: 0.5)
#' @param max_iter maximum number of iterations of the LBFGS optimizer (default: 200)
#' @param line_search_fn method used for line search in LBFGS (default: "strong_wolfe")
#' @param free_iter number of initial iterations for which the model in M-step is fully reseted (default: 3)
#' @return Object of type FLXMC for flexmix estimation.
#' @rdname PKBDNN_clust
#' @import flexmix
#' @import torch
#' @importFrom methods new
PKBDNN_clust <- function(formula = .~. , EPOCHS = 1, LR = 0.1, max_iter = 200, line_search_fn = "strong_wolfe", free_iter = 3){
  retval <- new ("FLXMC" , weighted = TRUE , formula = formula , dist = " PKBD " ,
                 name = " PKBD - based clustering using neural networks")
  retval@defineComponent <- function (para, df) {
    NNmodel = para$NNmodel
    logLik <- function (x , y ) {
      X = torch_tensor(x)
      Y = torch_tensor(y)
      
      NNmodel$eval()
      with_no_grad({ 
        para_new <- NNmodel(X)
      })
      pkbd_log_likelihood(mu = para_new$mu , rho = para_new$rho, Y)
      #(mu = para$mu , rho = para$rho, Y)
    }
    predict <- function ( x ) {
      X = torch_tensor(x)
      NNmodel$eval()
      with_no_grad({ 
        para_new <- NNmodel(X)
      })
      para_new$mu <- torch::as_array(para_new$mu)
      para_new$rho <- torch::as_array(para_new$rho)
      para_new
    }
    new ("FLXcomponent" , parameters = list(mu = torch::as_array(para$mu) , rho = torch::as_array(para$rho), model = NNmodel),
         df = para$df , logLik = logLik , predict = predict)
  }
  
  
  retval@fit <- function (x , y , w , component) {
    
    iteration <- eval(quote(get("iter")),parent.frame(n=8))
    n <- nrow(y)
    d <- ncol(y)
    input_dim = ncol(x)
    output_dim = d
    EPOCHS = EPOCHS
    LR = LR
    Y = torch_tensor(y)
    X = torch_tensor(x)
    W = torch_tensor(matrix(w/sum(w), ncol = 1))
    
    
    if(iteration <= free_iter){
      component$model = Spherical(input_dim, output_dim)
    } 
    NNmodel = component$model
    optimizer = optim_lbfgs(NNmodel$parameters, lr = LR, max_iter = max_iter, line_search_fn = line_search_fn)
    NNmodel$train()
    
    calc_loss <- function() {
      optimizer$zero_grad()
      res = NNmodel(X)
      loss = pkbd_weighted_neg_log_likelihood(res$mu, res$rho, Y, W)
      loss$backward()
      loss
    }
    for(epoch in seq_len(EPOCHS)){
      optimizer$step(calc_loss)
    }
    
    NNmodel$eval()
    with_no_grad({ 
      para <- NNmodel(X)
    })
    df <- (d+1)
    retval@defineComponent(c(para , df = df, NNmodel = NNmodel))
  }
  retval
}



#' Flexible Mixture Component for PKBD Clustering
#'
#' @description
#' Provides a flexible mixture component for PKBD (Projected Kent-Based Distribution) clustering,
#' offering different estimation methods for use with the flexmix package.
#'
#' @param formula An R formula object specifying the model. Default is `.~.`.
#' @param method A character string specifying the estimation method. 
#'        Options are "NN" for a neural network-based approach and "direct" for standard estimation.
#' @param ... Additional arguments to fine-tune the clustering process.
#'
#' @return An object of class "FLXMC" for use in flexmix estimation.
#'
#' @details
#' FLXMCpkbd offers two different PKBD clustering methods:
#' \itemize{
#'   \item "NN": A neural network-based approach, which may offer more flexibility for complex data structures.
#'   \item "direct": A standard estimation method, typically faster and suitable for many common scenarios.
#' }
#' 
#' The choice of method can significantly affect clustering results and computation time.
#' Users are encouraged to experiment with both methods to determine which works best for their specific data and requirements.
#'
#' @examples
#' \dontrun{
#' library(flexmix)
#' 
#' # Using neural network method
#' model_nn <- flexmix(~., data=your_data, k=2, model=FLXMCpkbd(method="NN"))
#'
#' # Using direct method
#' model_direct <- flexmix(~., data=your_data, k=2, model=FLXMCpkbd(method="direct"))
#'
#' # Compare results
#' summary(model_nn)
#' summary(model_direct)
#' }
#'
#' @seealso
#' \code{\link[flexmix]{flexmix}} for the main clustering function in the flexmix package.
#'
#' @export
FLXMCpkbd <- function(formula = .~. , method = "direct", ...){
  if (method == "NN"){
    return(PKBDNN_clust_adam(formula = formula, ...))
  } else if (method == "direct"){
    return(PKBD_clust(formula = formula, ...))
  } else {
    stop("Invalid method specified. Please use 'NN' or 'direct'.")
  }
}


#' Flexible Mixture Component for Spherical Cauchy Clustering
#'
#' @description
#' Provides a flexible mixture component for spherical Cauchy distribution clustering,
#' offering different estimation methods for use with the flexmix package.
#'
#' @param formula An R formula object specifying the model. Default is `.~.`.
#' @param method A character string specifying the estimation method. 
#'        Options are "NN" for a neural network-based approach and "direct" for standard estimation.
#' @param ... Additional arguments to fine-tune the clustering process.
#'
#' @return An object of class "FLXMC" for use in flexmix estimation.
#'
#' @details
#' FLXMCsCauchy offers two different spherical Cauchy clustering methods:
#' \itemize{
#'   \item "NN": A neural network-based approach, which may offer more flexibility for complex data structures.
#'   \item "direct": A standard estimation method, typically faster and suitable for many common scenarios.
#' }
#' 
#' The choice of method can significantly affect clustering results and computation time.
#' Users are encouraged to experiment with both methods to determine which works best for their specific data and requirements.
#'
#' @examples
#' \dontrun{
#' library(flexmix)
#' 
#' # Using neural network method
#' model_nn <- flexmix(~., data=your_data, k=2, model=FLXMCsCauchy(method="NN"))
#'
#' # Using direct method
#' model_direct <- flexmix(~., data=your_data, k=2, model=FLXMCsCauchy(method="direct"))
#'
#' # Compare results
#' summary(model_nn)
#' summary(model_direct)
#' }
#'
#' @seealso
#' \code{\link[flexmix]{flexmix}} for the main clustering function in the flexmix package.
#' @export
FLXMCsCauchy <- function(formula = .~., method = c("NN", "direct"), ...) {
  method <- match.arg(method)
  
  if (method == "NN") {
    return(SCauchyNN_clust_adam(formula, ...))
  } else if (method == "direct") {
    return(SCauchy_clust(formula, ...))
  } else {
    stop("Invalid method specified. Use 'NN' or 'direct'.")
  }
}

# 
# /*
#   // [[Rcpp::export]]
# void M_step(const arma::mat &data, const arma::mat &beta_matrix, 
#             arma::mat mu_matrix, arma::vec rho_vector,
#             int k, int n, int d, double tol, int maxiter){ 
#   arma::rowvec alpha = sum(beta_matrix)/n;
#   
#   
#   arma::mat crossmat = data*mu_matrix;
#   arma::mat rho_mat(n, k);
#   rho_mat.each_row() = rho_vector.t();
#   arma::mat wscale_mat =  1 + pow(rho_mat, 2) - 2*rho_mat%crossmat;
#   arma::mat scaled_weight_matrix = beta_matrix/wscale_mat;
#   arma::mat mu = scaled_weight_matrix.t() * data; 
#   arma::vec mu_norms = arma::vecnorm(mu, 2, 1);
#   mu_matrix = mu.each_col()/mu_norms;
#   Rcout << "mu_matrix : " << mu_matrix << "\n";
#   arma::rowvec sums_scaled_weight_matrix = sum(scaled_weight_matrix, 0);
#   //standardize each
#   double c1, c2, c3;
#   for(int i = 0; i < k; i++){
#     c1 = d*sums_scaled_weight_matrix(i);
#     c2 = 2*n*alpha(i);
#     c3 = d*mu_norms(i);
#     rho_vector(i) = hybridnewton(c1,c2,c3,tol,maxiter); 
#   }
#   Rcout << "rho_vector : " << rho_vector << "\n";
# }
# */
  



# //old code
# /*
#   
#   
#   // [[Rcpp::export]]
# void M_step_sCauchy2(const arma::mat &data, const arma::mat &beta_matrix,
#                      int k, int n, int d, double tol = 1e-6, int maxiter = 100){ 
#   arma::rowvec sums = sum(beta_matrix);
#   arma::rowvec alpha = sums/n;
#   arma::mat weights =  beta_matrix.each_row()/sums;
#   arma::mat weighted_means = data.t() * weights;
#   
#   // we iterate over all clusters k, calculate method of moments estimate and use it for MLE estimate
#   int niter;
#   double norm, rho0;
#   arma::vec mu0, psi, psiold, w, results_rho(k);
#   arma::mat weighted_trans_data(n, d), results_mu(d, k);
#   for(int i = 0; i < k; i++){
#     niter = 0;
#     mu0 = weighted_means.col(i);
#     w = weights.col(i);
#     psiold = 2*mu0;
#     norm = arma::norm(mu0, 2);
#     mu0 = mu0/norm;
#     rho0 = hybridnewton(d, norm, tol = tol, maxiter = maxiter);
#     psi = rho0*mu0;
#     Rcout << "psi0 : " << psi << "\n";
#     
#     while(arma::norm(psi-psiold, 2) > tol && niter < maxiter){
#       psiold = psi;
#       weighted_trans_data = Moebius_S(data, - mu0, rho0).t() * w;
#       psi = psiold + ((d+1)*(1-rho0*rho0)/(2*d))*weighted_trans_data; 
#       Rcout << "psi : " << psi << "\n";
#       rho0 = arma::norm(psi, 2);
#       mu0 = psi/rho0;
#       niter += 1;
#     }
#     results_mu.col(i) = mu0;
#     results_rho(i) = rho0;
#   }
#   Rcout << "rho_vector : " << results_rho << "\n";
#   Rcout << "results_mu : " << results_mu << "\n";
# } 
# 
# 
# 
# double n1d(const double d, const double x){
#   return 1 + (1-x)*(1-x)*(1-hyper2F1(d/2, d, 4*x/((1+x)*(1+x))))/(2*x);
# }  
# 
# 
# 
# double hyperg_2F1_series(const double a, const double b, const double c, const double x){
#   double sum = 1.0;
#   double del = 1.0;
#   double del_prev;
#   double k = 0.0;
#   
#   while(fabs(del/sum) > 1e-7){
#     
#     del_prev = del;
#     del *= (a+k)*(b+k) * x / ((c+k) * (k+1.0));  
#     sum +=  del;
#     Rcout << "del : " << del << "\n";
#     Rcout << "sum : " << sum << "\n";
#     if (fabs(del_prev/sum) < 1e-7 && fabs(del/sum) < 1e-7) break;
#     k += 1.0;
#   }
#   
#   return sum;
# } 
# 
# 
# 
# 
# double n1d_deriv(const double d, const double x){
#   double z = 4*x/((1+x)*(1+x));
#   double F = hyper2F1(d/2, d, z);
#   return (1/2-1/(2*x*x))*(1-F) - (1-x)*(1-x)*((d/2+1)*hyper2F1(d/2+2, d+1, z) + d*(1-F)/z )/(4*x);
# } 
# 
# 
# 
# 
# */
#   
