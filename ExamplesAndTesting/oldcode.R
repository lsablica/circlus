
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