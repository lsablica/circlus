#' @title Spherical Cauchy routine for flexmix
#' @description  \code{SCauchy_clust} offers a flexmix routine for spherical Cauchy distribution. 
#' @param formula formula
#' @return  Object of type FLXcomponent for flexmix estimation.
#' @rdname SCauchy_clust
#' @import flexmix
#' @importFrom methods new
#' @export
SCauchy_clust <- function (formula = .~.){
  retval <- new ("FLXMC", weighted = TRUE ,
                 formula = formula , dist = " Scauchy " ,
                 name = " Spherical Cauchy - based clustering ")
  retval@defineComponent <- function (para, df) {
    logLik <- function (x, y){
      logLik_sCauchy(y , mu_vec = para$mu , rho = para$rho)
    }
    predict <- function(x){
      para$mu
    }
    new ("FLXcomponent" , parameters = list(mu = para$mu, rho = para$rho),
         df = para$df , logLik = logLik , predict = predict)
  }
  retval@fit <- function (x , y , w , ...) {
    n <- nrow(y)
    d <- ncol(y)
    para <- M_step_sCauchy(y, w, n, d)
    df <- (d+1)
    retval@defineComponent(c(para, df = df))
  }
  retval
}

#################################################################################################

#' @title PKBD routine for flexmix
#' @description  \code{PKBD_clust} offers a flexmix routine for PKBD distribution. 
#' @param formula formula
#' @return Object of type FLXcomponent for flexmix estimation.
#' @rdname PKBD_clust
#' @import flexmix
#' @importFrom methods new
#' @importFrom stats runif
#' @export
PKBD_clust <- function (formula = .~.){
  retval <- new ("FLXMC" , weighted = TRUE ,
                 formula = formula , dist = " PKBD " ,
                 name = " PKBD - based clustering ")
  retval@defineComponent <- function (para, df) {
    logLik <- function (x , y ) {
      logLik_PKBD(y , mu_vec = para$mu , rho = para$rho)
    }
    predict <- function ( x ) {
      para$mu
    }
    new ("FLXcomponent" , parameters = list( mu = para$mu , rho = para$rho ),
         df = para$df , logLik = logLik , predict = predict )
  }
  retval@fit <- function (x , y , w , component) {
    n <- nrow(y)
    d <- ncol(y)
    if(length(component)==0){
      component <- list(mu = rep(0,d), rho = runif(1,0.7,0.95)) 
    } 
    para <- M_step_PKBD(y, w, component$mu, component$rho, n, d)
    df <- (d+1)
    retval@defineComponent(c( para , df = df))
  }
  retval
}

#################################################################################################


Spherical <- nn_module(
  "Spherical",
  initialize = function(input_dim, output_dim) {
    self$fc = nn_linear(input_dim, output_dim, bias = FALSE)
    self$output_dim = output_dim
  },
  forward = function(x) {
    mu = self$fc(x)
    normm = mu$norm(dim=-1, keepdim=TRUE)
    rho = normm / (1 + normm)
    mu = mu / normm
    list(mu = mu, rho = rho)
  }
)

scauchy_log_likelihood <- function(mu, rho, Y){
  d = Y$shape[2]
  term1 = (1-rho^2)$log()
  term2 = 1 + rho^2 - 2* rho*((mu$unsqueeze(2)$matmul(Y$unsqueeze(3)))$squeeze(3))
  log_likelihood = (d-1)*term1 - (d-1)*term2$log()
  return(as_array(log_likelihood))
}

scauchy_weighted_neg_log_likelihood <- function(mu, rho, Y, W){
  d = Y$shape[2]
  term1 = (1-rho^2)$log()
  term2 = 1 + rho^2 - 2* rho*((mu$unsqueeze(2)$matmul(Y$unsqueeze(3)))$squeeze(3))
  neg_log_likelihood = (d-1)*term2$log() - (d-1)*term1
  result = (neg_log_likelihood * W)$sum()
  return(result)
}

#' @title Spherical Cauchy routine using neural networks for flexmix
#' @description  \code{SCauchyNN_clust} offers a flexmix routine for spherical Cauchy distribution using neural networks. 
#' @param formula formula
#' @param EPOCHS number of epochs in the M-step estimation (default: 1)
#' @param LR learning rate used in the M-steo estimation (default: 0.5)
#' @param max_iter maximum number of iterations of the LBFGS optimizer (default: 200)
#' @param line_search_fn method used for line search in LBFGS (default: "strong_wolfe")
#' @param free_iter number of initial iterations for which the model in M-step is fully reseted (default: 3)
#' @return Object of type FLXcomponent for flexmix estimation.
#' @rdname SCauchyNN_clust
#' @import flexmix
#' @import torch
#' @importFrom methods new
#' @export
SCauchyNN_clust <- function(formula = .~. , EPOCHS = 1, LR = 0.5, max_iter = 200, line_search_fn = "strong_wolfe"){
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
      scauchy_log_likelihood(mu = para_new$mu , rho = para_new$rho, Y)
      #(mu = para$mu , rho = para$rho, Y)
    }
    predict <- function ( x ) {
      X = torch_tensor(x)
      NNmodel$eval()
      with_no_grad({ 
        para_new <- NNmodel(X)
      })
      para_new$mu <- torch::as_array(para$mu)
      para_new$rho <- torch::as_array(para$rho)
      para_new
    }
    new ("FLXcomponent" , parameters = list(mu = torch::as_array(para$mu) , rho = torch::as_array(para$rho), model = NNmodel),
         df = para$df , logLik = logLik , predict = predict )
  }
  retval@fit <- function (x , y , w , ...) {
    
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


#################################################################################################

pkbd_log_likelihood <- function(mu, rho, Y){
  d = Y$shape[2]
  term1 = (1-rho^2)$log()
  term2 = 1 + rho^2 - 2* rho*((mu$unsqueeze(2)$matmul(Y$unsqueeze(3)))$squeeze(3))
  log_likelihood = term1 - (d/2)*term2$log()
  return(as_array(log_likelihood))
}

pkbd_weighted_neg_log_likelihood <- function(mu, rho, Y, W){
  d = Y$shape[2]
  term1 = (1-rho^2)$log()
  term2 = 1 + rho^2 - 2* rho*((mu$unsqueeze(2)$matmul(Y$unsqueeze(3)))$squeeze(3))
  neg_log_likelihood = (d/2)*term2$log() - term1
  result = (neg_log_likelihood * W)$sum()
  return(result)
}

#' @title PKBD routine using neural networks for flexmix
#' @description  \code{PKBDNN_clust} offers a flexmix routine for PKBD distribution using neural networks. 
#' @param formula formula
#' @param EPOCHS number of epochs in the M-step estimation (default: 1)
#' @param LR learning rate used in the M-steo estimation (default: 0.5)
#' @param max_iter maximum number of iterations of the LBFGS optimizer (default: 200)
#' @param line_search_fn method used for line search in LBFGS (default: "strong_wolfe")
#' @param free_iter number of initial iterations for which the model in M-step is fully reseted (default: 3)
#' @return Object of type FLXcomponent for flexmix estimation.
#' @rdname PKBDNN_clust
#' @import flexmix
#' @import torch
#' @importFrom methods new
#' @export
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
      para_new$mu <- torch::as_array(para$mu)
      para_new$rho <- torch::as_array(para$rho)
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

#######################################################################################


#' @title PKBD routine using neural networks for flexmix
#' @description  \code{PKBDNN_clust_adam} offers a flexmix routine for PKBD distribution using neural networks and adam optimizer. 
#' @param formula formula
#' @param EPOCHS number of epochs in the M-step estimation (default: 1)
#' @param LR learning rate used in the M-steo estimation (default: 0.5)
#' @param free_iter number of initial iterations for which the model in M-step is fully reseted (default: 3)
#' @return Object of type FLXcomponent for flexmix estimation.
#' @rdname PKBDNN_clust_adam
#' @export
PKBDNN_clust_adam <- function(formula = .~. , EPOCHS = 100, LR = 0.1, max_iter = 200, 
                              line_search_fn = "strong_wolfe", free_iter = 5){
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
      para_new$mu <- torch::as_array(para$mu)
      para_new$rho <- torch::as_array(para$rho)
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
      print("adam")
      NNmodel = Spherical(input_dim, output_dim)
      optimizer = optim_adam(NNmodel$parameters, lr = LR)
      NNmodel$train()
      for(epoch in seq_len(EPOCHS)){
        optimizer$zero_grad()
        res = NNmodel(X)
        loss = pkbd_weighted_neg_log_likelihood(res$mu, res$rho, Y, W)
        loss$backward()
        optimizer$step()
      }
      para <- res  
    } else{
      print("lbfgs")
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
      optimizer$step(calc_loss)
      
      NNmodel$eval()
      with_no_grad({ 
        para <- NNmodel(X)
      })
      
    }
    df <- (d+1)
    retval@defineComponent(c(para , df = df, NNmodel = NNmodel))
  }
  retval
}

#######################################################################################



#' @title PKBD routine using neural networks for flexmix
#' @description  \code{PKBDNN_clust_adam} offers a flexmix routine for PKBD distribution using neural networks and adam optimizer. 
#' @param formula formula
#' @param EPOCHS number of epochs in the M-step estimation (default: 1)
#' @param LR learning rate used in the M-steo estimation (default: 0.5)
#' @param free_iter number of initial iterations for which the model in M-step is fully reseted (default: 3)
#' @return Object of type FLXcomponent for flexmix estimation.
#' @rdname PKBDNN_clust_adam
#' @export
PKBDNN_clust_adam <- function(formula = .~. , EPOCHS = 100, LR = 0.1, max_iter = 200, 
                              line_search_fn = "strong_wolfe", free_iter = 5){
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
      para_new$mu <- torch::as_array(para$mu)
      para_new$rho <- torch::as_array(para$rho)
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
      print("adam")
      NNmodel = Spherical(input_dim, output_dim)
      optimizer = optim_adam(NNmodel$parameters, lr = LR)
      NNmodel$train()
      for(epoch in seq_len(EPOCHS)){
        optimizer$zero_grad()
        res = NNmodel(X)
        loss = pkbd_weighted_neg_log_likelihood(res$mu, res$rho, Y, W)
        loss$backward()
        optimizer$step()
      }
      para <- res  
    } else{
      print("lbfgs")
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
      optimizer$step(calc_loss)
      
      NNmodel$eval()
      with_no_grad({ 
        para <- NNmodel(X)
      })
      
    }
    df <- (d+1)
    retval@defineComponent(c(para , df = df, NNmodel = NNmodel))
  }
  retval
}

#######################################################################################