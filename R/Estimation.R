#' @title Spherical Cauchy Driver for FlexMix
#' @description This model driver for flexmix implements model-based
#'     clustering of spherical Cauchy distributions.
#' @param formula A formula.
#' @return Returns an object of class `FLXMC`.
#' @rdname FLXMCspcauchy
#' @import flexmix
#' @importFrom methods new
#' @examples
#' mix <- rbind(rpkbd(30, 0.95, c(1, 0, 0)), rpkbd(30, 0.9, c(-1, 0, 0)))
#' m1 <- flexmix::flexmix(mix ~ 1, k = 2, model = FLXMCspcauchy())
#' @export
FLXMCspcauchy <- function(formula = .~.) {
  retval <- new("FLXMC", weighted = TRUE,
                formula = formula, dist = "SCauchy",
                name = "Spherical Cauchy-based clustering")
  retval@defineComponent <- function(para) {
    logLik <- function(x, y) {
      logLik_sCauchy(y, mu_vec = para$mu, rho = para$rho)
    }
    predict <- function(x) {
      para$mu
    }
    new("FLXcomponent", parameters = list(mu = para$mu, rho = para$rho),
        df = para$df, logLik = logLik, predict = predict)
  }
  retval@preproc.y <- function(y){
    norms <- rowSums(y^2)
    return(y/sqrt(norms))
  }
  
  retval@preproc.x <- function(x){
    if (ncol(x) > 1) 
      stop(paste("for the FLXMCspcauchy x must be univariate, use FLXMRspcauchy for problems with covariates"))
    x
  }
  retval@fit <- function(x, y, w, ...) {
    n <- nrow(y)
    d <- ncol(y)
    para <- M_step_sCauchy(y, w, n, d)
    para$df <- d
    retval@defineComponent(para)
  }
  retval
}

#################################################################################################

#' @title PKBD Driver for FlexMix
#' @description This model driver for flexmix implements model-based
#'     clustering of PKBD distributions.
#' @param formula A formula.
#' @return Returns an object of class `FLXMC`.
#' @rdname FLXMCpkbd
#' @import flexmix
#' @importFrom methods new
#' @importFrom stats runif
#' @examples
#' mix <- rbind(rpkbd(30, 0.95, c(1, 0, 0)), rpkbd(30, 0.9, c(-1, 0, 0)))
#' m1 <- flexmix::flexmix(mix ~ 1, k = 2, model = FLXMCpkbd())
#' @export
FLXMCpkbd <- function(formula = .~.){
  retval <- new("FLXMC", weighted = TRUE,
                formula = formula, dist = "PKBD",
                name = "PKBD-based clustering")
  retval@defineComponent <- function(para) {
    logLik <- function(x, y) {
      logLik_PKBD(y, mu_vec = para$mu, rho = para$rho)
    }
    predict <- function(x) {
      para$mu
    }
    new("FLXcomponent", parameters = list(mu = para$mu, rho = para$rho),
         df = para$df, logLik = logLik, predict = predict)
  }
  retval@preproc.y <- function(y){
    norms <- rowSums(y^2)
    return(y/sqrt(norms))
  }
  
  retval@preproc.x <- function(x){
    if (ncol(x) > 1) 
      stop(paste("for the FLXMCpkbd x must be univariate, use FLXMRpkbd for problems with covariates"))
    x
  }
  retval@fit <- function(x, y, w, component) {
    n <- nrow(y)
    d <- ncol(y)
    if (length(component) == 0) {
      component <- list(mu = rep(0, d), rho = runif(1, 0.7, 0.95)) 
    } 
    para <- M_step_PKBD(y, w, component$mu, component$rho, n, d)
    para$df <- d
    retval@defineComponent(para)
  }
  retval
}

#################################################################################################

Spherical <- nn_module(
  "Spherical",
  initialize = function(input_dim, output_dim) {
    self$fc <- nn_linear(input_dim, output_dim, bias = FALSE)
    self$output_dim <- output_dim
  },
  forward = function(x) {
    mu <- self$fc(x)
    normm <- mu$norm(dim = -1, keepdim = TRUE)
    rho <- normm /(1 + normm)
    mu <- mu / normm
    list(mu = mu, rho = rho)
  }
)

scauchy_log_likelihood <- function(mu, rho, Y) {
  d <- Y$shape[2]
  term1 <- (1-rho^2)$log()
  term2 <- 1 + rho^2 - 2 * rho * ((mu$unsqueeze(2)$matmul(Y$unsqueeze(3)))$squeeze(3))
  log_likelihood <- (d-1) * term1 - (d-1) * term2$log()
  return(as_array(log_likelihood))
}

scauchy_weighted_neg_log_likelihood <- function(mu, rho, Y, W){
  d <- Y$shape[2]
  term1 <- (1-rho^2)$log()
  term2 <- 1 + rho^2 - 2 * rho * ((mu$unsqueeze(2)$matmul(Y$unsqueeze(3)))$squeeze(3))
  neg_log_likelihood <- (d-1) * term2$log() - (d-1) * term1
  result <- (neg_log_likelihood * W)$sum()
  return(result)
}

#' @title Spherical Cauchy Driver for FlexMix Using Neural Networks
#' @description This model driver for flexmix implements model-based
#'     clustering of spherical Cauchy distributions using neural
#'     networks in the M-step.
#' @param formula A formula.
#' @param EPOCHS EPOCHS The number of epochs in the M-step estimation (default: 100).
#' @param LR The learning rate used in the M-steo estimation (default: 0.1).
#' @param max_iter The maximum number of iterations of the LBFGS optimizer (default: 200).
#' @param adam_iter The number of iteration for which the adam optimizer is used before the algorithm switches to L-BFGS (default: 5).
#' @param free_iter The number of initial iterations for which the model in M-step is fully reseted (default: adam_iter).
#' @param line_search_fn The method used for line search in LBFGS (default: "strong_wolfe").
#' @return Returns an object of class `FLXMC`.
#' @examples
#' \donttest{
#' if(torch::torch_is_installed()){
#' mix <- rbind(rpkbd(30, 0.95, c(1, 0, 0)), rpkbd(30, 0.9, c(-1, 0, 0)))
#' m1 <- flexmix::flexmix(mix ~ 1, k = 2, model = FLXMRspcauchy())
#' }
#' }
#' @rdname FLXMRspcauchy
#' @import flexmix
#' @import torch
#' @importFrom methods new
#' @export
FLXMRspcauchy <- function(formula = .~., EPOCHS = 100, LR = 0.1, max_iter = 200, 
                                 adam_iter = 5, free_iter = adam_iter, line_search_fn = "strong_wolfe"){
  retval <- new ("FLXMC", weighted = TRUE, formula = formula, dist = "PKBD",
                 name = "Spherical Cauchy-based clustering using neural networks")
  retval@defineComponent <- function(para, df) {
    NNmodel = para$NNmodel
    logLik <- function(x, y) {
      X = torch_tensor(x)
      Y = torch_tensor(y)
      
      NNmodel$eval()
      with_no_grad({ 
        para_new <- NNmodel(X)
      })
      scauchy_log_likelihood(mu = para_new$mu, rho = para_new$rho, Y)
      #(mu = para$mu, rho = para$rho, Y)
    }
    predict <- function(x) {
      X = torch_tensor(x)
      NNmodel$eval()
      with_no_grad({ 
        para_new <- NNmodel(X)
      })
      para_new$mu <- torch::as_array(para_new$mu)
      para_new$rho <- torch::as_array(para_new$rho)
      return(para_new)
    }
    new("FLXcomponent", parameters = list(mu = torch::as_array(para$mu),
                                          rho = torch::as_array(para$rho),
                                          model = NNmodel),
         df = para$df, logLik = logLik, predict = predict)
  }
  
  retval@preproc.y <- function(y){
    norms <- rowSums(y^2)
    return(y/sqrt(norms))
  }
  
  retval@fit <- function(x, y, w, component) {
    iteration <- eval(quote(get("iter")), parent.frame(n = 8))
    n <- nrow(y)
    d <- ncol(y)
    input_dim <- ncol(x)
    output_dim <- d
    EPOCHS <- EPOCHS
    LR <- LR
    Y <- torch_tensor(y)
    X <- torch_tensor(x)
    W <- torch_tensor(matrix(w/sum(w), ncol = 1))
    
    
    if (iteration <= adam_iter) {
      if (iteration <= free_iter) {
        component$model <- Spherical(input_dim, output_dim)
      } 
      NNmodel <- component$model
      optimizer <- optim_adam(NNmodel$parameters, lr = LR)
      NNmodel$train()
      for (epoch in seq_len(EPOCHS)) {
        optimizer$zero_grad()
        res <- NNmodel(X)
        loss <- scauchy_weighted_neg_log_likelihood(res$mu, res$rho, Y, W)
        loss$backward()
        optimizer$step()
      }
      para <- res  
    } else{
      if(iteration <= free_iter){
        component$model <- Spherical(input_dim, output_dim)
      } 
      NNmodel <- component$model
      optimizer <- optim_lbfgs(NNmodel$parameters, lr = LR,
                               max_iter = max_iter,
                               line_search_fn = line_search_fn)
      NNmodel$train()
      
      calc_loss <- function() {
        optimizer$zero_grad()
        res <- NNmodel(X)
        loss <- scauchy_weighted_neg_log_likelihood(res$mu, res$rho, Y, W)
        loss$backward()
        loss
      }
      optimizer$step(calc_loss)
      
      NNmodel$eval()
      with_no_grad({ 
        para <- NNmodel(X)
      })
      
    }
    para$df <- d*input_dim 
    para$NNmodel <- NNmodel
    retval@defineComponent(para)
  }
  retval
}


#################################################################################################

pkbd_log_likelihood <- function(mu, rho, Y) {
  d <- Y$shape[2]
  term1 <- (1-rho^2)$log()
  term2 <- 1 + rho^2 - 2 * rho * ((mu$unsqueeze(2)$matmul(Y$unsqueeze(3)))$squeeze(3))
  log_likelihood <- term1 - (d/2) * term2$log()
  return(as_array(log_likelihood))
}

pkbd_weighted_neg_log_likelihood <- function(mu, rho, Y, W) {
  d <- Y$shape[2]
  term1 <- (1-rho^2)$log()
  term2 <- 1 + rho^2 - 2 * rho * ((mu$unsqueeze(2)$matmul(Y$unsqueeze(3)))$squeeze(3))
  neg_log_likelihood <- (d/2) * term2$log() - term1
  result <- (neg_log_likelihood * W)$sum()
  return(result)
}


#######################################################################################


#' @title PKBD Driver for FlexMix Using Neural Networks
#' @description This model driver for flexmix implements model-based
#'     clustering of PKBD distributions using neural network in the M-step.
#' @param formula A formula.
#' @param EPOCHS The number of epochs in the M-step estimation (default: 100).
#' @param LR The learning rate used in the M-step estimation (default: 0.1).
#' @param max_iter The maximum number of iterations of the LBFGS optimizer (default: 200).
#' @param adam_iter The number of iteration for which the adam optimizer is used before the algorithm switches to L-BFGS (default: 5).
#' @param free_iter The number of initial iterations for which the model in M-step is fully reseted (default: adam_iter).
#' @param line_search_fn The method used for line search in LBFGS (default: "strong_wolfe").
#' @return Returns an object of class `FLXMC`.
#' @examples
#' \donttest{
#' if(torch::torch_is_installed()){
#' mix <- rbind(rpkbd(30, 0.95, c(1, 0, 0)), rpkbd(30, 0.9, c(-1, 0, 0)))
#' m1 <- flexmix::flexmix(mix ~ 1, k = 2, model = FLXMRpkbd())
#' }
#' }
#' @rdname FLXMRpkbd
#' @import flexmix
#' @import torch
#' @importFrom methods new
#' @export
FLXMRpkbd <- function(formula = .~., EPOCHS = 100, LR = 0.1, max_iter = 200, 
                              adam_iter = 5, free_iter = adam_iter, line_search_fn = "strong_wolfe"){
  retval <- new ("FLXMC", weighted = TRUE, formula = formula, dist = "PKBD",
                 name = " PKBD-based clustering using neural networks")
  retval@defineComponent <- function(para) {
    NNmodel <- para$NNmodel
    logLik <- function(x, y) {
      X <- torch_tensor(x)
      Y <- torch_tensor(y)
      
      NNmodel$eval()
      with_no_grad({ 
        para_new <- NNmodel(X)
      })
      pkbd_log_likelihood(mu = para_new$mu, rho = para_new$rho, Y)
    }
    predict <- function(x) {
      X <- torch_tensor(x)
      NNmodel$eval()
      with_no_grad({ 
        para_new <- NNmodel(X)
      })
      para_new$mu <- torch::as_array(para_new$mu)
      para_new$rho <- torch::as_array(para_new$rho)
      para_new
    }
    new("FLXcomponent", parameters = list(mu = torch::as_array(para$mu),
                                          rho = torch::as_array(para$rho),
                                          model = NNmodel),
         df = para$df, logLik = logLik, predict = predict)
  }
  
  retval@preproc.y <- function(y){
    norms <- rowSums(y^2)
    return(y/sqrt(norms))
  }
  
  retval@fit <- function(x, y, w, component) {
    
    iteration <- eval(quote(get("iter")), parent.frame(n = 8))
    n <- nrow(y)
    d <- ncol(y)
    input_dim <- ncol(x)
    output_dim <- d
    EPOCHS <- EPOCHS
    LR <- LR
    Y <- torch_tensor(y)
    X <- torch_tensor(x)
    W <- torch_tensor(matrix(w/sum(w), ncol = 1))
    
    
    if (iteration <= adam_iter) {
      if (iteration <= free_iter) {
        component$model <- Spherical(input_dim, output_dim)
      } 
      NNmodel <- component$model
      optimizer <- optim_adam(NNmodel$parameters, lr = LR)
      NNmodel$train()
      for (epoch in seq_len(EPOCHS)) {
        optimizer$zero_grad()
        res <- NNmodel(X)
        loss <- pkbd_weighted_neg_log_likelihood(res$mu, res$rho, Y, W)
        loss$backward()
        optimizer$step()
      }
      para <- res  
    } else {
      if (iteration <= free_iter) {
        component$model <- Spherical(input_dim, output_dim)
      } 
      NNmodel <- component$model
      optimizer <- optim_lbfgs(NNmodel$parameters, lr = LR, max_iter = max_iter,
                               line_search_fn = line_search_fn)
      NNmodel$train()
      
      calc_loss <- function() {
        optimizer$zero_grad()
        res <- NNmodel(X)
        loss <- pkbd_weighted_neg_log_likelihood(res$mu, res$rho, Y, W)
        loss$backward()
        loss
      }
      optimizer$step(calc_loss)
      
      NNmodel$eval()
      with_no_grad({ 
        para <- NNmodel(X)
      })
      
    }
    para$df <- d*input_dim 
    para$NNmodel <- NNmodel
    retval@defineComponent(para)
  }
  retval
}

#######################################################################################

#' Density Function for PKBD
#'
#' @description
#' Calculates the density of the PKBD for given data points.
#'
#' @param y A matrix or data frame where each row represents a data point on the unit hypersphere.
#' @param mu A vector or matrix representing the mean direction parameter(s). 
#'        If a vector, it must be normalized (unit length) and is applied to all data points.
#'        If a matrix, it must have the same number of rows as \code{y}, and each row must be normalized.
#' @param rho A scalar or a vector representing the concentration parameter. 
#'        If a vector, its length must match the number of rows in \code{y}. 
#'        Each \code{rho[i]} is used to evaluate the density for \code{y[i, ]}.
#'        Must be between 0 (inclusive) and 1 (exclusive).
#' @param log Logical; if TRUE, the log-density is returned. Default is FALSE.
#'
#' @return A vector of density values (or log-density if log = TRUE) for each row in y.
#'
#' @details
#' This function calculates the density of the PKBD for each data point in y, given the
#' parameters mu and rho. 
#' @examples
#' y <- matrix(c(1, 0, 0, 0, 0, 1), ncol = 3, byrow = TRUE)
#' mu <- c(1, 0, 0)
#' rho <- 0.5
#' dpkbd(y, mu, rho)
#' @export
dpkbd <- function(y, mu, rho, log = FALSE) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  n <- nrow(y)
  p <- ncol(y)
  
  if (is.vector(mu)) {
    if (length(mu) != p) {
      stop("mu vector length must match the number of columns in y")
    }
    # Check if mu is normalized
    if (abs(sum(mu^2) - 1) > 1e-5) {
      stop("mu must be normalized (unit length)")
    }
    mu_is_matrix <- FALSE
  } else if (is.matrix(mu)) {
    if (nrow(mu) != n) {
      stop("mu matrix must have the same number of rows as y")
    }
    if (ncol(mu) != p) {
      stop("mu matrix must have the same number of columns as y")
    }
    # Check if each row of mu is normalized
    mu_norms <- sqrt(rowSums(mu^2))
    if (any(abs(mu_norms - 1) > 1e-5)) {
      stop("All rows of mu must be normalized (unit length)")
    }
    mu_is_matrix <- TRUE
  } else {
    stop("mu must be either a vector or a matrix")
  }
  
  rho <- c(rho)
  if (!(is.numeric(rho) && (length(rho) == 1 || length(rho) == n))) {
    stop("rho must be a numeric scalar or a vector with the same length as the number of rows in y")
  }

  # Check if all rho values are between 0 (inclusive) and 1 (exclusive)
  if (any(rho < 0) || any(rho >= 1)) {
    stop("All rho values must be between 0 (inclusive) and 1 (exclusive)")
  }
  if(mu_is_matrix) rho_vector <- rep(rho, length.out = nrow(mu))
  # Check if y is normalized
  y_norms <- sqrt(rowSums(y^2))
  if (any(abs(y_norms - 1) > 1e-5)) {
    stop("All data points in y must be normalized (unit length)")
  }
  
  if(mu_is_matrix){
    density <- mapply(function(row, mu_row, r){
      logLik_PKBD(matrix(row, nrow = 1), mu_vec = mu_row, rho = r)},
      split(y, seq(n)), split(mu, seq(n)), rho_vector)
    names(density) <- NULL
  } else{
    density <- c(logLik_PKBD(y, mu_vec = mu, rho = rho))
  }
  
  if (log) {
    return(density)
  } else {
    return(exp(density))
  }
}


#' Density Function for Spherical Cauchy Distribution
#'
#' @description
#' Calculates the density of the spherical Cauchy distribution for given data points.
#'
#' @param y A matrix or data frame where each row represents a data point on the unit hypersphere.
#' @param mu A vector or matrix representing the mean direction parameter(s). 
#'        If a vector, it must be normalized (unit length) and is applied to all data points.
#'        If a matrix, it must have the same number of rows as \code{y}, and each row must be normalized.
#' @param rho A scalar or a vector representing the concentration parameter. 
#'        If a vector, its length must match the number of rows in \code{y}. 
#'        Each \code{rho[i]} is used to evaluate the density for \code{y[i, ]}.
#'        Must be between 0 (inclusive) and 1 (exclusive).
#' @param log Logical; if TRUE, the log-density is returned. Default is FALSE.
#'
#' @return A vector of density values (or log-density if log = TRUE) for each row in y.
#'
#' @details
#' This function calculates the density of the spherical Cauchy distribution for each data point in y, given the
#' parameters mu and rho. 
#' @examples
#' y <- matrix(c(1, 0, 0, 0, 0, 1), ncol = 3, byrow = TRUE)
#' mu <- c(1, 0, 0)
#' rho <- 0.5
#' dspcauchy(y, mu, rho)
#' @export
dspcauchy <- function(y, mu, rho, log = FALSE) {
  # Ensure y is a matrix
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  n <- nrow(y)
  p <- ncol(y)
  
  if (is.vector(mu)) {
    if (length(mu) != p) {
      stop("mu vector length must match the number of columns in y")
    }
    # Check if mu is normalized
    if (abs(sum(mu^2) - 1) > 1e-5) {
      stop("mu must be normalized (unit length)")
    }
    mu_is_matrix <- FALSE
  } else if (is.matrix(mu)) {
    if (nrow(mu) != n) {
      stop("mu matrix must have the same number of rows as y")
    }
    if (ncol(mu) != p) {
      stop("mu matrix must have the same number of columns as y")
    }
    # Check if each row of mu is normalized
    mu_norms <- sqrt(rowSums(mu^2))
    if (any(abs(mu_norms - 1) > 1e-5)) {
      stop("All rows of mu must be normalized (unit length)")
    }
    mu_is_matrix <- TRUE
  } else {
    stop("mu must be either a vector or a matrix")
  }
  
  rho <- c(rho)
  if (!(is.numeric(rho) && (length(rho) == 1 || length(rho) == n))) {
    stop("rho must be a numeric scalar or a vector with the same length as the number of rows in y")
  }
  
  # Check if all rho values are between 0 (inclusive) and 1 (exclusive)
  if (any(rho < 0) || any(rho >= 1)) {
    stop("All rho values must be between 0 (inclusive) and 1 (exclusive)")
  }
  if(mu_is_matrix) rho_vector <- rep(rho, length.out = nrow(mu))
  
  # Check if y is normalized
  y_norms <- sqrt(rowSums(y^2))
  if (any(abs(y_norms - 1) > 1e-5)) {
    stop("All data points in y must be normalized (unit length)")
  }
  
  if(mu_is_matrix){
    density <- mapply(function(row, mu_row, r){
      logLik_sCauchy(matrix(row, nrow = 1), mu_vec = mu_row, rho = r)},
      split(y, seq(n)), split(mu, seq(n)), rho_vector)
    names(density) <- NULL
  } else{
    density <- c(logLik_sCauchy(y, mu_vec = mu, rho = rho))
  }
  
  if (log) {
    return(density)
  } else {
    return(exp(density))
  }
}



