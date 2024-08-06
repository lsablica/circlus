library(flexmix)
library(Rcpp) 
sourceCpp("~/Documents/GitHub/PKBD---code/src/pkbd.cpp")
sourceCpp("~/Documents/GitHub/PKBD---code/src/rpkbd.cpp")
sourceCpp("~/Documents/GitHub/PKBD---code/src/sCauchy.cpp")

SCauchy_clust <- function ( formula = .~. , diagonal = TRUE ){
  retval <- new ("FLXMC", weighted = TRUE ,
                 formula = formula , dist = " Scauchy " ,
                 name = " Spherical Cauchy - based clustering ")
  retval@defineComponent <- function (para, df) {
    logLik <- function (x, y){
      #print("new iteration")
      #print(para$center)
      print(para$mu)
      print(para$rho)
      logLik_sCauchy(y , mu_vec = para$mu , rho = para$rho)
    }
    predict <- function(x){
      print(x)
    }
    new ("FLXcomponent" , parameters = list(mu = para$mu, rho = para$rho),
         df = para$df , logLik = logLik , predict = predict)
  }
  retval@fit <- function (x , y , w , ...) {
    #print(w)
    n <- nrow(y)
    d <- ncol(y)
    para <- M_step_sCauchy(y, w, n, d)
    #print(para)
    df <- (d+1)
    retval@defineComponent(c(para, df = df))
  }
  retval
}

#################################################################################################

PKBD_clust <- function ( formula = .~. , diagonal = TRUE ){
  retval <- new ("FLXMC" , weighted = TRUE ,
                 formula = formula , dist = " PKBD " ,
                 name = " PKBD - based clustering ")
  retval@defineComponent <- function (para, df) {
    logLik <- function (x , y ) {
      #print("new iteration")
      #print(para$mu)
      #print(para$rho)
      logLik_PKBD(y , mu_vec = para$mu , rho = para$rho)
    }
    predict <- function ( x ) {
      print("predict here")
      para$mu
      #print(x)
    }
    new ("FLXcomponent" , parameters = list( mu = para$mu , rho = para$rho ),
         df = para$df , logLik = logLik , predict = predict )
  }
  retval@fit <- function (x , y , w , component) {
    #print(w)
    n <- nrow(y)
    d <- ncol(y)
    if(length(component)==0){
      component <- list(mu = rep(0,d), rho = runif(1,0.7,0.95)) 
      print("sup")
    } 
    print(component)
    para <- M_step_PKBD(y, w, component$mu, component$rho, n, d)
    #print(para)
    df <- (d+1)
    retval@defineComponent(c( para , df = df))
  }
  retval
}

#################################################################################################

SCauchyNN_clust <- function ( formula = .~. , diagonal = TRUE ){
  retval <- new ("FLXMC", weighted = TRUE ,
                 formula = formula , dist = " Scauchy " ,
                 name = " Spherical Cauchy - based clustering using neural networks")
  retval@defineComponent <- function (para, df) {
    logLik <- function (x, y){
      #print("new iteration")
      #print(para$center)
      print(para$mu)
      print(para$rho)
      logLik_sCauchy(y , mu_vec = para$mu , rho = para$rho)
    }
    predict <- function(x){
      print(x)
    }
    new ("FLXcomponent" , parameters = list(mu = para$mu, rho = para$rho),
         df = para$df , logLik = logLik , predict = predict)
  }
  retval@fit <- function (x , y , w , ...) {
    #print(w)
    n <- nrow(y)
    d <- ncol(y)
    para <- M_step_sCauchy(y, w, n, d)
    #print(para)
    df <- (d+1)
    retval@defineComponent(c(para, df = df))
  }
  retval
}


#################################################################################################
library(torch)

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

pkbd_log_likelihood <- function(mu, rho, Y){
  d = Y$shape[2]
  term1 = (1-rho^2)$log()
  term2 = 1 + rho^2 - 2* rho*((mu$unsqueeze(2)$matmul(Y$unsqueeze(3)))$squeeze(3))
  neg_log_likelihood = term1 - (d/2)*term2$log()
  return(as_array(neg_log_likelihood))
}

pkbd_weighted_neg_log_likelihood <- function(mu, rho, Y, W){
  d = Y$shape[2]
  term1 = (1-rho^2)$log()
  term2 = 1 + rho^2 - 2* rho*((mu$unsqueeze(2)$matmul(Y$unsqueeze(3)))$squeeze(3))
  neg_log_likelihood = (d/2)*term2$log() - term1
  result = (neg_log_likelihood * W)$sum()
  return(result)
}

PKBDNN_clust <- function ( formula = .~. , diagonal = TRUE ){
  retval <- new ("FLXMC" , weighted = TRUE ,
                 formula = formula , dist = " PKBD " ,
                 name = " PKBD - based clustering using neural networks")
  retval@defineComponent <- function (para, df) {
    logLik <- function (x , y ) {
      #print(para$mu)
      #print(para$rho)
      Y = torch_tensor(y)
      pkbd_log_likelihood(mu = para$mu , rho = para$rho, Y)
    }
    predict <- function ( x ) {
      print("predict here")
      para$mu
      #print(x)
    }
    new ("FLXcomponent" , parameters = list( mu = para$mu , rho = para$rho, NNmodel = para$NNmodel),
         df = para$df , logLik = logLik , predict = predict )
  }
  retval@fit <- function (x , y , w , ...) {
    #print(w)
    n <- nrow(y)
    d <- ncol(y)
    input_dim = ncol(x)
    output_dim = d
    EPOCHS = 1
    LR = 0.5
    Y = torch_tensor(y)
    X = torch_tensor(x)
    W = torch_tensor(matrix(w/sum(w), ncol = 1))
    print(w)
    print(W)
    #print(Y)
    #print(X)
    
    NNmodel = Spherical(input_dim, output_dim)
    optimizer = optim_lbfgs(NNmodel$parameters, lr = LR, max_iter = 200, line_search_fn = "strong_wolfe")
    
    calc_loss <- function() {
      optimizer$zero_grad()
      res = NNmodel(X)
      #print(res$rho[1]$item())
      loss = pkbd_weighted_neg_log_likelihood(res$mu, res$rho, Y, W)
      #print(loss$item())
      loss$backward()
      loss
    }
    for(epoch in seq_len(EPOCHS)){
      optimizer$step(calc_loss)
    }
    
    para <- NNmodel(X)
    print(para)
    #print(para)
    df <- (d+1)
    retval@defineComponent(c(para , df = df, NNmodel = NNmodel))
  }
  retval
}



mix <- rbind(rPKBD_ACG(30, 0.95, c(1,0,0)), rPKBD_ACG(30, 0.9, c(-1,0,0)))
m1 <- flexmix(mix ~ 1, k = 2, model = PKBDNN_clust(), control = list(iter.max = 5))



p = PKBDNN_clust()



input_dim = 2
output_dim = 3
EPOCHS = 2
LR = 0.5
batch_size = 32

Y = rbind(rPKBD_ACG(1000, 0.95, c(1,0,0)), rPKBD_ACG(2000, 0.7, c(0,0,1)))
Y = torch_tensor(Y, dtype = torch_float32())
X = torch_tensor(cbind(rep(1,3000), c(rep(1,1000),rep(0,2000))),  dtype = torch_float32()) 
W = torch_ones(Y$shape[1], 1, dtype = torch_float32())/3000




model = PKBD(input_dim, output_dim)
optimizer = optim_lbfgs(model$parameters, lr = LR, max_iter = 200, line_search_fn = "strong_wolfe")


calc_loss <- function() {
  optimizer$zero_grad()
  res = model(X)
  print(res$rho[1]$item())
  loss = pkbd_log_likelihood(res$mu, res$rho, Y, W)
  print(loss$item())
  loss$backward()
  loss
}
for(epoch in seq_len(1)){
  cat("\nIteration: ", epoch, "\n")
  optimizer$step(calc_loss)
}


mix <- rbind(rPKBD_ACG(10, 0.95, c(1,0,0)), rPKBD_ACG(10, 0.9, c(-1,0,0)))
m1 <- flexmix(mix ~ 1, k = 2, model = PKBD_clust())
