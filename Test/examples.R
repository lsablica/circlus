numData = 500
alpha_current_h = 0.2
numVar = 50
mu_denom_h = 2
sum_h_weightMat = 10

root_func <- function(x) x*numVar*sum_h_weightMat + (2*numData*alpha_current_h)*x/(1 - x^2) - numVar*mu_denom_h 


rho_current <- (uniroot(
  f = root_func,
  interval = c(0,1),
  f.lower = root_func(0),
  f.upper = root_func(1),
  tol = 1e-9
))$root


c1 = numVar*sum_h_weightMat
c2 = (2*numData*alpha_current_h)
c3 = numVar*mu_denom_h

library(Rcpp)
library(microbenchmark)
sourceCpp("src/pkbd.cpp")
hybridnewton(c1, c2,c3, tol = 1e-9, maxiter = 100)

microbenchmark((uniroot(
  f = root_func,
  interval = c(0,1),
  f.lower = root_func(0),
  f.upper = root_func(1),
  tol = 1e-9
))$root, hybridnewton(c1, c2,c3, tol = 1e-9, maxiter = 100), times = 1000)




dat <- matrix(rnorm(50*5), nrow = 50, ncol = 5)
dat <- dat/sqrt(rowSums(dat^2))
numClust <- 3 

numVar <- ncol(dat)
numData <- nrow(dat)

alpha_current <- rep(1/numClust, numClust)
rho_current <- rep(0.5,numClust)
mu_current <- dat[sample(1:numData, size=numClust, replace=FALSE) ,]


mu_matrix = t(mu_current)
rho_vector = rho_current

currentIter <- 1
membCurrent <- rep.int(0, times = numData)
loglikCurrent <- -Inf

# Begin EM iterations
{
  v_mat <- dat %*% t(mu_current)
  alpha_mat_current <- matrix(alpha_current,nrow = numData,ncol = numClust,byrow = TRUE)
  rho_mat_current <- matrix(rho_current,nrow = numData,ncol = numClust,byrow = TRUE)
  log_probMat_denom <- log(1 + rho_mat_current^2 - 2*rho_mat_current*v_mat)
  # eq (11) in ICMLD16 paper
  log_probMat <- log(1 - (rho_mat_current^2)) - (numVar / 2) * log_probMat_denom
  ######### E step done#############################################
  ########## M step now#############################################
  # Denominator of eq (18) in ICMLD16 paper
  probSum <- matrix(exp(log_probMat) %*% alpha_current,nrow = numData,ncol = numClust)
  # eq (18) in ICMLD16 paper
  log_normProbMat_current <- log(alpha_mat_current) + log_probMat - log(probSum)
  beta_matrix = exp(log_normProbMat_current)
  # denominator of eq (20) in ICMLD16 paper
  log_weightMat <- log_normProbMat_current - log_probMat_denom
  
  # eq (19) in ICMLD16 paper
  alpha_current <- colSums(exp(log_normProbMat_current)) / numData
  # numerator of fraction in eq (21) of ICMLD16 paper
  mu_num_sum_MAT <- t(exp(log_weightMat)) %*% dat
  norm_vec <- function(x) sqrt(sum(x^2))
  mu_denom <- apply(mu_num_sum_MAT, 1, norm_vec)
  # eq (21) of ICMLD16 paper without sign function
  mu_current <- mu_num_sum_MAT / mu_denom
  for(h in 1:numClust){
    # update rho params
    sum_h_weightMat <- sum(exp(log_weightMat[,h]))
    alpha_current_h <- alpha_current[h]
    mu_denom_h <- mu_denom[h]
    root_func <- function(x) {
      (-2*numData*x*alpha_current_h)/(1 - x^2) + numVar*mu_denom_h -numVar*x*sum_h_weightMat
    }
    rho_current[h] <- (uniroot(
      f = root_func,
      interval = c(0,1),
      f.lower = root_func(0),
      f.upper = root_func(1),
      tol = 0.001
    ))$root
  }
  # Update counter to NEXT iteration number.
  currentIter <- currentIter + 1
}

mu_current ; rho_current
M_step(dat, beta_matrix, mu_matrix, rho_vector, k =3, n = 50, d = 5, tol = 1e-6, maxiter = 100)
M_step_PKBD(dat, beta_matrix[,1], mu_matrix[,1], rho_vector[1], n = 50, d = 5, tol = 1e-6, maxiter = 100)


library(Rcpp) 
library(microbenchmark)
sourceCpp("src/sCauchy.cpp")



gsl::hyperg_2F1(2, 5, 9, 0.7 )
5*hyper2F1( 6, 9, 0.7 )-4*hyper2F1( 5, 9, 0.7 )
5*hyper2F1( 6, 9, 0.7 )+(8/0.7)*(1-hyper2F1( 4, 8, 0.7 ))


n1d(8, 0.7)
n1d2(8, 0.7)

n1d_deriv(8, 0.7)
n1d_deriv2(8, 0.7)


hybridnewton(8,0.8)
n1d(8, 0.539277)

n1d(8, hybridnewton(8,0.99))

X = movMF::rmovMF(4, 0.7*c(1,0,0))


rho = 0.5
mu = c(1,1,1)/sqrt(sum(c(1,1,1)^2))
psi = rho * mu
V = (1-rho^2)*(X + matrix(psi, byrow = TRUE, ncol = 3, nrow = 4))/matrix(1+rho^2+2*rho*X%*%mu, nrow = 4, ncol = 3) + matrix(psi, byrow = TRUE, ncol = 3, nrow = 4)

Moebius_S(X, mu , rho)

k=2
n=4
d=3

beta_matrix = matrix(abs(rnorm(n*k)), nrow = n, ncol = k)
beta_matrix = beta_matrix/matrix(rowSums(beta_matrix), nrow = n, ncol = k)

weights_matrix = beta_matrix/matrix(colSums(beta_matrix), nrow = n, ncol = k, byrow = TRUE)
weighted_means = t(X)%*%weights_matrix

for(i in 1:k){
  means = weighted_means[,i]
  weightss = weights_matrix[,i]
  mu0 = means/sqrt(sum(means^2))
  rho0 = hybridnewton(d, sqrt(sum(means^2)));
  psi = rho0 * mu0
  psiold = psi +10
  print(psi)
  
  while(norm(psi-psiold, type = "2") > 1e-6){
    psiold = psi
    trans_data_weighted = t(Moebius_S(X, -mu0 , rho0)) %*% weightss
    psi = psiold + ((d+1)*(1-rho0^2)/(2*d))*trans_data_weighted
    rho0 = sqrt(sum(psi^2))
    mu0 = psi/rho0
  }
  print(mu0)
  print(rho0)
}


k=2
n=100
d=3
vv = c(rep(c(0.99,0.99,0.98,0.98,0.97),10), rep(c(0.01,0.01,0.03,0.02,0.01),10))

X = rbind(rsCauchy(50, 0.7, c(1,0,0)), rsCauchy(50, 0.2, c(0,1,0)))
beta_matrix = matrix(c(vv, 1-vv), nrow = n, ncol = k)
beta_matrix = beta_matrix/matrix(rowSums(beta_matrix), nrow = n, ncol = k)


M_step_sCauchy2(X, beta_matrix, k =2, n = 4, d = 3, tol = 1e-6, maxiter = 100)
M_step_sCauchy(X, beta_matrix[,1], n = 4, d = 3, tol = 1e-6, maxiter = 100)


library(flexmix)
mymclust <- function ( formula = .~. , diagonal = TRUE ){
  retval <- new ("FLXMC" , weighted = TRUE ,
                 formula = formula , dist = " mvnorm " ,
                 name = " my model - based clustering ")
  retval@defineComponent <- function ( para ) {
     logLik <- function (x , y ) {
       #print("new iteration")
       #print(para$center)
       #print(para$cov)
       #print(mvtnorm::dmvnorm(y , mean = para$center ,sigma = para$cov , log = TRUE ))
       mvtnorm::dmvnorm(y , mean = para$center ,
                        sigma = para$cov , log = TRUE )
     }
     predict <- function ( x ) {
       print(x)
       matrix ( para$center , nrow = nrow ( x ) ,
               ncol = length ( para$center ) , byrow = TRUE )
     }
     new ("FLXcomponent" , parameters =
         list ( center = para$center , cov = para$cov ) ,
         df = para$df , logLik = logLik , predict = predict )
  }
  retval@fit <- function (x , y , w , ...) {
    print(w)
    para <- cov.wt (y , wt = w )[ c ("center" , "cov")]
    df <- (3 * ncol ( y ) + ncol ( y )^2)/2
    if (diagonal) {
      para$cov <- diag ( diag ( para$cov ))
      df <- 2 * ncol ( y )
    }
    retval@defineComponent ( c ( para , df = df ))
  }
  retval
}

mix <- rbind(mvtnorm::rmvnorm(10, mean = c(1,1), sigma = diag(c(1,1))), mvtnorm::rmvnorm(10, mean = c(10,10), sigma = diag(c(1,1))))

m1 <- flexmix(mix ~ 1, k = 2, model = mymclust())


mvtnorm::dmvnorm(mix , mean = c(6.667134, 7.083449),
                 sigma = matrix(c(17.2475,0,0,21.54636),2,2) , log = TRUE )


predict(m1, as.data.frame(mvtnorm::rmvnorm(10, mean = c(1,1))))


library(Rcpp) 
library(microbenchmark)
sourceCpp("src/sCauchy.cpp")

mymclust <- function ( formula = .~. , diagonal = TRUE ){
  retval <- new ("FLXMC" , weighted = TRUE ,
                 formula = formula , dist = " Scauchy " ,
                 name = " my model - based clustering ")
  retval@defineComponent <- function ( para ) {
    logLik <- function (x , y ) {
      #print("new iteration")
      #print(para$center)
      print(para$mu)
      print(para$rho)
      logLik_sCauchy(y , mu_vec = para$mu , rho = para$rho)
    }
    predict <- function ( x ) {
      print(x)
    }
    new ("FLXcomponent" , parameters =
           list ( mu = para$mu , rho = para$rho ),
         df = para$df , logLik = logLik , predict = predict )
  }
  retval@fit <- function (x , y , w , ...) {
    #print(w)
    n <- nrow(y)
    d <- ncol(y)
    para <- M_step_sCauchy(y, w, n, d)
    #print(para)
    df <- (d+1)
    retval@defineComponent ( c ( para , df = df ))
  }
  retval
}

mix <- rbind(rsCauchy(10, 0.95, c(1,0,0)), rsCauchy(10, 0.9, c(-1,0,0)))
m1 <- flexmix(mix ~ 1, k = 2, model = mymclust())


sourceCpp("src/pkbd.cpp")
sourceCpp("src/rpkbd.cpp")


mymclust2 <- function ( formula = .~. , diagonal = TRUE ){
  retval <- new ("FLXMC" , weighted = TRUE ,
                 formula = formula , dist = " PKBD " ,
                 name = " my model - based clustering ")
  retval@defineComponent <- function ( para ) {
    logLik <- function (x , y ) {
      #print("new iteration")
      #print(para$center)
      #print(para$mu)
      #print(para$rho)
      logLik_PKBD(y , mu_vec = para$mu , rho = para$rho)
    }
    predict <- function ( x ) {
      print(x)
    }
    new ("FLXcomponent" , parameters =
           list ( mu = para$mu , rho = para$rho ),
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
    retval@defineComponent ( c ( para , df = df ))
  }
  retval
}

mix <- rbind(rPKBD_ACG(10, 0.95, c(1,0,0)), rPKBD_ACG(10, 0.9, c(-1,0,0)))
m1 <- flexmix(mix ~ 1, k = 2, model = mymclust2())


library(reticulate)
np <- import("numpy")
X <- rbind(rPKBD_ACG(1000, 0.95, c(1,0,0)), rPKBD_ACG(1000, 0.9, c(-1,0,0)))
X <- np$asarray(X)
np$save("X.npy",r_to_py(X))


