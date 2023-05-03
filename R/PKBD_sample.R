library(Tinflex)

Null <- function(M) {
    tmp <- qr(M)
    set <- if(tmp$rank == 0) 1 : ncol(M) else - (1 : tmp$rank)
    qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
  }

PKBD_generator <- function(lambda = 0.95, d = 5){
  lpdf <- function(x) (-d*log(1-lambda*x)+(d-3)*log(1-x^2))/2
  dlpdf <- function(x) (((d-6)*lambda*x^2-2*(d-3)*x+d*lambda)/((1-lambda*x)*(1-x^2)))/2
  d2lpdf <- function(x)   ((d*lambda*lambda*(1-x^2)^2-2*(d-3)*(1+x^2)*(1-lambda*x)^2)/((1-lambda*x)^2*(1-x^2)^2))/2
  theta <- 1-10*.Machine$double.neg.eps
  r <- Re(sort(polyroot(c(d*lambda*lambda-2*d+6,4*(d-3)*lambda,(6-4*d)*lambda*lambda-2*d+6, 4*(d-3)*lambda, (6-d)*lambda*lambda))))
  x <- mean(r[-1<=r & r<=1])
  v <- crossprod(c(d*lambda*lambda-2*d+6,4*(d-3)*lambda,(6-4*d)*lambda*lambda-2*d+6, 4*(d-3)*lambda, (6-d)*lambda*lambda),c(1,x,x^2,x^3,x^4))
  ib <- if(d <=3 || v < 0) c(-theta, theta) else c(-theta, x, theta)
  Tinflex::Tinflex.setup(lpdf, dlpdf, d2lpdf, ib)
}

rPKBD <-function(n, lambda, mu){
    d <- length(mu)
    mu <- mu/sqrt(sum(mu^2))

    if(lambda == 0) {
      y <- matrix(rnorm(n * d), n, d)
      y <- y / sqrt(rowSums(y ^ 2))
    }
    else {
      gen <- PKBD_generator(lambda, d)
      w <- Tinflex.sample(gen, n)
      v <- matrix(rnorm(n * (d - 1)), n, d - 1)
      v <- v / sqrt(rowSums(v ^ 2))
      y <- tcrossprod(cbind(sqrt(1 - w ^ 2) * v, w),
                      cbind(Null(mu), mu))
    }
    
    y
}

PKBD_generatorC <- function(lambda = 0.95, d = 5){
  lpdf <- function(x) (-d*log(1-lambda*x)+(d-3)*log(1-x^2))/2
  dlpdf <- function(x) (((d-6)*lambda*x^2-2*(d-3)*x+d*lambda)/((1-lambda*x)*(1-x^2)))/2
  d2lpdf <- function(x)   ((d*lambda*lambda*(1-x^2)^2-2*(d-3)*(1+x^2)*(1-lambda*x)^2)/((1-lambda*x)^2*(1-x^2)^2))/2
  theta <- 1-.Machine$double.neg.eps
  r <- Re(sort(polyroot(c(d*lambda*lambda-2*d+6,4*(d-3)*lambda,(6-4*d)*lambda*lambda-2*d+6, 4*(d-3)*lambda, (6-d)*lambda*lambda))))
  x <- mean(r[-1<=r & r<=1])
  v <- crossprod(c(d*lambda*lambda-2*d+6,4*(d-3)*lambda,(6-4*d)*lambda*lambda-2*d+6, 4*(d-3)*lambda, (6-d)*lambda*lambda),c(1,x,x^2,x^3,x^4))
  ib <- if(d <=3 || v < 0) c(-theta, theta) else c(-theta, x, theta)
  Tinflex::Tinflex.setup.C(lpdf, dlpdf, d2lpdf, ib)
}

rPKBD_newC <-function(n, lambda, mu){
  d <- length(mu)
  mu <- mu/sqrt(sum(mu^2))
  if(lambda == 0) {
    y <- matrix(rnorm(n * d), n, d)
    y <- y / sqrt(rowSums(y ^ 2))
  }
  else{
    gen <- PKBD_generatorC(lambda, d)
    w <- Tinflex::Tinflex.sample.C(gen, n)
    Z <- matrix(rnorm(n*d), n, d)
    Y <- matrix(1, nrow=n, ncol=d)
    v <- Z - tcrossprod(Z%*%mu, mu)
    v <- v / sqrt(rowSums(v ^ 2))
    y <- tcrossprod(w,mu) + sqrt(1 - w^2) * v
  }
  y
}


V <- rPKBD(1000, 0.99, c(0,0,1))

library(Rcpp)
library(RcppArmadillo)
sourceCpp("~/Dropbox/Watson/rpkbd.cpp")

V <- rpkbd(1000, 0.99, c(0,0,1))

V <- rpkbd3(1000, 0.99, c(0,0,1))





library(rgl)
rgl.open()
view3d(zoom =0.5)
rgl.bg(color = "white")
rgl.points(V[,1],V[,2],V[,3], ylim=c(-1,1),col = "red" ,xlim=c(-1,1), zlim = c(-1,1), xlab = "x", ylab = "y", zlab = "z", alpha = 0.8)
spheres3d(x = 0, y = 0, z = 0, radius = 0.98, col = "green", alpha = 0.6, back = "lines")
rgl.lines(c(-1.5,1.5), c(0, 0), c(0, 0), color = "black")
rgl.lines(c(0, 0), c(-1.5,1.5), c(0, 0), color = "blue")
rgl.lines(c(0, 0), c(0, 0), c(-1.5,1.5), color = "red")


library(microbenchmark)
for(n in c(10, 100, 1000, 10000)){
  rhos <- c(0.01, 0.1, 0.25, 0.4, 0.6, 0.75, 0.9, 0.99)
  ds <- c( 3, 5, 10, 20, 50, 100, 200, 1000)
  
  set.seed(1)
  i=1
  results <- as.data.frame(matrix(numeric(length(rhos)*length(ds)*4), ncol = 4))
  colnames(results) <- c("d", "rho", "Tinflex", "ACG")
  for(d in ds){
    for(rho in rhos){
      lambda <- 2*rho/(1+rho^2)
      mu = rep(1, times = d)
      mu = mu/sqrt(sum(mu^2))
      m = microbenchmark(rPKBD_newC(n, lambda, mu), rpkbd3(n, lambda, mu) , times = 300L, unit = "ms")
      mm=summary(m)
      or = order(mm[,"expr"])
      results[i,] <- c(d, rho, mm[,"mean"][or])
      print(c(d, rho))
      i = i+1
    }
  }
  save(results, file = paste0("timesPKBD",n,".Rda"))
}

A <- matrix(results[,3], length(ds), length(rhos), byrow = TRUE)
colnames(A) <- rhos
rownames(A) <- ds

B <- matrix(results[,4], length(ds), length(rhos), byrow = TRUE)
colnames(B) <- rhos
rownames(B) <- ds

G = (B-A)/A
rownames(G) <- ds

A <- matrix(results[,3], length(ds), length(rhos), byrow = TRUE)
colnames(A) <- rhos
rownames(A) <- ds

B <- matrix(results[,4], length(ds), length(rhos), byrow = TRUE)
colnames(B) <- rhos
rownames(B) <- ds

G2 = (B-A)/A
rownames(G2) <- ds

par(mfrow=c(1,2))
corrplot2(G,method = "color", cl.align.text="l",  cl.offset = 0.2 , title = "n = 10", addCoef.col = "black")
text(-0.7, 4.5, bquote(d == .("")), col = 2 )
text(4.5, 9.4, bquote(rho == .("")), col = 2 )
corrplot2(G2,method = "color", cl.align.text="l",  cl.offset = 0.2 , title = "n = 100", addCoef.col = "black")
text(-0.7, 4.5, bquote(d == .("")), col = 2 )
text(4.5, 9.4, bquote(rho == .("")), col = 2 )
par(mfrow=c(1,1))

corrplot2 <- function(...){
  cb <- function(corrPlot, ..., rectArgs = list() ){
    lst <- list(...)
    n <- ncol(corrPlot)
    nms <- colnames(corrPlot)
    colnames(corrPlot) <- if(is.null(nms)) 1:ncol(corrPlot) else nms
    
    xleft <- match(lst$x, colnames(corrPlot)) - 0.5
    ybottom <- n - match(lst$y, rownames(corrPlot)) + 0.5
    
    lst <- list(xleft=xleft, ybottom=ybottom, xright=xleft+1, ytop=ybottom+1)
    do.call(rect, c(lst, rectArgs))
  }
  RR <- corrplot(...)
  fun <- function(d, kappa) ifelse(d>3, kappa >= (d-3)/2 , FALSE)*1
  YY <- outer(ds,rhos, FUN = fun)
  colnames(YY) <- rhos
  rownames(YY) <- ds
  y = rep(ds, length(rhos))[c(YY)==1]
  x = rep(rhos, each = length(ds))[c(YY)==1]
  cb(RR, x=x, y=y, rectArgs=list(border="black", lwd=3))
}

# R <- function(b,lambda){
#   t = 1/lambda - sqrt(1/lambda^2-1/b)
#   (1-b*t^2)/(1-lambda*t)
# }
# 
# RPKBD <- function(n, lambda = 0.9){
#   d = length(mu)
#   mu <- c(0,0,1)
#   
#   b = sort(Re(polyroot(c(-d^2*lambda^2, 2*d*(d-2)*lambda^2, 4*d-lambda*lambda*(d-2)^2, -4*(d-1)))))[2]
#   A = matrix(0, nrow = n, ncol = d)
#   num <- 0
#   nun <- 0
#   while(num < n){
#     m = MASS::mvrnorm(1, mu = rep(0,d), Sigma = diag(c(rep(1,d-1), 1/(1-b))) )
#     m = m/sqrt(sum(m^2))
#     #R <- (2*sqrt(1-lambda^2)/(1+sqrt(1-lambda^2)))*(1/sqrt(1-b))*((1+sqrt(1-lambda^2))/(1+sqrt(1-lambda^2/b)))^(d/2)
#     #ratio <- (1-lambda^2)*((1+sqrt(1-lambda^2))/(2))^(d/2-1)*(1/(sqrt(1-b)))*(1-lambda*crossprod(mu,m))^(-d/2)/((1-b*crossprod(mu,m)^2)^(-d/2)*R)
#     ratio <- (1-lambda*crossprod(mu,m))^(-d/2)/((1-b*crossprod(mu,m)^2)^(-d/2)*(R(b,lambda))^(d/2))
#     #print(ratio)
#     u <- runif(1)
#     if(u<ratio){
#       num = num + 1
#       A[num, ] <- m
#     }
#     nun <- nun + 1
#   }
#   print(num/nun)
#   A
# }