PKBD_generatorC <- function(lambda = 0.95, d = 5){
  lpdf <- function(x) (-d*log(1-lambda*x)+(d-3)*log(1-x^2))/2
  dlpdf <- function(x) (((d-6)*lambda*x^2-2*(d-3)*x+d*lambda)/((1-lambda*x)*(1-x^2)))/2
  d2lpdf <- function(x) ((d*lambda*lambda*(1-x^2)^2-2*(d-3)*(1+x^2)*(1-lambda*x)^2)/((1-lambda*x)^2*(1-x^2)^2))/2
  theta <- 1-.Machine$double.neg.eps
  r <- Re(sort(polyroot(c(d*lambda*lambda-2*d+6,
                          4*(d-3)*lambda,(6-4*d)*lambda*lambda-2*d+6,
                          4*(d-3)*lambda, (6-d)*lambda*lambda))))
  x <- mean(r[-1 <= r & r <= 1])
  v <- crossprod(c(d*lambda*lambda-2*d+6,
                   4*(d-3)*lambda,
                   (6-4*d)*lambda*lambda-2*d+6,
                   4*(d-3)*lambda,
                   (6-d)*lambda*lambda),
                 c(1, x, x^2, x^3, x^4))
  ib <- if (d <= 3 || v < 0) c(-theta, theta) else c(-theta, x, theta)
  Tinflex::Tinflex.setup.C(lpdf, dlpdf, d2lpdf, ib)
}

#' @title Random Sampling from PKB Distributions using Projected Saw Distribution
#' @description Generates a random sample from PKB distributions. 
#' @param n The number of random draws.
#' @param rho A numeric giving the concentration parameter.
#' @param mu A numeric vector giving the mean direction parameter.
#' @return A matrix with the generated values.
#' @rdname rPKBD_Saw
#' @keywords internal
#' @import Tinflex
#' @importFrom stats rnorm
rPKBD_Saw <-function(n, rho, mu) {
  lambda <- 2*rho/(1+rho*rho)
  d <- length(mu)
  mu <- mu/sqrt(sum(mu^2))
  if (lambda == 0) {
    y <- matrix(stats::rnorm(n * d), n, d)
    y <- y / sqrt(rowSums(y ^ 2))
  }
  else {
    gen <- PKBD_generatorC(lambda, d)
    w <- Tinflex::Tinflex.sample.C(gen, n)
    Z <- matrix(stats::rnorm(n * d), n, d)
    v <- Z - tcrossprod(Z%*%mu, mu)
    v <- v / sqrt(rowSums(v ^ 2))
    y <- tcrossprod(w, mu) + sqrt(1 - w^2) * v
  }
  y
}

#' @title Random Sampling from PKB Distributions using ACG or
#'     Projected Saw Distributions
#' @description Generates a random sample from PKB distributions.
#' @param n The number of random draws.
#' @param rho A numeric giving the concentration parameter.
#' @param mu A numeric vector giving the mean direction parameter.
#' @param method A character indicating the method to use, `"ACG"` for
#'     angular central Gaussian distribution envelopes and `"Saw"` for
#'     the use of projected Saw distributions.
#' @return A matrix with the generated values.
#' @rdname rpkb
#' @examples
#' rpkb(10, 0.95, c(1, 0, 0))
#' @export
rpkb <- function(n, rho, mu, method = "ACG") {
  if (method == "ACG") {
    return(rPKBD_ACG(n, rho, mu)) 
  } else if (method == "Saw") {
    return(rPKBD_Saw(n, rho, mu))
  } else {
    stop("No such method is available.")
  }
}
