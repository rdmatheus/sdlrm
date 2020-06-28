#' @name sdlrm
#'
#' @title Sdl Regression for Integer-Valued Data
#'
#' @description Fit of the Sdl regression model via maximum
#' likelihood for a new parameterization of this distribution
#' that is indexed by the mean and a dispersion parameter.
#'
#' @param formula description of the model, of type \code{y ~ x} for covariates in mean only, or \code{y ~ x | z} to enter covariates in the dispersion parameter.
#' @param data an optional data frame containing the variables in the formula. By default the variables are taken from environment(formula).
#' @param link character specification of the link function for the dispersion parameter. The links \code{"log"} (default) \code{"sqrt"} and \code{"identity"} can be used.
#' @param y a integer-valued vector of the response variable. This argument is used in the \code{fit_sdl} function to get the estimates and required quantities that will return in the object resulting from the sdlrm function.
#' @param X,Z model matrices associated with the mean and the dispersion parameter, respectively, which are used in the \code{fit_sdl} function.
#' @param disp_test logical, if TRUE, the function \code{sdlrm} returns the test for constant dispersion.
#' @param control a list of control arguments specified via \code{sdl_control}.
#' @param ... arguments passed to \code{sdl_control}.
#'
#' @return  The \code{sdlrm} function returns an object of class "sdlrm", which consists of a list with the following components.
#' \describe{
#'   \item{coefficients}{ a list containing the elements "mean" and
#'         "dispersion" that consist of the estimates of the
#'          coefficients associated with the mean and the dispersion
#'          parameter, respectively,}
#'   \item{link}{the link function used for the dispersion parameter model,}
#'   \item{vcov}{asymptotic covariance matrix of the maximum likelihood estimator of the model parameters vector,}
#'   \item{logLik}{log0likelihood of the ditted model,}
#'   \item{AIC}{model Akaike information criteria,}
#'   \item{BIC}{model bayesian information criteria,}
#'   \item{n.obs, p, k}{Sample size, number of coefficients in the mean model, and number of coefficients in the dispersion parameter model, respectively,}
#'   \item{feasible}{Logical. If \code{TRUE}, the estimates obtained belong to the parametric space,}
#'   \item{pearson.residuals}{a vector with the Pearson residuals,}
#'   \item{residuals}{a vector with the randomized quantile residuals,}
#'   \item{fitted.values}{a vector with the fitted means,}
#'   \item{phi.hat}{a vector with the fitted dispersion parameters,}
#'   \item{eta}{a vector with the fitted linear predictor for the dispersion parameters,}
#'   \item{response}{the vector of the response,}
#'   \item{Z, X}{model matrices associated with the mean and the dispersion parameter, respectively,}
#'   \item{call}{the function call,}
#'   \item{formula}{the formula used to specify the model in \code{sdlrm}.}
#'  }
#'
#'
#' @references Medeiros, R. M. R., Bourguignon, M. (2020). Mean and
#'     dispersion regression model for paired count data with
#'     application to quality of life in prisons.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples #
NULL

###################################################################
# Link function                                                   #
##################################################################
g <- function(link){

  switch(link,

         identity = {
           fun <- function(theta) theta
           inv <- function(eta) eta
           deriv. <- function(theta) rep.int(1, length(theta))
           deriv.. <- function(theta) rep.int(0, length(theta))
           valideta <- function(eta) TRUE
         },

         log = {
           fun <- function(theta) log(theta)
           inv <- function(eta) pmax(exp(eta), .Machine$double.eps)
           deriv. <- function(theta) 1 / theta
           deriv.. <- function(theta) -1 / (theta ^ 2)
           valideta <- function(eta) TRUE
         },

         sqrt = {
           fun <- function(theta) sqrt(theta)
           inv <- function(eta) eta^2
           deriv. <- function(theta) 1 / (2 * sqrt(theta))
           deriv.. <- function(theta) -1 / (4 * (theta ^ (3 / 2)))
           valideta <- function(eta) all(is.finite(eta)) && all(eta > 0)
         },

         stop(gettextf("link %s not available", sQuote(link)), domain = NA))

  environment(fun) <- environment(inv) <- environment(deriv.) <-
    environment(deriv..) <- environment(valideta) <- asNamespace("stats")

  structure(list(fun = fun, inv = inv, deriv. = deriv.,
                 deriv.. = deriv.., valideta = valideta,
                 name = link), class = "link-sdlrm")
}

###################################################################
# Log-lokelihood function                                         #
##################################################################
ll_sdl <- function(theta, y, X, Z, link = "log"){

  # Link function
  inv <- g(link)$inv

  # Design matrices and necessary quantities
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- NCOL(X); k <- NCOL(Z)

  # Relations
  mu <- X%*%theta[1 : p]
  phi <- inv(Z%*%theta[(p + 1):(p + k)])

  l1 = -log(1 + phi[y < 0]) - y[y < 0] *
        (log(phi[y < 0] - mu[y < 0]) -
        log(2 + phi[y < 0] - mu[y < 0]))

  l2 = -log(1 + phi[y >= 0]) + y[y >= 0] *
        (log(phi[y >= 0] + mu[y >= 0]) -
        log(2 + phi[y >= 0] + mu[y >= 0]))

  return(sum(c(l1, l2)))
}

###################################################################
# Score function                                                  #
##################################################################
U_sdl <- function(par, y, X, Z, link = "log"){

  # Link function
  inv <- g(link)$inv
  deriv. <- g(link)$deriv.

  # Design matrices and necessary quantities
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- NCOL(X); k <- NCOL(Z)
  delta <- (as.numeric(y >= 0))

  # Relations
  mu <- X%*%par[1 : p]
  phi <- inv(Z%*%par[(p + 1):(p + k)])

  ci <- di <- rep(0, length(y))
  ci[delta == 1] = 1 / ((phi[delta == 1] + mu[delta == 1]) *
                          (2 + phi[delta == 1] + mu[delta == 1]))
  di[delta == 0] = 1 / ((phi[delta == 0] - mu[delta == 0]) *
                          (2 + phi[delta == 0] - mu[delta == 0]))

  u1 <- 2*y*ci + 2*y*di
  u2 <- -1/(1+phi) + 2*y*ci - 2*y*di

  D <- diag(as.numeric(1 / deriv.(phi)))

  Ub <- t(X)%*%u1
  Ug = t(Z)%*%D%*%u2

  U <- c(Ub, Ug)

  return(U)
}

###################################################################
# Fisher information matrix                                       #
##################################################################
K_sdl <- function(par, X, Z, link = "log"){

  # Link function
  inv <- g(link)$inv
  deriv. <- g(link)$deriv.
  deriv.. <- g(link)$deriv..

  # Design matrices and necessary quantities
  X <- as.matrix(X); Z <- as.matrix(Z)
  p <- NCOL(X); k <- NCOL(Z)

  # Relations
  mu <- X%*%par[1 : p]
  phi <- inv(Z%*%par[(p + 1):(p + k)])

  ci <- 1 / ((phi + mu) * (2 + phi + mu))
  di <- 1 / ((phi - mu) * (2 + phi - mu))

  D <- diag(as.numeric(1 / deriv.(phi)))

  K <- matrix(rep(NA, (p + k) * (p + k)), p + k, p + k)

  W1 <- diag(as.numeric((ci * (1 + phi + mu) / (1 + phi) + di *
                           (1 + phi - mu) / (1 + phi))))

  W2 <- diag(as.numeric((ci * (1 + phi + mu) / (1 + phi) - di *
                           (1 + phi - mu) /(1+phi))))

  W3 <- diag(as.numeric((- (1 + phi)^(-2) + ci * (1 + phi + mu) / (1 + phi) +
                           di * (1 + phi - mu) / (1 + phi)) * (1 / deriv.(phi))))

  Kbb <- t(X)%*%W1%*%X
  Kbg <- t(Z)%*%W2%*%D%*%X
  Kgb <- t(Kbg)
  Kgg <- t(Z)%*%W3%*%D%*%Z

  K[1:p,1:p] <- Kbb
  K[1:p, (p + 1):(p + k)] <- Kbg
  K[(p + 1):(p + k), (1:p)] <- Kgb
  K[(p + 1):(p + k), (p + 1):(p + k)] <- Kgg

  return(K)
}

###################################################################
# Control optimization function in nloptr                         #
##################################################################
sdl_control <- function(start = NULL,
                        start2 = NULL,
                        constant = 1e-8,
                        error = 1e-8,
                        optimizer = "nloptr",
                        algorithm = "NLOPT_LD_SLSQP", ...){

  rval <- list(start = start,
               start2 = start2,
               constant = constant,
               error = error,
               optimizer = optimizer,
               algorithm = algorithm)

  rval <- c(rval, list(...))
  return(rval)
}

#' @rdname sdlrm
#' @export
sdlrm <- function(formula, data, link = NULL, disp_test = FALSE, control = sdl_control(...), ...)
{
  cl <- match.call()
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  oformula <- stats::as.formula(formula)
  formula <- Formula::as.Formula(formula)

  if (length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  }else {
    if (length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }

  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- stats::terms(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)
  mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2L))
  y <- stats::model.response(mf, "numeric")
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)
  p <- NCOL(X); k <- NCOL(Z)

  if (length(y) < 1) stop("empty model")

  if (is.null(link)){
    link <- if (k == 1) "identity" else "log"
  }

  out <- fit_sdl(y, X, Z, link = link, control = control)
  out$call <- cl
  out$formula <- formula
  out$names.mean <- c("(Intercept)", colnames(X)[2:p])

  if (k > 1) {
    out$names.dispersion <- c("(Intercept)", colnames(Z)[2:k])

    if (disp_test == TRUE){
      start2 <- control$start2
      out$test <- round(disp_test(y, X, Z, cols = 2:k, link = link, start = start2), 5)
    }

  }else{
    if (out$link == "identity"){
      out$names.dispersion <- "phi"
    }else{
      out$names.dispersion <- "g(phi)"
    }
    out$test <- NULL
  }

  est <- c((out$coe)$mean, (out$coe)$disp)
  vcov <- out$vcov

  # Confidence interval (95%)
  #out$CI <- cbind(est - 1.96*sqrt(diag(vcov)),
  #                est + 1.96*sqrt(diag(vcov)))
  #colnames(out$CI) <- c("Lower","Upper")
  #rownames(out$CI) <- c(out$names.dispersion, out$names.mean)

  class(out) <- "sdlrm"
  out
}

###################################################################
# Fit function                                                    #
##################################################################
#' @rdname sdlrm
#' @export
fit_sdl <- function (y, X, Z = NULL, link = "log", control = sdl_control())
{

  n <- length(y); delta <- (as.numeric(y == 0))
  if (is.null(Z)) Z <- matrix(1, nrow = n, ncol = 1)

  # Vector lenght
  p <- NCOL(X); k <- NCOL(Z)

  # Link function
  g.inv <- g(link)$inv

  mle <- mle_sdl(y, X, Z, link = link, control = control)

  # Estimates
  est <- mle$est

  # Coefficients
  beta <- est[1:p]; gama <- est[(p + 1):(p + k)]

  # Feasible?
  feasible <- mle$feasible

  # Fit
  mu  <- X%*%beta
  phi <- g.inv(Z%*%gama)

  # Assymptotic variance and covariance matrix of the coefficients
  vcov <- solve(K_sdl(est, X, Z, link = link))

  # Loglikelihood and informations criteria
  logLik <- mle$logLik
  AIC <- - 2 * logLik + 2 * (p + k)
  BIC <- - 2 * logLik + log(n) * (p + k)

  # Pearson residuals
  rp = as.numeric((y - mu) / sqrt((phi * (2 + phi) + mu ^ 2)/2))

  # Randomized quantile residuals
  rq = rqr_sdl(y, mu, phi)

  rval <- list(coefficients = list(mean = beta, dispersion = gama),
               link = link, vcov = vcov, logLik = logLik, AIC = AIC,
               BIC = BIC, n.obs = n, p = p, k = k, feasible = feasible,
               pearson.residuals = rp, residuals = rq,
               fitted.values = structure(mu, .Names = names(y)),
               eta = g(link)$fun(phi),
               phi.hat = phi, response=y, Z=Z, X=X)
  rval
}
