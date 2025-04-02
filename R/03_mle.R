#' Optimization Control Parameters Passed to optim
#'
#' Optimization parameters passed to \code{\link[stats]{optim}} for the fit of an modified skew
#'     discrete Laplace (SDL) regression model via \code{\link{sdlrm}}. This function acts in the
#'     same spirit as \code{\link[betareg]{betareg.control}} from the \code{betareg} package. Its
#'     primary purpose is to gather all the optimization control arguments in a single function.
#'
#' @param method the method to be used. See "Details" in \code{\link[stats]{optim}}. The default
#'     method (\code{"BFGS"}) is a quasi-Newton method (also known as a variable metric algorithm),
#'     specifically that published simultaneously in 1970 by Broyden, Fletcher, Goldfarb and Shanno.
#' @param maxit the maximum number of iterations of the algorithm. Defaults to \code{2000}.
#' @param hessian logical. Should a numerically differentiated Hessian matrix be returned?
#' @param start an optional vector with starting values for all parameters for fitting an SDL
#'     regression model. It must be passed in the order: \code{(beta, gamma)}, where
#'     \code{beta} and \code{gamma} are regression coefficients associated with the mean and dispersion
#'     regression submodels, respectively.
#' @param reltol relative convergence tolerance. The algorithm stops if it is unable to reduce the
#'     value by a factor of reltol * (abs(val) + reltol) at a step. Defaults to \code{1e-10}.
#' @param ... further arguments to be passed to \code{\link[stats]{optim}}.
#'
#' @references Cribari-Neto, F., and Zeileis, A. (2010). Beta regression in R.
#'     \emph{Journal of statistical software}, 34, 1-24.
#'
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @return A list with the arguments specified.
#'
#' @examples
#' # Data set: pss (for description run ?pss)
#' barplot(table(pss$difference), xlab = "PSS index difference", ylab = "Frequency")
#' boxplot(pss$difference ~ pss$group, xlab = "Group", ylab = "PSS index difference")
#'
#' ## Fit of the model using the Fisher information matrix to obtain the covariance
#' ## matrix of the coefficients
#' fit1 <- sdlrm(difference ~ group, data = pss, xi = 1)
#'
#' ## Fit of the model using the numerical Hessian matrix provided by optim
#' fit2 <- sdlrm(difference ~ group, data = pss, xi = 1, hessian = TRUE)
#'
#' ## Compare the reported standard errors
#' summary(fit1)
#' summary(fit2)
#' @export
sdl_control <- function(method = "BFGS", maxit = 8000, hessian = FALSE,
                        start = NULL, reltol = 1e-10, ...) {

  rval <- list(method = method, maxit = maxit, hessian = hessian,
               start = start, reltol = reltol)

  rval <- c(rval, list(...))

  if (!is.null(rval$fnscale))
    warning("fnscale must not be modified")

  rval$fnscale <- -1

  rval
}


sdl_mle <- function(y, X = NULL, Z = NULL, link = "log", xi, control = sdl_control(...), ...){

  ### Control fit list
  method <- control$method
  maxit <- control$maxit
  hessian <- control$hessian
  start  <- control$start

  control$start <- control$hessian <- NULL

  ### Data setting
  n <- length(y)

  if (is.null(X)) X <- matrix(1, nrow = n)
  if (is.null(Z)) Z <- matrix(1, nrow = n)

  p <- ncol(X)
  k <- ncol(Z)

  ### Initial values
  if(is.null(start)){
    b <- solve(t(X)%*%X)%*%t(X)%*%y
    m. <- X%*%b
    # p. <- sqrt(1 - (m. - xi)^2 + 2 * stats::var(y)) - 1 + 2 * abs(m. - xi)
    # a <- solve(t(Z)%*%Z)%*%t(Z)%*%g(link)$fun(p.)
    #p. <- stats::var(y) + abs(m. - xi)
    a <- c(max(abs(m. - xi) + 1), rep(0L, k - 1))
    start <- c(b, a)
  }

  # Log-likelihood
  ll <- function(par) ll_sdl(par, y, X, Z, link, xi)

  # Score function
  U <- function(par) U_sdl(par, y, X, Z, link, xi)

  opt <- suppressWarnings(stats::optim(par = start,
                                       fn = ll,
                                       gr = U,
                                       method = method,
                                       control = control,
                                       hessian = hessian))

  # Convergence status
  if (opt$convergence > 0)
    warning("optimization failed to converge")

  opt

}
