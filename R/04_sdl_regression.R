#' @name sdlrm
#'
#' @title Modified SDL Regression for Integer-Valued Data
#'
#' @description Fit of the modified SDL regression model via maximum likelihood for a
#'     parameterization of this distribution that is indexed by the mean, a dispersion parameter,
#'     and the mode (\code{xi}).
#'
#' @param formula a symbolic description of the model, of type \code{y ~ x} for covariates in mean
#'     only, or \code{y ~ x | z} to enter covariates in the dispersion parameter.
#' @param data an optional data frame containing the variables in the formula. By default the
#'     variables are taken from environment(formula).
#' @param subset an optional vector specifying a subset of observations to be used in the fitting
#'     process. (See additional details about how this argument interacts with data-dependent bases
#'     in the ‘Details’ section of the \link{model.frame} documentation.)
#' @param na.action a function which indicates what should happen when the data contain \code{NAs}.
#'     The default is set by the \code{na.action} setting of \link{options}, and is \link{na.fail}
#'     if that is unset. The ‘factory-fresh’ default is \link{na.omit}. Another possible value is
#'     \code{NULL}, no action. Value \link{na.exclude} can be useful.
#' @param phi.link character specification of the link function for the dispersion parameter. The links
#'     \code{"log"} (default) \code{"sqrt"} and \code{"identity"} can be used.
#' @param xi the mode of the distribution, an integer value.
#' @param control a list of control arguments specified via \code{\link{sdl_control}}.
#' @param ... arguments passed to \code{\link[stats]{optim}} via \code{\link{sdl_control}}.
#'
#' @return  The \code{sdlrm} function returns an object of class "\code{sdlrm}", which consists of a
#'     list with the following components:
#' \describe{
#'   \item{coefficients}{a list containing the elements "\code{mean}" and "\code{dispersion}" that
#'          consist of the estimates of the coefficients associated with the mean and the dispersion
#'          parameter, respectively.}
#'   \item{fitted.values}{a vector with the fitted means.}
#'   \item{phi}{a vector with the fitted dispersion parameters.}
#'   \item{phi.link}{the link function used for the dispersion parameter model.}
#'   \item{xi}{the specified mode for the model.}
#'   \item{logLik}{log-likelihood of the ditted model.}
#'   \item{vcov}{asymptotic covariance matrix of the maximum likelihood estimator of the model parameters vector.}
#'   \item{nobs}{Sample size.}
#'   \item{y}{the response vector.}
#'   \item{x}{a list with elements "\code{mean}" and "\code{dispersion}" containing the
#'       model matrices from the respective models.}
#'   \item{optim.pars}{object returned by \code{optim} function in the \code{sdlrm} function.}
#'   \item{call}{the function call.}
#'   \item{formula}{the formula used to specify the model in \code{sdlrm}.}
#'   \item{terms}{ a list with elements "mean", "dispersion" and "full"
#'          containing the terms objects for the respective models.}
#'  }
#'
#'
#' @references Medeiros, R. M. R., and Bourguignon, M. (2025). Modified skew discrete Laplace
#'     regression models for integer valued data with applications to paired samples
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @examples
#' # Data set: pss (for description run ?pss)
#' barplot(table(pss$difference), xlab = "PSS index difference", ylab = "Frequency")
#' boxplot(pss$difference ~ pss$group, xlab = "Group", ylab = "PSS index difference")
#'
#' # Fit with a model only for the mean (mode = 1)
#' fit0 <- sdlrm(difference ~ group, data = pss, xi = 1)
#' summary(fit0)
#'
#' # Fit a double model (mean and dispersion)
#' fit <- sdlrm(difference ~ group | group, data = pss, xi = 1)
#' summary(fit)
NULL

#' @rdname sdlrm
#' @export
sdlrm <- function(formula, data, subset, na.action, phi.link = "log", xi = 0,
                  control = sdl_control(...), ...)
{

  # Function call
  cl <- match.call()

  # Model frame
  if (missing(data)) data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  formula <- Formula::as.Formula(formula)

  if (length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
  }else {
    if (length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
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
  p <- NCOL(X)
  k <- NCOL(Z)

  if (length(y) < 1) stop("empty model")

  opt <- sdl_mle(y, X, Z, phi.link, xi, control)

  ### out
  beta <- opt$par[1:p]
  gamma <- opt$par[1:k + p]

  names(beta) <- colnames(X)
  names(gamma) <- colnames(Z)

  if (!is.null(opt$hessian)) {
    vcov <- solve(-opt$hessian)
  } else {
    vcov <- K_sdl(opt$par, X, Z, phi.link, xi, inverse = TRUE)
  }
  rownames(vcov) <- colnames(vcov) <- c(colnames(X), colnames(Z))

  out <- list(coefficients = list(mean = beta, dispersion = gamma),
              fitted.values = c(X%*%beta),
              phi = c(g(phi.link)$inv(Z%*%gamma)),
              phi.link = phi.link,
              xi = xi,
              logLik = opt$value,
              vcov = vcov,
              nobs = length(y),
              y = y,
              x = list(mean = X, dispersion = Z),
              optim.pars = opt,
              call = cl,
              formula = formula,
              terms = list(mean = mtX, dispersion = mtZ, full = mt))

  class(out) <- "sdlrm"
  out
}

