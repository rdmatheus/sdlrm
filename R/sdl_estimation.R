#' @name sdlrm_estimation
#'
#' @title Maximum likelihood estimates for the Sdl regression coefficients
#'
#' @description Maximum likelihood estimates for the coefficients of a fit
#' of the Sdl regression model.
#'
#' @param y a numeric integer-valued vector of the response variable.
#' @param X,Z model matrices associated with the mean and the dispersion parameter, respectively.
#' @param link character specification of the link function for the dispersion parameter. The links \code{"log"} (default) \code{"sqrt"} and \code{"identity"} can be used.
#' @param control a list of control arguments specified via \code{sdl_control}.
#' @param eq_constraint a function to evaluate equality constraints under the null hypothesis, in hypothesis testing contexts.
#' @param eq_constraint_jac a function to evaluate the jacobian of the equality constraints, , which were passed by the \code{eq_constraint} argument.
#' @param ... arguments passed to \code{sdl_control}.
#'
#' @return A list with the values of the maximum likelihood estimates (\code{est});
#' the value of the log likelihood at estimates (\code{logLik}); and a logical value,
#' where \code{TRUE} indicates that the point obtained is feasible (\code{feasible}).
#'
#' @references Medeiros, R. M. R., Bourguignon, M. (2020). Mean and
#'     dispersion regression model for paired count data with
#'     application to quality of life in prisons.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
mle_sdl <- function(y, X, Z, link = "log", control = sdl_control(...),
                    eq_constraint = NULL,
                    eq_constraint_jac = NULL,...){

  # Control list
  start     <- control$start
  constant  <- control$constant
  error     <- control$error
  optimizer <- control$optimizer
  algorithm <- control$algorithm

  # Link funtions

  # Dispersion link function
  g.fun <- g(link)$fun
  g.inv <- g(link)$inv
  g. <- g(link)$deriv.

  # Design matrices and necessary quantities
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  p <- NCOL(X)
  k <- NCOL(Z)
  n <- as.numeric(length(y))

  if(is.null(start)){
    b = solve(t(X)%*%X)%*%t(X)%*%y
    m. = X%*%b
    p. = -1 + sqrt(1 - m. ^ 2 + 2 * stats::var(y)) + 2 * abs(m.)
    a = solve(t(Z)%*%Z)%*%t(Z)%*%g.fun(p.)
    start = c(b, a)
  }

  # Log-likelihood
  ll <- function(par) -ll_sdl(par, y, X, Z, link)

  # Score function
  U <- function(par) -U_sdl(par, y, X, Z, link)

  # Constraints and its Jacobian
  hj <- function(par){

    # Relations
    mu  <- X%*%par[1:p]
    phi <- g.inv(Z%*%par[(p + 1):(p + k)])

    # Function h (h (theta) < 0 if theta is admissible)
    hj <- rbind(- phi - mu + constant, - phi + mu + constant)

    return(hj)
  }

  Jh <- function(par){

    # Relations
    mu  <- X%*%par[1:p]
    phi <- g.inv(Z%*%par[(p + 1):(p + k)])

    # Diagonal matrix
    D <- diag(as.numeric(1 / g.(phi)))

    # Jacobian matrix
    Jh = matrix(NA, 2 * n, p + k)
    Jh[1:n, 1:p] = - X
    Jh[1:n, (p + 1):(p + k)] = - D%*%Z
    Jh[(n + 1):(2 * n), 1:p] = X
    Jh[(n + 1):(2 * n),(p + 1):(p + k)] = - D%*%Z

    return(Jh)
  }


  if (optimizer == "nloptr"){

    if (!is.null(eq_constraint)){
      eq_constraint <- match.fun(eq_constraint)
      eq_constraint_jac <- match.fun(eq_constraint_jac)
    }

    est <- nloptr::nloptr(x0 = start,
                          eval_f = ll,
                          eval_grad_f = U,
                          eval_g_ineq = hj,
                          eval_jac_g_ineq = Jh,
                          eval_g_eq = eq_constraint,
                          eval_jac_g_eq = eq_constraint_jac,
                          opts = list("algorithm" = algorithm,
                                      "xtol_rel"= error))$solution

    logLik <- -ll(est)

    feasible <- all(hj(est) < 0)
  }

  if (optimizer == "optim"){

    est <- suppressWarnings(stats::optim(par = start,
                                         fn = ll,
                                         gr = U,
                                         method = "BFGS")$par)
    logLik <- -ll(est)

    feasible <- all(hj(est) < 0)
  }

  return(list(est = est, logLik = logLik, feasible = feasible))
}
