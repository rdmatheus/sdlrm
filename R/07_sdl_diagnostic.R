# #' @name expected_sdl
# #'
# #' @title Expected Frequencies by the SDL Regression
# #'
# #' @description Provides the expected frequencies for a fit of
# #' the Sdl regression model according to the definition given in
# #' Kleiber and Zeileis (2016).
# #'
# #' @param y vector of the observed values of the response variable.
# #' @param mu vector of the fitted values.
# #' @param phi vector of the fitted values for the dispersion parameter.
# #'
# #' @return Expected frequencies for a fit of the Sdl regression model
# #' to a dataset.
# #'
# #'
# #' @references  Kleiber, C., & Zeileis, A. (2016). Visualizing count
# #'     data regressions using rootograms. The American Statistician,
# #'     70, 296-303
# #'
# #' @references Medeiros, R. M. R., Bourguignon, M. (2020). Mean and
# #'     dispersion regression model for paired count data with
# #'     application to quality of life in prisons.
# #'
# #' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
# #'
# #' @export
# #'
# #' @examples
# #'
# #' # Sample size
# #' n <- 100
# #
# #' # Covariates
# #' X <- cbind(rep(1,n), runif(n, 0,1), runif(n,0,1))
# #' Z <- cbind(rep(1,n), runif(n, 0,1), runif(n,0,1))
# #'
# #' # Parameters and relations
# #' beta <- c(1, 1.2, 0.2)
# #' gama <- c(2, 1.5, 1.2)
# #'
# #' mu <- exp(X%*%beta)
# #' phi <- exp(Z%*%gama)
# #'
# #' y <- rsdl(n, mu, phi, xi = 0)
# #'
# #' expect_sdl(y, mu, phi, xi = 0)
# #'
expect_sdl <- function(y, mu, phi, xi){

  x <- sort(unique(y))
  n <- length(y)
  s <- rep(0, length(x))

  if (length(phi) == 1) phi = as.matrix(rep(phi, n))

  for(i in 1:n){
    s = s + dsdl(x, mu[i], phi[i], xi)
  }

  return(s)
}

#' @name envelope
#'
#' @title Envelope Plot for the Residuals of a Modified SDL Regression Fit
#'
#' @description Provides the normal probability plot with simulated
#'     envelope of Pearson residuals and randomized quantile residuals
#'     resulting from the modified SDL regression fit.
#'
#' @param object,x an object of class \code{"\link{sdlrm}"}.
#' @param type character; specifies which residual should be produced in the
#'     envelope plot. The available options are \code{"quantile"} (default) and
#'     \code{"pearson"} ((y - mean) / sd).
#' @param nsim the number of replicates. The default is \code{nsim = 99}.
#' @param level level of the sample quantile of the residual used in the
#'     construction of confidence bands.
#' @param plot logical; if \code{TRUE}, the envelope plot of the residuals is displayed.
#' @param progressBar logical; if \code{TRUE}, a progress bar is displayed giving
#'     the progress of making the graph. It can slow down the function
#'     considerably in applications with a large sample size.
#' @param ... further arguments passed to or from other methods.
#'
#' @references
#' Medeiros, R. M. R., and Bourguignon, M. (2025). Modified skew discrete Laplace regression models
#'     for integer valued data with applications to paired samples
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @return \code{envelope} returns an \code{"sdlrm_envel"} object which consists of
#'     a list with the following components:
#' \describe{
#'   \item{residuals}{a list with the quantile and pearson residuals resulting from the
#'       fit of the SDL regression model.}
#'   \item{simulation}{a list whose components are matrices containing the
#'       ordered quantile and pearson residuals of the simulation for the
#'       plot envelope.}
#'
#'  The method \code{plot} makes the envelope plot.
#'  }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Data set: pss (for description run ?pss)
#' barplot(table(pss$difference), xlab = "PSS index difference", ylab = "Frequency")
#' boxplot(pss$difference ~ pss$group, xlab = "Group", ylab = "PSS index difference")
#'
#' ## Fit with a model only for the mean (mode = 1)
#' fit <- sdlrm(difference ~ group, data = pss, xi = 1)
#'
#' ## Building the envelope plot
#' envel <- envelope(fit, plot = FALSE)
#'
#' # Class
#' class(envel)
#' envel
#'
#' # Plot for the randomized quantile residuals (default)
#' plot(envel)
#'
#' # Plot for the Pearson residuals
#' plot(envel, type = "pearson")
#' }
envelope <- function(object, nsim = 99, progressBar = TRUE, plot = TRUE, ...)
{

  # Model specifications
  y <- object$y
  X <- object$x$mean
  Z <- object$x$dispersion
  p <- NCOL(X)
  k <- NCOL(Z)
  n <- length(y)

  phi.link <- object$phi.link

  # Link functions
  g.inv <- g(phi.link)$inv
  g. <- g(phi.link)$deriv.

  # Estimates
  mu <- stats::fitted.values(object)
  phi <- object$phi
  xi <- object$xi

  # Randomized quantile residuals
  rq <- stats::residuals(object, "quantile")

  # Pearson residuals
  rp <- stats::residuals(object, "pearson")

  # Envelope plot
  rs_q <- matrix(0, n, nsim)
  rs_p <- matrix(0, n, nsim)

  rqr_sdl <- function(y, mu, phi){

    n <- length(y)
    a <- vector()
    b <- vector()
    u <- vector()
    for(i in 1:n){
      a[i] <- psdl(y[i] - 1, mu[i], phi[i], xi)
      b[i] <- psdl(y[i], mu[i], phi[i], xi)
      u[i] <- stats::runif(1, a[i], b[i])
    }

    stats::qnorm(u)
  }

  pb <- utils::txtProgressBar(1, nsim, style = 3)
  for(i in 1:nsim){

    # Simulated sample
    y.tilde <- rsdl(n, mu, phi, xi)

    # Estimates
    est.tilde <- sdl_mle(y.tilde, X, Z, phi.link, xi, start = object$optim.pars$par)$par

    mu.tilde  <-  X%*%est.tilde[1:p]
    phi.tilde <- g.inv(Z%*%est.tilde[(p + 1):(p + k)])

    # Empirical residuals
    rs_p[, i] <- sort(as.numeric((y.tilde - mu.tilde) /
                             sqrt(0.5 * (phi.tilde^2 + 2 * phi.tilde + (mu.tilde - xi)^2))))

    rs_q[, i] <- sort(rqr_sdl(y.tilde, mu.tilde, phi.tilde))

    if(progressBar){ utils::setTxtProgressBar(pb, i)}
  }

  out <- list(residuals = list(quantile = rq, pearson = rp),
              simulation = list(quantile = rs_q, pearson = rs_p))

  class(out) <- "envelope"

  if (plot)
    plot(out)

  out

}

#' @rdname envelope
#' @export
print.envelope <- function(x, ...){
  cat("\nAn 'envelope' object\n\n")
  utils::str(x)
}

#' @rdname envelope
#' @export
plot.envelope <- function(x, type = c("quantile", "pearson"), level = 0.95, ...){

  alpha <- (1 - level)/2
  type <- match.arg(type)

  if(type == "quantile"){

    res <- x$residuals$quantile
    n <- length(res)
    rs <- x$simulation$quantile

    # alpha and 1 - alpha quantiles
    lower <- apply(rs, 1, stats::quantile, probs = alpha)
    upper <- apply(rs, 1, stats::quantile, probs = 1 - alpha)

    # Median
    md <- apply(rs, 1, stats::quantile, probs = 0.5)

    # Theoretical quantiles for the normal distribution
    qq <- stats::qqnorm(1:n, plot.it = FALSE)$x

    gp <- cbind(qq, sort(res), md, lower, upper)
    graphics::plot(gp[,1], gp[,2],
                   xlab = "Normal quantiles", ylab = "Randomized quantile residuals", type = "n", ...)
    graphics::polygon (c(gp[,1], sort(gp[,1], decreasing = T)),
                       c(gp[,4], sort(gp[,5], decreasing = T)),
                       col = "lightgray", border=NA)
    graphics::lines(gp[,1], gp[,3], lty = 2)
    graphics::points(gp[,1], gp[,2], pch = "+")
    graphics::box()
  }

  if(type == "pearson"){

    res <- x$residuals$pearson
    n <- length(res)
    rs <- x$simulation$pearson

    # alpha and 1 - alpha quantiles
    lower <- apply(rs, 1, stats::quantile, probs = alpha)
    upper <- apply(rs, 1, stats::quantile, probs = 1 - alpha)

    # Median
    md <- apply(rs, 1, stats::quantile, probs = 0.5)

    # Theoretical quantiles for the normal distribution
    qq <- stats::qqnorm(1:n, plot.it = FALSE)$x

    gp <- cbind(qq, sort(res), md, lower, upper)
    graphics::plot(gp[,1], gp[,2],
                   xlab = "Normal quantiles", ylab = "Pearson residuals", type = "n", ...)
    graphics::polygon (c(gp[,1], sort(gp[,1], decreasing = T)),
                       c(gp[,4], sort(gp[,5], decreasing = T)),
                       col = "lightgray", border=NA)
    graphics::lines(gp[,1], gp[,3], lty = 2)
    graphics::points(gp[,1], gp[,2], pch = "+")
    graphics::box()
  }

  invisible(x)

}
