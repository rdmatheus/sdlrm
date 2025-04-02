# Additional methods --------------------------------------------------------------------------
#' @name sdlrm-methods
#'
#' @title Extract Information From a Modified Skew Discrete Laplace Regression Fit
#'
#' @description Additional methods for \code{"sdlrm"} objects.
#'
#' @param object an object of class \code{"sdlrm"}, a result of a call to \code{\link{sdlrm}}.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC.
#' @param formula a model \link{formula} or \link{terms} object or an \code{"sdlrm"} object.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @returns
#' \itemize{
#' \item \code{model.frame} returns a \code{data.frame} containing the variables required
#'     by \code{formula} and any additional arguments provided via \code{...}.
#' \item \code{model.matrix} returns the design matrix used in the regression structure,
#'     as specified by the \code{parm} argument.
#' \item \code{coef} returns a numeric vector of estimated regression coefficients, based
#'     on the \code{parm} argument. If \code{parm = "full"}, it returns a list with the
#'     components \code{"mean"} and \code{"dispersion"}, each containing the corresponding
#'     coefficient estimates.
#' \item \code{vcov} returns the asymptotic covariance matrix of the regression coefficients,
#'     based on the \code{parm} argument.
#' \item \code{logLik} returns the log-likelihood value of the fitted model.
#' \item \code{AIC} returns a numeric value representing the Akaike Information Criterion
#'     (AIC), Bayesian Information Criterion, or another criterion, depending on \code{k}.
#' }
#'
#' @examples
#' # Data set: pss (for description run ?pss)
#' barplot(table(pss$difference), xlab = "PSS index difference", ylab = "Frequency")
#' boxplot(pss$difference ~ pss$group, xlab = "Group", ylab = "PSS index difference")
#'
#' # Fit a double model (mode = 1)
#' fit <- sdlrm(difference ~ group | group, data = pss, xi = 1)
#'
#' # Coef
#' coef(fit)
#' coef(fit, parm = "dispersion")
#' coef(fit, parm = "full")
#'
#' # vcov
#' vcov(fit)
#' vcov(fit, parm = "dispersion")
#' vcov(fit, parm = "full")
#'
#' # Log-likelihood value
#' logLik(fit)
#'
#' # AIC and BIC
#' AIC(fit)
#' AIC(fit, k = log(fit$nobs))
#'
#' # Model matrices
#' model.matrix(fit)
#' model.matrix(fit, "dispersion")
NULL

## Model frame
#' @export
#' @rdname sdlrm-methods
model.frame.sdlrm <- function(formula, ...) {
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
}

## Model matrix
#' @export
#' @rdname sdlrm-methods
model.matrix.sdlrm <- function(object,
                               parm = c("mean", "dispersion"), ...) {
  parm <- match.arg(parm)
  rval <- if(!is.null(object$x[[parm]])) object$x[[parm]]
  else stats::model.matrix(object$terms[[parm]], stats::model.frame(object))
  rval
}

# Parameter estimates
#' @rdname sdlrm-methods
#' @export
#' @param parm a character indicating which group of parameters is to be considered in the function.
#'     The options are \code{"mean"} (default) to obtain the output in relation to the coefficients
#'     associated with the mean, \code{"dispersion"} to obtain the output in relation to
#'     the coefficients associated with the dispersion, or \code{"full"} (when applicable)
#'     to consider all the parameters,
coef.sdlrm <- function(object,
                       parm = c("mean", "dispersion", "full"), ...) {

  parm <- match.arg(parm)
  switch(parm,
         "full"        = list(
           mean       = (object$coef)$mean,
           dispersion = (object$coef)$disp),
         "mean"       = (object$coef)$mean,
         "dispersion" = (object$coef)$disp)
}


#  Variance-covariance matrix
#' @rdname sdlrm-methods
#' @export
vcov.sdlrm <- function(object,
                       parm = c("mean", "dispersion", "full"), ...) {

  parm <- match.arg(parm)
  covm <- object$vcov

  X <- object$x$mean
  Z <- object$x$dispersion
  p <- ncol(X)
  k <- ncol(Z)


  switch(parm,
         "mean" = {
           covm <- covm[seq.int(length.out = p),
                        seq.int(length.out = p), drop = FALSE]
           covm
         },
         "dispersion" = {
           covm <- covm[seq.int(length.out = k) + p,
                        seq.int(length.out = k) + p, drop = FALSE]
           covm
         },
         "full" = {
           covm
         })

}


# Log-likelihood
#' @rdname sdlrm-methods
#' @export
logLik.sdlrm <- function(object, ...) {

  p <- length(stats::coefficients(object, "mean"))
  k <- length(stats::coefficients(object, "dispersion"))

  structure(object$logLik,
            df = p + k + 1,
            class = "logLik")
}


# AIC
#' @export
#' @rdname sdlrm-methods
AIC.sdlrm <- function(object, ..., k = 2) {

  p <- length(stats::coefficients(object, "mean"))
  k2 <- length(stats::coefficients(object, "dispersion"))

  AIC <- - 2 * object$logLik + k * (p + k2 + 1)

  class(AIC) <- "AIC"
  return(AIC)
}


# Residuals -----------------------------------------------------------------------------------
#' @name residuals.sdlrm
#' @title Extract Model Residuals for a Modified Skew Discrete Laplace Regression Fit
#'
#' @description Residuals resulting from fitting a modified Laplace discrete skew regression.
#'
#' @param object an object of class \code{"sdlrm"}, a result of a call to \code{\link{sdlrm}}.
#' @param type character; specifies which residual should be extracted.
#'     The available arguments are \code{"quantile"} (randomized quantile
#'     residuals; default), \code{"pearson"} (Pearson residuals, i.e., (y - mean) / sd),
#'     and \code{"response"} (raw residuals, i.e., y - mean).
#' @param ... further arguments passed to or from other methods.
#'
#' @returns A vector with the required residuals.
#'
#' @export
#'
#' @examples
#' ## Data set: pss (for description run ?pss)
#' barplot(table(pss$difference), xlab = "PSS index difference", ylab = "Frequency")
#' boxplot(pss$difference ~ pss$group, xlab = "Group", ylab = "PSS index difference")
#'
#' ## Fit with a model only for the mean (mode = 1)
#' fit <- sdlrm(difference ~ group, data = pss, xi = 1)
#'
#' ## Randomized quantile residuals
#' rq <- residuals(fit)
#'
#' ## Pearson residuals
#' rp <- residuals(fit, type = "pearson")
#'
#' ## Raw response residuals
#' rr <- residuals(fit, type = "response")
#'
#' cbind(quantile = rq, pearson = rp, raw = rr)
residuals.sdlrm <- function(object,
                            type = c("quantile", "pearson", "response"), ...)
{

  ## extract fit information

  y <- object$y
  mu <- stats::fitted.values(object)
  phi <- object$phi
  xi <- object$xi

  ## raw response residuals and desired type
  res <- y - mu

  type <- match.arg(type, c("quantile", "pearson", "response"))
  if(type == "response") return(res)

  rqr_sdl <- function(y){

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

  res <- switch(type,
                "pearson" = {
                  as.numeric(res / sqrt(0.5 * (phi^2 + 2 * phi + (mu - xi)^2)))
                },


                "quantile" = {
                  rqr_sdl(y)
                })

  res
}


# Summary method ------------------------------------------------------------------------------
#'
#' @name summary.sdlrm
#'
#' @title Summarizing a Modified Skew Discrete Laplace Regression Fit
#'
#' @description \code{summary} method for class \code{"sdlrm"}.
#'
#' @param object an object of class \code{"sdlrm"}, a result of a call to \code{\link{sdlrm}}.
#' @param x an object of class \code{"summary.sdlrm"}, a result of a call to \code{summary.sdlrm}.
#' @param digits a non-null value for digits specifies the minimum number of significant digits to
#'     be printed in values.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The function \code{summary.sdlrm} returns an object of class \code{"summary.sdlrm"},
#'     which consists of a list with the following components:
#'  \describe{
#'     \item{call}{the original function call, given in \code{object}.}
#'     \item{mean}{summary statistics for the mean regression structure.}
#'     \item{dispersion}{summary statistics for the dispersion regression structure.}
#'     \item{xi}{the specified mode for the model.}
#'     \item{phi.link}{the link function used for the dispersion parameter model.}
#'     \item{residuals}{the randomized quantile residuals.}
#'     \item{pR2}{the pseudo-R2 for integer-valued regression models, as introduced by
#'         Medeiros and Bourguignon (2025).}
#'     \item{logLik}{log-likelihood value of the fitted model.}
#'     \item{AIC, BIC}{Akaike and Bayesian information criteria.}
#'  }
#'
#' @export
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @references Medeiros, R. M. R., and Bourguignon, M. (2025). Modified skew discrete Laplace
#'     regression models for integer valued data with applications to paired samples.
#'     \emph{Manuscript submitted for publication.}
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
summary.sdlrm <- function(object, ...)
{

  ## Link function
  phi.link <- object$phi.link

  ## Summary for quantile residuals
  res <- stats::residuals(object, type = "quantile")

  # Summary for the mean coefficients
  est.beta <- stats::coef(object, "mean")
  se.beta <- sqrt(diag(stats::vcov(object, "mean")))
  zval.beta <- est.beta / se.beta
  pval.beta <- 2 * stats::pnorm(abs(zval.beta), lower.tail = FALSE)

  TAB.beta <- cbind(Estimate = est.beta,
                    `Std. Error` = se.beta,
                    `z value` = zval.beta,
                    `Pr(>|z|)` = pval.beta)

  # Summary for sigma
  est.gamma <- stats::coef(object, "dispersion")
  se.gamma <- sqrt(diag(stats::vcov(object, "dispersion")))
  zval.gamma <- est.gamma / se.gamma
  pval.gamma <- 2 * stats::pnorm(abs(zval.gamma), lower.tail = FALSE)

  TAB.gamma <- cbind(Estimate = est.gamma,
                     `Std. Error` = se.gamma,
                     `z value` = zval.gamma,
                     `Pr(>|z|)` = pval.gamma)

  y <- object$y
  mu <- object$fitted.values
  phi <- object$phi

  obs <- as.numeric(table(y))
  exp <- expect_sdl(y, mu, phi, object$xi)

  pR2 <- stats::cor(obs, exp)^2

  out <- list(call = object$call,
              mean = TAB.beta,
              dispersion = TAB.gamma,
              xi = object$xi,
              phi.link = object$phi.link,
              residuals = res,
              pR2 = pR2,
              logLik = stats::logLik(object),
              AIC = stats::AIC(object),
              BIC = stats::AIC(object, k = log(object$nobs)))



  class(out) <- "summary.sdlrm"
  out
}

# Print summary
#' @rdname summary.sdlrm
#' @export
print.summary.sdlrm <- function(x, digits = getOption("digits"), ...)
{
  cat("Call:\n")
  print(x$call)

  cat("\n\nSummary for randomized quantile residuals:\n")
  res <- x$residuals
  skewness <- mean((res - mean(res))^3) / (stats::sd(res)^3)
  kurtosis <- mean((res - mean(res))^4) / (stats::sd(res)^4)
  TAB.residuals <- round(cbind(mean(res), stats::sd(res),
                               skewness, kurtosis), 6)
  colnames(TAB.residuals) <- c("Mean", "Std. dev.", "Skewness", "Kurtosis")
  rownames(TAB.residuals) <- " "
  print(round(TAB.residuals, digits))

    cat("\n\nMean coefficients:\n")
  stats::printCoefmat(round(x$mean, digits))

  cat("\n\nDispersion coefficients with", x$phi.link, "link:\n")
  stats::printCoefmat(round(x$dispersion, digits))

  cat("\n---",
      "\nLog-lik value: ", round(x$logLik, digits),
      "\nMode: ", x$xi,
      "\nPseudo-R2:", round(x$pR2, digits),
      "\nAIC: ", round(x$AIC, digits),
      " and BIC: ", round(x$BIC, digits), "\n", sep = "")

  invisible(x)
}




# Predict -------------------------------------------------------------------------------------
#' Predict Method for a Modified Skew Discrete Laplace Regression Fit
#'
#' Obtains predictions from a fitted modified skew discrete Laplace regression object.
#'
#' @param object an object of class \code{"sdlrm"}, a result of a call to \code{\link{sdlrm}}.
#' @param newdata optionally, a data frame in which to look for variables
#'     with which to predict. If omitted, the fitted linear predictors are
#'     used.
#' @param type the type of prediction required. The default is on the scale of
#'     the response variable \code{("response")}, that is, the fitted values
#'     (fitted means). The alternative \code{"dispersion"} provides the fitted
#'     dispersion parameters, while \code{"variance"} provides the fitted variances. Finally,
#'     the option \code{"quantile"} gives the fitted quantiles in the order
#'     specified via \code{at}.
#' @param at the order of the quantile to be predicted if
#'     \code{type = "quantile"}. The default is to predict the median,
#'     that is, \code{at = 0.5}.
#' @param na.action function determining what should be done with missing
#'     values in \code{newdata}. The default is to predict \code{NA}.
#' @param ...  arguments passed to or from other methods.
#'
#' @return A vector with the required predictions.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @references Medeiros, R. M. R., and Bourguignon, M. (2025). Modified skew discrete Laplace
#'     regression models for integer valued data with applications to paired samples.
#'     \emph{Manuscript submitted for publication.}
#'
#' @export
#'
#' @examples
#' ## Data set: pss (for description run ?pss)
#' barplot(table(pss$difference), xlab = "PSS index difference", ylab = "Frequency")
#' boxplot(pss$difference ~ pss$group, xlab = "Group", ylab = "PSS index difference")
#'
#' ## Fit a double model (mode = 1)
#' fit <- sdlrm(difference ~ group | group, data = pss, xi = 1)
#'
#' ## Fitted values (fitted means)
#' means <- predict(fit)
#' means
#'
#' ## Fitted dispersion parameter
#' phi <- predict(fit, type = "dispersion")
#' phi
#'
#' ## Fitted variances
#' vars <- predict(fit, type = "variance")
#' vars
#'
#' ## Fitted medians
#' medians <- predict(fit, type = "quantile")
#' medians
#'
#' ## Fitted third quartiles
#' quantiles <- predict(fit, type = "quantile", at = 0.75)
#' quantiles
predict.sdlrm <- function(object, newdata = NULL,
                           type = c("response", "dispersion",
                                    "variance", "quantile"),
                           at = 0.5,
                           na.action = stats::na.pass, ...)
{

  type <- match.arg(type)
  xi <- object$xi

  qfun <- function(at, mu, phi) {

    rval <- sapply(at, function(p) qsdl(rep(p, length(mu)), mu, phi, xi))

    if(length(at) > 1L) {
      if(NCOL(rval) == 1L)
        rval <- matrix(rval, ncol = length(at),
                       dimnames = list(unique(names(rval)), NULL))

      colnames(rval) <- paste("q_", at, sep = "")
    } else {
      rval <- drop(rval)
    }
    rval
  }


  if(missing(newdata)) {

    mu <- object$fitted.values
    phi <- object$phi

    rval <- switch(type,
                   "response" = {
                     mu
                   },

                   "dispersion" = {
                     phi
                   },

                   "variance" = {
                     0.5 * (phi^2 + 2 * phi + (mu - xi)^2)
                   },
                   "quantile" = {
                     qfun(at, mu, phi)
                   })

    return(rval)

  } else {

    tnam <- switch(type,
                   "response" = "mean",
                   "dispersion" = "dispersion",
                   "variance" = "full",
                   "quantile" = "full")

    mf <- stats::model.frame(stats::delete.response(object$terms[[tnam]]),
                             newdata, na.action = na.action)
    newdata <- newdata[rownames(mf), , drop = FALSE]

    if(type %in% c("response", "variance", "quantile"))
      X <- stats::model.matrix(stats::delete.response(object$terms$mean), mf)

    if(type %in% c("dispersion", "variance", "quantile"))
      Z <- stats::model.matrix(object$terms$dispersion, mf)

    mu <- as.numeric(drop(X %*% object$coefficients$mean))
    phi <- g(object$phi.link)$inv(drop(Z %*% object$coefficients$dispersion))

    rval <- switch(type,
                   "response" = {
                     mu
                   },
                   "dispersion" = {
                     phi
                   },
                   "variance" = {
                     0.5 * (phi^2 + 2 * phi + (mu - xi)^2)

                   },
                   "quantile" = {
                     qfun(at, mu, phi)
                   }
    )

    return(rval)

  }
}



# Plot ----------------------------------------------------------------------------------------
#' Diagnostic Plots for the Modified Skew Discrete Laplace Regression
#'
#' This function provides plots for diagnostic analysis of a modified skew discrete
#'     Laplace regression fit.
#'
#' @param x an object of class \code{"sdlrm"}, a result of a call to \code{\link{sdlrm}}.
#' @param which numeric; if a subset of the plots is required, specify a subset
#'     of the numbers \code{1:6}.
#' @param type character; specifies which residual should be produced in the
#'     envelope plot. The available options are \code{"quantile"} (default),
#'     \code{"pearson"}, and \code{"response"} (raw residuals, y - mu).
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param pch,lty,... graphical parameters (see \code{\link[graphics]{par}})
#'
#' @details The \code{plot} method for \code{"\link{sdlrm}"} objects provides six types
#'     of diagnostic plots in the following order:
#'     \describe{
#'         \item{Residuals vs fitted values}{a plot of the residuals
#'             against fitted values.}
#'         \item{Residuals vs observation indices.}{an index plot of the residuals
#'             against  observation indices.}
#'         \item{Normal probability plot}{a normal probability plot of the residuals.}
#'         \item{Fitted vs observed frequencies}{a bar plot with comparisons of the observed
#'             and fitted frequencies.}
#'         \item{Sample autocorrelation plot}{sample autocorrelation function plot of the residuals.}
#'         \item{Sample partial autocorrelation plot}{sample partial autocorrelation function plot of the residuals.}
#'      }
#'
#'      The \code{which} argument can be used to select a subset of the implemented plots.
#'      Default is \code{which = 1:4}.
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @return \code{plot} method for \code{"\link{sdlrm}"} objects returns six types
#'     of diagnostic plots.
#'
#' @export
#'
#' @examples
#' ## Data set: pss (for description run ?pss)
#' barplot(table(pss$difference), xlab = "PSS index difference", ylab = "Frequency")
#' boxplot(pss$difference ~ pss$group, xlab = "Group", ylab = "PSS index difference")
#'
#' ## Fit with a model only for the mean (mode = 1)
#' fit <- sdlrm(difference ~ group, data = pss, xi = 1)
#'
#' ## Available plots (using the randomized quantile residuals):
#' # Residuals versus fitted values
#' plot(fit, which = 1)
#'
#' # Residuals versus observation indices
#' plot(fit, which = 2)
#'
#' # Normal Q-Q plot
#' plot(fit, which = 3)
#'
#' # Observed versus fitted frequencies
#' plot(fit, which = 4)
#'
#' # Sample autocorelation function of residuals
#' plot(fit, which = 5)
#'
#' # Sample partial autocorelation of residuals
#' plot(fit, which = 6)
plot.sdlrm <- function(x, which = 1:4,
                       type = c("quantile", "pearson", "response"),
                        ask = prod(graphics::par("mfcol")) < length(which) &&
                          grDevices::dev.interactive(),
                       pch = "+", lty = 2, ...)
{

  if(!is.numeric(which) || any(which < 1) || any(which > 6))
    stop("`which' must be in 1:6")

  ## Reading
  type <- match.arg(type, c("quantile", "pearson", "response"))
  res <- stats::residuals(x, type = type)

  ## Legends
  types <- c("quantile", "pearson", "response")
  Types <- c("Quantile residuals", "Pearson residuals", "Raw response residuals")
  Type <- Types[type == types]

  ## Graphical parameters setting
  if (ask) {
    op <- graphics::par(ask = TRUE)
    on.exit(graphics::par(op))
  }

  ## Plots to shown
  show <- rep(FALSE, 6)
  show[which] <- TRUE

  ## Residuals versus Fitted values
  if (show[1]){
    graphics::plot(stats::fitted(x), res,
                   xlab = "Fitted values", ylab = Type, pch = pch, ...)
    graphics::abline(h = 0, col = "royalblue", lty = lty)
  }

  ## Residuals versus index observation
  if (show[2]){
    n <- x$nobs
    graphics::plot(1:n, res, xlab = "Index", ylab = Type, pch = pch, ...)
    graphics::abline(h = 0, col= "royalblue", lty = lty)
  }

  ## Normal probability plot
  if(show[3]) {
    stats::qqnorm(res, main = "Normal Q-Q Plot",
                  xlab = "Theoretical Quantiles",
                  ylab = "Sample Quantiles",
                  plot.it = TRUE,
                  frame.plot = TRUE, pch =  pch, ...)
    graphics::abline(0, 1, col = "royalblue", lty = lty)
  }

  ## Expected frequencies
  if(show[4]) {

    y <- x$y
    mu <- x$fitted.values
    phi <- x$phi

    obs <- as.numeric(table(y))
    exp <- expect_sdl(x$y, x$fitted.values, x$phi, x$xi)

    xcoord <- graphics::barplot(table(y), xlab = "y", ylab = "Frequency",
                                ylim = c(0, max(max(obs), max(exp)) + 0.5),
                                col = "white", ...)
    graphics::points(xcoord, exp, col = "royalblue", type = "b", pch = 16)
    graphics::legend("topleft", "Fitted frequencies", col = "royalblue", lty = 1, pch = 16,
                     bty = "n")
  }

  ## ACF of residuals
  if(show[5]) {
    stats::acf(res, main = " ", xlab = "Lags",
               ylab = paste("Sample ACF of", type, "residuals"), ...)
  }

  ## PACF of residuals
  if(show[6]) {
    stats::pacf(res, main = " ", xlab = "Lags",
                ylab = paste("Sample PACF of", type, "residuals"), ...)
  }

  invisible(x)

}
