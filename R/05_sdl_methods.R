#' @name sdlrm-methods
#' @title Methods for \code{"sdlrm"} objects.
#' @param x,object an object of class \code{"sdlrm"}.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC.
#' @param formula a model \link{formula} or \link{terms} object or an \code{"sdlrm"} object.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
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


# Residuals
#' @name residuals.sdlrm
#' @title Extract Model Residuals for a Modified SDL Regression
#'
#' @param object an \code{"sdlrm"} object.
#' @param type character; specifies which residual should be extracted.
#'     The available arguments are "quantile" (default), "pearson",
#'     and "response" (raw residuals, y - mu).
#' @param ... further arguments passed to or from other methods.
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


# Print
#' @rdname sdlrm-methods
#' @param digits a non-null value for digits specifies the minimum number of significant digits to
#'     be printed in values. The default, \code{getOption("digits")}.
#' @export
print.sdlrm <- function(x, digits = getOption("digits"), ...)
{

  p <- length(x$coefficients$mean)
  k <- length(x$coefficients$dispersion)
  n <- x$nobs

  cat("\nCall:\n")
  print(x$call)
  if(x$optim.pars$convergence != 0) {

   cat("\nmodel did not converge\n")

  }else{

    cat("\nMean Coefficients:\n")
    print(round(stats::coefficients(x, "mean"), digits))

    cat("\nDispersion Coefficients:\n")
    print(round(stats::coefficients(x, "dispersion"), digits))

    cat("\n---",
        "\nLog-lik value: ", round(stats::logLik(x), digits),
        "\nMode: ", x$xi,
        "\nAIC: ", round(stats::AIC(x), digits),
        " and BIC: ", round(stats::AIC(x, k = log(n)), digits), "\n", sep = "")

  }

  invisible(x)
}

# Summary
#' @rdname sdlrm-methods
#' @export
summary.sdlrm <- function(object, ...)
{

  ## Link function
  phi.link <- object$phi.link

  ## Summary for quantile residuals
  res <- stats::residuals(object, type = "quantile")
  skewness <- mean((res - mean(res))^3) / (stats::sd(res)^3)
  kurtosis <- mean((res - mean(res))^4) / (stats::sd(res)^4)
  TAB.residuals <- round(cbind(mean(res), stats::sd(res),
                               skewness, kurtosis), 6)
  colnames(TAB.residuals) <- c("Mean", "Std. dev.", "Skewness", "Kurtosis")
  rownames(TAB.residuals) <- " "

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
              residuals = TAB.residuals,
              phi.link = object$phi.link,
              mean = TAB.beta,
              dispersion = TAB.gamma,
              xi = object$xi,
              pR2 = pR2,
              logLik = stats::logLik(object),
              AIC = stats::AIC(object),
              BIC = stats::AIC(object, k = log(object$nobs)))



  class(out) <- "summary.sdlrm"
  out
}

# Print summary
#' @rdname sdlrm-methods
#' @export
print.summary.sdlrm <- function(x, digits = getOption("digits"), ...)
{
  cat("Call:\n")
  print(x$call)

  cat("\n\nSummary for quantile residuals:\n")
  print(round(x$residuals, digits))

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




#' Predict Method for a Modified SDL Fit
#'
#' Obtains predictions from a fitted modified SDL regression object.
#'
#' @param object an \code{"sdlrm"} object.
#' @param newdata optionally, a data frame in which to look for variables
#'     with which to predict. If omitted, the fitted linear predictors are
#'     used.
#' @param type the type of prediction required. The default is on the scale of
#'     the response variable \code{("response")}, that is, the fitted values
#'     (fitted means). The alternative \code{"dispersion"} provides the fitted
#'     dispersion, while \code{"variance"} provides the fitted variances. Finally,
#'     the option \code{"quantile"} gives the fitted quantiles in the order
#'     specified via \code{at}.
#' @param at the order of the quantile to be predicted if
#'     \code{type = "quantile"}. The default is to predict the median,
#'     that is, \code{at = 0.5}.
#' @param na.action function determining what should be done with missing
#'     values in \code{newdata}. The default is to predict \code{NA}.
#' @param ...  arguments passed to or from other methods.
#'
#' @return A vector of predictions.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
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


# Plot
#' Diagnostic Plots for the Modified SDL Regression
#'
#' Six plots (selectable by \code{which}) are currently available:
#' a plot of residuals against fitted values, a plot of residuals against
#' the observation indices, a Normal Q-Q plot, a barplot with comparisons of the
#' observed and fitted frequencies, a plot of the sample autocorelations
#' of the residuals, and a plot of the sample partial autocorelations
#' of the residuals.
#'
#' @param x an object of class \code{"sdlrm"}.
#' @param which numeric; if a subset of the plots is required, specify a subset
#'     of the numbers \code{1:6}.
#' @param type character; specifies which residual should be produced in the
#'     envelope plot. The available options are \code{"quantile"} (default),
#'     \code{"pearson"}, and \code{"response"} (raw residuals, y - mu).
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
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
plot.sdlrm <- function(x, which = 1:2,
                       type = c("quantile", "pearson", "response"),
                        ask = prod(graphics::par("mfcol")) < length(which) &&
                          grDevices::dev.interactive(),
                        ...)
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
                   xlab = "Fitted values", ylab = Type, pch = "+", ...)
    graphics::abline(h = 0, col = "royalblue", lty = 3)
  }

  ## Residuals versus index observation
  if (show[2]){
    n <- x$nobs
    graphics::plot(1:n, res, xlab = "Index", ylab = Type, pch = "+", ...)
    graphics::abline(h = 0, col= "royalblue", lty = 3)
  }

  ## Normal probability plot
  if(show[3]) {
    stats::qqnorm(res, main = "Normal Q-Q Plot",
                  xlab = "Theoretical Quantiles",
                  ylab = "Sample Quantiles",
                  plot.it = TRUE,
                  frame.plot = TRUE, pch =  "+", ...)
    graphics::abline(0, 1, col = "royalblue", lty = 2)
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
