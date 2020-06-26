#' @name sdlrm-methods
#' @title Methods for 'sdlrm' objects
#' @param object an object of class \code{sdlrm}.
#'
#' @return .
#'
#' @references Medeiros, R. M. R., Bourguignon, M. (2020). Mean and
#'     dispersion regression model for paired count data with
#'     application to quality of life in prisons.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
NULL

sdlrm <- function(x) UseMethod("sdlrm")

# Print
#' @export
print.sdlrm <- function(object)
{
  cat("Call:\n")
  print(object$call)
  cat("\nmean Coefficients:\n")
  print((object$coefficients)$mean)
  cat("\ndispersion Coefficients:\n")
  print((object$coefficients)$dispersion)
}

# Summary
#' @export
summary.sdlrm <- function(object)
{
  n <- object$n.obs
  p <- object$p
  k <- object$k

  # Summary for residuals
  rq <- object$residuals
  TAB.residuals <- round(cbind(mean(rq), stats::sd(rq), skewness(rq), kurtosis(rq) + 3), 6)
  colnames(TAB.residuals) <- c("Mean", "Sd", "Skewness", "Kurtosis")
  rownames(TAB.residuals) <- " "

  # Summary for mu
  est.beta <- stats::coef(object, what = "mean")
  se.beta <- sqrt(diag(stats::vcov(object))[1:p])
  tval.beta <- est.beta/se.beta
  pval.beta <- 2*stats::pnorm(abs(tval.beta), lower.tail = FALSE)

  TAB.mu <- cbind(Estimate = est.beta,
                  `Std. Error` = se.beta,
                  `t value` = tval.beta,
                  `Pr(>|t|)` = pval.beta)
  rownames(TAB.mu) <- object$names.mean

  # Summary for phi
  est.gamma <- stats::coef(object, what = "dispersion")
  se.gamma <- sqrt(diag(stats::vcov(object))[(p + 1):(p + k)])
  tval.gamma <- est.gamma/se.gamma
  pval.gamma <- 2*stats::pnorm(abs(tval.gamma), lower.tail = FALSE)

  TAB.phi <- cbind(Estimate = est.gamma,
                   `Std. Error` = se.gamma,
                   `t value` = tval.gamma,
                   `Pr(>|t|)` = pval.gamma)
  rownames(TAB.phi) <- object$names.dispersion

  Log.lik <- stats::logLik(object)#; names(Log.lik) <- " "
  AIC <- as.numeric(object$AIC); names(AIC) <- " "
  BIC <- as.numeric(object$BIC); names(BIC) <- " "

  if (!is.null(object$test)){
    p.val <- pchisq(object$test, k-1, lower.tail = FALSE)
    TAB.test <- rbind(object$test,p.val)
    colnames(TAB.test) <- c("S","W","LR","G")
    rownames(TAB.test) <- c("Value", "P value")
  }else{
    TAB.test <- NULL
  }

  out <- list(call = object$call, residuals = TAB.residuals,link = object$link,
              link.phi = object$link.phi, mean = TAB.mu, dispersion = TAB.phi, test = TAB.test,
              Log.lik = Log.lik, AIC = AIC, BIC = BIC)

  class(out) <- "summary.sdlrm"
  out
}

# Print summary
#' @export
print.summary.sdlrm <- function(object)
{
  cat("Call:\n")
  print(object$call)
  cat("\nSummary for residuals:\n")
  print(object$residuals)
  cat("\n----------------------------------------------------------------\n")
  cat("Mean:\n")
  cat("Coefficients:\n")
  stats::printCoefmat(object$mean)
  cat("\n----------------------------------------------------------------\n")
  cat("Dispersion:\n")
  cat("\nLink function:",object$link,"\n")
  cat("Coefficients:\n")
  stats::printCoefmat(object$dispersion)
  cat("\n----------------------------------------------------------------")
  if (!is.null(object$test)){
    cat("\n\nTest for constant dispersion:\n")
    print(object$test)
  }
  cat("\nIn addition, Log-lik value:",object$Log.lik,
      "\nAIC:",object$AIC,"and BIC:",object$BIC)
}


# Plot
#' @export
#' @rdname  sdlrm-methods
plot.sdlrm <- function(object)
{
  y <- object$response
  rq <- object$residuals
  mu.h <- object$fitted.values
  phi.h <- object$phi.h
  n <- object$n.obs

  graphics::par(mfrow=c(2,2))
  graphics::plot(mu.h, rq, xlab = "Fitted values", ylab = "Residuals", pch = "+")
  graphics::plot(1:n, rq, xlab = "Index", ylab = "Residuals", pch = "+")
  graphics::plot(stats::density(rq), xlab = "Residuals", ylab = "Density", main = " ", ylim = c(0, stats::dnorm(0)))
  graphics::curve(stats::dnorm(x), lty = 2, col = 2, add = T)
  stats::qqnorm(rq, xlab = "Theoretical quantile", ylab = "Residuals", pch = "+", main = " ")
  graphics::abline(0, 1, lty = 2)
  graphics::par(mfrow=c(1,1))

  ob <- sort(unique(y))
  obs <- table(y)
  esp <-expect_sdl(y, mu.h, phi.h)

  # Observed vs expected
  graphics::barplot(t(cbind(sqrt(obs), sqrt(esp))), beside = T, col = c("gray40", "lightgray"), ylab = "Frequency", xlab = "y")
  graphics::legend("topright", c("Observed", "Expected"), col = c("gray40", "lightgray"),bty = "n", pch = 15)

  # Rootogram
  lab = (- length(ob[ob < 0])):(length(ob[ob >= 0]) - 1)
  graphics::plot(lab, sqrt(esp), type = "n", ylim = c(min(sqrt(esp) - sqrt(obs)), max(sqrt(esp))),
                 xlab = "Observed", ylab = "Fit", xaxt = "n", xlim = c(min(lab) - 0.4, max(lab) + 0.4))
  graphics::axis(1, at = (- length(ob[ob < 0])):(length(ob[ob >= 0]) - 1), labels = ob)
  graphics::rect(lab - 0.47, sqrt(esp) - sqrt(obs), lab + 0.47, sqrt(esp), col = "lightgray")
  graphics::abline(h = 0)
  graphics::points(as.numeric(lab), sqrt(esp), col = "red4", pch = 16, type = 'b', lwd = 1.5)
}

# Log-likelihood
#' @export
#' @rdname sdlrm-methods
logLik.sdlrm <- function(object) {
  ll <- object$logLik
  class(ll) <- "logLik"
  return(ll)
}

# AIC
#' @export
#' @rdname sdlrm-methods
AIC.sdlrm <- function(object, k = 2) {
  AIC <- - 2 * object$logLik + k * (object$p + object$k)
  class(AIC) <- "AIC"
  return(AIC)
}

# BIC
#' @export
#' @param numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @rdname sdlrm-methods
BIC.sdlrm <- function(object) {
  n <- object$n.obs
  BIC <- - 2 * object$logLik + log(n) * (object$p + object$k)
  class(BIC) <- "BIC"
  return(BIC)
}

# Parameter estimates
#' @rdname sdlrm-methods
#' @export
#' @param what a character indicating which parameter coefficients are
#'   required, parameters for the \code{"mean"} or for the
#'   \code{"dispersion"} model. If \code{"all"} (default), a list with
#'   coefficients for the \code{mean} and for the \code{dispersion}
#'   model is returned.
coef.sdlrm <- function(object,
                       what = c("all", "mean", "dispersion")) {
  what <- match.arg(what)
  out <- switch(what,
                "all"        = list(
                  mean       = (object$coef)$mean,
                  dispersion = (object$coef)$disp),
                "mean"       = (object$coef)$mean,
                "dispersion" = (object$coef)$disp)
  return(out)
}

#  Variance-covariance matrix
#' @rdname sdlrm-methods
#' @export
vcov.sdlrm <- function(object) {
  return(object$vcov)
}

# Design matrices
#' @rdname sdlrm-methods
#' @export
#' @param matrix a character indicating which model matrix is
#'   required, the model matrix for the mean (\code{"mean"}) or for the
#'   dispersion parameter (\code{"dispersion"}). If \code{"all"} (default), a list with
#'   with both matrices are returned.
model.matrix.sdlrm <- function (object, matrix = c("all", "mean", "dispersion")) {

    what <- match.arg(what)
    out <- switch(what,
                  "all" = list(mean = object$X,
                               dispersion = object$Z),
                  "mean" = object$X,
                  "dispersion" = object$Z)
    return(out)
}
