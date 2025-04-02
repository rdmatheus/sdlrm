#' @name choose_mode
#'
#' @title Choose the Mode of the Modified Skew Discrete Laplace Regression
#'
#' @description Estimation of the mode in a modified skew discrete Laplace (SDL) regression fit
#'     via profile log-likelihood.
#'
#' @param object an object of class \code{"sdlrm"}, a result of a call to \code{\link{sdlrm}}.
#' @param grid grid of values that will be used to evaluate the profile log-likelihood function.
#' @param trace logical; if \code{TRUE}, a summary with the profile log-likelihood value, the AIC,
#'     and the BIC of the fit is displayed.
#' @param plot logical; if \code{TRUE}, a graph of the profiled log-likelihood evaluated in the
#'     considered grid of values is shown.
#' @param control a list of control arguments specified via \code{\link{sdl_control}}.
#' @param ... further arguments passed to \code{\link{sdl_control}}.
#'
#' @return An object of class \code{"choose_mode"}. More specifically, it returns a list in which
#'     each element consists of the fit of the modified SDL regression with each value of the mode
#'     specified in \code{grid}. In addition, it has the elements \code{“logLik”} with the vector
#'     of log-likelihood values for each adjustment and \code{“grid”} with the specified grid of
#'     values.
#'
#'     The \code{print} function summarizes the fits by displaying, for each value in \code{grid},
#'     the log-likelihood value and the Akaike (AIC) and Bayesian (BIC) information criteria. The
#'     \code{plot} function returns a graph of the profiled likelihood of the mode, highlighting its
#'     maximum.
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
#' # Fit with a model only for the mean with xi = 0 (default)
#' fit0 <- sdlrm(difference ~ group, data = pss)
#'
#' # Choosing the mode on the grid {-10, -9, ..., 0, ..., 9, 10}
#' fit <- choose_mode(fit0, grid = -10:10)
#'
#' # Class
#' class(fit)
#'
#' # It is possible to recovery the plot:
#' plot(fit)
#'
#' # and the trace:
#' fit
#'
#' # Fit with xi = 1
#' fit[[1]]
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @export
#'
choose_mode <- function(object, grid = -5:5, trace = TRUE, plot = TRUE,
                        control = sdl_control(...), ...) {

  n <- object$nobs
  xi_obj <- object$xi

  fit_update <- lapply(grid, function(xi) {

    if (xi != xi_obj) {
      init <- Sys.time()
      opt <- try(stats::update(object, xi = xi, control = control), silent = TRUE)
      end <- Sys.time()
    } else {
      opt <- object
    }

    if (trace) {
      cat(
        "\nmode:", xi,
        "|",
        "logLik:", if (unique(grepl("Error", opt))) NA else round(stats::logLik(opt), 3),
        "|",
        "AIC:", if (unique(grepl("Error", opt))) NA else round(stats::AIC(opt), 3),
        "|",
        "BIC:", if (unique(grepl("Error", opt))) NA else round(stats::AIC(opt, k = log(n)), 3)
      )
    }

    opt
  })

  ll <- vector()
  for (i in 1:length(grid)) {
    ll[i] <- if (unique(grepl("Error", fit_update[i]))) NA else stats::logLik(fit_update[[i]])
  }

  if (plot) {

    plot(grid, ll, type = "o", pch = 16, cex = 0.6,
         xlab = "Mode", ylab = "Profile log-likelihood")
    graphics::abline(v = grid[which.max(ll)], lty = 3, col = "grey", lwd = 2)
    graphics::points(c(grid[which.max(ll)], grid[which.max(ll)]), c(ll[which.max(ll)], ll[which.max(ll)]),
                     col = c("#56B1F7", 1), pch = c(16, 1))


  }

  fit_update <- stats::setNames(fit_update, grid)
  fit_update$logLik <- ll
  fit_update$grid <- grid
  class(fit_update) <- "choose_mode"
  fit_update

}


# Print
#' @rdname choose_mode
#' @param x an object of class \code{"choose_mode"}.
#'
#' @export
print.choose_mode <- function(x, ...) {

  cat("SDL regression fit with mode", x$xi, "\n")

  grid <- x$grid
  i <- 1
  ll <- AIC <- BIC <- vector("numeric", length(grid))
  for (xi in grid) {
    ll[i] <- if (unique(grepl("Error", x[[i]]))) NA else round(as.numeric(stats::logLik(x[[i]])), 3)
    AIC[i] <- if (unique(grepl("Error", x[[i]]))) NA else round(stats::AIC(x[[i]]), 3)
    BIC[i] <- if (unique(grepl("Error", x[[i]]))) NA else round(stats::AIC(x[[i]], k = log(x[[i]]$nobs)), 3)

    cat(
      "mode:", xi,
      "|",
      "logLik:", ll[i],
      "|",
      "AIC:", AIC[i],
      "|",
      "BIC:", BIC[i], "\n"
    )

    i <- i + 1
  }

  cat("\n\nBest value for the mode according to AIC:", grid[which.min(AIC)],
      "BIC:", grid[which.min(BIC)],
      "and logLik:", grid[which.max(ll)])

  invisible(x)
}


# Plot
#' @rdname choose_mode
#' @param x an object of class \code{"choose_mode"}.
#'
#' @export
plot.choose_mode <- function(x, ...) {
  grid <- x$grid

  ll <- vector()
  for (i in 1:length(grid)) {
    ll[i] <- if (unique(grepl("Error", x[[i]]))) NA else stats::logLik(x[[i]])
  }

  plot(grid, ll, type = "o", pch = 16, cex = 0.6,
       xlab = "Mode", ylab = "Profile log-likelihood")
  graphics::abline(v = grid[which.max(ll)], lty = 3, col = "grey", lwd = 2)
  graphics::points(c(grid[which.max(ll)], grid[which.max(ll)]), c(ll[which.max(ll)], ll[which.max(ll)]),
                   col = c("#56B1F7", 1), pch = c(16, 1))

  invisible(x)
}
