#' @name expected_sdl
#'
#' @title Expected Frequencies by the Sdl Regression Model
#'
#' @description Provides the expected frequencies for a fit of
#' the Sdl regression model according to the definition given in
#' Kleiber and Zeileis (2016).
#'
#' @param y vector of the observed values of the response variable.
#' @param mu vector of the fitted values.
#' @param phi vector of the fitted values for the dispersion parameter.
#'
#' @return Expected frequencies for a fit of the Sdl regression model
#' to a dataset.
#'
# #' @details For a count y with possible outcomes y = ..., -2, -1, 0, 1, 2, \ldots, the
# #' expected frequencies were defined in Kleiber and Zeileis (2016).
#'
#' @references  Kleiber, C., & Zeileis, A. (2016). Visualizing count
#'     data regressions using rootograms. The American Statistician,
#'     70, 296-303
#'
#' @references Medeiros, R. M. R., Bourguignon, M. (2020). Mean and
#'     dispersion regression model for paired count data with
#'     application to quality of life in prisons.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
#' @examples
#'
#' # Sample size
#' n <- 100
#
#' # Covariates
#' X <- cbind(rep(1,n), runif(n, 0,1), runif(n,0,1))
#' Z <- cbind(rep(1,n), runif(n, 0,1), runif(n,0,1))
#'
#' # Parameters and relations
#' beta <- c(1, 1.2, 0.2)
#' gama <- c(2, 1.5, 1.2)
#'
#' mu <- exp(X%*%beta)
#' phi <- exp(Z%*%gama)
#'
#' y <- rsdl(n, mu, phi)
#'
#' expect_sdl(y, mu, phi)
#'
expect_sdl <- function(y, mu, phi){

  x <- sort(unique(y))
  n <- length(y)
  s <- rep(0, length(x))

  if (length(phi) == 1) phi = as.matrix(rep(phi, n))

  for(i in 1:n){
    s = s + dsdl(x, mu[i], phi[i])
  }

  return(s)
}

skewness = function (x, na.rm = FALSE, type = 3)
{
  if (any(ina <- is.na(x))) {
    if (na.rm)
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3)))
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  y <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2))
  if (type == 2) {
    if (n < 3)
      stop("Need at least 3 complete observations.")
    y <- y * sqrt(n * (n - 1))/(n - 2)
  }
  else if (type == 3)
    y <- y * ((1 - 1/n))^(3/2)
  y
}

kurtosis = function (x, na.rm = FALSE, type = 3)
{
  if (any(ina <- is.na(x))) {
    if (na.rm)
      x <- x[!ina]
    else return(NA)
  }
  if (!(type %in% (1:3)))
    stop("Invalid 'type' argument.")
  n <- length(x)
  x <- x - mean(x)
  r <- n * sum(x^4)/(sum(x^2)^2)
  y <- if (type == 1)
    r - 3
  else if (type == 2) {
    if (n < 4)
      stop("Need at least 4 complete observations.")
    ((n + 1) * (r - 3) + 6) * (n - 1)/((n - 2) * (n - 3))
  }
  else r * (1 - 1/n)^2 - 3
  y
}

#' @name rqr_sdl
#'
#' @title Randomized quantile residuals in the Sdl regression model
#'
#' @description Randomized quantile residuals resulting from the fit of
#' the Sdl regression model.
#'
#' @param y vector; observed values of the response variable.
#' @param mu fitted values.
#' @param phi fitted values for the dispersion parameter.
#'
#' @return Randomized quantile residuals resulting from the fit of
#' the Sdl regression model.
#'
#' @export
#'
#' @references
#' Medeiros, R. & Bourguignon, M. (2020). \cr
#' Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics, 5, 236-244.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
rqr_sdl <- function(y, mu, phi){

  n <- length(y)
  a <- vector()
  b <- vector()
  u <- vector()

  for(i in 1:n){
    a[i] <- psdl(y[i] - 1, mu[i], phi[i])
    b[i] <- psdl(y[i], mu[i], phi[i])
    u[i] <- stats::runif(1, a[i], b[i])
  }

  return(stats::qnorm(u))
}

#' @name envel_sdl
#'
#' @title Envelope Graph of the Sdl Regression Model
#'
#' @description Provides the normal probability plot with simulated
#' envelope of Pearson residuals and randomized quantile residuals
#' resulting from the Sdl regression model fit.
#'
#' @param object object of class 'sdlrm'.
#' @param R number of replicates.
#' @param control a list of control arguments specified via \code{sdl_control}.
#'
#' @references
#' Medeiros, R. & Bourguignon, M. (2020).
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @return Returns the normal envelope simulated probability plot of
#' Pearson residuals and random quantile residuals, respectively.
#'
#' @export
envel_sdl = function(object, R = 99, control = sdl_control())
{
  # Model specifications
  y <- object$response
  X <- object$X
  Z <- object$Z
  p <- NCOL(X)
  k <- NCOL(Z)
  n <- length(y)
  delta <- as.numeric(y == 0)

  link <- object$link

  # Link functions
  g.inv <- g(link)$inv
  g. <- g(link)$deriv.

  # Estimatives
  mu <- stats::fitted.values(object)
  phi <- object$phi.hat

  ######################
  #     Residuals      #
  #####################

  # Randomized quantile residuals
  rq <- object$residuals

  # Pearson residuals
  rp <- object$pearson.residuals

  ###################
  #     Envelope    #
  ###################
  rs_1 <- matrix(0, n, R)
  rs_2 <- matrix(0, n, R)

  for(i in 1:R){

    # Simulated sample
    y.tilde <- rsdl(n, mu, phi)

    # Estimatives
    est.tilde <- mle_sdl(y.tilde, X, Z, link = link,
                         control = sdl_control(start = c(stats::coef(object, what = "mean"),
                                                         stats::coef(object, what = "dispersion"))))$est

    mu.tilde  <-  X%*%est.tilde[1:p]
    phi.tilde <- g.inv(Z%*%est.tilde[(p + 1):(p + k)])

    # Empirical residuals
    rs_1[,i] <- as.numeric((y.tilde - mu.tilde) /
                             sqrt((phi.tilde * (2 + phi.tilde) + mu.tilde ^ 2) / 2))

    rs_2[,i] <- rqr_sdl(y.tilde, mu.tilde, phi.tilde)
  }

  # Sort
  rs_1 = apply(rs_1, 2, sort)
  rs_2 = apply(rs_2, 2, sort)

  # Min and max
  mint = apply(rs_1, 1, min);  minq = apply(rs_2, 1, min)
  Maxt = apply(rs_1, 1, max);  Maxq = apply(rs_2, 1, max)

  # 0.5 and 99.5 quantiles for the envelope
  mmt = apply(rs_1, 1, stats::quantile, probs = 0.005);   mmq = apply(rs_2, 1, stats::quantile, probs = 0.005)
  MMt = apply(rs_1, 1, stats::quantile, probs = 0.995); MMq = apply(rs_2, 1, stats::quantile, probs = 0.995)

  # 2.5, 5 e 97.5 quantiles for the envelope
  mt = apply(rs_1, 1, stats::quantile, probs = 0.025);   mq = apply(rs_2, 1, stats::quantile, probs = 0.025)
  Mt = apply(rs_1, 1, stats::quantile, probs = 0.975);  Mq = apply(rs_2, 1, stats::quantile, probs = 0.975)

  # Median
  at = apply(rs_1, 1, stats::quantile, probs = 0.5);aq = apply(rs_2, 1, stats::quantile, probs = 0.5)

  # Theoretical quantiles for the normal distribution
  qq = stats::qqnorm(1:n, axes = FALSE, xlab = " ", ylab = " ",
                     type = "l", lty = 1, plot.it = FALSE)$x

  ######################################
  # Envelope                           #
  ######################################

  # Pearson residual
  stats::qqnorm(rp, main = " ", xlab = "Theoretical quantile",
                ylab = "Pearson residual", type = "n", ylim = c(-5, 5))
  graphics::polygon (c(qq, sort(qq, decreasing = T)),
                     c(mint, sort(Maxt, decreasing = T)), col = "lightgray", border=NA)
  graphics::polygon(c(qq, sort(qq, decreasing = T)),
                    c(mmt, sort(MMt, decreasing = T)), col = "gray", border = NA)
  graphics::polygon(c(qq, sort(qq, decreasing = T)),
                    c(mt, sort(Mt, decreasing = T)), col = "darkgray", border = NA)

  graphics::points(qq, sort(rp), pch = "+")
  graphics::points(qq, at, type="l", lty=2)
  graphics::box()


  # Randomized quantile residual
  stats::qqnorm(sort(rq), main = " ", xlab = "Theoretical quantile",
                ylab = "Randomized quantile residuals", type = "n", ylim = c(-3.2, 3.2))
  graphics::polygon(c(qq, sort(qq, decreasing = T)),
                    c(minq, sort(Maxq, decreasing = T)), col = "lightgray", border = NA)
  graphics::polygon(c(qq, sort(qq, decreasing = T)),
                    c(mmq, sort(MMq, decreasing = T)),col = "gray", border = NA)
  graphics::polygon(c(qq, sort(qq, decreasing = T)),
                    c(mq, sort(Mq, decreasing = T)), col = "darkgray", border = NA)

  graphics::points(qq, sort(rq), pch = "+")
  graphics::points(qq, aq, type = "l", lty = 2)
  graphics::box()

}
