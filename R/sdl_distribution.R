#' @name sdl
#'
#' @title The Skew Discrete Laplace (Sdl) Distribution
#'
#' @description Probability mass function, distribution function,
#' quantile function and random generation for the Sdl distribution
#' with mean \code{mu} and dispersion parameter \code{phi}.
#'
#' @param x vector of integer quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of means.
#' @param phi vector of dispersion indexes (greather than \code{abs(mu)}).
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P(X <= x)}, otherwise, \code{P(X > x)}.
#'
#' @return \code{dsdl} returns the probability function, \code{psdl}
#' gives the distribution function, \code{qsdl} gives the quantile function,
#' and \code{rsdl} generates random observations.
#'
#' @details The Sdl distribution was proposed by Kozubowski and Inusah (2006). This
#' probability distribution has some attractive features as the support
#' in integers \code{Z = {..., -2, -1, 0, 1, 2, ...}}; mathematical expressions in closed-form
#' for the cumulative distribution function, quantile function, characteristic function,
#' ordinal and absolute moments; assumes a skewed-left and skewed-right distribution,
#' and also symmetry; in addition, it can be represented stochastically as the
#' difference between two independent random variables, both with geometric distribution.
#'
#' This set of functions represents the probability function, the cumulative distribution
#' function, quantile function and a random number generator for the Sdl distribution
#' parameterized in terms of the mean and a dispersion parameter. This new
#' parameterization was proposed by Medeiros and Bourguignin (2020).
#'
# #' Let X be a discrete random variable following a Sdl distribution
# #' with mean \code{mu} and dispersion parameter \code{phi}. The probability
# #' mass function of X is
#'
#' @references Medeiros, R. and Bourguignon, M. (2020).
#'
#' @references Kozubowski, T. J., and Inusah, S. (2006). A skew Laplace distribution on integers. Annals of the Institute of Statistical Mathematics, 58, 555--571.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#'
#' ### Probability function ###
#'
#' # Parameters
#' mu <- c(-2, 0, 2); phi <- 3
#
#' # Skewed-left distribution (mu < 0)
#' plot(-15:15, dsdl(-15:15, mu[1], phi),
#'       type = "h", xlab = "x", ylab = "Pmf")
#'
#' # Simmetric distribution (mu = 0)
#' plot(-15:15, dsdl(-15:15, mu[2], phi),
#'       type = "h", xlab = "x", ylab = "Pmf")
#'
#' # Skewed-right distribution (mu > 0)
#' plot(-15:15, dsdl(-15:15, mu[3], phi),
#'       type = "h", xlab = "x", ylab = "Pmf")
#'
#' ### Difference between paired samples of non-negative observations  ###
#'
#' # Parameters
#' mu <- 2; phi <- 4
#'
#' # Paired samples of a pre-post treatment experimental study
#' before <- rgeom(1000, 2 / (2 + phi - mu))
#' after <- rgeom(1000, 2 / (2 + phi + mu))
#'
#' # Response variable
#' y <- after - before
#'
#' # Barplot
#' obj <- barplot(prop.table(table(y)),
#'                xlab = "Response",
#'                ylab = "Proportion",
#'                col = "white",
#'                ylim = c(0, dsdl(0, mu, phi)))
#'
#' # Sdl model for the differences
#' points(obj, dsdl(sort(unique(y)), mu, phi), col = "red", pch = 16)
#'
NULL

#' @rdname sdl
#' @export
dsdl <- function(x, mu, phi){
  if (any(phi <= 0))
    stop("The dispersion parameter must be positive")
  if (any(phi < abs(mu)))
    warning("Constraints are not satisfied")

  prob.minus = (1 / (1 + phi)) * (((phi - mu) / (2 + phi - mu)) ^ abs(x[x < 0]))
  prob.plus = (1 / (1 + phi)) * (((mu + phi) / (2 + mu + phi)) ^ x[x >= 0])

  p=c(prob.minus, prob.plus)
  index=c(which(x < 0), which(x >= 0))

  return(p[sort(index,index.return=T)$ix])
}

#' @rdname sdl
#' @export
psdl <- function(q, mu, phi, lower.tail = TRUE){
  if (any(phi <= 0))
    stop("The dispersion parameter must be positive")
  if (any(phi < abs(mu)))
    warning("Constraints are not satisfied")

  p1 <- (floor(q[q < 0]) + 1) * log(2 + phi - mu) - log(2 * (1 + phi)) -
    floor(q[q < 0]) * log(phi - mu)

  p2 <- log(mu + phi) + floor(q[q >= 0]) * log((mu + phi) / (2 + phi + mu)) - log(2 * (1 + phi))

  prob <- c(exp(p1), 1 - exp(p2))
  index <- c(which(q < 0), which(q >= 0))

  prob[sort(index,index.return=T)$ix]
  ifelse(lower.tail == TRUE, prob, 1 - prob)

  return(prob)
}

#' @rdname sdl
#' @export
qsdl <- function(p, mu, phi, lower.tail = TRUE){
  if ((any(p < 0)) || (any(p > 1)))
    stop("p must be in the unit interval: (0, 1)")
  if (any(phi <= 0))
    stop("The dispersion parameter must be positive")
  if (any(phi < abs(mu)))
    warning("Constraints are not satisfied")

  if(lower.tail == FALSE)
    p <- 1 - p

  p0 <- (2 + phi - mu) / (2 * (1 + phi))

  ind1 <- which(p < p0)
  ind2 <- which(p >= p0)

  ifelse(length(p) > 1, p.1 <- p[ind1], p.1 <- p)
  ifelse(length(mu) > 1, mu.1 <- mu[ind1], mu.1 <- mu)
  ifelse(length(phi) > 1, phi.1 <- phi[ind1], phi.1 <- phi)

  q1 <- ceiling(round(log((2 * (1 + phi.1) * p.1) / (2 + phi.1 - mu.1))/
                        log((2 + phi.1 - mu.1) / (phi.1 - mu.1)), 8))

  ifelse(length(p) > 1, p.2 <- p[ind2], p.2 <- p)
  ifelse(length(mu) > 1, mu.2 <- mu[ind2], mu.2 <- mu)
  ifelse(length(phi) > 1, phi.2 <- phi[ind2], phi.2 <- phi)

  q2 <- ceiling(round(log((2 * (1 + phi.2) * (1 - p.2)) / (mu.2 + phi.2)) /
                        log((mu.2 + phi.2) / (2 + phi.2 + mu.2)),8))

  quanti <- c(q1, q2)

  return(quanti[sort(c(ind1, ind2), index.return = TRUE)$ix])
}

#' @rdname sdl
#' @export
rsdl <- function(n, mu, phi){
  if (any(phi <= 0))
    stop("The dispersion parameter must be positive")
  if (any(phi < abs(mu)))
    warning("Constraints are not satisfied")

  u <- stats::runif(n)

  return(qsdl(u, mu, phi))
}
