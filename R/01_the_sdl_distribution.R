#' @name sdl
#'
#' @title The Modified Skew Discrete Laplace Distribution
#'
#' @description Probability mass function, distribution function, quantile function, and a random
#'     generation for the modified skew discrete Laplace (SDL) distribution with mean \code{mu},
#'     dispersion parameter \code{phi}, and mode \code{xi}.
#'
#' @param x,q vector of integer quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of means.
#' @param phi vector of dispersion parameters (greater than \code{abs(mu - xi)}).
#' @param xi the mode of the distribution, an integer value.
#' @param log logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \code{P(X <= x)}, otherwise, \code{P(X > x)}.
#'
#' @return \code{dsdl} returns the probability mass function, \code{psdl}
#'     gives the distribution function, \code{qsdl} gives the quantile function,
#'     and \code{rsdl} generates random observations.
#'
#' @details The SDL distribution was introduced by Kozubowski and Inusah (2006) as the discrete part
#'     of the continuous skew Laplace distribution centered at zero  (Kotz et al., 2001, Ch. 3).
#'     Although the SDL distribution has attractive properties, the discrete version of the
#'     zero-centered skew Laplace distribution induces that the mode of the resulting model
#'     is always equal to zero.
#'
#'     To overcome this limitation, Medeiros and Bourguignon (2025) proposed to obtain the discrete
#'     version of the Laplace skew distribution without setting its location parameter to zero,
#'     defining a new probability model that generalizes the SDL distribution.
#'
#'     This set of functions represents the probability mass function, the cumulative distribution
#'     function, the quantile function, and a random number generator for the modified SDL
#'     distribution parameterized in terms of \code{mu} (mean), \code{phi} (a dispersion parameter),
#'     and \code{xi} (the mode of the distribution).
#'
#'     Let \eqn{X} be a discrete random variable following a SDL distribution with mean \eqn{\mu},
#'     dispersion parameter \eqn{\phi}, and mode \eqn{\xi}. The probability mass function of X is
#'
#'     \eqn{
#'     \textrm{P}(X = x) =  \left\{\begin{array}{ll}
#'     \dfrac{1}{1 + \phi}\left(\dfrac{\phi - \mu + \xi}{2+ \phi - \mu + \xi}\right)^{-(x - \xi)}, & x \in \{\xi - 1, \xi - 2, \ldots\}, \\ \\
#'     \dfrac{1}{1 + \phi}\left(\dfrac{\phi + \mu - \xi}{2+ \phi + \mu - \xi}\right)^{x - \xi}, & x \in \{\xi, \xi + 1, \xi + 2, \ldots\}.
#'     \end{array}\right.
#'     }
#'
#'     The parametric space of this parameterization satisfies the constraint \eqn{\mu \in \mathbb{R}},
#'     \eqn{\phi > |\mu - \xi| }, and \eqn{\xi \in \mathbb{Z}}. Additionally, the expected value and
#'     the variance of \eqn{X} are given, respectively, by
#'
#'     \eqn{
#'     \textrm{E}(Y) = \mu \quad \mbox{ and } \quad \textrm{Var}(Y) = \dfrac{\phi(\phi + 2) + (\mu - \xi)^2}{2}.
#'     }
#'
#' @references Kotz, S., Kozubowski, T. J., and Podgórski, K. (2001). \emph{The Laplace Distribution
#'     and Generalizations: A Revisit with Applications to Communications, Economics, Engineering,
#'     and Finance}. Birkhauser, Boston.
#'
#'     Kozubowski, T. J., and Inusah, S. (2006). A skew Laplace distribution on integers.
#'     \emph{Annals of the Institute of Statistical Mathematics}, \bold{58}, 555---571.
#'
#'     Medeiros, R. M. R., and Bourguigon, M. (2025). Modified skew discrete Laplace regression
#'         models for integer valued data with applications to paired samples
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @examples
#' ### Probability function ###
#'
#' # Parameters
#' mu <- c(-4, 2, 4)
#' phi <- 6.5
#' xi <- 2
#'
#' xvals <- -30:30
#'
#' # Skewed-left distribution (mu < xi)
#' plot(xvals, dsdl(xvals, mu[1], phi, xi),
#'      type = "h", xlab = "x", ylab = "Pmf")
#'
#' # Symmetric distribution (mu = xi)
#' plot(xvals, dsdl(xvals, mu[2], phi, xi),
#'      type = "h", xlab = "x", ylab = "Pmf")
#'
#' # Skewed-right distribution (mu > 0)
#' plot(xvals, dsdl(xvals, mu[3], phi, xi),
#'      type = "h", xlab = "x", ylab = "Pmf")
#'
#' ### Difference between paired samples of non-negative observations  ###
#'
#' # Parameters
#' mu <- 3
#' phi <- 4
#' xi <- 0
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
#'                ylim = c(0, mean(y == 0) + 0.01))
#'
#' # Sdl model for the differences
#' points(obj, dsdl(sort(unique(y)), mu, phi, xi), col = "red", pch = 16)
NULL

#' @rdname sdl
#' @export
dsdl <- function(x, mu, phi, xi = 0, log = FALSE){

  if (any(phi <= 0))
    warning("The dispersion parameter must be positive")

  if (any(phi <= abs(mu - xi)))
    warning("Constraints are not satisfied")

  if (is.vector(x))
    x <- matrix(x, nrow = length(x))

  n <- nrow(x)
  mu <- matrix(mu, nrow = n)
  phi <- matrix(phi, nrow = n)

  lpmf <- matrix(-Inf, nrow = n)

  # NaN index
  lpmf[which(phi <= abs(mu - xi), arr.ind = TRUE)] <- NaN

  # Pmf index
  id1 <- which(x < xi & !is.nan(lpmf), arr.ind = TRUE)
  id2 <- which(x >= xi & !is.nan(lpmf), arr.ind = TRUE)

  lpmf[id1] <- (-x[id1] + xi) * (log(phi[id1] - mu[id1] + xi) - log(2 + phi[id1] - mu[id1] + xi)) -
    log(1 + phi[id1])

  lpmf[id2] <- (x[id2] - xi) * (log(phi[id2] + mu[id2] - xi) - log(2 + phi[id2] + mu[id2] - xi)) -
    log(1 + phi[id2])


  #   (1 / (1 + phi[id1])) * (((phi[id1] - mu[id1] + xi) /
  #                                        (2 + phi[id1] - mu[id1] + xi)) ^ (-x[id1] + xi))
  # pmf[id2] <- (1 / (1 + phi[id2])) * (((phi[id2] + mu[id2] - xi) /
  #                                        (2 + phi[id2] + mu[id2] + xi)) ^ (x[id12] - xi))

  if(n == 1L)
    lpmf <- as.vector(lpmf)

  if(log) lpmf else exp(lpmf)
}

#' @rdname sdl
#' @export
psdl <- function(q, mu, phi, xi = 0, lower.tail = TRUE){

  if (any(phi <= 0))
    stop("The dispersion parameter must be positive")

  if (any(phi < abs(mu - xi)))
    warning("Constraints are not satisfied")

  if (is.vector(q))
    q <- matrix(q, nrow = length(q))

  n <- nrow(q)
  mu <- matrix(mu, nrow = n)
  phi <- matrix(phi, nrow = n)

  cdf <- matrix(NaN, nrow = n, ncol = ncol(q))


  # Cdf index
  id1 <- which(q < xi & phi >= abs(mu - xi), arr.ind = TRUE)
  id2 <- which(q >= xi & phi >= abs(mu - xi), arr.ind = TRUE)

  cdf[id1] <- (2 + phi[id1] - mu[id1] + xi)^(floor(q[id1]) - xi + 1) /
    (2 * (1 + phi[id1]) * (phi[id1] - mu[id1] + xi)^(floor(q[id1] - xi)))
  cdf[id2] <- 1 - (phi[id2] + mu[id2] - xi) * ((phi[id2] + mu[id2] - xi) / (2 + phi[id2] + mu[id2] - xi))^(floor(q[id2]) - xi) /(2 * (1 + phi[id2]))

  if(nrow(cdf) == 1L)
    cdf <- as.vector(cdf)

  if(lower.tail) cdf else 1 - cdf

}

#' @rdname sdl
#' @export
qsdl <- function(p, mu, phi, xi = 0, lower.tail = TRUE){

  if ((any(p < 0)) || (any(p > 1)))
    stop("p must be in the unit interval: (0, 1)")

  if (any(phi <= 0))
    warning("The dispersion parameter must be positive")

  if (any(phi <= abs(mu - xi)))
    warning("Constraints are not satisfied")

  if(lower.tail == FALSE)
    p <- 1 - p

  if (is.vector(p))
    p <- matrix(p, nrow = length(p))

  n <- nrow(p)
  mu <- matrix(mu, nrow = n)
  phi <- matrix(phi, nrow = n)

  q <- matrix(NaN, nrow = n, ncol = ncol(p))

  p0 <- (2 + phi - mu + xi) / (2 * (1 + phi))

  id1 <- which(p < p0 & phi >= abs(mu - xi), arr.ind = TRUE)
  id2 <- which(p >= p0 & phi >= abs(mu - xi), arr.ind = TRUE)

  q[id1] <- ceiling(log(2 * (1 + phi[id1]) * p[id1] / (2 + phi[id1] - mu[id1] + xi)) /
                      log((2 + phi[id1] - mu[id1] + xi) / (phi[id1] - mu[id1] + xi)))

  q[id2] <- ceiling(round(log(2 * (1 + phi[id2]) * (1 - p[id2]) / (phi[id2] + mu[id2] - xi)) /
                      log((phi[id2] + mu[id2] - xi) / (2 + phi[id2] + mu[id2] - xi)), 8))

  q + xi
}

#' @rdname sdl
#' @export
rsdl <- function(n, mu, phi, xi = 0){

  if (any(phi <= 0))
    stop("The dispersion parameter must be positive")
  if (any(phi < abs(mu - xi)))
    warning("Constraints are not satisfied")

  stats::rgeom(n, 2 / (2 + phi + mu - xi)) -
    stats::rgeom(n, 2 / (2 + phi - mu + xi)) + xi

}
