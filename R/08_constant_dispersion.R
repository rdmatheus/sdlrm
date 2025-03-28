#' Test for Constant Dispersion in the Modified SDL Regression
#'
#' Hypothesis test on constant dispersion in the modified SDL regression.
#'
#' @param object an \code{"sdlrm"} object.
#'
#' @returns the function \code{disp_test} returns the values and corresponding asymptotic
#'     \emph{p}-values of the score, Wald, likelihood ratio, and gradient test statistics
#'
#' @export
#'
#' @examples
#' ## Data set: pss (for description run ?pss)
#' barplot(table(pss$difference), xlab = "PSS index difference", ylab = "Frequency")
#' boxplot(pss$difference ~ pss$group, xlab = "Group", ylab = "PSS index difference")
#'
#' ## Fit a double model (mode = 1)
#' fit0 <- sdlrm(difference ~ group | group, data = pss, xi = 1)
#'
#' ## Constant dispersion test
#' disp_test(fit0)
#'
disp_test <- function(object){

  # Necessary quantities
  y <- object$y
  X <- object$x$mean
  Z <- object$x$dispersion

  p <- NCOL(X)
  k <- NCOL(Z);

  phi.link <- object$phi.link
  xi <- object$xi

  n <- length(y)
  q <- k - 1

  gamma0 <- as.matrix(rep(0, q))

  # Partition in Z
  Z1 <- as.matrix(Z[, 1])
  Z2 <- as.matrix(Z[, 2:k])

  # Unrestricted mle
  theta_hat <- object$optim.pars$par

  # Restricted mle
  theta_tilde <- c(sdl_mle(y, X, Z1, phi.link, xi)$par, gamma0)

  # Par_id
  id <- (p + k - q + 1):(p + k)

  # Tests statistics
  Ug2 <- U_sdl(theta_tilde, y, X, Z, phi.link, xi)[id]
  Kg2 <- chol2inv(Rfast::cholesky(K_sdl(theta_tilde, X, Z, phi.link, xi)))[id, id]

  # Score
  S <- t(Ug2)%*%Kg2%*%Ug2

  # Wald
  W <- t(theta_hat[id])%*%chol2inv(Rfast::cholesky(stats::vcov(object, "dispersion")[-1, -1]))%*%theta_hat[id]


  # Likelihood ratio
  LR <- 2 * (ll_sdl(theta_hat, y, X, Z, phi.link, xi) -
               ll_sdl(theta_tilde, y, X, Z, phi.link, xi))

  # Gradient
  G <- t(Ug2)%*%(theta_hat[id] - gamma0)

  out <- matrix(c(S, W, LR, G,
                  stats::pchisq(S, q, lower.tail = FALSE),
                  stats::pchisq(W, q, lower.tail = FALSE),
                  stats::pchisq(LR, q, lower.tail = FALSE),
                  stats::pchisq(G, q, lower.tail = FALSE)), ncol = 4, byrow = TRUE)
  colnames(out) <- c("Score","Wald","Lik. Ratio","Gradient")
  rownames(out) <- c("Statistic", "p-Value")

  round(out, 4)

}
