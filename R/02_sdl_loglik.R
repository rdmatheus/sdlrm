# Link function
g <- function(link){

  switch(link,

         identity = {
           fun <- function(theta) theta
           inv <- function(eta) eta
           deriv. <- function(theta) rep.int(1, length(theta))
           deriv.. <- function(theta) rep.int(0, length(theta))
           valideta <- function(eta) TRUE
         },

         log = {
           fun <- function(theta) log(theta)
           inv <- function(eta) pmax(exp(eta), .Machine$double.eps)
           deriv. <- function(theta) 1 / theta
           deriv.. <- function(theta) -1 / (theta ^ 2)
           valideta <- function(eta) TRUE
         },

         sqrt = {
           fun <- function(theta) sqrt(theta)
           inv <- function(eta) eta^2
           deriv. <- function(theta) 1 / (2 * sqrt(theta))
           deriv.. <- function(theta) -1 / (4 * (theta ^ (3 / 2)))
           valideta <- function(eta) all(is.finite(eta)) && all(eta > 0)
         },

         stop(gettextf("link %s not available", sQuote(link)), domain = NA))

  environment(fun) <- environment(inv) <- environment(deriv.) <-
    environment(deriv..) <- environment(valideta) <- asNamespace("stats")

  structure(list(fun = fun, inv = inv, deriv. = deriv.,
                 deriv.. = deriv.., valideta = valideta,
                 name = link), class = "link-sdlrm")
}

# Log-likelihood function
ll_sdl <- function(theta, y, X, Z, link = "log", xi){

  # Inverse link function
  inv <- g(link)$inv

  # Design matrices
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  p <- NCOL(X)
  k <- NCOL(Z)

  # Parameters
  mu <- X%*%theta[1 : p]
  phi <- inv(Z%*%theta[1:k + p])

  # LogLik
  if (any(!is.finite(mu)) | any(!is.finite(phi)) | any(phi <= abs(mu - xi)))
    NaN
  else {
    ll <- suppressWarnings(dsdl(y, mu, phi, xi, log = TRUE))
    if (any(!is.finite(ll)))
      NaN
    else sum(ll)
  }


}

# Score function
U_sdl <- function(par, y, X, Z, link = "log", xi){

  # Link function
  inv <- g(link)$inv
  deriv. <- g(link)$deriv.

  # Design matrices
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  p <- NCOL(X)
  k <- NCOL(Z)

  # Parameters
  mu <- c(X%*%par[1 : p])
  phi <- c(inv(Z%*%par[(p + 1):(p + k)]))

  # Matrix representation
  delta <- (y >= xi)
  ci <- 1 / ((phi + mu - xi) * (2 + phi + mu - xi))
  di <- 1 / ((phi - mu + xi) * (2 + phi - mu + xi))

  u1 <- u2 <- rep(NA, length(y))

  u1[delta] <- 2 * (y[delta] - xi) * ci[delta]
  u1[!delta] <- 2 * (y[!delta] - xi) * di[!delta]

  u2[delta] <- - 1 / (1 + phi[delta]) + 2 * (y[delta] - xi) * ci[delta]
  u2[!delta] <- - 1 / (1 + phi[!delta]) - 2 * (y[!delta] - xi) * di[!delta]

  D <- diag(as.numeric(1 / deriv.(phi)))

  Ub <- t(X)%*%u1
  Ug <- t(Z)%*%D%*%u2

  c(Ub, Ug)

}

# Hessian matrix
J_sdl <- function(par, y, X, Z, link = "log", xi, inverse = FALSE) {

  # Link function
  inv <- g(link)$inv
  deriv. <- g(link)$deriv.
  deriv.. <- g(link)$deriv..

  # Design matrices
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  p <- NCOL(X)
  k <- NCOL(Z)

  # Relations
  mu <- c(X %*% par[1:p])
  phi <- c(inv(Z %*% par[(p + 1):(p + k)]))

  # Pré-computação para evitar repetições
  one_plus_phi <- 1 + phi
  deriv_phi <- deriv.(phi)
  deriv_inv <- 1 / deriv_phi
  deriv_double <- deriv..(phi)

  # Indicador delta
  delta <- as.numeric(y >= xi)

  # Cálculo de ci e di evitando valores muito pequenos
  denom1 <- (phi + mu - xi) * (2 + phi + mu - xi)
  denom2 <- (phi - mu + xi) * (2 + phi - mu + xi)

  ci <- pmax(1 / denom1, .Machine$double.eps^(1/2))
  di <- pmax(1 / denom2, .Machine$double.eps^(1/2))

  # Matriz diagonal D
  D <- diag(deriv_inv)

  # Matrizes V
  common_term1 <- -4 * delta * (y - xi) * ci^2 * (one_plus_phi + mu - xi)
  common_term2 <- 4 * (1 - delta) * (y - xi) * di^2 * (one_plus_phi - mu + xi)

  V1 <- diag(common_term1 + common_term2)
  V2 <- diag(common_term1 - common_term2)

  V3 <- diag(((1 / one_plus_phi^2) + common_term1 + common_term2) * deriv_inv +
               ((1 / one_plus_phi) - 2 * delta * (y - xi) * ci +
                  2 * (1 - delta) * (y - xi) * di) * deriv_double / deriv_phi^2)

  # Cálculo das submatrizes
  Jbb <- crossprod(X, V1 %*% X)
  Jbg <- crossprod(X, V2 %*% D %*% Z)
  Jgg <- crossprod(Z, V3 %*% D %*% Z)

  # Construção da matriz completa
  J <- matrix(0, p + k, p + k)
  J[1:p, 1:p] <- Jbb
  J[1:p, (p + 1):(p + k)] <- Jbg
  J[(p + 1):(p + k), 1:p] <- t(Jbg)  # Aproveita a transposta já calculada
  J[(p + 1):(p + k), (p + 1):(p + k)] <- Jgg

  # Se inverse = TRUE, calcular inversa de J
  if (inverse) {

    # Uso da decomposição de Cholesky para maior eficiência
    Jbb_inv <- solve(Jbb)
    Jgg_inv <- solve(Jgg)

    temp1 <- Jbb_inv %*% Jbg
    temp2 <- Jgg_inv %*% t(Jbg)

    J_inv <- matrix(0, p + k, p + k)
    J_inv[1:p, 1:p] <- solve(Jbb - Jbg %*% Jgg_inv %*% t(Jbg))
    J_inv[1:p, (p + 1):(p + k)] <- -temp1 %*% solve(Jgg - t(Jbg) %*% temp1)
    J_inv[(p + 1):(p + k), 1:p] <- -temp2 %*% solve(Jbb - Jbg %*% Jgg_inv %*% t(Jbg))
    J_inv[(p + 1):(p + k), (p + 1):(p + k)] <- solve(Jgg - t(Jbg) %*% temp1)

    return(J_inv)
  }

  return(J)
}


J_sdl2 <- function(par, y, X, Z, link = "log", xi){

  # Link function
  inv <- g(link)$inv
  deriv. <- g(link)$deriv.
  deriv.. <- g(link)$deriv..

  # Design matrices
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  p <- NCOL(X)
  k <- NCOL(Z)

  # Relations
  mu <- c(X%*%par[1 : p])
  phi <- c(inv(Z%*%par[1:k + p]))

  # Matrix representation
  delta <- (as.numeric(y >= xi))
  ci <- pmax(1 / ((phi + mu - xi) * (2 + phi + mu - xi)), .Machine$double.eps^(1/2))
  di <- pmax(1 / ((phi - mu + xi) * (2 + phi - mu + xi)), .Machine$double.eps^(1/2))

  D <- diag(as.numeric(1 / deriv.(phi)))

  J <- matrix(rep(NA, (p + k) * (p + k)), p + k, p + k)

  V1 <- diag(as.numeric(

        -4 * delta * (y - xi) * ci^2 * (1 + phi + mu - xi) +
          4 * (1 - delta) * (y - xi) * di^2 * (1 + phi - mu + xi)

    ))

  V2 <- diag(as.numeric(

        -4 * delta * (y - xi) * ci^2 * (1 + phi + mu - xi) -
          4 * (1 - delta) * (y - xi) * di^2 * (1 + phi - mu + xi)

    ))

  V3 <- diag(as.numeric(

        (1 / (1 + phi)^2 - 4 * delta * (y - xi) * ci^2 * (1 + phi + mu - xi) +
           4 * (1 - delta) * (y - xi) * di^2 * (1 + phi - mu + xi)) / deriv.(phi) +

        (1 / (1 + phi) - 2 * delta * (y - xi) * ci +
           2 * (1 - delta) * (y - xi) * di) * deriv..(phi) / deriv.(phi)^2
  ))

  Jbb <- t(X)%*%V1%*%X
  #Jbg <- t(Z)%*%V2%*%D%*%X
  Jbg <- t(X)%*%V2%*%D%*%Z
  Jgb <- t(Jbg)
  Jgg <- t(Z)%*%V3%*%D%*%Z

  J[1:p,1:p] <- Jbb
  J[1:p, (p + 1):(p + k)] <- Jbg
  J[(p + 1):(p + k), (1:p)] <- Jgb
  J[(p + 1):(p + k), (p + 1):(p + k)] <- Jgg

  J
}


# Fisher information matrix
K_sdl <- function(par, X, Z, link = "log", xi, inverse = FALSE) {

  # Link function
  inv <- g(link)$inv
  deriv. <- g(link)$deriv.

  # Design matrices
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  p <- NCOL(X)
  k <- NCOL(Z)

  # Relations
  mu <- c(X %*% par[1:p])
  phi <- c(inv(Z %*% par[(p + 1):(p + k)]))

  # Pré-computação para evitar repetições
  one_plus_phi <- 1 + phi
  deriv_phi <- deriv.(phi)
  deriv_inv <- 1 / deriv_phi

  # Cálculo de ci e di evitando valores muito pequenos
  denom1 <- (phi + mu - xi) * (2 + phi + mu - xi)
  denom2 <- (phi - mu + xi) * (2 + phi - mu + xi)

  ci <- pmax(1 / denom1, .Machine$double.eps^(1/2))
  di <- pmax(1 / denom2, .Machine$double.eps^(1/2))

  # Matriz diagonal D
  D <- diag(deriv_inv)

  # Matrizes W
  common_term1 <- ci * (one_plus_phi + mu - xi) / one_plus_phi
  common_term2 <- di * (one_plus_phi - mu + xi) / one_plus_phi

  W1 <- diag(common_term1 + common_term2)
  W2 <- diag(common_term1 - common_term2)
  W3 <- diag((-(one_plus_phi)^(-2) + common_term1 + common_term2) * deriv_inv)

  # Cálculo das submatrizes
  Kbb <- crossprod(X, W1 %*% X)
  Kbg <- crossprod(X, W2 %*% D %*% Z)
  Kgg <- crossprod(Z, W3 %*% D %*% Z)

  # Construção da matriz completa
  K <- matrix(0, p + k, p + k)
  K[1:p, 1:p] <- Kbb
  K[1:p, (p + 1):(p + k)] <- Kbg
  K[(p + 1):(p + k), 1:p] <- t(Kbg)  # Aproveita a transposta já calculada
  K[(p + 1):(p + k), (p + 1):(p + k)] <- Kgg

  # Se inverse = TRUE, calcular inversa de K
  if (inverse) {

    # Uso da decomposição de Cholesky para maior eficiência
    Kbb_inv <- solve(Kbb)
    Kgg_inv <- solve(Kgg)

    temp1 <- Kbb_inv %*% Kbg
    temp2 <- Kgg_inv %*% t(Kbg)

    K_inv <- matrix(0, p + k, p + k)
    K_inv[1:p, 1:p] <- solve(Kbb - Kbg %*% Kgg_inv %*% t(Kbg))
    K_inv[1:p, (p + 1):(p + k)] <- -temp1 %*% solve(Kgg - t(Kbg) %*% temp1)
    K_inv[(p + 1):(p + k), 1:p] <- -temp2 %*% solve(Kbb - Kbg %*% Kgg_inv %*% t(Kbg))
    K_inv[(p + 1):(p + k), (p + 1):(p + k)] <- solve(Kgg - t(Kbg) %*% temp1)

    return(K_inv)
  }

  return(K)
}

K_sdl2 <- function(par, X, Z, link = "log", xi, inverse = FALSE){

  # Link function
  inv <- g(link)$inv
  deriv. <- g(link)$deriv.
  deriv.. <- g(link)$deriv..

  # Design matrices
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  p <- NCOL(X)
  k <- NCOL(Z)

  # Relations
  mu <- X%*%par[1 : p]
  phi <- inv(Z%*%par[1:k + p])

  # Matrix representation
  ci <- pmax(1 / ((phi + mu - xi) * (2 + phi + mu - xi)), .Machine$double.eps^(1/2))
  di <- pmax(1 / ((phi - mu + xi) * (2 + phi - mu + xi)), .Machine$double.eps^(1/2))

  D <- diag(as.numeric(1 / deriv.(phi)))

  K <- matrix(NA, p + k, p + k)

  W1 <- diag(as.numeric((ci * (1 + phi + mu - xi) / (1 + phi) +
                           di * (1 + phi - mu + xi) / (1 + phi))))

  W2 <- diag(as.numeric((ci * (1 + phi + mu - xi) / (1 + phi) -
                           di * (1 + phi - mu + xi) /(1 + phi))))

  W3 <- diag(as.numeric((- (1 + phi)^(-2) + ci * (1 + phi + mu - xi) / (1 + phi) +
                           di * (1 + phi - mu + xi) / (1 + phi)) * (1 / deriv.(phi))))

  Kbb <- t(X)%*%W1%*%X
  #Kbg <- t(Z)%*%W2%*%D%*%X
  Kbg <- t(X)%*%W2%*%D%*%Z
  Kgb <- t(Kbg)
  Kgg <- t(Z)%*%W3%*%D%*%Z

  K[1:p,1:p] <- Kbb
  K[1:p, (p + 1):(p + k)] <- Kbg
  K[(p + 1):(p + k), (1:p)] <- Kgb
  K[(p + 1):(p + k), (p + 1):(p + k)] <- Kgg

  if (inverse) {

    K_inv <- matrix(NA, p + k, p + k)

    Kbb_inv <- solve(Kbb)
    Kgg_inv <- solve(Kgg)

    K_inv[1:p,1:p] <- solve(Kbb - Kbg%*%Kgg_inv%*%Kgb)
    K_inv[1:p, (p + 1):(p + k)] <- -Kbb_inv%*%Kbg%*%solve(Kgg - Kgb%*%Kbb_inv%*%Kbg)
    K_inv[(p + 1):(p + k), (1:p)] <- -Kgg_inv%*%Kgb%*%solve(Kbb - Kbg%*%Kgg_inv%*%Kgb)
    K_inv[(p + 1):(p + k), (p + 1):(p + k)] <- solve(Kgg - Kgb%*%Kbb_inv%*%Kbg)

    return(K_inv)

  } else {

    return(K)

  }

}
