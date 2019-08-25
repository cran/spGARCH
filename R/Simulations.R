##################################################
#
# Simulation of spatial ARCH and GARCH models
#
##################################################


# spatial ARCH process

sim.spARCH <- function(n = dim(W)[1], rho, alpha, W, b = 2, type = "spARCH", control = list()){

  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
    runif(1) # initialize the RNG
  }

  con <- list(seed = runif(1) * 10^9, silent = FALSE, triangular = FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  }

  R.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))

  if (length(n) > 1){
    n <- length(n)
  }
  if (!is.numeric(n)){
    stop("non-numeric argument n. \n")
  }
  if (n < 1){
    stop("n must be greater than 1 \n")
  }
  if (dim(W)[1] != dim(W)[2]){
    stop("Spatial weight matrix W must be quadratic \n")
  }
  if (n != dim(W)[1]){
    stop("Spatial weight matrix W must have dimension n x n \n")
  }
  if (class(W) != "dgCMatrix"){
    W <- .asdgCMatrix(rho * as.matrix(W))
  } else {
    W <- rho * W
  }

  set.seed(con$seed)

  if (type == "spARCH"){
    if(!con$triangular){
      eigen_values <- eigen(W, only.values = TRUE)$values
      if(any(eigen_values != 0)){
        a <- 1 / sqrt(sqrt(norm(W %*% W, "1")))
        eps <- rtruncnorm(n, a = -a, b = a)
      } else {
        eps <- rnorm(n)
      }
    } else {
      eps <- rnorm(n)
    }
  } else {
    eps <- rnorm(n)
  }

  if (con$silent == FALSE){
    cat("Current random seed: ", con$seed, " (To reproduce the results, add argument 'seed = ", con$seed, "' to control argument) \n", sep = "")
  }
  if (type == "log-spARCH"){
    return(.Simulation_CPP_E(eps, W, alpha, b))
  } else {
    h <- .Simulation_CPP(eps, W, alpha)
  }
  if (type == "spARCH"){
    return(as.numeric(sqrt(h) * eps))
  }
  if (type == "log-spARCH"){
    return(as.numeric(exp(h) * eps))
  }
  if (type == "complex-spARCH"){
    return(as.complex(sqrt(as.complex(h)) * eps))
  }
  stop("invalid type: ", type)
}


sim.spGARCH <- function(n = dim(W1)[1], rho, lambda, alpha, W1, W2, b = 2, zeta = 0.5, theta = 0.5, type = "spGARCH", control = list()){

  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
    runif(1) # initialize the RNG
  }

  con <- list(seed = runif(1) * 10^9, silent = FALSE, triangular = FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  }

  R.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))

  if (length(n) > 1){
    n <- length(n)
  }
  if (!is.numeric(n)){
    stop("non-numeric argument n. \n")
  }
  if(!is.element(type, c("spGARCH", "e-spGARCH", "log-spGARCH", "h-spGARCH", "complex-spGARCH"))){
    stop("invalid type: ", type)
  }
  if (n < 1){
    stop("n must be greater than 1 \n")
  }
  if (dim(W1)[1] != dim(W1)[2]){
    stop("Spatial weight matrix W1 must be quadratic \n")
  }
  if (n != dim(W1)[1]){
    stop("Spatial weight matrix W1 must have dimension n x n \n")
  }
  if (dim(W2)[1] != dim(W2)[2]){
    stop("Spatial weight matrix W2 must be quadratic \n")
  }
  if (n != dim(W2)[1]){
    stop("Spatial weight matrix W2 must have dimension n x n \n")
  }
  if (dim(W2)[1] != dim(W1)[1] | dim(W2)[2] != dim(W1)[2]){
    stop("Spatial weight matrices W1 and W2 must have equal dimensions \n")
  }
  if (class(W1) != "dgCMatrix"){
    W1 <- .asdgCMatrix(rho * as.matrix(W1))
  } else {
    W1 <- rho * W1
  }
  if (dim(W2)[1] != dim(W2)[2]){
    stop("Spatial weight matrix W2 must be quadratic \n")
  }
  if (n != dim(W2)[1]){
    stop("Spatial weight matrix W2 must have dimension n x n \n")
  }
  if (class(W2) != "dgCMatrix"){
    W2 <- .asdgCMatrix(lambda * as.matrix(W2))
  } else {
    W2 <- lambda * W2
  }

  set.seed(con$seed)

  functions <- .choose_functions(type)
  f_inv     <- functions$f_inv
  tau_eps   <- functions$tau_eps

  # simulation

  alpha     <- alpha * rep(1, n)

  if (type == "spARCH"){
    if(!con$triangular){
      eigen_values <- eigen(W1, only.values = TRUE)$values
      if(any(eigen_values != 0)){
        a <- 1 / sqrt(sqrt(norm(W1 %*% W1, "1")))
        eps <- rtruncnorm(n, a = -a, b = a)
      } else {
        eps <- rnorm(n)
      }
    } else {
      eps <- rnorm(n)
    }
  } else {
    eps <- rnorm(n)
  }

  tau   <- tau_eps(eps, W1, W2, alpha, theta, zeta, b)
  S     <- solve(diag(n) - W2)
  X     <- S %*% (alpha + W1 %*% tau)

  if(type == "complex-spGARCH"){
    Y     <- diag(as.complex(as.vector(f_inv(X))))^(0.5) %*% eps
  } else {
    Y     <- diag(as.vector(f_inv(X)))^(0.5) %*% eps
  }

  return(as.vector(Y))

}



##################################################
#
# Choosing Functions for Unified Approach
#
##################################################



.choose_functions <- function(model){

  if (model == "spGARCH" | model == "complex-spGARCH")
  {

    f_inv <- function(x){
      return(x)
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      eps_2 <- eps^2
      n     <- length(eps)
      return(diag(eps_2) %*% solve(diag(n) - Matrix(W_1) %*% diag(eps_2) - Matrix(W_2)) %*% alpha)
    }
    tau_y <- function(y){
      return(y^2)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      n   <- length(y)
      X   <- solve(diag(n) - Matrix(lambdaW_2)) %*% (alpha + Matrix(rhoW_1) %*% tau_y(y))
      eps <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    # sourceCpp("Cpp_Functions.cpp")

    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2){
      n <- length(eps)
      aux <- solve(diag(n) - rhoW_1 %*% diag(as.vector(eps^2)) - lambdaW_2)
      result <- array(, dim = c(n,n))
      for(j in 1:n){
        aux_rhoW_1 <- rhoW_1
        aux_rhoW_1[,-j] <- array(0, dim = c(n, n-1))
        eps_j <- eps[j]
        result[, j] <- 2 * eps_j * aux %*% aux_rhoW_1 %*% aux %*% alpha
      }
      return(result)
    }

  }
  else if (model == "log-spGARCH")
  {

    f_inv <- function(x){
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      return(log(abs(eps)^b))
    }
    tau_y <- function(y){
      return(NULL)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      # is only valid if g(eps) = (ln(|eps_1|^b), ..., ln(|eps_n|^b))
      n     <- length(y)
      X     <- solve(diag(n) + 0.5 * b * Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% (alpha + b * Matrix(rhoW_1) %*% log(abs(y)))
      eps   <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2){
      tau_eps_prime <- function(eps){
        return(b / eps) # Achtung evtl. zu ändern, hard coded!
      }
      n     <- length(eps)
      aux   <- solve(diag(n) - Matrix(lambdaW_2)) %*% Matrix(rhoW_1)
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i   <- array(rep(h, n), dim = c(n, n))
      return(aux * tau_eps_prime(eps_j) * h_i)
    }

  }
  else if (model == "e-spGARCH")
  {

    f_inv <- function(x){
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      return(theta * eps)
    }
    tau_y <- function(y){
      return(NULL)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      n     <- length(y)
      function_y_eps <- function(x, y, alpha, rhoW_1, lambdaW_2, theta, zeta){
        n <- length(x)
        eps <- x
        return(    (diag(n) - lambdaW_2) %*% log(abs(y)^2) - alpha     -     ( rhoW_1 %*% (theta * eps) + (diag(n) - lambdaW_2) %*% log(abs(eps)^2) )  )
      }
      out <- nleqslv(x = y, fn = function_y_eps, alpha = alpha, rhoW_1 = rhoW_1, lambdaW_2 = lambdaW_2, theta = theta, zeta = zeta, y = y, control = list(ftol = 1e-10))
      eps <- out$x
      X <- log(abs(y)^2) - log(abs(eps)^2)
      return(list(eps = eps, h =  f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2){
      tau_eps_prime <- function(eps, theta){
        return(theta) # Achtung evtl. zu ändern, hard coded!
      }
      n     <- length(eps)
      aux   <- solve(diag(n) - Matrix(lambdaW_2)) %*% Matrix(rhoW_1)
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i   <- array(rep(h, n), dim = c(n, n))
      return(aux * tau_eps_prime(eps_j, theta) * h_i)
    }

  }
  else if (model == "h-spGARCH")
  {

    f_inv <- function(x){
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      eps_2 <- log(eps^2)
      n     <- length(eps)
      return( solve(diag(n) - Matrix(W_1) - Matrix(W_2)) %*% alpha + (diag(n) + solve(diag(n) - Matrix(W_1) - Matrix(W_2)) %*% Matrix(W_1)) %*% eps_2  )
    }
    tau_y <- function(y){
      log(y^2)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      n   <- length(y)
      X   <- solve(diag(n) - Matrix(lambdaW_2)) %*% (alpha + rhoW_1 %*% tau_y(y))
      eps <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2){
      n     <- length(eps)
      aux   <- solve(diag(n) - Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% Matrix(rhoW_1)
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i   <- array(rep(h, n), dim = c(n, n))
      return(2 * aux / eps_j * h_i)
    }
  }

  return(list(f_inv = f_inv, tau_eps = tau_eps, tau_y = tau_y, g = g, d_h_d_eps = d_h_d_eps))
}


