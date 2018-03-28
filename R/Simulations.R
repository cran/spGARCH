##################################################
#
# Simulation of spatial ARCH and GARCH models
#
##################################################


# spatial ARCH process

sim.spARCH <- function(n = dim(W)[1], rho, alpha, W, b = 2, type = "gaussian", control = list()){

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
    stop("Spatial weight matrix must be quadratic \n")
  }
  if (n != dim(W)[1]){
    stop("Spatial weight matrix must have dimension n x n \n")
  }
  if (class(W) != "dgCMatrix"){
    W <- .asdgCMatrix(rho * as.matrix(W))
  } else {
    W <- rho * W
  }

    set.seed(con$seed)

  if (type == "gaussian"){
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
  if (type == "exp"){
    return(.Simulation_CPP_E(eps, W, alpha, b))
  } else {
    h <- .Simulation_CPP(eps, W, alpha)
  }
  if (type == "gaussian"){
    return(as.numeric(sqrt(h) * eps))
  }
  if (type == "exp"){
    return(as.numeric(exp(h) * eps))
  }
  if (type == "complex"){
    return(as.complex(sqrt(as.complex(h)) * eps))
  }
  stop("invalid type: ", type)
}




