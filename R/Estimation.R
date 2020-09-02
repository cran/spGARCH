##################################################
#
# Estimation of spatial ARCH and GARCH models
#
##################################################


# spatial ARCH process


qml.spARCH <- function(formula, W, type = "spARCH", data = NULL, b = 2, start = NULL, control = list()){
  cl <- match.call()
  con <- list(trace = FALSE, rho = 1, outer.iter = 400, inner.iter = 800, delta = 1.0e-7, tol = 1.0e-8)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  if(!is.element(type, c("spARCH", "log-spARCH"))){
    stop("invalid type: ", type)
  }

  model_frame <- model.frame(formula, data)
  model_terms <- attr(model_frame, "terms")
  y <- model.response(model_frame)
  X <- model.matrix(model_terms, model_frame)

  if (length(dim(y)) > 0 | !inherits(y, "numeric") | mode(y) != "numeric"){
    stop("y must be a numeric vector \n")
  }

  if (dim(W)[1] != dim(W)[2]){
    stop("Spatial weight matrix must be quadratic \n")
  }
  if (length(y) != dim(W)[1]){
    stop(paste("Spatial weight matrix must have dimension", length(y), "x", length(y), "(length of y) \n"))
  }
  if (!inherits(W, "dgCMatrix")){
    W <- .asdgCMatrix(as.matrix(W))
  }

  if(is.null(names(y))){
    names(y) <- as.character(1:length(y))
  }

  if (type == "spARCH"){

    if(dim(X)[2] == 0){

      param <- list(y, W, isSymmetric(W)) # parameters and regressors

      LB <- c(1e-6, 0)
      UB <- c(Inf, Inf)
      default <- c(0.4, 1)
      start <- .teststart(start, LB, UB, default)

      snp <- solnp(pars = start, fun = .LL_spARCH, param = param, # .LL_spARCH_oriented_r
                   LB = LB,
                   UB = UB,
                   control =  con)

    } else if (dim(X)[2] > 0) {

      param <- list(y, W, X, isSymmetric(W)) # parameters and regressors

      LB <- c(1e-6, 0, rep(-Inf, dim(X)[2]))
      UB <- c(Inf, Inf, rep(Inf, dim(X)[2]))
      default <- c(0.4, 1, rep(0, dim(X)[2]))
      start <- .teststart(start, LB, UB, default)

      snp <- solnp(pars = start, fun = .LL_spARCHX, param = param,
                   LB = LB,
                   UB = UB,
                   control =  con)

    } else {
      stop("X must be a matrix")
    }

  }

  if (type == "log-spARCH"){

    if(dim(X)[2] == 0){

      param <- list(y, W, b) # parameters and regressors

      LB <- c(1e-6, 0)
      UB <- c(Inf, Inf)
      default <- c(0.4, 1)
      start <- .teststart(start, LB, UB, default)

      snp <- solnp(pars = start, fun = .LL_EspARCH, param = param,
                   LB = LB,
                   UB = UB,
                   control =  con)

    } else if (dim(X)[2] > 0) {

      param <- list(y, W, b, X) # parameters and regressors

      LB <- c(1e-6, 0, rep(-Inf, dim(X)[2]))
      UB <- c(Inf, Inf, rep(Inf, dim(X)[2]))
      default <- c(0.4, 1, rep(0, dim(X)[2]))
      start <- .teststart(start, LB, UB, default)

      snp <- solnp(pars = start, fun = .LL_EspARCHX, param = param,
                   LB = LB,
                   UB = UB,
                   control =  con)

    } else {
      stop("X must be a matrix")
    }

  }

  if (snp$convergence != 0){
    stop("Iterative maximization did not converge. \n")
  }

  stderr <- tryCatch(sqrt(diag(solve(snp$hessian))),
                     error = function(e) {warning("Hessian matrix is computationally singular (no standard errors)\n"); return(rep(NA, length(snp$pars)))})

  if(dim(X)[2] == 0){

    h <- as.vector(snp$pars[1] + snp$pars[2] * W %*% (y^2))
    residuals <- as.vector(y/sqrt(h))
    names(residuals) <- names(y)
    out <- list(coefficients = snp$pars, residuals = residuals, fitted.values = y - residuals,
                terms = model_terms, df.residual = length(residuals) - length(snp$pars),
                stderr = stderr, hessian = snp$hessian,
                LL = -snp$values[length(snp$values)],
                y = y,  h = h, type = type, W = W, call = cl, regressors = FALSE, AR = FALSE)
  } else {

    xi <- y - X %*% snp$pars[3:length(snp$pars)]
    h <- as.vector(snp$pars[1] + snp$pars[2] * W %*% (xi^2))
    residuals <- as.vector(xi/sqrt(h))
    names(residuals) <- names(y)
    out <- list(coefficients = snp$pars, residuals = xi, fitted.values = y - xi,
                terms = model_terms, df.residual = length(xi) - length(snp$pars),
                stderr = stderr, hessian = snp$hessian,
                LL = -snp$values[length(snp$values)],
                y = y,  h = h, type = type, W = W, call = cl, regressors = TRUE, X = X, AR = FALSE)

  }


  return(structure(out, class = c("spARCH")))
}


# qml.spGARCH <- function(formula, W1, W2, type = "spGARCH", data = NULL, b = 2, start = NULL, control = list()){
#
#   # to-do
#   # add possibility for fixed parameters (eGARCH, logGARCH)
#   # parameter space
#   # X-models
#   # spGARCH object, summary
#   # documentation
#   # cpp functions
#
#
#   cl <- match.call()
#   con <- list(trace = FALSE, rho = 1, outer.iter = 400, inner.iter = 800, delta = 1.0e-7, tol = 1.0e-8)
#   nmsC <- names(con)
#   con[(namc <- names(control))] <- control
#   if (length(noNms <- namc[!namc %in% nmsC]))
#     warning("unknown names in control: ", paste(noNms, collapse = ", "))
#
#   if(!is.element(type, c("spGARCH", "log-spGARCH", "e-spGARCH"))){
#     stop("invalid type: ", type)
#   }
#
#   model_frame <- model.frame(formula, data)
#   model_terms <- attr(model_frame, "terms")
#   y <- model.response(model_frame)
#   X <- model.matrix(model_terms, model_frame)
#
#   if (length(dim(y)) > 0 | !is(y, "numeric") | mode(y) != "numeric"){
#     stop("y must be a numeric vector \n")
#   }
#
#   if (dim(W1)[1] != dim(W1)[2] | dim(W2)[1] != dim(W2)[2]){
#     stop("Spatial weight matrices must be quadratic \n")
#   }
#
#   if (length(y) != dim(W1)[1] | length(y) != dim(W2)[1]){
#     stop(paste("Spatial weight matrices must have dimension", length(y), "x", length(y), "(length of y) \n"))
#   }
#   if (!is(W1, "dgCMatrix")){
#     W1 <- .asdgCMatrix(as.matrix(W))
#   }
#   if (!is(W2, "dgCMatrix")){
#     W2 <- .asdgCMatrix(as.matrix(W))
#   }
#
#   if(is.null(names(y))){
#     names(y) <- as.character(1:length(y))
#   }
#
#
#   functions <- .choose_functions(type)
#
#   f_inv     <- functions$f_inv
#   tau_y     <- functions$tau_y
#   g         <- functions$g
#   d_h_d_eps <- functions$d_h_d_eps
#
#
#
#   if(dim(X)[2] == 0){
#
#     param <- list(y = y, W_1 = W_1, W_2 = W_2,
#                   f_inv = f_inv, tau_y = tau_y, g = g, d_h_d_eps = d_h_d_eps,
#                   model = type) # parameters and regressors
#
#     if (type == "spGARCH"){
#       LB <- c(1e-6, 0  , 0  )
#       UB <- c(Inf , Inf, Inf)
#       default <- c(1, 0.5, 0.3)
#       start <- .teststart(start, LB, UB, default)
#     }
#     if (type == "e-spGARCH"){
#       LB <- c(1e-6, 0  , 0  , 0  , -Inf )
#       UB <- c(Inf , Inf, Inf, Inf,  Inf )
#       default <- c(1, 0.5, 0.3, 1, 0)
#       start <- .teststart(start, LB, UB, default)
#     }
#     if (type == "log-spGARCH"){
#       LB <- c(1e-6, 0  , 0  , 0  )
#       UB <- c(Inf , Inf, Inf, Inf)
#       default <- c(1, 0.5, 0.3, 1)
#       start <- .teststart(start, LB, UB, default)
#     }
#
#     snp <- solnp(pars = start, fun = .LL_unified, param = param, # .LL_spARCH_oriented_r
#                  LB = LB,
#                  UB = UB,
#                  control =  con)
#   } else if (dim(X)[2] > 0) {
#
#     # param <- list(y, W, X, isSymmetric(W)) # parameters and regressors
#     #
#     # LB <- c(1e-6, 0, rep(-Inf, dim(X)[2]))
#     # UB <- c(Inf, Inf, rep(Inf, dim(X)[2]))
#     # default <- c(0.4, 1, rep(0, dim(X)[2]))
#     # start <- .teststart(start, LB, UB, default)
#     #
#     # snp <- solnp(pars = start, fun = .LL_spARCHX, param = param,
#     #              LB = LB,
#     #              UB = UB,
#     #              control =  con)
#
#   } else {
#     stop("X must be a matrix")
#   }
#
#
#   if (snp$convergence != 0){
#     stop("Iterative maximization did not converge. \n")
#   }
#
#   stderr <- tryCatch(sqrt(diag(solve(snp$hessian))),
#                      error = function(e) {warning("Hessian matrix is computationally singular (no standard errors)\n"); return(rep(NA, length(snp$pars)))})
#
#   if(dim(X)[2] == 0){
#
#     h <- as.vector(snp$pars[1] + snp$pars[2] * W %*% (y^2))
#     residuals <- as.vector(y/sqrt(h))
#     names(residuals) <- names(y)
#     out <- list(coefficients = snp$pars, residuals = residuals, fitted.values = y - residuals,
#                 terms = model_terms, df.residual = length(residuals) - length(snp$pars),
#                 stderr = stderr, hessian = snp$hessian,
#                 LL = -snp$values[length(snp$values)],
#                 y = y,  h = h, type = type, W = W, call = cl, regressors = FALSE, AR = FALSE)
#   } else {
#
#     xi <- y - X %*% snp$pars[3:length(snp$pars)]
#     h <- as.vector(snp$pars[1] + snp$pars[2] * W %*% (xi^2))
#     residuals <- as.vector(xi/sqrt(h))
#     names(residuals) <- names(y)
#     out <- list(coefficients = snp$pars, residuals = xi, fitted.values = y - xi,
#                 terms = model_terms, df.residual = length(xi) - length(snp$pars),
#                 stderr = stderr, hessian = snp$hessian,
#                 LL = -snp$values[length(snp$values)],
#                 y = y,  h = h, type = type, W = W, call = cl, regressors = TRUE, X = X, AR = FALSE)
#
#   }
#
#
#   return(structure(out, class = c("spGARCH")))
# }





# SARspARCH process

qml.SARspARCH <- function(formula, B, W, type = "spARCH", data = NULL, b = 2, start = NULL, eigen_v = NULL, control = list()){

  cl <- match.call()
  con <- list(trace = FALSE, rho = 1, outer.iter = 400, inner.iter = 800, delta = 1.0e-7, tol = 1.0e-8)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  if(!is.element(type, c("spARCH", "log-spARCH"))){
    stop("invalid type: ", type)
  }

  model_frame <- model.frame(formula, data)
  model_terms <- attr(model_frame, "terms")
  y <- model.response(model_frame)
  X <- model.matrix(model_terms, model_frame)

  if (length(dim(y)) > 0 | !inherits(y, "numeric") | mode(y) != "numeric"){
    stop("'y' must be a numeric vector \n")
  }

  if (dim(B)[1] != dim(B)[2] | dim(W)[1] != dim(W)[2]){
    stop("Spatial weighting matrices must be quadratic \n")
  }
  if (length(y) != dim(B)[1] | length(y) != dim(W)[1]){
    stop(paste("Spatial weight matrix must have dimension", length(y), "x", length(y), "(length of y) \n"))
  }
  if (!inherits(B, "dgCMatrix")){
    B <- .asdgCMatrix(as.matrix(B))
  }
  if (!inherits(W, "dgCMatrix")){
    W <- .asdgCMatrix(as.matrix(W))
  }

  if(is.null(eigen_v)){
    eigen <- eigen(B, only.values = TRUE)
    B_eigen <- eigen$values
  } else {
    eigen <- list(values = eigen_v, vector = NULL)
    B_eigen <- eigen_v
  }
  if(min(eigen$values) == 0){
    min_lambda <- -Inf
  } else {
    min_lambda <- 1 / min(eigen$values) + 1e-6
  }
  if(max(eigen$values) == 0){
    max_lambda <- Inf
  } else {
    max_lambda <- 1 / max(eigen$values) - 1e-6
  }

  if(is.null(names(y))){
    names(y) <- as.character(1:length(y))
  }

  if (type == "spARCH"){

    if(dim(X)[2] == 0){

      param <- list(y, B, W, B_eigen, isSymmetric(W)) # parameters and regressors

      LB <- c(1e-6, 0, min_lambda)
      UB <- c(Inf, Inf, max_lambda)
      st_lambda <- ifelse(any(!is.finite(c(min_lambda, max_lambda))), 0, (UB[3] + LB[3])/2)
      default <- c(0.4, 0.4, st_lambda)
      start <- .teststart(start, LB, UB, default)

      snp <- solnp(pars = start, fun = .LL_SARspARCH, param = param,
                   LB = LB,
                   UB = UB,
                   control =  con)

    } else if (dim(X)[2] > 0) {

      param <- list(y, B, W, X, B_eigen, isSymmetric(W)) # parameters and regressors

      LB <- c(1e-6, 0, min_lambda, rep(-Inf, dim(X)[2]))
      UB <- c(Inf, Inf, max_lambda, rep(Inf, dim(X)[2]))
      st_lambda <- ifelse(any(!is.finite(c(min_lambda, max_lambda))), 0, (UB[3] + LB[3])/2)
      default <- c(0.4, 0.4, st_lambda, rep(0, dim(X)[2]))
      start <- .teststart(start, LB, UB, default)

      snp <- solnp(pars = start, fun = .LL_SARspARCHX, param = param,
                   LB = LB,
                   UB = UB,
                   control =  con)

    } else {
      stop("X must be a matrix")
    }

  }

  if (type == "log-spARCH"){

    if(dim(X)[2] == 0){

      param <- list(y, B, W, b, B_eigen) # parameters and regressors

      LB <- c(1e-6, 0, min_lambda)
      UB <- c(Inf, Inf, max_lambda)
      st_lambda <- ifelse(any(!is.finite(c(min_lambda, max_lambda))), 0, (UB[3] + LB[3])/2)
      default <- c(0.4, 0.4, st_lambda)
      start <- .teststart(start, LB, UB, default)

      snp <- solnp(pars = start, fun = .LL_SAREspARCH, param = param,
                   LB = LB,
                   UB = UB,
                   control =  con)

    } else if (dim(X)[2] > 0) {

      param <- list(y, B, W, b, X, B_eigen) # parameters and regressors

      LB <- c(1e-6, 0, min_lambda, rep(-Inf, dim(X)[2]))
      UB <- c(Inf, Inf, max_lambda, rep(Inf, dim(X)[2]))
      st_lambda <- ifelse(any(!is.finite(c(min_lambda, max_lambda))), 0, (UB[3] + LB[3])/2)
      default <- c(0.4, 0.4, st_lambda, rep(0, dim(X)[2]))
      start <- .teststart(start, LB, UB, default)

      snp <- solnp(pars = start, fun = .LL_SAREspARCHX, param = param,
                   LB = LB,
                   UB = UB,
                   control =  con)

    } else {
      stop("X must be a matrix")
    }

  }

  if (snp$convergence != 0){
    stop("Iterative maximization did not converge. \n")
  }

  stderr <- tryCatch(sqrt(diag(solve(snp$hessian))),
                     error = function(e) {warning("Hessian matrix is computationally singular (no standard errors)\n"); return(rep(NA, length(snp$pars)))})

  if(dim(X)[2] == 0){

    xi <- as.vector(y - snp$pars[3] * B %*% y)
    h <- as.vector(snp$pars[1] + snp$pars[2] * W %*% (xi^2))
    residuals <- as.vector(xi/sqrt(h))
    names(residuals) <- names(y)
    out <- list(coefficients = snp$pars, residuals = xi, fitted.values = y - xi,
                stderr = stderr, hessian = snp$hessian,
                LL = -snp$values[length(snp$values)],
                y = y,  h = h, type = type, W = W, B = B, call = cl, regressors = FALSE, AR = TRUE)
  } else {

    xi <- as.vector(y - snp$pars[3] * B %*% y - X %*% snp$pars[4:length(snp$pars)])
    h <- as.vector(snp$pars[1] + snp$pars[2] * W %*% (xi^2))
    residuals <- as.vector(xi/sqrt(h))
    names(residuals) <- names(y)
    out <- list(coefficients = snp$pars, residuals = xi, fitted.values = y - xi,
                terms = model_terms, df.residual = length(xi) - length(snp$pars),
                stderr = stderr, hessian = snp$hessian,
                LL = -snp$values[length(snp$values)],
                y = y,  h = h, type = type, W = W, B = B, call = cl, regressors = TRUE, X = X, AR = TRUE)

  }

  return(structure(out, class = c("spARCH")))
}

##################################################
#
# Generic Functions for spARCH objects
#
##################################################

spARCH <- setClass("spARCH", slots = c(coefficients = "numeric", residuals = "numeric", fitted.values = "numeric",
                                       terms = "formula", df.residual = "numeric",
                                       stderr = "numeric", hessian = "matrix",
                                       LL = "numeric",
                                       y = "numeric", h = "numeric", type = "character", W = "dgCMatrix", B = "dgCMatrix", call = "call",
                                       regressors = "logical", X = "matrix", AR = "logical"))




summary.spARCH <- function(object, ...){

  # Coefficients

  coef <- array(, dim = c(length(object$coefficients), 4))
  coef[,1] <- object$coefficients
  coef[,2] <- object$stderr
  coef[,3] <- coef[,1] / coef[,2]
  coef[,4] <- 2 * (1 - pnorm(abs(coef[,3])))
  colnames(coef) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  if(object$AR){
    par_names <- c("alpha (spARCH)", "rho (spARCH)", "lambda (SAR)")
  } else {
    par_names <- c("alpha", "rho")
  }
  if(object$regressors){
    rownames(coef) <- c(par_names, colnames(object$X))
  } else {
    rownames(coef) <- c(par_names)
  }


  object$coef <- coef

  # Information criteria

  object$AIC <- 2 * length(object$coefficients) - 2 * object$LL
  object$BIC <- log(length(residuals)) * length(object$coefficients) - 2 * object$LL

  # Spatial autocorrelation of the residuals

  if(object$AR){
    object$moran_res <- moran.test(object$residuals, mat2listw(object$B), zero.policy = TRUE, alternative = "two.sided")
  } else {
    object$moran_res <- moran.test(object$residuals, mat2listw(object$W), zero.policy = TRUE, alternative = "two.sided")
  }
  object$moran_sq_res <- moran.test(object$residuals^2, mat2listw(object$W), zero.policy = TRUE, alternative = "two.sided")

  # return

  return(structure(object, class = c("summary.spARCH", class(object))))
}

print.spARCH <- function(x, ...){
  cat("Use the summary function to print a detailed summary of the estimation results. \n ?summary.spARCH \n")
  return(str(x))
}

print.summary.spARCH <- function(x, digits = max(5, .Options$digits - 3), signif.stars = TRUE, ...){

  # Call
  cat("\n Call: \n")
  print(x$call)

  # Residuals
  cat("\n Residuals: \n")
  print(summary(as.vector(x$residuals)), digits = digits)

  # Moran's I of residuals and squared residuals - TO_DO
  # Coefficients
  cat("\n Coefficients: \n")
  printCoefmat(x$coef, signif.stars = signif.stars, digits = digits, na.print = "NA")
  cat("\n Note: all p-values (incl. significance stars) are asymptotic \n")
  # Further stuff
  cat("\n AIC: ", format(signif(x$AIC, digits)), ", BIC: ", format(signif(x$BIC, digits)), " (Log-Likelihood: ", format(signif(x$LL, digits)), ") \n", sep = "")
  cat("\n Moran's I (residuals): ",  format(signif(x$moran_res$estimate[1], digits)), ", p-value: ", format(signif(x$moran_res$p.value, digits)), "\n", sep = "")
  cat("\n Moran's I (squared residuals): ",  format(signif(x$moran_sq_res$estimate[1], digits)), ", p-value: ",  format(signif(x$moran_sq_res$p.value, digits)), "\n", sep = "")

  cat("\n")
  invisible(x)
}

logLik.spARCH <- function(object, ...){

  LL <- c(object$LL)
  class(LL) <- "logLik"
  N <- length(object$residuals)
  attr(LL, "nall") <- N
  attr(LL, "nobs") <- N
  attr(LL, "df") <- length(object$coefficients)

  return(LL)
}

extractAIC.spARCH <- function(fit, scale = 0, k = 2, ...){

  n <- length(fit$residuals)
  edf <- n - fit$df.residual
  LL <- logLik.spARCH(fit)

  return(c(edf, - 2 * LL + k * edf))
}

fitted.spARCH <- function(object, ...){
  return(object$fitted.values)
}

residuals.spARCH <- function(object, ...){
  return(object$residuals)
}

plot.spARCH <- function (x, which = c(1:3),
            ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...,
            qqline = TRUE){

  exp_dots <- match.call(expand.dots = FALSE)$`...`

  if (!inherits(x, "spARCH")){
    stop("use only with \"spARCH\" objects")
  }
  if(!is.numeric(which) || any(which < 1) || any(which > 3)){
    stop("'which' must be between 1 and 3")
  }
  show <- rep(FALSE, 3)
  show[which] <- TRUE

  one.fig <- prod(par("mfcol")) == 1
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  xlab <- c("Residuals", "Squared Residuals")
  ylab <- c("Spatially Lagged Residuals", "Spatially Lagged Squared Residuals", "Standardized Residuals")
  zero <- TRUE

  if(show[1]){
    dev.hold()

    eps <- residuals(x)
    W <- as.matrix(x$W)
    moran.plot(eps, mat2listw(W), zero.policy = zero, xlab = xlab[1], ylab = ylab[1], ...)

    cat("Reproduce the results as follows: \n")
    cat("\t eps <- residuals(x) \n\t W <- as.matrix(x$W) \n\t moran.plot(eps, mat2listw(W), zero.policy = TRUE, xlab = \"Residuals\", ylab = \"Spatially Lagged Residuals\") \n")

    dev.flush()
  }
  if(show[2]){
    dev.hold()

    eps <- residuals(x)
    W <- as.matrix(x$W)
    moran.plot(eps^2, mat2listw(W), zero.policy = zero, xlab = xlab[2], ylab = ylab[2], ...)

    cat("Reproduce the results as follows: \n")
    cat("\t eps <- residuals(x) \n\t W <- as.matrix(x$W) \n\t moran.plot(eps^2, mat2listw(W), zero.policy = TRUE, xlab = \"Residuals\", ylab = \"Spatially Lagged Residuals\") \n")

    dev.flush()
  }
  if(show[3]){
    dev.hold()

    eps <- residuals(x)
    std_eps <- (eps - mean(eps))/sd(eps)
    qqnorm(eps, ylab = ylab[3], ...)

    cat("Reproduce the results as follows: \n")
    cat("\t eps <- residuals(x) \n\t std_eps <- (eps - mean(eps))/sd(eps) \n\t qqnorm(eps, ylab = \"Standardized Residuals\") \n")

    if(qqline){
      qqline(eps)
      cat("\t qqline(eps) \n")
    }


    dev.flush()
  }

  invisible()
}

##################################################
#
# Other stuff
#
##################################################

.teststart <- function(start, LB, UB, default){
  if(!is.null(start)){
    if(length(start) != length(LB)){
      stop(paste("For the chosen model, the vector of starting values must have length ", length(UB), ". Current size: ", length(start), sep = ""))
    }
    if(any(start > UB | start < LB)){
      if(any(start < LB) & all(start < UB)){
        stop(paste("The ", .ord(which(start < LB)), " starting value(s) is/are smaller than the lower bound (", LB[which(start < LB)],")\n", sep = ""))
      }
      if(any(start > UB) & all(start > LB)){
        stop(paste("The ", .ord(which(start > UB)), " starting value(s) is/are larger than the upper bound (", UB[which(start > UB)],")\n", sep = ""))
      }
      if(any(start > UB) & any(start < LB)){
        stop(paste("The ", .ord(which(start > UB)), " starting value(s) is/are larger than the upper bound (", UB[which(start > UB)],"). ",
                   "Moreover, the ", .ord(which(start < LB)), " starting value(s) is/are smaller than the lower bound (", LB[which(start < LB)],")\n", sep = ""))
      }
    }
    return(start)
  } else {
    return(default)
  }
}

.ord <- function(numbers){
  rules_inner <- function(x){
    if(x%%10 == 1 & x != 11){
      return(paste(as.character(x), "st", sep = ""))
    }
    if(x%%10 == 2 & x != 12){
      return(paste(as.character(x), "nd", sep = ""))
    }
    if(x%%10 == 3 & x != 13){
      return(paste(as.character(x), "rd", sep = ""))
    } else {
      return(paste(as.character(x), "th", sep = ""))
    }
  }
  return(unlist(lapply(numbers, rules_inner)))
}


##################################################
#
# Log-Likelihood Functions
#
##################################################

.LL_unified <- function(pars, param){

  model  <- param$model

  if(is.element(model, c("spGARCH", "h-spGARCH"))){
    alpha  <- pars[1]
    rho    <- pars[2]
    lambda <- pars[3]
    theta  <- 1
    zeta   <- 1
    b      <- 1
  } else if(is.element(model, c("log-spGARCH"))) {
    alpha  <- pars[1]
    rho    <- pars[2]
    lambda <- pars[3]
    theta  <- NULL
    zeta   <- NULL
    b      <- pars[4]
  } else if(is.element(model, c("e-spGARCH"))) {
    alpha  <- pars[1]
    rho    <- pars[2]
    lambda <- pars[3]
    theta  <- pars[4]
    zeta   <- pars[5]
    b      <- NULL
  }


  y         <- param$y
  f_inv     <- param$f_inv
  tau_y     <- param$tau_y
  W_1       <- param$W_1
  W_2       <- param$W_2
  g         <- param$g
  d_h_d_eps <- param$d_h_d_eps

  n            <- length(y)
  rhoW_1       <- rho * W_1
  lambdaW_2    <- lambda * W_2
  alpha        <- alpha * rep(1, n)

  result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
  eps          <- as.vector(result_g[[1]])
  h            <- as.vector(result_g[[2]])

  J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))

  return((-1) * (sum(log(dnorm(eps))) - log(det(J))))

}




#
# .LL_spARCHX_r <- function(pars, param){
#   # spARCH coefficients
#   alpha <- pars[1]
#   rho <- pars[2]
#   # regression coefficients
#   beta <- pars[3:length(pars)]
#   # data and regressors
#   Y <- param[[1]] # observations Y
#   X <- param[[3]] # matrix of regressors
#   W <- param[[2]] # weighting matrix W_tilde (non-triangular) (spARCH term)
#   n <- length(Y) # number of observations
#   # computation of the log likelihood
#   Xbeta <- X %*% beta
#   xi <- Y - Xbeta
#   h <- as.vector(alpha + rho * W %*% (xi^2))
#   eps <- as.vector(xi/sqrt(h))
#   # if W is not symmetric, the following line must be changed: symmetric = FALSE
#   eigen_v <- eigen(diag(as.vector(h/xi^2)) - t(rho * W), only.values = TRUE, symmetric = FALSE)
#   log_abs_det_partial <- sum( log(xi^2 / h^(3/2)))  + sum(log(abs(eigen_v$values))) # log determinant of the spARCH term (Jacobian)
#   sum_f_residuals <- sum(log(dnorm(eps))) # sum of the logarithmic residuals
#   # negative log-likelihood function
#   # (due to the minimization algorithm 'solnp' the LL is multiplied by (-1))
#   return( -1 * (sum_f_residuals + log_abs_det_partial) )
# }
#
# .LL_spEARCH_r <- function(pars, param){
#   # spARCH coefficients
#   alpha <- pars[1]
#   rho <- pars[2]
#   # data and regressors
#   Y <- param[[1]] # observations Y
#   W <- param[[2]] # weighting matrix W_tilde (non-triangular) (spARCH term)
#   b <- param[[3]] # parameter b
#   n <- length(Y) # number of observations
#   # computation of the log likelihood
#   S <- solve(diag(n) + 0.5 * b * rho * W)
#   lnh <- S %*% (alpha + b * rho * W %*% log(abs(Y)))
#   h <- exp(lnh)
#   eps <- as.vector(Y/sqrt(h))
#   # if W is not symmetric, the following line must be changed: symmetric = FALSE
#   eigen_v <- eigen(diag(as.vector(2 * h / b)) - t(rho * W * S), only.values = TRUE, symmetric = FALSE)
#   log_abs_det_partial <- sum( log(b / (2*h^(3/2))))  + sum(log(abs(eigen_v$values))) # log determinant of the spARCH term (Jacobian)
#   sum_f_residuals <- sum(log(dnorm(eps))) # sum of the logarithmic residuals
#   # negative log-likelihood function
#   # (due to the minimization algorithm 'solnp' the LL is multiplied by (-1))
#   # return(log_abs_det_partial)
#   return( -1 * (sum_f_residuals + log_abs_det_partial) )
# }
#
#
# .LL_spARCH_oriented_r <- function(theta, param){ #param y, By
#   alpha <- theta[1]
#   rho <- theta[2]
#
#   y <- param[[1]]
#   W <- param[[2]]
#   Wy2 <- W %*% (y^2)
#   N <- length(y)
#   h <- as.vector(alpha + rho * Wy2)
#   xi <- y/sqrt(h)
#   return( -1 * ( -N/2 * log(2 * pi) - 0.5*sum(log(h)) - 1/2*sum(xi^2) )    ) #
# }
#
# .LL_SARspARCH_triangular <- function(pars, param){
#   lambda <- pars[1]
#   alpha <- pars[2]
#   rho <- pars[3]
#   beta <- pars[4:length(pars)]
#   Y <- param[[1]]
#   BY <- param[[2]]
#   B <- param[[3]]
#   X <- param[[4]]
#   W <- param[[5]]
#   n <- length(Y)
#   Xbeta <- X %*% beta
#   eps <- Y - lambda * BY - Xbeta
#   h <- as.vector(alpha + rho * W %*% (eps^2))
#   xi <- eps/sqrt(h)
#   return( -1 * ( -n/2 * log(2 * pi) - 0.5*sum(log(h)) + log(det(diag(n) - lambda * B)) - 1/2*sum(xi^2) )    ) #
# }
