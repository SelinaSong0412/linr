library("Rcpp")
library("RcppArmadillo")

#' Function to Fit Linear Regression
#'
#' "linr" is used to fit a simple linear model. In this function, linear regression can be done by three matrix decomposition methods, which are the QR decomposition, Cholesky decomposition and the singular value decomposition (SVD).
#' The defalt fitting method used is the QR method which is used in "lm" function in "stats" package. All three methods can fit linear model with high efficiency.
#'
#' @param formula a symbolic description of the model to be fitted with specific pattern (i.e. y ~ x1 + x2 where y is the responding variable; x1 and x2 are the covariates for the model)
#' @param data an optional data frame, list or environment containing the variables in the model. If not found in data, the variables are taken from formula.
#' @param method an optional character string specifying the fitting method. It must be one of the strings "qr", "cholesky", or "svd".
#' @param
#'
#'
#' @return A list containing the following elements
#'         \itemize{
#'            \item call - return the fitted linear model fomula and the corresponding data.
#'            \item coefficients - a table with the estimated values for each covariates and the intercept as well as their corresponding standard error, t-statistics and (two-sided) p-value.
#'            \item residuals - summary of the usual residuals with min, 1st quantile, median, 3rd quantile, max being computed.
#'            \item fitted.values - the fitted mean values.
#'            \item Residual standard error - the square root of the estimated variance of the random error and a corresponding degrees of freedom will also be computed.
#'            \item R.squared - R^2, the 'fraction of variance explained by the model', SSR (variation in fitted values about the overall mean) / SSY (total variation in Y about its overall mean).
#'            \item Adj.R.squared - adjusted version of R^2, penalized based on number of covariates.
#'            \item SE.beta - the corresponding standard error for the estimated coefficient.
#'            \item t value - the t-statistic for corresponding variable.
#'            \item Pr(>|t|) - the corresponding (two-sided) p-value for the t-statistic
#'            \item F.statistic - the test statistic for F-tests. Variation Between Sample Means / Variation Within the Samples.
#'         }
#'
#' @export
#'
#'
#'

linr <- function(formula, data, method = "cholesky") {

  # Exclude NA in data if exist. This function uses complete data to fit regression model.
  # data = na.omit(data)

  cl = match.call()

  if (hasArg(data)) {mf = model.frame(formula, na.omit(data))}
  else {mf = model.frame(formula)}

  # Defining the outcome matrix Y and the Observation matrix X. Y's elements should be all numeric.
  Y = mf[, 1]
  X = mf[, -1]
  if (!is.numeric(Y)) {stop("The outcomes should be numeric")}

  # convert vector to matrix
  if (is.null(attributes(Y))) {attr(Y, "dim") = c(length(Y), 1)}
  if (is.null(attributes(X))) {attr(X, "dim") = c(length(X), 1)}

  # Checking dimension of X and Y. If wrong, return error.
  if (nrow(Y) != nrow(X)) {
    stop("Number of the outcomes and observations do not match.")
  } else if (nrow(X) < ncol(X)) {
    stop("Number of observations is less than predictors.")
  }

  n = nrow(X)
  p = ncol(X)

# test
# --------------
#   obs = rownames(mf)
#   var = colnames(mf)[-1]
#
#   if (is.null(var)) {
#     var = c("(Intercept)", paste0("x", 1:p))
#   } else {
#     var = c("(Intercept)", var)
#   }
#
#   if (is.null(obs)) {
#     obs = as.character(1:n)
#   }
#
# # ------------------


  # NOTE: should contain value of betahat, MSE, SE.betahat, fit.val, residuals, R.squared, Adj.R.squared, T.stat, p.val

  if (is.null(n)) stop("'x' must be a matrix")
  if (n == 0L) stop("0 (non-NA) cases")
  if (p == 0L) { # The Null model
    return(list(coefficients = numeric(), residuals = Y,
                fitted.values = 0 * Y, rank = 0))
  }

  # High Efficiency Simple Linear Regression.
  if (p == 1L) {
    x = as.vector(X) - mean(X)
    y = as.vector(Y) - mean(Y)
    SSX = sum(x * x)
    SSY = sum(y * y)
    SSXY = sum(x * y)
    R.XY = SSXY / sqrt(SSY * SSX)                           # Corr(X, Y)
    betahat1 = SSXY / SSX
    betahat0 = mean(Y) - betahat1 * mean(X)
    betahat = c(betahat0, betahat1)                         # Coefficient estimates
    fit.val = betahat0 + betahat1 * X                       # Fitted values
    residuals = Y - fit.val                                 # Residuals
    SSE = sum(residuals ^ 2)
    MSE = SSE / (n - 2)                                     # Mean Square Error
    SE.betahat1 = sqrt(MSE / SSX)
    SE.betahat0 = sqrt(MSE * (1/ n + mean(X)^2 / SSX))
    SE.betahat = c(SE.betahat0, SE.betahat1)                # Coefficient Standard Error
    p = p + 1
  }

  else {

    X = cbind(rep(1, n), X)
    p = p + 1

    if (method != "svd" & method != "qr" & method != "cholesky") {
      warning(gettextf("method = '%s' is not supported. Using 'qr'", method), domain = NA)
    }

    else if (method == "svd") { # if using SVD decomposition method
      svd_decom = svd(X, LINPACK = TRUE)
      tuy = crossprod(svd_decom$u, Y)
      betahat = crossprod(t(svd_decom$v), (tuy / svd_decom$d))
    }

    else if (method == "qr") { # if using QR decomposition method

      for (j in 1:p) { # Householder transformation
        id = j:n
        sigma = sum(X[id,j]^2)
        s = sqrt(sigma)
        diag_j = X[j,j]
        gamma = 1 / (sigma + abs(s * diag_j))
        kappa = ifelse(diag_j < 0, s, -s)
        X[j,j] = X[j,j] - kappa
        if (j < p) {
          for (k in (j + 1):p) {
            Y.prime = sum(X[id, j] * X[id, k]) * gamma
            X[id, k] = X[id, k] - X[id, j] * Y.prime
          }
        }
        Y.prime = sum(X[id, j] * Y[id]) * gamma
        Y[id] = Y[id] - X[id, j] * Y.prime
        X[j,j] = kappa
      }

      betahat = rep(NA, p)
      for (j in p:1) { # Using Backward solvin to obtain betahat
        betahat[j] <- Y[j]
        if (j < p)
          for (i in (j + 1):p)
            betahat[j] = betahat[j] - X[j, i] * betahat[i]
          betahat[j] = betahat[j] / X[j, j]
      }
    }

    else { # Else, do Cholesky decomposition as defult.
      tXX = crossprod(X)
      tXY = crossprod(X,Y)
      R = chol(tXX)
      z = forwardsolve(R, tXY, upper.tri = TRUE, transpose = TRUE)
      betahat = backsolve(R, z)
    }

    fit.val = crossprod(t(X), betahat)                      # Fitted value
    residuals = Y - fit.val                                 # Residuals
    SSE = sum(residuals ^ 2)                                # residual sum of square
    SSY = sum((Y - mean(Y)) ^ 2)                            # total sum of square
    MSE = SSE / (n - p)                                     # mean square error
    varmat = solve(crossprod(X), diag(MSE, p, p), tol = 0, transpose = FALSE)
    SE.betahat = sqrt(diag(varmat))
  }

  R.square = 1 - SSE / SSY                                  # Coefficient determination R^2
  Adj.R.square = 1 - (n - 1) * MSE / SSY                    # Adjusted coefficient determination R_adj^2
  T.stat = betahat / SE.betahat                             # T statistics
  p.val.T = pt(abs(T.stat), n - p, lower.tail = FALSE) * 2  # p value of T test
  F.stat = (SSY - SSE) / (p - 1) / MSE                      # F statistics
  p.val.F = pf(F.stat, (p - 1), (n - p), lower.tail = FALSE)# p value of F test

  # Give names to outputs
  observ = rownames(mf)                                     # Observation names
  variab = colnames(mf)[-1]                                 # Variable names
  if (is.null(variab)) {
    variab = c("(Intercept)", paste0("X", 1:p-1))
  } else {
    variab = c("(Intercept)", variab)
  }
  if (is.null(observ)) {
    observ = as.character(1:n)
  }
  names(betahat) = variab
  names(fit.val) = observ
  names(residuals) = observ
  names(SE.betahat) = variab
  names(T.stat) = variab
  names(p.val.T) = variab

  return(list(coefficients = betahat, fitted.values = fit.val,
              residuals = residuals, MSE = MSE, std.error = SE.betahat,
              R.square = R.square, Adj.R.square = Adj.R.square,
              T_statistic = T.stat, p_value.T = p.val.T,
              F_statistic = F.stat, p_value.F = p.val.F))
}



X = matrix(c(1, 5, 3, 4, 5, 6), 3, 2)
Y = matrix(c(1, 2, 3), 3, 1)
betahat = c(1, 2, 3)


# Example for Simple linear regression
y = rnorm(5000000)
x = rnorm(5000000)

# Tests:
round(linr(y~x)$coefficients[0], 10) == round(lm(y~x)$coefficients[0], 10)
round(linr(y~x)$coefficients[1], 10) == round(lm(y~x)$coefficients[1], 10)


# Example for Multiple linear regression
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  n = 10000L
  p = 800L
  q = 5L
  rho = 0.99
} else {
  n = as.integer(args[1])
  p = as.integer(args[2])
  q = as.integer(args[3])
  rho = as.numeric(args[4])
}
X = matrix(rnorm(p * n, sd = sqrt(1 - rho)), nrow = n, ncol = p) + matrix(rnorm(n, sd = sqrt(rho)), nrow = n,ncol = p)
beta = c(rep(c(1, -1), length = q), rep(0, length = p - q))
eps = rnorm(n, sd = 1)
Y = X %*% beta + eps

# Tests:
round(linr(Y~X)$coefficients[1], 10) == round(lm(Y~X)$coefficients[1], 10)
round(linr(Y~X)$coefficients[84], 10) == round(lm(Y~X)$coefficients[84], 10)
round(linr(Y~X)$fitted.value[84], 10) == round(lm(Y~X)$fitted.value[84], 10)
round(linr(Y~X)$fitted.value[999], 10) == round(lm(Y~X)$fitted.value[999], 10)
round(linr(Y~X)$residuals[9898], 10) == round(lm(Y~X)$residuals[9898], 10)










