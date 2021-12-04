#' @title Function to Fit Linear Models with High Efficiency
#'
#' @description "linr" is used to fit a linear model. In this function, linear regression can be done by three matrix decomposition methods, which are the QR decomposition, Cholesky decomposition and the singular value decomposition (SVD). The defalt fitting method used is the Cholesky decomposition method. All three decomposition methods can fit linear model with high efficiency.
#'
#' @details See {\code{vignette("Intro_to_linr", package = "linr")}} for an overview of the package. Also see {\code{vignette("Efficiency_tests", package = "linr")}} for efficiency testing of linr.
#'
#' @param formula an object of class 'formula' (or one that can be coerced to that class): a symbolic description of the model to be fitted. (e.g. Y ~ X + Z, Y is the outcome, X and Z are predictors)
#' @param data an optional data frame, list or environment containing the variables in the model. If not found in data, the variables are defultly taken from formula (environment), typically the environment from which lm is called.
#' @param method an optional character string specifying the fitting method of matrix decomposition. It must be one of the strings "qr", "cholesky", or "svd".
#'
#' @return linr returns an object of class "linr".
#' The functions linr.summary is used to obtain and print a summary and analysis of variance table of the results. The generic accessor functions coefficients, effects, fitted.values and residuals extract various useful features of the value returned by linr.
#' An object of class "lm" is a list containing at least the following components:
#'         \itemize{
#'           \item {Call} - {The fitted linear model fomula and the corresponding data.}
#'           \item {coefficients} - {A vector of coefficients estimates. Containing the estimated regression parameters for intercept and each covariates.}
#'           \item {fitted.values} - {The fitted mean values.}
#'           \item {residuals} - {A vector of the difference between the observed value and the fitted mean values for that observation}
#'           \item {MSE} - {The residual standard error, the square root of the residual sum of squares divided by the residual degrees of freedom. It is a measure used to assess how well a linear regression model fits the data.}
#'           \item {R.squared} - {The coefficient determination, which indicates fraction of variance explained by the fitted model.}
#'           \item {Adj.R.squared} - {A modified version of R-squared that has been adjusted for the number of predictors in the model.}
#'           \item {std.error} - {A vector of standatd error corresponds to each estimated coefficient.}
#'           \item {T_statistic} - {A vector of T-statistic corresponds to each estimated coefficient.}
#'           \item {p_value.T} - {The p-value (two-sided) for the T-statistic}
#'           \item {F_statistic} - {F-statistic, The ratio of the mean regression sum of squares divided by the mean error sum of squares.}
#'           \item {p_value.F} - {The p-value for the F-statistic}
#'         }
#' @seealso  Useful links:
#'        \itemize{
#'          \item {Github page} {\url{https://github.com/SelinaSong0412/linr}}
#'          \item {Report bug at} {\url{https://github.com/SelinaSong0412/linr/issues}}
#'         }
#'
#' @examples
#' Linear_model.lm <- lm(dist~speed,data=cars)
#' print(Linear_model.lm)
#' Linear_model.linr <- linr(dist~speed,data=cars)
#' print(Linear_model.linr$Call)
#' print(Linear_model.linr$Coefficients)
#'
#' @importFrom methods hasArg
#' @importFrom stats na.omit
#' @importFrom stats pf
#' @importFrom stats pt
#' @importFrom stats model.frame
#'
#' @docType package
#' @name linr
#' @export

linr <- function(formula, data, method = "cholesky") {

  cl = match.call()
  if (hasArg(data)) {
    mf = model.frame(formula, na.omit(data))
  } else {
    mf = model.frame(formula)
  }
  # Defining the outcome matrix Y and the Observation matrix X.
  Y = mf[, 1]
  X = as.matrix(mf[, -1])
  colnames(X) = NULL
  rownames(X) = NULL
  if (!is.numeric(Y)) {stop("The outcomes should be numeric")}

  # convert vector to matrix
  if (is.null(attributes(Y))) {attr(Y, "dim") = c(length(Y), 1)}
  if (is.null(attributes(X))) {attr(X, "dim") = c(length(X), 1)}

  # # Checking dimension of X and Y. If wrong, return error.
  # if (nrow(Y) != nrow(X)) {
  #   stop("Number of the outcomes and observations do not match.")
  # } else if (nrow(X) < ncol(X)) {
  #   stop("Number of observations is less than predictors.")
  # }

  n = nrow(Y)
  p = ncol(X) + 1

  # if (n == 0L) stop("0 (non-NA) cases")
  # if (p == 1L) { # The Null model
  #   return(list(coefficients = numeric(), residuals = Y,
  #               fitted.values = 0 * Y, rank = 0))
  # }

  if (p == 2L) { # High Efficient Simple Linear Regression.
    x = as.vector(X) - mean(X)
    y = as.vector(Y) - mean(Y)
    SSX = sum(x * x)
    SSY = sum(y * y)
    SSXY = sum(x * y)
    R.XY = SSXY / sqrt(SSY * SSX)                           # Corr(X, Y)
    betahat1 = SSXY / SSX
    betahat0 = mean(Y) - betahat1 * mean(X)
    betahat = c(betahat0, betahat1)                         # Coefficient estimates
    fit.val = betahat0 + betahat1 * as.vector(X)            # Fitted values
    residuals = as.vector(Y) - fit.val                      # Residuals
    SSE = sum(residuals ^ 2)
    MSE = SSE / (n - 2)                                     # Mean Square Error
    SE.betahat1 = sqrt(MSE / SSX)
    SE.betahat0 = sqrt(MSE * (1 / n + mean(X) ^ 2 / SSX))
    SE.betahat = c(SE.betahat0, SE.betahat1)                # Coefficient Standard Error
  }

  else { # High Efficient Multiple Linear Regression partly wrote with Rcpp

    X = cbind(rep(1, n), X)

    if (method != "svd" & method != "qr" & method != "cholesky") {
      warning(gettextf("method = '%s' is not supported. Using 'qr', 'svd' or 'cholesky'", method), domain = NA)
    }

    else if (method == "svd") { # if using SVD decomposition method
      svd_decom = svd(X)
      tuy = crossprod(svd_decom$u, Y)
      betahat = crossprod(t(svd_decom$v), (tuy / svd_decom$d))
    }

    else if (method == "qr") { # if using QR decomposition method
      res_qr <- qr(X, lapack=TRUE)
      betahat = qr.coef(res_qr, Y)
    }

    else { # Else, do Cholesky decomposition as defult.
      XtX = crossprod(X)
      XtY = crossprod(X,Y)
      R = chol(XtX)
      z = forwardsolve(R, XtY, upper.tri = TRUE, transpose = TRUE)
      betahat = backsolve(R, z)
    }

    fit.val = as.vector(crossprod(t(X), betahat))           # Fitted value
    residuals = as.vector(Y) - fit.val                      # Residuals
    SSE = sum(residuals ^ 2)                                # residual sum of square
    SSY = sum((Y - mean(Y)) ^ 2)                            # total sum of square
    MSE = SSE / (n - p)                                     # mean square error
    varmat = solve(crossprod(X), diag(MSE, p, p), tol = 0, transpose = FALSE)
    SE.betahat = sqrt(diag(varmat))                         # standard error of estimates
    betahat = as.vector(betahat)
  }

  R.square = 1 - SSE / SSY                                  # Coefficient determination R^2
  Adj.R.square = 1 - (n - 1) * MSE / SSY                    # Adjusted coefficient determination R_adj^2
  T.stat = betahat / SE.betahat                             # T statistics
  p.val.T = pt(abs(T.stat), n - p, lower.tail = FALSE) * 2  # p value of T test
  F.stat = (SSY - SSE) / (p - 1) / MSE                      # F statistics
  p.val.F = pf(F.stat, (p - 1), (n - p), lower.tail = FALSE)# p value of F test

  # Give names to variables and observations.
  vars = colnames(mf)[-1]
  if (hasArg(data)) {                                       # If have input data names
    variab = c("(Intercept)", vars)
    observ = rownames(mf)
  } else {                                                  # If not have input data names
    if (p > 2) {
      variab = c("(Intercept)", paste0(vars, 1:(p - 1)))
    } else {
      variab = c("(Intercept)", vars)
    }
    observ = 1:n
  }

  names(betahat) = variab
  names(fit.val) = observ
  names(residuals) = observ
  names(SE.betahat) = variab
  names(T.stat) = variab
  names(p.val.T) = variab

  fit = list(Call = cl, coefficients = betahat, fitted.values = fit.val,
              residuals = residuals, MSE = MSE, std.error = SE.betahat,
              R.square = R.square, Adj.R.square = Adj.R.square,
              T_statistic = T.stat, p_value.T = p.val.T,
              F_statistic = F.stat, p_value.F = p.val.F)
  class(fit) = "linr"
  return(fit)
}









