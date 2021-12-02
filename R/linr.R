#' Function to Fit Linear Regression
#'
#' "linr" is used to fit a simple linear model. In this function, linear regression can be done by three matrix decomposition methods, which are the QR decomposition, Cholesky decomposition and the singular value decomposition (SVD).
#' The defalt fitting method used is the QR method which is used in "lm" function in "stats" package. All three methods can fit linear model with high efficiency.
#'
#' @param formula a symbolic description of the model to be fitted with specific pattern (i.e. y ~ x1 + x2 where y is the responding variable; x1 and x2 are the covariates for the model)
#' @param data an optional data frame, list or environment containing the variables in the model. If not found in data, the variables are taken from formula.
#' @param method an optional character string specifying the fitting method. It must be one of the strings "qr", "cholesky", or "svd".
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
#' @examples
#' Linear_model.lm <- lm(dist~speed,data=cars)
#' print(Linear_model.lm)
#' Linear_model.linr <- linr(dist~speed,data=cars)
#' print(Linear_model.linr$Call)
#' print(Linear_model.linr$Coefficients)
#'
#'
#' @seealso
#' @export
##


library(Rcpp)
library(RcppArmadillo)

linr <- function(formula, data, method = "cholesky") {

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

  if (is.null(n)) stop("'x' must be a matrix")
  if (n == 0L) stop("0 (non-NA) cases")
  if (p == 0L) { # The Null model
    return(list(coefficients = numeric(), residuals = Y,
                fitted.values = 0 * Y, rank = 0))
  }

  if (p == 1L) { # High Efficient Simple Linear Regression.
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
      warning(gettextf("method = '%s' is not supported. Using 'qr', 'svd' or 'cholesky'", method), domain = NA)
    }

    else if (method == "svd") { # if using SVD decomposition method
      svd_decom = svd(X, LINPACK = TRUE)
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
    variab = c("(Intercept)", paste0("X", 1:p - 1))
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

  fit = list(Call = cl, Coefficients = betahat, fitted.values = fit.val,
              residuals = residuals, MSE = MSE, std.error = SE.betahat,
              R.square = R.square, Adj.R.square = Adj.R.square,
              T_statistic = T.stat, p_value.T = p.val.T,
              F_statistic = F.stat, p_value.F = p.val.F)
  class(fit) = "linr"
  return(fit)
}



#' @title Summarizing Linear Model Fits
#'
#' @description "linr.summary" is a generic function used to produce result summaries of the results from linear regression done by "linr" function.
#'
#' @param object the fitting results of "linr".
#' @param correlation logical; if TRUE, the correlation matrix of the estimated parameters is returned and printed.
#'
#' @return A list containing the following elements
#'         \itemize{
#'            \item call - the fitted linear model fomula.
#'            \item Residuals - summary of the usual residuals with min, 1st quantile, median, 3rd quantile, max being computed.
#'            \item Coefficients - a table with the estimated values for each covariates and the intercept as well as their corresponding standard error, t-statistics and (two-sided) p-value.
#'            \item Residual standard erro - the square root of the estimated variance of the random error and a corresponding degrees of freedom will also be computed.
#'            \item r.squared - R^2, the 'fraction of variance explained by the model', SSR (variation in fitted values about the overall mean) / SSY (total variation in Y about its overall mean).
#'            \item adj.r.squared - adjusted version of R^2, penalized based on number of covariates.
#'            \item cov.unscaled - a table of (unscaled) covariances of the coeficients.
#'            \item correlation - the correlation table corresponding to the above cov.unscaled table, if correlation = TRUE is specified.
#'         }
#'
#' @export


# linr.summary = function(fit) { # test
#
#   summ = fit
#   var.name = names( fit$coefficients )
#   summ$aliased = drop( is.na( coef(fit) ) )
#   summ$coefficients = cbind( fit$coefficients, fit$sd.beta, fit$`t value`, fit$`Pr(>|t|)` )
#   dimnames(summ$coefficients) = list( var.name, c("Estimate", "Std. Error", "t value", "Pr(>|t|)") )
#
#   fdf = length(fit$coefficients) - 1
#   summ$fstatistic = c(fit$fstatistic, fdf, fit$df.residual)
#   names(summ$fstatistic) = c("value", "numdf", "dendf")
#   summ$fpval = pf(fit$fstatistic, fdf, fit$df.residual, lower.tail = F)
#
#   if ( prt ) {
#     cat("\nCall:\n")
#     print( summ$call )
#     cat("\nResiduals:\n")
#     print( round(summary(fit$residuals)[-4], 5) )
#     cat("\nCoefficients:\n")
#     print( round(summ$coefficients, 5) )
#     cat( sprintf("\nResidual standard erro:%.2f on %d degrees of freedom\n", fit$sigma, fit$df.residual) )
#     cat( sprintf("Multiple R-squared:\t%.4f,\tAdjusted R-squared:\t%.4f\n", fit$r.squared, fit$adj.r.squared) )
#     cat( sprintf("F-statistic: %.2f on %d and %d DF, p-value: %e\n", fit$fstatistic, fdf, fit$df.residual, summ$fpval) )
#   }
#
#   if ( correlation ) {
#     summ$correlation = (fit$cov.unscaled * fit$sigma^2) / outer(fit$sd.beta, fit$sd.beta)
#     dimnames( summ$correlation ) = list( var.name, var.name )
#     cor.tri = summ$correlation
#
#     if ( prt ) {
#       cat("\nCorrelation of Coefficients:\n")
#       cor.tri = round(cor.tri, 2)
#       cor.tri[upper.tri(cor.tri, diag = TRUE)] = ""
#       print( cor.tri[-1, -ncol(cor.tri)], quote = FALSE )
#     }
#
#   }
#   summ[["sd.beta"]] = NULL
#   summ[["t value"]] = NULL
#   summ[["Pr(>|t|)"]] = NULL
#   return( summ )
#
# }






