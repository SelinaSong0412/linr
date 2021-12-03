#' @title Summarizing Linear Model Fits
#'
#' @description "linr.summary" is a generic function used to produce result summaries of the results from linear regression done by "linr" function.
#'
#' @param object the fitting results returned from "linr".
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
#   var.name = names(fit$coefficients)
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
