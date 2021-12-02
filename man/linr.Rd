\name{linr}
\alias{linr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
linr(formula, data, method = "cholesky")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, data, method = "cholesky") 
{
    cl = match.call()
    if (hasArg(data)) {
        mf = model.frame(formula, na.omit(data))
    }
    else {
        mf = model.frame(formula)
    }
    Y = mf[, 1]
    X = mf[, -1]
    if (!is.numeric(Y)) {
        stop("The outcomes should be numeric")
    }
    if (is.null(attributes(Y))) {
        attr(Y, "dim") = c(length(Y), 1)
    }
    if (is.null(attributes(X))) {
        attr(X, "dim") = c(length(X), 1)
    }
    if (nrow(Y) != nrow(X)) {
        stop("Number of the outcomes and observations do not match.")
    }
    else if (nrow(X) < ncol(X)) {
        stop("Number of observations is less than predictors.")
    }
    n = nrow(X)
    p = ncol(X)
    if (is.null(n)) 
        stop("'x' must be a matrix")
    if (n == 0L) 
        stop("0 (non-NA) cases")
    if (p == 0L) {
        return(list(coefficients = numeric(), residuals = Y, 
            fitted.values = 0 * Y, rank = 0))
    }
    if (p == 1L) {
        x = as.vector(X) - mean(X)
        y = as.vector(Y) - mean(Y)
        SSX = sum(x * x)
        SSY = sum(y * y)
        SSXY = sum(x * y)
        R.XY = SSXY/sqrt(SSY * SSX)
        betahat1 = SSXY/SSX
        betahat0 = mean(Y) - betahat1 * mean(X)
        betahat = c(betahat0, betahat1)
        fit.val = betahat0 + betahat1 * X
        residuals = Y - fit.val
        SSE = sum(residuals^2)
        MSE = SSE/(n - 2)
        SE.betahat1 = sqrt(MSE/SSX)
        SE.betahat0 = sqrt(MSE * (1/n + mean(X)^2/SSX))
        SE.betahat = c(SE.betahat0, SE.betahat1)
        p = p + 1
    }
    else {
        X = cbind(rep(1, n), X)
        p = p + 1
        if (method != "svd" & method != "qr" & method != "cholesky") {
            warning(gettextf("method = '\%s' is not supported. Using 'qr', 'svd' or 'cholesky'", 
                method), domain = NA)
        }
        else if (method == "svd") {
            svd_decom = svd(X, LINPACK = TRUE)
            tuy = crossprod(svd_decom$u, Y)
            betahat = crossprod(t(svd_decom$v), (tuy/svd_decom$d))
        }
        else if (method == "qr") {
            res_qr <- qr(X, lapack = TRUE)
            betahat = qr.coef(res_qr, Y)
        }
        else {
            XtX = crossprod(X)
            XtY = crossprod(X, Y)
            R = chol(XtX)
            z = forwardsolve(R, XtY, upper.tri = TRUE, transpose = TRUE)
            betahat = backsolve(R, z)
        }
        fit.val = crossprod(t(X), betahat)
        residuals = Y - fit.val
        SSE = sum(residuals^2)
        SSY = sum((Y - mean(Y))^2)
        MSE = SSE/(n - p)
        varmat = solve(crossprod(X), diag(MSE, p, p), tol = 0, 
            transpose = FALSE)
        SE.betahat = sqrt(diag(varmat))
    }
    R.square = 1 - SSE/SSY
    Adj.R.square = 1 - (n - 1) * MSE/SSY
    T.stat = betahat/SE.betahat
    p.val.T = pt(abs(T.stat), n - p, lower.tail = FALSE) * 2
    F.stat = (SSY - SSE)/(p - 1)/MSE
    p.val.F = pf(F.stat, (p - 1), (n - p), lower.tail = FALSE)
    observ = rownames(mf)
    variab = colnames(mf)[-1]
    if (is.null(variab)) {
        variab = c("(Intercept)", paste0("X", 1:p - 1))
    }
    else {
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
        R.square = R.square, Adj.R.square = Adj.R.square, T_statistic = T.stat, 
        p_value.T = p.val.T, F_statistic = F.stat, p_value.F = p.val.F)
    class(fit) = "linr"
    return(fit)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
