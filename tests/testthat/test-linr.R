

test_that("Whether 'linr' function have the same output values with the 'lm' function", {
  # Test model 1
  y = rnorm(300)
  x = rnorm(300)
  model.linr <- linr(y ~ x)
  model.lm <- lm(y ~ x)
  sum = summary(model.lm)
  expect_equal(model.linr$coefficients, model.lm$coefficients)
  expect_equal(model.linr$fitted.values, model.lm$fitted.values)
  expect_equal(model.linr$residuals, model.lm$residuals)
  expect_equal(model.linr$F_statistic, sum$fstatistic[[1]])
  expect_equal(model.linr$R.square, sum$r.squared)
  expect_equal(model.linr$Adj.R.square, sum$adj.r.squared)
  expect_equal(model.linr$std.error, sum$coefficients[,2])
  expect_equal(model.linr$T_statistic, sum$coefficients[,3])
  expect_equal(model.linr$p_value.T, sum$coefficients[,4])
})

test_that("Whether 'linr' function have the same output values with the 'lm' function", {
  # Test model 2
  Y = rnorm(300)
  X = matrix(rnorm(6000), nrow = 300, ncol = 20)
  model.linr.mul <- linr(Y ~ X, method = "cholesky")
  model.lm.mul <- lm(Y ~ X)
  sum = summary(model.lm.mul)
  expect_equal(model.linr.mul$coefficients, model.lm.mul$coefficients)
  expect_equal(model.linr.mul$fitted.values, model.lm.mul$fitted.values)
  expect_equal(model.linr.mul$residuals, model.lm.mul$residuals)
  expect_equal(model.linr.mul$F_statistic, sum$fstatistic[[1]])
  expect_equal(model.linr.mul$R.square, sum$r.squared)
  expect_equal(model.linr.mul$Adj.R.square, sum$adj.r.squared)
  expect_equal(model.linr.mul$std.error, sum$coefficients[,2])
  expect_equal(model.linr.mul$T_statistic, sum$coefficients[,3])
  expect_equal(model.linr.mul$p_value.T, sum$coefficients[,4])
  })


test_that("Whether 'linr' function have the same output values with the 'lm' function", {
  # Test model 3
  model.linr.cars <- linr(dist ~ speed, data = cars)
  model.lm.cars <- lm(dist ~ speed, data = cars)
  expect_equal(model.linr.cars$coefficients, model.lm.cars$coefficients)
  expect_equal(model.linr.cars$fitted.values, model.lm.cars$fitted.values)
  expect_equal(model.linr.cars$residuals, model.lm.cars$residuals)
  })

test_that("Whether 'linr' function have the same output values with the 'lm' function", {
  # Test model 4
  model.linr.mtcars <- linr(disp ~ mpg + wt + carb, data = mtcars, method = 'qr')
  model.lm.mtcars <- lm(disp ~ mpg + wt + carb, data = mtcars)
  expect_equal(model.linr.mtcars$coefficients, model.lm.mtcars$coefficients)
  expect_equal(model.linr.mtcars$fitted.values, model.lm.mtcars$fitted.values)
  expect_equal(model.linr.mtcars$residuals, model.lm.mtcars$residuals)
  })

test_that("Whether 'linr' function have the same output values with the 'lm' function", {
  # Test model 5
  model.linr.mtcars <- linr(disp ~ mpg + wt + carb, data = mtcars, method = 'svd')
  model.lm.mtcars <- lm(disp ~ mpg + wt + carb, data = mtcars)
  expect_equal(model.linr.mtcars$coefficients, model.lm.mtcars$coefficients)
  expect_equal(model.linr.mtcars$fitted.values, model.lm.mtcars$fitted.values)
  expect_equal(model.linr.mtcars$residuals, model.lm.mtcars$residuals)
})

test_that("Output of Error", {
  # Test error outputs:
  Y.nonnumeric = rep("1", 300)
  x = rnorm(300)
  expect_error(linr(Y.nonnumeric~x), "The outcomes should be numeric")
})

# test_that("Output of Error", {
#   # Test error outputs:
#   Y.nonnumeric = rep("1", 300)
#   x = rnorm(300)
#   expect_error(linr(Y.nonnumeric~x), "The outcomes should be numeric")
# })
#
# test_that("Output of Error", {
#   # Test error outputs:
#   y = rnorm(300)
#   x = rnorm(301)
#   expect_error(linr(y~x), "Number of the outcomes and observations do not match.")
# })
#
# test_that("Output of Error", {
#   # Test error outputs:
#   Y = rnorm(300)
#   X = matrix(rnorm(6000), nrow = 30, ncol = 200)
#   expect_error(linr(Y~X), "Number of observations is less than predictors.")
# })
#
# test_that("Output of Error is.null(n)", {
#   # Test error outputs:
#   Y = rnorm(300)
#   X = c()
#   expect_error(linr(Y~X), "Number of observations is less than predictors.")
# })






