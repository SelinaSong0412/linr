y = rnorm(300)
x = rnorm(300)
model.linr <- linr(y ~ x)
model.lm <- lm(y ~ x)

Y = rnorm(300)
X = matrix(rnorm(6000), nrow = 300, ncol = 20)
model.linr.mul <- linr(Y ~ X, method = "cholesky")
model.lm.mul <- linr(Y ~ X)

model.linr.cars <- linr(dist ~ speed, data = cars)
model.lm.cars <- lm(dist ~ speed, data = cars)

model.linr.mtcars <- linr(disp ~ mpg + wt + carb, data = mtcars)
model.lm.mtcars <- lm(disp ~ mpg + wt + carb, data = mtcars)


test_that("usage", {
  expect_equal(model.linr$coefficients, model.lm$coefficients)
  expect_equal(model.linr$fitted.values, model.lm$fitted.values)
  expect_equal(model.linr$residuals, model.lm$residuals)

  expect_equal(model.linr.mul$coefficients, model.lm.mul$coefficients)
  expect_equal(model.linr.mul$fitted.values, model.lm.mul$fitted.values)
  expect_equal(model.linr.mul$residuals, model.lm.mul$residuals)

  expect_equal(model.linr.cars$coefficients, model.lm.cars$coefficients)
  expect_equal(model.linr.cars$fitted.values, model.lm.cars$fitted.values)
  expect_equal(model.linr.cars$residuals, model.lm.cars$residuals)

  expect_equal(model.linr.mtcars$coefficients, model.lm.mtcars$coefficients)
  expect_equal(model.linr.mtcars$fitted.values, model.lm.mtcars$fitted.values)
  expect_equal(model.linr.mtcars$residuals, model.lm.mtcars$residuals)
})



