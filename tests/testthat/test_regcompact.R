# the generated G explains 50% of liabilty
n <- 10000
# This results in Y~N(0,1)
vr <- 1/sqrt(2)
set.seed(123)
x <- rnorm(n, sd=vr)
er <- rnorm(n, sd=vr)

prevalence <- 0.1
prevalenceSample <- 0.5
# We set the threshold corresponding to the prevalence above
th <- -qnorm(prevalence)
# Create the outcome with an error
dta <- data.frame(y = x + er, x=x)
# Create cases
dta$yth <- as.integer(dta$y > th)
library(fmsb, quietly=T)
test_that('Testing regCompact', {
  dta$x <- stdVar(dta$x)
  m <- glm(yth~x, data=dta, family=binomial)
  m0 <- glm(yth~1, data=dta, family=binomial)
  r2 <- NagelkerkeR2(m)$R2 - NagelkerkeR2(m0)$R2
  mc <- regCompact('yth', 'x', '', T, dta)
  expect_equal(mc$beta, coef(m)[['x']])
  expect_equal(mc$r2Nagelkerke, r2)
  print(mc$R2)
})
