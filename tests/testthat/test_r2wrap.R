source('helper.R')
# "prevalence" is set in helper.R
# Sequences are tested by creating an incremented vector
prevalenceVector <- c(prevalence, prevalence + seq(0.01,0.03, 0.01))

# Compare with DescTools::PseudoR2 and r2redux::cc_trf
test_that('Testing linear model on a y=1/0', {
  m <- glm(yth~x, data=dtaSamp)
  mP <- R2Wrap(m, cases=nCas, controls=nCont)
  expect_equal(class(mP)[1], 'glm-gaussian')
  expect_equal(class(mP@model), c('glm', 'lm'))
  expect_equal(mP@cases, nCas)
  expect_equal(mP@controls, nCont)
  expect_false(hasComparison(mP))
  # Nagelkerke
  r2 <- NagelkerkeR2(m)$R2
  expect_equal(R2.Nagelkerke(mP), r2)

  # Scale liability
  # We calculate the observed scale R2 manually as this is not done by r2redux::cc_trf to get liability scale R2
  r2 <- var(predict(m))/(prevalenceSample*(1-prevalenceSample))
  expect_equal(R2.scaleObserved(mP), r2)
  # This function takes a SD of R2 aswell, but this is not relevant for the calculation of libability scale R2
  r2 <- cc_trf(r2, 0.01, K=prevalence, P=prevalenceSample)$R2l
  expect_equal(R2.scaleLiability(mP, prevalence), r2)

  # Extra test of scaleObserved on the full population
  # This should give the simulated R2 of about 0.5
  m <- glm(y~x, data=dta)
  mP <- R2Wrap(m)
  r2 <- round(R2.scaleObserved(mP)*100)
  expect_equal(r2, 50, tolerance=1)

  # Test error with missing values
  dtaSampBad <- dtaSamp
  dtaSampBad$x[1] <- NA
  # Default with glm is to na.omit, ie rows with NAs will be ignored
  m <- glm(yth~x, data=dtaSampBad, na.action=na.exclude)
  mP <- R2Wrap(m)
  expect_error(R2.scaleObserved(mP))
})

test_that('Testing binomial model on a y=1/0', {
  m <- glm(yth~x, data=dtaSamp, family=binomial)
  mP <- R2Wrap(m)
  expect_equal(class(mP)[1], 'glm-binomial')
  expect_equal(class(mP@model), c('glm', 'lm'))
  expect_equal(mP@cases, nCas)
  expect_equal(mP@controls, nCont)

  ### Nagelkerke
  r2 <- NagelkerkeR2(m)$R2
  expect_equal(R2.Nagelkerke(mP), r2)

  ### Scale liability
  # We calculate the observed scale R2 manually as this is not done by r2redux::cc_trf to get liability scale R2
  r2 <- var(predict(m, type='response'))/(prevalenceSample*(1-prevalenceSample))
  expect_equal(R2.scaleObserved(mP), r2)
  # This function takes a SD of R2 aswell, but this is not relevant for the calculation of libability scale point estimate R2
  r2 <- cc_trf(r2, 0.01, K=prevalence, P=prevalenceSample)$R2l
  r2Liab <- R2.scaleLiability(mP, prevalence)
  expect_equal(r2Liab, r2)
  # Test using a vector
  r2LiabVec <- R2.scaleLiability(mP, prevalenceVector)
  expect_length(r2LiabVec, length(prevalenceVector))
  expect_equal(r2LiabVec[1], r2)
  # This is just testing using a different coefficient (the error term)
  m <- glm(yth~x + z, data=dtaSamp, family=binomial)
  mP <- R2Wrap(m, comparison=~.-z)
  expect_equal(round(10*R2.scaleLiability(mP, prevalence)), 0)

  # Test what happens with missingness
  dtaSampBad <- dtaSamp
  dtaSampBad$x[1] <- NA
  # Default with glm is to na.omit, ie rows with NAs will be ignored
  m <- glm(yth~x, data=dtaSampBad, family=binomial, na.action = na.exclude)
  mP <- R2Wrap(m)
  expect_error(R2.scaleLiability(mP, prevalence))
})

test_that('Testing errors', {
  m <- glm(yth~x, data=dtaSamp, family=poisson)
  expect_error(R2Wrap(m))
})
