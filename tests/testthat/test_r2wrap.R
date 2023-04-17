# the generated G explains 50% of liabilty
n <- 10000
# This results in Y~N(0,1)
vr <- 1/sqrt(2)
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

# Calculate number of cases in the population that should roughly be prevalence*n
nCas <- nrow(subset(dta, yth==1))
# N in sample will thus be
nSample <- nCas/prevalenceSample
nCont <- abs(nSample - nCas)

# Now we draw a random sample so that we get the sample prevalence defined above
set.seed(123)
dtaSamp <- rbind(subset(dta, yth==1),
                 subset(dta, yth==0)[sample(1:(n-nCas), nCas),])

#ml <- glm(yth~x, data=dtaSamp, family=binomial)
#library(fmsb)
#NagelkerkeR2(m)

# Desctools Nagelkerke: This package does not comply with fmsb
# (1 - exp((D.full - D.base)/n))/(1 - exp(-D.base/n))
# (1 - exp((rr$dev - rr$null)/n))/(1 - exp(-rr$null/n)) #fmsb

summary(lm(y~x, data=dta))
library(DescTools, quietly=T)
library(fmsb, quietly=T)
library(r2redux, quietly=T)
# Compare with DescTools::PseudoR2 and r2redux::cc_trf
test_that('Testing linear model on a y=1/0', {
  m <- glm(yth~x, data=dtaSamp)
  mP <- R2Wrap(m, cases=nCas, nCont)
  expect_equal(class(mP)[1], 'glm-gaussian')
  expect_equal(class(mP@model), c('glm', 'lm'))
  expect_equal(mP@cases, nCas)
  expect_equal(mP@controls, nCont)
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

  #likelihoodR2(mP)
  #nagelkerkeR2(mP)
  #CoxSnellR2(mP)
})

test_that('Testing binomial model on a y=1/0', {
  m <- glm(yth~x, data=dtaSamp, family=binomial)
  mP <- R2Wrap(m)
  expect_equal(class(mP)[1], 'glm-binomial')
  expect_equal(class(mP@model), c('glm', 'lm'))
  expect_equal(mP@cases, nCas)
  expect_equal(mP@controls, nCont)

  # Nagelkerke
  r2 <- NagelkerkeR2(m)$R2
  expect_equal(R2.Nagelkerke(mP), r2)

  # Scale liability
  # We calculate the observed scale R2 manually as this is not done by r2redux::cc_trf to get liability scale R2
  r2 <- var(predict(m, type='response'))/(prevalenceSample*(1-prevalenceSample))
  expect_equal(R2.scaleObserved(mP), r2)
  # This function takes a SD of R2 aswell, but this is not relevant for the calculation of libability scale R2
  r2 <- cc_trf(r2, 0.01, K=prevalence, P=prevalenceSample)$R2l
  expect_equal(R2.scaleLiability(mP, prevalence), r2)
})
