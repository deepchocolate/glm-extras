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
# N controls
nCont <- abs(nSample - nCas)

# Now we draw a random sample so that we get the sample prevalence defined above
set.seed(123)
dtaSamp <- rbind(subset(dta, yth==1),
                 subset(dta, yth==0)[sample(1:(n-nCas), nCas),])
# Add an error term
dtaSamp$z <- rnorm(nrow(dtaSamp))
# Desctools Nagelkerke: This package does not comply with fmsb
# (1 - exp((D.full - D.base)/n))/(1 - exp(-D.base/n))
# (1 - exp((rr$dev - rr$null)/n))/(1 - exp(-rr$null/n)) #fmsb

#library(DescTools, quietly=T) # Install fail
library(fmsb, quietly=T)
library(r2redux, quietly=T)
