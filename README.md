# GLM Extras

Just some convenience functions for running regressions with R's `glm` and
calculating R-squared.

## Installation

If the R-package `devtools` is installed, this package can be installed by issuing:
```
devtools::install_github('https://github.com/deepchocolate/glm-extras')
```

Otherwise, download one of the release .tar.gz-files and issue:
```
install.packages('glm-extras-0.0.0.tar.gz')

```

## Usage: `regCompact`

```
library(glmExtras)
### Simulate data
# Sample size
n <- 10000
# This results in Y~N(0,1)
vr <- 1/sqrt(2)
set.seed(123)
x <- rnorm(n, sd=vr)
er <- rnorm(n, sd=vr)

prevalence <- 0.1
# We set the threshold corresponding to the prevalence above
th <- -qnorm(prevalence)
# Create the outcome with an error
# x explains ~50% of y
dta <- data.frame(y = x + er, x=x)
# Create cases (about 10% will be y=1)
dta$yth <- as.integer(dta$y > th)


phenotypes <- list('yBin'="Binary y", 
                   'yCont'="Continuous y")
binaries <- c('yBin')
covs <- ""

results <- NULL
for (pheno in names(phenotypes)) {
  m <- regCompact(phenotype=pheno, exposure='x', covs='', dta=dta, binary=pheno %in% binaries)
  results <- rbind(results, m)
}
```

## Usage: `R2Wrap`

This functionality is useful for calculating various R-squared statistics.

```
m <- glm(yth~x, data=dta, family=binomial)
mR2 <- R2Wrap(m)
R2.Nagelkerke(mR2)
R2.scaleLiability(mR2, prevalence=0.1)
R2.scaleLiability(mR2, prevalence=seq(0.05, 0.15, 0.01))
```

R2.scaleLiability calculates R-squared according to Lee (2012). Briefly, this method
converts the observed scale R-squared (i.e., variance of probability predictions of y=1),
and converts this back to the liability scale. It depends on a given prevalence of a binary 
phenotype, which is most certainly an uncertain parameter.

To my knowledge, it is appropriate when studying the R2 in models where case status is the
phenotype is the outcome being modeled. Whether it is appropriate when the outcome phenotype 
differs from the phenotype defining a case is another question...

If using a linear (gaussian) model to analyse a binary phenotype it is possible
to pass the number of cases (in the meaning of that `y=1`) and controls to override
the observed variance calculation in `R2.scaleObserved`.
```
mR2 <- R2Wrap(m, cases=#, controls#)
```

# References
Lee, S.H., Goddard, M.e., Wray, N.R. and Visscher, P.M. (2012), A Better Coefficient of Determination for Genetic Profile Analysis. Genet. Epidemiol., 36: 214-224. https://doi.org/10.1002/gepi.21614
