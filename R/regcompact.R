#' Run a series of regressions.
#'
#' Run a regression phenotype~exposure + covariate 1 + covariate 2 and extract
#' parameters in a tabular format.
#' For binary=F, both phenotype and exposure will be standardized to mean 0 and
#' unit variance. For binary=T, exposure will be standardized.
#'
#' `regCompact` calculates incremental R-squared for `exposure` by default.
#'
#' For clustered data (eg samples from families) the clustering argument can be
#' used to calculate cluster robust confidence intervals (see sandwhich::vcovCL)
#'
#' @param phenotype Variable to analyze.
#' @param exposure Exposure variable
#' @param covariates Variables to adjust the exposure for, eg "COV1 + COV2"
#' @param binary Pass TRUE to run a binomial/logistic model.
#' @param dta A data.frame containing all necessary variables
#' @param r2s A vector of R2 measures to calculate in addition to the default.
#' @param clustering Pass the name of the clustering variable in dta to add cluster robust confidence intervals.
#' @returns `regCompact` returns a data.frame with the following statistics.
#'
#' | `phenotype` | Phenotype variable name |
#' | --- | --- |
#' `exposure` | Exposure variable name
#' `beta` | Beta coefficient
#' `stderr` | Standard error of `beta`
#' `z` | Z-statistic of `beta`
#' `p` | P-value of `beta`
#' `l95` | Lower bound of 95% confidence interval
#' `u95` | Upper bound of 95% confidence interval
#' `l95r` | Robust version of `l95` if clustering is used
#' `u95r` | Robust version of `u95` if clustering is used
#' `r2` | R-squared, for a continuous exposure and outcome, this is the square of `beta`
#' |  | For a binary phenotype, this is the observed scale R-squared
#' `family` | Regression family, either gaussian or binomial
#' `n` | Number of observations
#' `pheno_na` | Missingness in phenotype
#' `exp_na` | Missingness in exposure
#'
regCompact <- function (phenotype, exposure, covariates, binary, dta, clustering=F) {
  frm <- paste0(phenotype, '~', exposure)
  if (covariates != '') frm <- paste0(frm, '+', covariates)
  frm <-  as.formula(frm)
  fam <- ifelse(binary == T, 'binomial', 'gaussian')
  dta[, exposure] <- stdVar(dta[, exposure])
  if (binary == F) dta[, phenotype] <- stdVar(dta[, phenotype])
  m <- glm(frm, data=dta, family=fam)
  cfs <- coef(summary(m))[exposure, ]
  # Calculate incremental R2
  if (binary == F) {
    R2 <- cfs[1]^2
    R2N <- NA
  } else {
    R2.m1 <- R2Wrap(m)
    m2 <- update(m, as.formula(paste0('. ~ . -', exposure)))
    R2.m0 <- R2Wrap(m2)
    R2 <- R2.scaleObserved(R2.m1) - R2.scaleObserved(R2.m0)
    R2N <- R2.Nagelkerke(R2.m1) - R2.Nagelkerke(R2.m0)
  }
  confin <- confint.default(m)[exposure, ] # Use asymptotic std errs
  out <- data.frame(
    phenotype=phenotype,
    extract(m, exposure, F),
    l95r='', u95r='', # Cluster-robust CI if applicable
    r2=R2,
    r2Nagelkerke=R2N,
    family=fam,
    n=nrow(dta),
    pheno_na=sum(is.na(dta[, phenotype])), # NAs in phenotype
    exp_na=sum(is.na(dta[, exposure])), # NAs in exposure
    stringsAsFactors = F)
  if (is.character(clustering)) {
    require(sandwhich)
    vcm <- vcovCL(m, cluster=dta[, clustering])
    robErr <- qt(0.975, df=m$df.residual)*sqrt(vcm[exposure, exposure])
    out$u95r <- coef(m)[exposure] + robErr
    out$l95r <- coef(m)[exposure] - robErr
  }
  rownames(out) <- NULL
  return(out)
}
