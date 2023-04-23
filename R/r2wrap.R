#' Create a model container from a fitted glm model
#'
#' if passing a formula to argument `comparison`, R2.scaleLiability will provide
#' the incremental R2, i.e. the difference between model `mod` and `mod` with
#' the updated formula.
#'
#' The parameters `cases` and `controls` is relevant when analyzing a binary
#' phenotype with a linear model (`glm(..., family="gaussian")`), in which
#' `p*(1-p)` is used as a denominator for `R2.scaleObserved` instead of
#' `var(y)`. For binomial models these arguments can be ignored.
#'
#' @export
#' @param mod The fitted model object
#' @param comparison A formula to compare the model `mod` with, e.g. ~-exposure to remove `exposure` from `mod`
#' @param cases N of `y=1`
#' @param controls N of `y=0`
R2Wrap <- function (mod, comparison=~0, cases=NA_integer_, controls=NA_integer_) {
  determineFamily(mod, comparison, cases, controls)
}

#' An S4 class to contain a fitted model object
#'
#' @export
#' @slot model The type of model used
#' @slot validModels The accepted type of models
#' @slot cases The number of `y=1`
#' @slot controls The number of `y=0`
#' @slot N The sample size
setClass('R2Wrap', representation(
  model = 'glm',
  family = 'character',
  validModels = 'vector',
  cases = 'numeric',
  controls = 'numeric',
  comparison = 'formula',
  N = 'integer'
  ),
  prototype=list(
    validModels = c('glm-gaussian', 'glm-binomial')
  ))
setMethod('initialize', signature='R2Wrap',
          function (.Object, model, comparison, family, cases, controls) {
            .Object <- callNextMethod()
            .Object@model = model
            .Object@comparison = comparison
            .Object@family = family
            .Object@cases <- cases
            .Object@controls <- controls
            .Object@N = nobs(model)
            .Object
          })

setClass('glm-gaussian', contains='R2Wrap')
setMethod('initialize', signature('glm-gaussian'),
          function (.Object, ...) {
            .Object <- callNextMethod(.Object, ...)
            .Object
          })
setClass('glm-binomial', contains='R2Wrap')
setMethod('initialize', signature='glm-binomial',
          function (.Object, ...) {
            .Object <- callNextMethod(.Object, ...)
            .Object
          })

setGeneric('determineFamily', function (object, comparison, cases, controls) standardGeneric('determineFamily'))
setMethod('determineFamily', signature('glm', 'formula', 'numeric', 'numeric'),
          function (object, comparison, cases, controls) {
            fam <- family(object)$family
            if (fam == 'binomial') {
              cases <- sum(object$y == 1)
              controls <- sum(object$y == 0)
            }
            mod <- class(object)[1]
            cl <- paste0(mod, '-', fam)
            if (!fam %in% c('gaussian', 'binomial')) stop('Model of type ', mod, ' with family ', fam, ' is not supported')
            new(cl, object, comparison, fam, cases, controls)
          })

setGeneric('hasComparison', function (object) standardGeneric('hasComparison'))
setMethod('hasComparison', signature('R2Wrap'),
          function (object) {
            deparse(object@comparison) != '~0'
          })

#' @export
#' @docType methods
#' @rdname R2.scaleLiability-methods
setGeneric('R2.scaleLiability', function (object, prevalence, R2, prevalence.sample) standardGeneric('R2.scaleLiability'))

#' Liability scale R2 of a glm model
#'
#' Calculate liability scale R2 given a population prevalence or a set of prevalences
#' according to Lee et.al (2012).
#'
#' @rdname R2.scaleLiability-methods
#' @aliases R2.scaleLiability
#' @param object A R2Wrap object
#' @param prevalence Population prevalence or a vector of prevalence
#' @param R2 R-squared on the observed scale of the model
#' @param prevalence.sample The sample prevalence of Y=1
#' @references Lee, S.H., Goddard, M.e., Wray, N.R. and Visscher, P.M. (2012), A Better Coefficient of Determination for Genetic Profile Analysis. Genet. Epidemiol., 36: 214-224. https://doi.org/10.1002/gepi.21614
#' @details The liability scale R2 as described in Lee et al is particularly
#' appropriate for ascertained case-control studies in which cases have been sampled
#' from the general population based on a binary phenotype `y`. The motivation comes
#' from genetic epidemiology for studying genetic liability (eg, R2 of a polygenic score) of `y`.
#' Other pseudo R2's may be more biased in such situations as the sample prevalence differs from the population
#' prevalance (eg 1:1 matching). However, to my knowledge it is not known how this study
#' how this R2 performs when the phenotype differs from the phenotype cases were ascertained for.
#' @examples
#' # Fit a linear glm model
#' m <- glm(y~x, data=pop)
#'
#' # Create the R2Wrap object
#' modR2 <- R2Wrap(m, cases=100, controls=200)
#'
#' # If the fit is instead a binomial model
#' m <- glm(y~x, data=pop, family=binomial)
#' # the arguments cases and controls can be ignored
#'
#' modR2 <- R2Wrap(m)
#'
#' # Calculate the liability scale R2 for a population prevalence
#' scaleLiabilityR2(modR2, prevalence=0.1)
#'
#' # ... or for a sequence of prevalences
#' scaleLiabilityR2(modR2, prevalence=seq(0.01, 0.2, 0.01))
#' @return Liability scale R2 given provided prevalence

setMethod('R2.scaleLiability', signature('R2Wrap', 'numeric', 'numeric', 'numeric'),
          function (object, prevalence, R2, prevalence.sample) {
            threshold <- -qnorm(prevalence, 0, 1)
            zDens <- dnorm(threshold)
            liabCase <- zDens/prevalence
            liabControl <- -liabCase * prevalence/(1-prevalence)
            theta <- liabCase * (prevalence.sample - prevalence)/(1 - prevalence) * (liabCase*(prevalence.sample - prevalence)/(1 - prevalence) - threshold)
            bigC <- prevalence*(1-prevalence)/zDens^2 * prevalence*(1 - prevalence)/(prevalence.sample*(1 - prevalence.sample))
            R2liab <- R2 * bigC/(1 + R2 * theta * bigC)
            R2liab
          })
#' @rdname R2.scaleLiability-methods
#' @aliases R2.scaleLiability
setMethod('R2.scaleLiability', signature=c('glm-gaussian', 'numeric', 'missing', 'missing'),
          function (object, prevalence, R2, prevalence.sample) {
            R2 <- R2.scaleObserved(object)
            prevalenceSample <- object@cases/object@N
            R2.scaleLiability(object, prevalence, R2, prevalenceSample)
          })

#' @rdname R2.scaleLiability-methods
#' @aliases R2.scaleLiability
setMethod('R2.scaleLiability', signature=c('glm-gaussian', 'vector', 'missing', 'missing'),
          function (object, prevalence, R2, prevalence.sample) {
            sapply(prevalences, FUN=function (x) R2.scaleLiability(mP, x))
          })

#' @rdname R2.scaleLiability-methods
#' @aliases R2.scaleLiability
setMethod('R2.scaleLiability', signature=c('glm-binomial', 'vector', 'missing', 'missing'),
          function (object, prevalence) {
            sapply(prevalences, FUN=function (x) R2.scaleLiability(mP, x))
          })

#' @rdname R2.scaleLiability-methods
#' @aliases R2.scaleLiability
setMethod('R2.scaleLiability', signature=c('glm-binomial', 'numeric', 'missing', 'missing'),
          function (object, prevalence) {
            R2 <- R2.scaleObserved(object)
            prevalenceSample <- object@cases/object@N
            r2 <- R2.scaleLiability(object, prevalence, R2, prevalenceSample)
            if (hasComparison(object)) {
              o <- object
              o@model <- update(object=o@model, formula=o@comparison)
              o@comparison <- ~0
              r22 <- R2.scaleLiability(o, prevalence)
              r2 <- r2 - r22
              cat ('Main model:', deparse(object@model$formula), '\nComparison:', deparse(o@model$formula), '\n')
            }
            r2
          })

# Variance in predictions/Variance in Y Var(Ŷ)/Var(Y)
#' @export
setGeneric('R2.scaleObserved', function (object) standardGeneric('R2.scaleObserved'))
setMethod('R2.scaleObserved', signature='glm-binomial',
          function (object) {
            preds <- predict(object@model, type='response')
            r2 <- var(preds)/((object@cases/object@N) * (object@controls/object@N))
            r2
          })
setMethod('R2.scaleObserved', signature=c('glm-gaussian'),
          function (object) {
            preds <- predict(object@model)
            if (is.na(object@cases)) {
              den <- var(object@model$y)
            } else {
              den <- ((object@cases/object@N) * (object@controls/object@N))
            }
            var(preds)/den
          })

# Sometimes referred to as Cohen R2 (probably from a textbook)
#' @export
setGeneric('R2.likelihood', function (object) standardGeneric('R2.likelihood'))
setMethod('R2.likelihood', signature = 'R2Wrap',
          function (object) {
            with(object@model, 1 - deviance/null.deviance)
          })

#' @export
setGeneric('R2.Nagelkerke', function (object) standardGeneric('R2.Nagelkerke'))
setMethod('R2.Nagelkerke', signature='R2Wrap',
          function (object) {
            n <- object@N
            with(object@model, (1 - exp((deviance - null.deviance)/n))/(1 - exp(-null.deviance/n)))
          })

#' @export
setGeneric('R2.CoxSnell', function (object) standardGeneric('R2.CoxSnell'))
setMethod('R2.CoxSnell', signature='glm-gaussian',
          function (object) {
            n <- object@N
            resp <- all.vars(terms(object@model))[1]
            frm <- formula(paste0(resp, '~ 1'))
            m0 <- update(object@model, frm)
            l1 <- logLik(object@model)[1]
            l0 <- logLik(m0)[1]
            1 - exp((l0-l1)*(2/object@N))
          })
