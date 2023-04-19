#' Extract
#'
#' @param model A fitted glm model
#' @param coefficient A coefficient to extract statistics for
#' @param ci Whether to extract 95% confidence intervals
#' @export
setGeneric('extract.est', function (model, coefficient, ci, ci.clustered) standardGeneric('extract.est'))
setMethod('extract.est', signature('matrix', 'character', 'missing', 'missing'),
          function (model, coefficient, ci) {
            stats <- model[coefficient,]
            out <- data.frame(
              beta=stats[1],
              stderr=stats[2],
              z=stats[3],
              p=stats[4]
            )
            rownames(out) <- NULL
            out
          })
setMethod('extract.est', signature('glm', 'character', 'logical', 'missing'),
          function (model, coefficient, ci) {
            stats <- coef(summary(model))
            out <- extract.est(stats, coefficient)
            if (ci) {
              ci <- confint.default(model)[coefficient,]
              out$l95 <- ci[1]
              out$u95 <- ci[2]
            }
            out
          })

