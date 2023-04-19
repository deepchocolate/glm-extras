setGeneric('extract', function (model, coefficient, ci, ci.clustered) standardGeneric('extract'))
setMethod('extract', signature('matrix', 'character', 'missing', 'missing'),
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
setMethod('extract', signature('glm', 'character', 'logical', 'missing'),
          function (model, coefficient, ci) {
            stats <- coef(summary(model))
            out <- extract(stats, coefficient)
            if (ci) {
              ci <- confint.default(model)[coefficient,]
              out$l95 <- ci[1]
              out$u95 <- ci[2]
            }
            out
          })

