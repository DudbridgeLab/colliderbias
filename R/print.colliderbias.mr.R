# Print results of collider bias adjustment
#
# After regressing effects on the outcome on the effects on the exposure and the selection variable,
# print out the estimated causal effect and its confidence interval
#
#' @export
print.colliderbias.mr = function(x, ...) {
  print(paste("Coefficient", x$b))
  print(paste("Standard error", x$b.se))
  print(paste("95% CI", x$b.ci[1], x$b.ci[2]))
}
