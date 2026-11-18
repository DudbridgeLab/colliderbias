#' Adjust for collider bias in Mendelian randomisation
#'
#' Estimates the causal effect of an exposure X on an outcome Y, where Y is conditional on a selection variable D.
#' 
#' The \code{colliderbias.mr} function performs two-sample Mendelian randomisation, adjusting the instrumental variable associations
#' with Y for collider bias through D.
#'
#' Effect sizes are on a linear scale, so they could be linear regression coefficients, or log odds ratios, or log hazard ratios.
#' Effects on Y are jointly regressed on the effects on X and D, with the coefficient of X estimating the causal effect.
#' This approach is called generalised instrument effect regression.
#' Correction for weak instruments is implemented by the corrected weighted bivariate least squares formula.
#'
#' The regression should be performed on a subset of predictors that are independent.
#' In a genome-wide association study, these would be LD-pruned SNPs.
#' In terms of the input parameters, the generalised instrument effect regression command is
#' \code{lm(ybeta[prune] ~ xbeta[prune] + dbeta[prune], weights = yse[prune]^(-2))}.
#'
#' The effects in \code{xbeta}, \code{ybeta} and \code{dbeta} should be aligned for the same variables
#' and the same direction prior to calling \code{colliderbias.mr}.
#'

#' @param xbeta Vector of effects on the exposure variable.
#' @param xse Vector of standard errors of \code{xbeta}.
#' @param dbeta Vector of effects on the selection variable.
#' @param dse Vector of standard errors of \code{dbeta}.
#' @param ybeta Vector of effects on the outcome variable.
#' @param yse Vector of standard errors of \code{ybeta}.
#' @param prune Vector containing indices of an approximately independent subset of the predictors in \code{xbeta},
#' \code{ybeta} and \code{dbeta}. If unspecified, all predictors will be used.
#' @param method Method for correcting weak instrument bias.
#' "CWBLS"(default) applies the Corrected Weighted Bivariate Least Squares estimator (Cai & Dudbridge 2025).
#' "WBLS" applies weighted bivariate least squares without correction for weak instruments.
#' @param intercept if \code{TRUE} (default), includes an intercept in the generalised instrument effect regression.
#' @param weighted Weights used in the generalised instrument effect regression.
#' TRUE (default) uses inverse variance weights.
#' FALSE uses equal weights, giving ordinary least squares.
#' @param nboot Number of bootstrap replicates (default 100) to estimate standard error of CWBLS estimator.
#'
#' @import stats
#'
#'
#' @return An object of class "colliderbias.mr" which is a list containing:
#'   \itemize{
#'     \item{\code{b} Estimated causal effect of exposure on outcome}
#'     \item{\code{b.se} Standard error of \code{b}}
#'     \item{\code{b.ci} Lower and upper 95\% confidence limits for \code{b}}
#'     \item{\code{b.raw} Estimated causal effect, without correction for weak instruments}
#'     \item{\code{b.raw.se} Standard error of \code{b.raw}}
#'     \item{\code{b.raw.ci} Lower and upper 95\% confidence limits for \code{b.raw}}
#'   }
#'
#' @examples
#'
#' # Load the test dataset
#' data(testData)
#' attach(testData)
#'
#' # Estimate the causal effect of exposure on outcome
#' result <- colliderbias.MR(xbeta, xse, ybeta, yse, dbeta, dse, method = "CWBLS")
#' 
#' # True causal effect is 0.4
#' result
#'

#' @author Siyang Cai, Frank Dudbridge
#'
#' @references
#' Cai S, Dudbridge F (2025) Estimating causal effects on a disease progression trait Using bivariate Mendelian randomisation.
#' Genet Epidemiol 49:e22600. \href{https://doi.org/10.1002/gepi.22600}{https://doi.org/10.1002/gepi.22600}
#'
#' @export
colliderbias.mr = function(xbeta,
                           xse,
                           dbeta,
                           dse,
                           ybeta,
                           yse,
                           prune = NULL,
                           method = "CWBLS",
                           intercept = T,
                           weighted = T,
                           nboot = 100) {

  # check input lengths
  if (length(xbeta) != length(xse))
    stop("xbeta and xse have different lengths")
  if (length(ybeta) != length(yse))
    stop("ybeta and yse have different lengths")
  if (length(dbeta) != length(dse))
    stop("ybeta and yse have different lengths")
  if (length(xbeta) != length(ybeta))
    stop("xbeta and ybeta have different lengths")
  if (length(xbeta) != length(dbeta))
    stop("xbeta and dbeta have different lengths")
  
  # check input method
  if (!(tolower(method) %in% c("cwbls", "wbls")))
    stop(paste("method", method, "not implemented"))
  
  # use all SNPs if no pruned set given
  if (is.null(prune)) {
    prune = 1:length(xbeta)
  }
  
  # bootstrapped CWBLS function
  cwbls = function(data, i){
    data = data[i,]
    n_i = sum(weight * data$dse ^ 2) * sum(weight * data$xbeta * data$ybeta)
    d_i = sum(weight * data$xbeta^2) * sum(weight * data$dse^2) + 
      sum(weight * data$xse^2) * sum(weight * data$dbeta^2) - sum(weight * data$xse^2) * sum(weight * data$dse^2)
    b = (sum(weight * data$dbeta^2) * sum(weight * data$xbeta * data$ybeta) - sum(weight * data$dbeta * data$xbeta) * sum(weight * data$dbeta * data$ybeta) - n_i)/
      (sum(weight * data$xbeta^2) * sum(weight * data$dbeta^2) - (sum(weight * data$xbeta * data$dbeta)) ^ 2 - d_i)
    return(b)
  }
  
  # select pruned SNPs
  xbeta.prune = xbeta[prune]
  xse.prune = xse[prune]
  ybeta.prune = ybeta[prune]
  yse.prune = yse[prune]
  dbeta.prune = dbeta[prune]
  dse.prune = dse[prune]
  
  #regression weights
  if (weighted) {
    # inverse variance weights
    weight = 1 / yse.prune^2
  }
  else {
    # unweighted
    weight = rep(1, length(prune))
  }

    # uncorrected for weak instruments
  if (intercept) {
    fit = lm(ybeta.prune ~ xbeta.prune + dbeta.prune, weights = weight)
    fit.summary = summary(fit)
    b.raw = as.numeric(fit.summary$coef[2,1])
    b.raw.se = as.numeric(fit.summary$coef[2,2])
  }
  else {
    fit = lm(ybeta.prune ~ xbeta.prune + dbeta.prune - 1, weights = weight)
    fit.summary = summary(fit)
    b.raw = as.numeric(fit.summary$coef[1,1])
    b.raw.se = as.numeric(fit.summary$coef[1,2])
  }
  b.raw.ci = c(b.raw - qnorm(0.975) * b.raw.se, b.raw + qnorm(0.975) * b.raw.se)
  
  # CWBLS with bootstrap SE
  if (method == "CWBLS") {
    data = data.frame(xbeta=xbeta.prune, xse=xse.prune,
                      ybeta=ybeta.prune, yse=yse.prune,
                      dbeta=dbeta.prune, dse=dse.prune)
    boot_cwbls = boot::boot(data = data, statistic = cwbls, R = nboot, stype = "i")
    b = boot_cwbls$t0
    b.se = sd(boot_cwbls$t)
  }
  # no correction for weak instruments
  else {
    b = b.raw
    b.se = b.raw.se
  }
  b.ci = c(b - qnorm(0.975) * b.se, b + qnorm(0.975) * b.se)
  
  results = list(
    b = b,
    b.se = b.se,
    b.ci = b.ci,
    b.raw = b.raw,
    b.raw.se = b.raw.se,
    b.raw.ci = b.raw.ci
  )
  class(results) = "colliderbias.mr"
  results
}



