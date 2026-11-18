#' Adjust association statistics for collider bias
#'
#' The \code{colliderbias} function adjusts association statistics for a trait Y for collider bias
#' through a selection variable X.
#'
#' The model assumes that predictors have direct effects on Y, which are the effects of interest,
#' as well as indirect effects via X, which may be confounded with Y.
#' The observed \code{ybeta} are conditional on X, creating collider bias on the direct effects.
#' 
#' Effect sizes are on a linear scale, so they could be linear regression coefficients, or log odds ratios, or log hazard ratios.
#' Effects on Y are regressed on the effects on X.
#' The residuals are then used to obtain adjusted effect sizes and standard errors for Y.
#' This approach is called instrument effect regression.
#'
#' The regression should be performed on a subset of predictors that are independent.
#' In a genome-wide association study, these would be LD-pruned SNPs.
#' In terms of the input parameters, the instrument effect regression command is
#' \code{lm(ybeta[prune] ~ xbeta[prune], weights = yse[prune]^(-2))}.
#'
#' The effects in \code{xbeta} and \code{ybeta} should be aligned for the same variables
#' and the same direction prior to calling \code{colliderbias}.
#'

#' @param xbeta Vector of effects on the selection variable.
#' @param xse Vector of standard errors of \code{xbeta}.
#' @param ybeta Vector of effects on the trait of interest.
#' @param yse Vector of standard errors of \code{ybeta}.
#' @param prune Vector containing indices of an approximately independent subset of the predictors in \code{xbeta} and \code{ybeta}. If unspecified, all predictors will be used.
#' @param method Method for correcting weak instrument bias.
#' "CWLS"(default) applies the Corrected Weighted Least Squares estimator (Cai et al 2022).
#' "WLS" applies weighted least squares without correction for weak instruments.
#' "mr.raps" applies the MR-RAPS method (Zhao et al 2020).
#' @param intercept if \code{TRUE} (default), includes an intercept in the instrument effect regression.
#' @param weighted Weights used in the instrument effect regression.
#' TRUE (default) uses inverse variance weights.
#' FALSE uses equal weights, giving ordinary least squares.
#' @param over.dispersion parameter passed to mr.raps.
#'
#'
#' @import stats
#'
#' @return An object of class "colliderbias" which is a list containing:
#'   \itemize{
#'     \item{\code{ybeta.adj} Adjusted effects on the trait of interest}
#'     \item{\code{yse.adj} Adjusted standard errors of \code{ybeta.adj}}
#'     \item{\code{ychisq.adj} Chi-square statistics for \code{ybeta.adj}}
#'     \item{\code{yp.adj} P-values for \code{ychisq.adj} on 1df}
#'     \item{\code{b} Estimated slope of the instrument effect regression}
#'     \item{\code{b.se} Standard error of \code{b}}
#'     \item{\code{b.ci} Lower and upper 95\% confidence limits for \code{b}}
#'     \item{\code{b.raw} Estimated slope, without correction for weak instruments}
#'     \item{\code{b.raw.se} Standard error of \code{b.raw}}
#'     \item{\code{b.raw.ci} Lower and upper 95\% confidence limits for \code{b.raw}}
#'   }
#'
#' @examples
#' # Load the test dataset
#' data(testData)
#' attach(testData)
#' 
#' # Adjust for collider bias using instrument effect regression,
#' # correcting weak instrument bias using CWLS.
#' result <- colliderbias(dbeta, dse, ybeta, yse, method = "CWLS")
#' 
#' # slope without correcting weak instruments is -0.07050519
#' result$b.raw
#'
#' # slope with CLWS correction is -0.1825845
#' result$b
#'
#' @author Siyang Cai, Frank Dudbridge
#'
#' @references
#' Dudbridge F, Allen RJ, Sheehan NA et al. (2019) Adjustment for index event bias in genome-wide association studies of subsequent events. Nat Commun 10:1561. \href{https://doi.org/10.1038/s41467-019-09381-w}{https://doi.org/10.1038/s41467-019-09381-w}
#'
#' Cai S, Hartley A, Mahmoud O et al. (2022) Adjusting for collider bias in genetic association studies using instrumental variable methods. Genet Epidemiol 46:303–16. \href{https://doi.org/10.1002/gepi.22455}{https://doi.org/10.1002/gepi.22455}
#'
#' Zhao Q, Wang J, Hemani G et al. (2020) Statistical inference in two-sample summary-data Mendelian randomization using robust adjusted profile score. The Annals of Statistics 48:1742-1769. \href{https://doi.org/10.1214/19-AOS1866}{https://doi.org/10.1214/19-AOS1866}
#'
#' @export
colliderbias = function(xbeta,
                        xse,
                        ybeta,
                        yse,
                        prune = NULL,
                        method = "CWLS",
                        intercept = T,
                        weighted = T,
                        over.dispersion = FALSE) {
  
  # check input lengths
  if (length(xbeta) != length(xse))
    stop("xbeta and xse have different lengths")
  if (length(ybeta) != length(yse))
    stop("ybeta and yse have different lengths")
  if (length(xbeta) != length(ybeta))
    stop("xbeta and ybeta have different lengths")
  
  # check input method
  if (!(tolower(method) %in% c("cwls", "wls", "mr.raps")))
    stop(paste("method", method, "not implemented"))
  
  # use all SNPs if no pruned set given
    if (is.null(prune)) {
    prune = 1:length(xbeta)
  }

  # select pruned SNPs
  xbeta.prune = xbeta[prune]
  xse.prune = xse[prune]
  ybeta.prune = ybeta[prune]
  yse.prune = yse[prune]

  # regression weights
  if (weighted){
    # inverse variance weights
    weight = 1 / yse.prune ^ 2
  }
  else {
    # unweighted
    weight = rep(1, length(prune))
  }    

  # instrument effect regression
  if (intercept) {
    fit = lm(ybeta.prune ~ xbeta.prune, weights = weight)
    fit.summary = summary(fit)
    b.raw = as.numeric(fit.summary$coef[2,1])
    b.raw.se = as.numeric(fit.summary$coef[2,2])
  }
  else {
    fit = lm(ybeta.prune ~ xbeta.prune - 1 , weights = weight)
    fit.summary = summary(fit)
    b.raw = as.numeric(fit.summary$coef[1])
    b.raw.se = as.numeric(fit.summary$coef[2])
  }
  b = b.raw
  b.se = b.raw.se

  # Corrected weighted least squares adjustment for weak instruments
  if (tolower(method) == "cwls") {
    sumweight = sum(weight, na.rm=TRUE)
    sumxbeta = 0
    if (intercept) sumxbeta = sum(weight * xbeta.prune, na.rm=TRUE)
    cwls = (sumweight * sum(weight * xbeta.prune^2, na.rm = TRUE) - sumxbeta^2) /
      (sumweight * sum(weight * xbeta.prune^2, na.rm = TRUE) - sumxbeta^2 - sumweight * sum(weight * xse.prune^2, na.rm = TRUE))
    
    b = b.raw * cwls
    b.se = b.raw.se * abs(cwls)
  }
  
  # MR-RAPS adjustment for weak instruments
  if(tolower(method) == "mr.raps"){
    MRRAPS.results = mr.raps::mr.raps(xbeta.prune, ybeta.prune, xse.prune, yse.prune, 
                                      over.dispersion = over.dispersion)
    b = MRRAPS.results$beta.hat
    b.se = MRRAPS.results$beta.se
  }

  # confidence intervals
  b.ci = c(b - qnorm(0.975) * b.se, b + qnorm(0.975) * b.se)
  b.raw.ci = c(b.raw - qnorm(0.975) * b.raw.se, b.raw + qnorm(0.975) * b.raw.se)
  
  # adjusted summary statistics
  ybeta.adj = rep(0, length(ybeta))
  ybeta.adj = ybeta - b * xbeta
  yse.adj = sqrt(yse ^ 2 + b ^ 2 * xse ^ 2)
  ychisq.adj = (ybeta.adj / yse.adj) ^ 2
  yp.adj = stats::pchisq(ychisq.adj, 1, lower.tail = FALSE)
  
  results = list(
    ybeta.adj = ybeta.adj,
    yse.adj = yse.adj,
    ychisq.adj = ychisq.adj,
    yp.adj = yp.adj,
    b = b,
    b.se = b.se,
    b.ci = b.ci,
    b.raw = b.raw,
    b.raw.se = b.raw.se,
    b.raw.ci = b.raw.ci
  )
  class(results) = "colliderbias"
  results
}
