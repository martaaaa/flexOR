#' Akaike's Information Criterion corrected
#'
#' This function calculates the corrected Akaike Information Criterion (AICc)
#' for a given statistical model. AICc is a modification of the AIC that
#' provides a correction for small sample sizes, which can lead to overfitting
#' in model selection. AICc = AIC + 2*p(p+1)/(n-p-1)
#'
#' @importFrom stats logLik AIC
#' @param object This object has to be of class 'GAM' or 'gam'.
#'
#' @return Returns the AICc coefficient.
#'
#' @examples
#' library(mgcv)
#' set.seed(123)
#' n <- 100
#' x <- rnorm(n)
#' y <- sin(x) + rnorm(n, sd = 0.2)
#' gam_fit <- gam(y ~ s(x))
#' AICc(gam_fit)
#'
#' @keywords internal
#' @export
#'
AICc <- function(object){
  if(class(object)[1] != "Gam" & class(object)[1] != "gam") stop("'object' must be of class 'Gam' or 'gam'")
  ll <- logLik(object)
  d <- attributes(ll)$df
  if(class(object)[1]=="Gam") n <- attributes(ll)$nobs
  if(class(object)[1]=="gam") n <- object$df.null + 1
  aic <- -2*logLik(ll)[1]+2*attributes(ll)$df
  caic <- aic + 2*d*(d+1)/(n-d-1)
  caic[d+1 >= n] <- Inf
  attributes(caic)[c('df','nobs','class')] <- NULL
  caic
}

