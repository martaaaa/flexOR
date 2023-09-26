#' AICc Calculation for GAM Models
#'
#' Calculate AICc (Akaike Information Criterion corrected for small sample sizes)
#' for Generalized Additive Models (GAM).
#'
#' @param object An object of class "Gam" or "gam" representing a fitted GAM model.
#'
#' @details
#' This function calculates the AICc value (Akaike Information Criterion corrected
#' for small sample sizes) for a given GAM model. AICc is a measure of model fit
#' that penalizes the number of parameters in the model to avoid overfitting.
#'
#' @return
#' A numeric value representing the AICc for the GAM model.
#'
#' @examples
#' \dontrun{
#' # Load necessary libraries
#' library(mgcv)
#'
#' # Simulate data
#' set.seed(123)
#' data1 <- data.frame(
#'   x = rnorm(100),
#'   y = rnorm(100)
#' )
#'
#' # Fit a GAM model
#' mod <- gam(y ~ s(x), data = data1)
#'
#' # Calculate AICc
#' AICc_value <- AICc(mod)
#' }
#' @seealso
#' \code{\link{gam}}, \code{\link{logLik}}, \code{\link{AIC}}
#'
#' @keywords GAM AICc model-selection
#' @export

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

