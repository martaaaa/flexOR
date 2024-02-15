#' @title AICc: Calculation for GAM Models
#'
#' @description
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
#' library(mlbench)
#' 
#' # Load dataset
#' data(PimaIndiansDiabetes2)
#' 
#' # Fit a GAM model
#' fit <- gam(diabetes ~ s(age) + s(mass) + s(pedigree) + pressure + glucose, data = PimaIndiansDiabetes2, family = binomial)
#'
#' # Calculate AICc
#' AICc_value <- AICc(fit)
#' }
#' @seealso
#' \code{\link{gam}}, \code{\link{logLik}}, \code{\link{AIC}}
#'
#' @keywords GAM AICc model-selection
#' @export

AICc <- function(object) {
  if (!(inherits(object, "Gam") || inherits(object, "gam"))) {
    stop("'object' must be of class 'Gam' or 'gam'")
  }
  
  if (!is.null(findMethod("logLik", class(object)))) {
    ll <- logLik(object)
    d <- attributes(ll)$df
    
    if (inherits(object, "Gam")) {
      n <- attributes(ll)$nobs
    } else if (inherits(object, "gam")) {
      n <- object$df.null + 1
    }
    
    aic <- -2 * logLik(ll)[1] + 2 * attributes(ll)$df
    caic <- aic + 2 * d * (d + 1) / (n - d - 1)
    caic[d + 1 >= n] <- Inf
    attributes(caic)[c('df', 'nobs', 'class')] <- NULL
    return(caic)
  } else {
    stop("No applicable method for 'logLik' applied to an object of class ", class(object))
  }
}
