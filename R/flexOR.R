#' flexOR: Flexible Odds Ratio Computation for GAM Models
#'
#' Calculate odds ratios (ORs) and associated confidence intervals (CIs) for specified
#' predictors in generalized additive models (GAMs).
#'
#' @aliases flexOR
#'
#' @description
#' The `flexOR` function computes odds ratios and CIs for predictors in GAM models.
#' It provides flexibility in specifying predictors using either a data frame, a response
#' variable and a formula, or a pre-fitted GAM model. The function is useful for
#' understanding the impact of predictors on binary outcomes in GAMs.
#'
#' @usage
#' flexOR(data, response = NULL, formula = NULL, gamfit)
#'
#' @param data A data frame containing the variables.
#' @param response The response variable as a character string.
#' @param formula A formula specifying the model if not using a pre-fitted GAM.
#' @param gamfit A pre-fitted GAM model (class 'Gam').
#'
#' @details
#' The `flexOR` function calculates odds ratios and CIs for specified predictors in
#' a GAM model. It accepts three different ways of specifying the model: by providing
#' the data frame and response variable, by specifying the formula, or by providing
#' a pre-fitted GAM model.
#'
#' @return
#' A list containing the following components:
#' \itemize{
#'   \item \code{dataset}: The dataset used for the analysis.
#'   \item \code{gamfit}: The fitted GAM model.
#' }
#'
#' @references
#' Include relevant references here.
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
#'   y = rbinom(100, 1, 0.5)
#' )
#'
#' # Fit a GAM model
#' fit <- gam(y ~ s(x), data = data1, family = binomial)
#'
#' # Calculate odds ratios using flexOR
#' result <- flexOR(data = data1, response = "y", gamfit = fit)
#'
#' # Print the odds ratios and CIs
#' print(result)
#' }
#'
#' @keywords GAM odds-ratio binary-data confidence-interval


flexOR <- function(data, response=NULL, formula=NULL, gamfit) {
  modelfit <- "TRUE";
  mydata2 <- deparse( substitute(data) );
  if ( missing(gamfit) ) {modelfit <- "FALSE";}
  if ( !missing(gamfit) ) {
    if ( !inherits(gamfit, "Gam") ) {stop("Argument gamfit must be of class Gam");}
  }
  if ( !is.data.frame(data) ) {stop("data must be of class data.frame");}
  if (modelfit == "TRUE") {
    if ( !inherits(gamfit, "Gam") ) {stop("Object gamfit must be of class Gam");}
    if ( is.null(gamfit$x) ) {stop("The argumment x in the Gam object is missing");}
    fit <- gamfit;
  }
  #if ( !missing(formula) & !missing(gamfit) ) {stop("....only one is requested");}
  if ( missing(data) ) {stop("The argumment data is missing");}
  if (missing(response) & modelfit == "FALSE") {stop("The argumment response is missing");}
  if (missing(formula) & modelfit == "FALSE") {stop("The argumment formula is missing");}
  mydata <- data;
  if (modelfit == "FALSE") {
    p0 <- match(names(data),response, nomatch=0);
    p1 <- which(p0 == 1);
    ny <- data[,p1];
    if ( !missing(response) ) {response <- ny;}
    fmla <- attr(terms(formula), "term.labels");
    ncov <- length(fmla);
    colvar <- rep(0, ncov);
    for (k in 1:ncov) {
      if ( fmla[k] %in% names(data) ) {
        colvar[k] <- which(names(data) == fmla[k]);
      } else {
        for ( j in 1:ncol(data) ) {
          if ( any( grep(names(data)[j], fmla[k]) ) ) {colvar[k] <- j;}
        }
      }
    }
    if ( any(colvar == 0) ) {stop("'formula' must contain the right variables");}

    covar <- as.formula( paste( " ny ~ ", paste(fmla, collapse = "+") ) );
    fit <- gam(covar, data=data, x=TRUE, family=binomial);

  }
  a1 <- c();
  if (modelfit == "TRUE") {
    for(k in 1:dim(mydata)[2]) {
      if ( any( grep(names(mydata)[k], fit$call) ) ) {a1 <- c(a1,k);}
    }
    mydata <- mydata[,a1];
  }
  if (modelfit == "FALSE") {
    a1 <- p1;
    for (k in 1:dim(mydata)[2]) {
      if ( any( grep(names(mydata)[k], formula) ) ) {a1<-c(a1,k);}
    }
    mydata <- mydata[,a1];
  }
  if (modelfit == "TRUE") {
    if (as.name( toString(fit$call[[3]]) ) != mydata2) {cat("Warning: check if argument 'data' is the same used in gamfit", "\n");}
  }
  mydata <- na.omit(mydata);
  nv <- fit$nevent;
  object <- list(dataset=mydata, gamfit=fit);
  class(object) <- "OR";
  return(object);
}
