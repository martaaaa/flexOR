#' Smooth OR
#'
#' @notes fits a Generalized Additive Model (GAM) and creating a resulting object
#' containing the model outcomes, related to odds ratio analysis.
#' Here's a step-by-step description of what the function does:
#'
#' @param data A data frame containing the dataset for analysis.
#' @param response The response variable in the data frame. Required when fitting a new model (modelfit = "FALSE")
#' @param formula The formula describing the model. Required when fitting a new model (modelfit = "FALSE").
#' @param gamfit A pre-fitted GAM object. When provided, modelfit is set to "TRUE" and the function processes this object.
#'
#' @examples
#' mod1 <- smoothOR(data=data1, response=variable, formula=~covariate1+covariate2)
#' @export
#' @return An object of class "OR" containing the processed data and analysis results.
#'
#' The returned object includes relevant information from the fitted Generalized Additive Model (GAM) or
#' analysis of the provided GAM object. It provides insights into odds ratios and their analysis.
#' Users can further explore and interpret the results using appropriate methods for class "OR".
#'
#' @seealso \code{\link{gam}}, \code{\link{summary}}, \code{\link{plot}}
#' @keywords internal


smoothOR <- function(data, response=NULL, formula=NULL, gamfit) {
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
