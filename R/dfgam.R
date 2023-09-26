#' dfgam: Degrees of Freedom Selection for GAM Models
#'
#' Calculate the degrees of freedom for generalized additive models (GAMs) using
#' a selection method based on AIC, AICc, or BIC criteria.
#'
#' @aliases dfgam
#'
#' @description
#' The `dfgam` function calculates the degrees of freedom for specified non-linear
#' predictors in a GAM model. The user can choose between AIC (Akaike Information
#' Criterion), AICc (AIC corrected for small sample sizes), or BIC (Bayesian
#' Information Criterion) as the selection criteria. This function is useful for
#' determining the appropriate degrees of freedom for smoothing terms in GAMs.
#'
#' @param response The response variable as a formula.
#' @param nl.predictors A character vector specifying the non-linear predictors.
#' @param other.predictors A character vector specifying other predictors if needed.
#' @param smoother The type of smoothing term, currently only "s" is supported.
#' @param method The selection method, one of "AIC", "AICc", or "BIC".
#' @param data The data frame containing the variables.
#' @param step The step size for grid search when there are multiple non-linear predictors.
#'
#' @details
#' The `dfgam` function calculates the degrees of freedom for specified non-linear
#' predictors in a GAM model. It fits multiple GAMs with different degrees of freedom
#' and selects the best model based on the chosen selection method (AIC, AICc, or BIC).
#'
#' @return
#' A list containing the following components:
#' \itemize{
#'   \item \code{fit}: The fitted GAM model.
#'   \item \code{df}: A numeric vector of degrees of freedom for each non-linear predictor.
#'   \item \code{method}: The selection method used (AIC, AICc, or BIC).
#'   \item \code{nl.predictors}: The non-linear predictors used in the model.
#'   \item \code{other.predictors}: Other predictors used in the model if specified.
#' }
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
#' # Calculate degrees of freedom using AIC
#' result <- dfgam(response = y ~ x, nl.predictors = "x", data = data1)
#' }
#'
#' @keywords GAM degrees-of-freedom model-selection smoothing AIC AICc BIC
#' @export

dfgam <- function(response, nl.predictors, other.predictors=NULL, smoother="s", method = "AIC", data, step=NULL) {
  #options(warn=-1);
  if ( missing(data) ) {stop("The argument data is missing");}
  if ( missing(response) ) {stop("The argument response is missing");}
  if ( missing(nl.predictors) ) {stop("The argument 'nl.predictors' is missing");}
  if ( missing(smoother) ) {smoother <- "s";}
  if (smoother != "s") {stop("argument 'smoother' must be 's'");}
  if ( missing(method) ) {method <- "AIC";}
  if ( method != "AIC" & method != "AICc" & method != "BIC") {stop("The argument 'method' is not valid");}  # include other methods from mgcv?
  if (missing(step)) step <- 3
  if (step < 1 | step > 10) stop("'step' must be between 1 and 10")
  p0 <- match(names(data), response, nomatch=0);
  p1 <- which(p0 == 1);
  ny <- data[,p1];
  if (sum(p0) == 0) {stop("variable defined in argument 'response' is not in the dataset 'data'");}
  nnl <- length(nl.predictors);
  df <- rep(1, length(nl.predictors));
  fmla <- c();
  nl.fmla <- c();

  # df from simple regression models using mgcv::gam
  for (i in 1:nnl) {
    fmla[i] <- paste("s(", nl.predictors[i], ")", collapse="+")
    covar <- as.formula( paste( names(data)[p1], "~", paste(fmla[i], collapse = "+") ) );
    fit <- mgcv::gam(covar, data=data, family=binomial);    #melhor utilizar o modelo reg multipla logo a partida
    df[i] <- sum(fit$edf)-fit$nsdf
    fmla[i] <- c(paste("s(", nl.predictors[i], ",df=",df[i], ")", collapse=""))
    nl.fmla[i] <- paste("s(", nl.predictors[i], ",df=",df[i], ")", collapse="+")
  }

  if ( !missing(other.predictors) ) {
    nop <- length(other.predictors);
    for (i in 1:nop) {
      p2 <- match(names(data), other.predictors[i], nomatch=0);
      if (sum(p2) == 0) {stop("Check variables in argument 'other.predictors'");}
    }
  }


  if ( missing(other.predictors) & nnl > 1){
    fmla <- paste("s(", nl.predictors, ")", collapse="+");
    fmla4 <- as.formula( paste( names(data)[p1]," ~ ", fmla, collapse = "+") ) ;
    fit0 <- mgcv::gam(fmla4, data=data, family=binomial)
    df <- summary(fit0)$edf
  }

  if ( !missing(other.predictors) & nnl > 1){
    fmla <- paste("s(", nl.predictors, ")", collapse="+");
    fmla2 <- paste(other.predictors, collapse="+");
    fmla3 <- paste (c(fmla, fmla2), collapse="+");
    fmla4 <- as.formula( paste( names(data)[p1]," ~ ", fmla, collapse = "+") ) ;
    fit0 <- mgcv::gam(fmla4, data=data, family=binomial)
    df <- summary(fit0)$edf
  }


  if(nnl > 1){
    opfmla <- paste(other.predictors, collapse="+")
    covar2 <- c(paste( nl.fmla, collapse = "+"), opfmla)
    covar2 <- paste(covar2, collapse ="+")
    fmla2 <- as.formula( paste( names(data)[p1]," ~ ", covar2, collapse = "+") ) ;
    fit <- gam(fmla2, data=data, family=binomial)
  }


  mat <- array(NA, dim=c(41,4,nnl)) #df, AIC, BIC, cAIC     #Reduce grid to became faster!
  dfAIC <- df
  dfBIC <- df
  dfAICc <- df
  lmat <- dim(mat)[1]
  df.AIC <- df
  df.BIC <- df
  df.AICc <- df
  msg <- c()

  if(nnl == 1 & missing(other.predictors) ) {   #OK!
    covar <- as.formula( paste( names(data)[p1], "~", nl.predictors ) );
    fit <- gam(covar, data=data, family=binomial)
    lAIC <- AIC(fit)
    auxi <- round(df,1) + seq(-2,2,0.1)
    if (auxi[1] < 1) auxi <- 1.1 + auxi - auxi[1]
    mat[,1,1] <- auxi
    for(h in 1:lmat){
      ndf <- mat[h,1,1]
      aux <- paste("s(", nl.predictors, ",df=",ndf, ")", collapse="+")
      covar <- as.formula( paste( names(data)[p1], "~", paste(aux, collapse = "+") ) );
      fit <- gam(covar, data=data, family=binomial);
      mat[h,2,1] <- AIC(fit); mat[h,3,1] <- BIC(fit); mat[h,4,1] <- AICc(fit)
    }
    df.AIC <- mat[which.min(mat[,2,1]),1,1]
    df.BIC <- mat[which.min(mat[,3,1]),1,1]
    df.AICc <- mat[which.min(mat[,4,1]),1,1]
    if(lAIC < min(mat[,2,1])) msg <- paste("The effect of",nl.predictors,"is linear")
    if(method == "AIC") df <- df.AIC
    if(method == "BIC") df <- df.BIC
    if(method == "AICc") df <- df.AICc
  }

  if(nnl == 1 & !missing(other.predictors) ) {   #OK!

    auxop <- paste(other.predictors, collapse="+")
    auxnl <- paste(nl.predictors, collapse="+")
    aux0 <- paste(c(auxnl,auxop), collapse="+")
    fmla4 <- as.formula( paste( names(data)[p1]," ~ ", aux0, collapse = "+") ) ;
    fit0 <- gam(fmla4, data=data, family=binomial)
    lAIC <- AIC(fit0)

    auxi <- round(df,1) + seq(-2,2,0.1)
    if (auxi[1] < 1) auxi <- 1.1 + auxi - auxi[1]
    mat[,1,1] <- auxi
    for(h in 1:lmat){
      ndf <- mat[h,1,1]

      aux1 <- paste("s(", nl.predictors[1], ",df=",ndf, ")", collapse="+")
      aux2 <- paste(c(aux1,auxop), collapse="+")
      fmla3 <- as.formula( paste( names(data)[p1]," ~ ", aux2, collapse = "+") ) ;
      fit <- gam(fmla3, data=data, family=binomial)
      mat[h,2,1] <- AIC(fit); mat[h,3,1] <- BIC(fit); mat[h,4,1] <- AIC(fit)
    }
    df.AIC <- mat[which.min(mat[,2,1]),1,1]
    df.BIC <- mat[which.min(mat[,3,1]),1,1]
    df.AICC <- mat[which.min(mat[,4,1]),1,1]
    if(lAIC < min(mat[,2,1])) msg <- paste("The effect of",nl.predictors,"is linear")
    if(method == "AIC") df <- df.AIC
    if(method == "BIC") df <- df.BIC
    if(method == "AICc") df <- df.AICc
  }


  if(nnl > 1 & !missing(other.predictors) ) { #OK!
    for (m in 1:step){
      auxop <- paste(other.predictors, collapse="+")

      for(k in 1:nnl){

        auxi <- round(df[k],1) + seq(-2,2,0.1)
        if (auxi[1] < 1) auxi <- 1.1 + auxi - auxi[1]
        mat[,1,k] <- auxi
        for(h in 1:lmat){
          ndf <- df
          ndf[k] <- mat[h,1,k]

          aux1 <- paste("s(", nl.predictors, ",df=",ndf, ")", collapse="+")
          aux2 <- paste(c(aux1,auxop), collapse="+")
          fmla3 <- as.formula( paste( names(data)[p1]," ~ ", aux2, collapse = "+") ) ;
          fit <- gam(fmla3, data=data, family=binomial)
          mat[h,2,k] <- AIC(fit); mat[h,3,k] <- BIC(fit); mat[h,4,k] <- AICc(fit)
        }
        df.AIC[k] <- mat[which.min(mat[,2,k]),1,k]
        df.BIC[k] <- mat[which.min(mat[,3,k]),1,k]
        df.AICc[k] <- mat[which.min(mat[,4,k]),1,k]
        if(method == "AIC") df[k] <- df.AIC[k]
        if(method == "BIC") df[k] <- df.BIC[k]
        if(method == "AICc") df[k] <- df.AICc[k]
      }
    }
  }

  if(nnl > 1 & missing(other.predictors) ) { #ok
    for (m in 1:step){
      auxop <- paste(other.predictors, collapse="+")

      for(k in 1:nnl){

        auxi <- round(df[k],1) + seq(-2,2,0.1)
        if (auxi[1] < 1) auxi <- 1.1 + auxi - auxi[1]
        mat[,1,k] <- auxi
        for(h in 1:lmat){
          ndf <- df
          ndf[k] <- mat[h,1,k]

          aux1 <- paste("s(", nl.predictors, ",df=",ndf, ")", collapse="+")
          #aux2 <- paste(c(aux1,auxop), collapse="+")
          fmla3 <- as.formula( paste( names(data)[p1]," ~ ", aux1, collapse = "+") ) ;
          fit <- gam(fmla3, data=data, family=binomial)
          mat[h,2,k] <- AIC(fit); mat[h,3,k] <- BIC(fit); mat[h,4,k] <- AICc(fit)
        }
        df.AIC[k] <- mat[which.min(mat[,2,k]),1,k]
        df.BIC[k] <- mat[which.min(mat[,3,k]),1,k]
        df.AICc[k] <- mat[which.min(mat[,4,k]),1,k]
        if(method == "AIC") df[k] <- df.AIC[k]
        if(method == "BIC") df[k] <- df.BIC[k]
        if(method == "AICc") df[k] <- df.AICc[k]
      }
    }
  }

  #return
  res <- matrix(NA,ncol=1, nrow=nnl)
  colnames(res) <- "df"
  rownames(res) <- nl.predictors
  res[1:nnl] <- df[1:nnl]

  if(missing(other.predictors)){
    aux1 <- paste("s(", nl.predictors, ",df=",df, ")", collapse="+")
    fmla3 <- as.formula( paste( names(data)[p1]," ~ ", aux1, collapse = "+") ) ;
    fit <- gam(fmla3, data=data, family=binomial)
  }
  if(!missing(other.predictors)){
    auxop <- paste(other.predictors, collapse="+")
    aux1 <- paste("s(", nl.predictors, ",df=",ndf, ")", collapse="+")
    aux2 <- paste(c(aux1,auxop), collapse="+")
    fmla3 <- as.formula( paste( names(data)[p1]," ~ ", aux2, collapse = "+") ) ;
    fit <- gam(fmla3, data=data, family=binomial)
  }

  if(!is.null(msg)) print(msg)
  ob <- list(fit=fit, df=res, method=method, nl.predictors=nl.predictors, other.predictors=other.predictors)
  return(ob)
}
