#' @title Predict Method for OR Objects
#'
#' @description
#' Predicts values using a fitted OR model.
#'
#' @param object An object of class "OR."
#' @param predictor The predictor variable for which you want to make predictions.
#' @param prob Probability value for prediction. Use 0 for point prediction, 0.5 for median, or a custom value between 0 and 1.
#' @param pred.value Optional custom prediction value (use with prob=NULL).
#' @param conf.level Confidence level for prediction intervals (default is 0.95).
#' @param prediction.values Vector of specific prediction values to calculate.
#' @param round.x Number of decimal places to round the prediction values (default is 5).
#' @param ref.label Label for the predictor variable in the output (optional).
#' @param ... Additional arguments (not used in this function).
#'
#' @details
#' This function predicts values and prediction intervals using a fitted OR model.
#'
#' @return
#' A matrix with predicted values and prediction intervals.
#'
#' @examples
#' \dontrun{
#' # Example of how to use predict.OR
#' library(gam)
#' library(mlbench)
#' data(PimaIndiansDiabetes2)
#' str(PimaIndiansDiabetes2)
#'
#' fit <- gam(diabetes ~ s(age), data=PimaIndiansDiabetes2, family=binomial(), x=TRUE)
#' summary(fit)
#' plot(fit, se=T, ask=TRUE)
#'
#' pdval <- c(21, 25, 29, 30, 31, 34, 35, 38, 40, 41, 41, 45, 50, 60, 65, 80)
#' predict.OR(fit, predictor="age", pred.value=30, conf.level=0.95, prediction.values=pdval)
#' }
#'
#'
#' @keywords
#' prediction, odds ratio, confidence interval, OR, predict method
#' @export

predict.OR <- function(object, predictor, prob=NULL, pred.value=NULL, conf.level=0.95, prediction.values=NULL, round.x=NULL, ref.label=NULL, ...) {
  if ( missing(object) ) {stop("Missing object");}
  if ( missing(predictor) ) {stop("Missing predictor");}
  if ( missing(round.x) ) {round.x <- 5;}
  if ( !missing(pred.value) ) {prob <- 0.5; }
  if ( missing(prob) & missing(pred.value) ) {prob <- 0;}

  model <- object
  mydata <- model$data                          #retirar obs com NA restringindo a bd as variaveis em object
  p.vars <- match(all.vars(model$formula), names(mydata))
  mydata <- mydata[,p.vars]
  mydata <- na.omit(mydata)

  fit <- model
  if (class(fit$x)[1] == "list") fit <- update(fit, . ~ ., x = T)
  ctype <- "FALSE";
  qvalue <- (1+conf.level)/2;
  linear.predictor <- FALSE;
  k1 <- 9999;
  k <- which(names(mydata) == predictor);
  k <- c(k, k1);
  if (k[1] == 9999) {stop("predictor must be in data");}
  k <- k[1];
  if ( missing(prediction.values) ) {prediction.values <- c( min(mydata[,k]), quantile( mydata[,k], c(0.05, 0.25, 0.5, 0.75, 0.95) ), max(mydata[,k]) );}
  prediction.values <- sort(unique(c(prediction.values, pred.value)))
  # print(prediction.values)
  if ( min(prediction.values) < min(mydata[,k]) | max(prediction.values) > max(mydata[,k]) ) {stop("prediction.values must be between minimum and maximum of the predictor");}
  a <- mydata; # a is our dataset
  n.predictor <- names(a)[k];
  n <- dim(a)[1];
  if (prob == 0) {
    eta.no.ref <- predict(fit, type="terms");
    if ( inherits(eta.no.ref, "numeric") ) {kp <- 1; eta.no.ref <- cbind(eta.no.ref, eta.no.ref);}
    else {kp <- grep( predictor, colnames(eta.no.ref) );}
    eta.xref <- min(eta.no.ref[,kp]);
    ii <- which.min(eta.no.ref[,kp]);
    xref <- a[ii,k];
    eta.ref <- eta.no.ref[,kp]-eta.xref;
    indices <- grep(names(a)[k], dimnames(fit$x)[[2]]);
    submatriz.diseno <- fit$x[,indices];
    if (is.matrix(submatriz.diseno) == FALSE) {linear.predictor <- TRUE;}
    submatriz.var <- vcov(fit)[indices,indices]
    xref1 <- rep(fit$x[ii,indices], dim(fit$x)[1]);
    if (linear.predictor == FALSE) {
      xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=dim(submatriz.diseno)[2], byrow=TRUE);
    }
    if (linear.predictor == TRUE) {
      xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=1, byrow=TRUE);
    }
    eta.ref1 <- fit$x[,indices]-xref1;
    var.eta.ref1 <- rep(NA, n);
    for (i in 1:n) {var.eta.ref1[i] <- eta.ref1[i,]%*%vcov(fit)[indices,indices]%*%eta.ref1[i,];}
    se.eta.ref1 <- sqrt(var.eta.ref1);
  }
  if (prob > 0 & prob < 1) {
    eta.no.ref <- predict(fit, type="terms");
    if ( inherits(eta.no.ref, "numeric") ) {kp <- 1; eta.no.ref <- cbind(eta.no.ref, eta.no.ref);}
    else {kp <- grep( predictor, colnames(eta.no.ref) );}
    ord <- order(a[,k]);
    if ( !missing(pred.value) ) {
      pp <- seq(0, 1, len=200);
      app <- quantile(a[,k], pp);
      qq <- which(app <= pred.value);
      qq1 <- max(qq);
      prob <- qq1/200;
    }
    ind.prob <- trunc(prob*n);
    xref <- a[,k][ord[ind.prob]];
    #print(xref)
    eta.xref <- eta.no.ref[,kp][ord[ind.prob]];
    eta.ref <- eta.no.ref[,kp]-eta.xref;
    indices <- grep(names(a)[k], dimnames(fit$x)[[2]]);
    submatriz.diseno <- fit$x[,indices];
    if ( !is.matrix(submatriz.diseno) ) {linear.predictor <- TRUE;}
    submatriz.var <- vcov(fit)[indices,indices]
    #submatriz.var <- fit$var[indices,indices];      #bb <- predict(fit, type = "terms", se.fit = T)$se.fit
    xref1 <- rep(fit$x[ord[ind.prob],indices], dim(fit$x)[1]);
    if (linear.predictor == FALSE) {
      xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=dim(submatriz.diseno)[2], byrow=TRUE);
    }
    if (linear.predictor == TRUE) {
      xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=1, byrow=TRUE);
    }
    eta.ref1 <- fit$x[,indices]-xref1;
    var.eta.ref1 <- rep(NA, n);
    for (i in 1:n) {var.eta.ref1[i] <- eta.ref1[i,]%*%vcov(fit)[indices,indices]%*%eta.ref1[i,];}
    se.eta.ref1 <- sqrt(var.eta.ref1);
    #print(se.eta.ref1)
  }
  if (prob == 1) {
    eta.no.ref <- predict(fit, type="terms");
    if ( inherits(eta.no.ref, "numeric") ) {kp <- 1; eta.no.ref <- cbind(eta.no.ref, eta.no.ref);}
    else {kp <- grep( predictor, colnames(eta.no.ref) );}
    eta.xref <- max(eta.no.ref[,kp]);
    ii <- which.max(eta.no.ref[,kp]);
    xref <- a[ii,k];
    eta.ref <- eta.no.ref[,kp]-eta.xref;
    indices <- grep(names(a)[k], dimnames(fit$x)[[2]]);
    submatriz.diseno <- fit$x[,indices];
    if ( !is.matrix(submatriz.diseno) ) {linear.predictor <- TRUE;}
    submatriz.var <- vcov(fit)[indices,indices]
    xref1 <- rep(fit$x[ii,indices], dim(fit$x)[1]);
    if (linear.predictor == FALSE) {
      xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=dim(submatriz.diseno)[2], byrow=TRUE);
    }
    if (linear.predictor == TRUE) {
      xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=1, byrow=TRUE);
    }
    eta.ref1 <- fit$x[,indices]-xref1;
    var.eta.ref1 <- rep(NA, n);
    for (i in 1:n) {var.eta.ref1[i] <- eta.ref1[i,]%*%vcov(fit)[indices,indices]%*%eta.ref1[i,];}
    se.eta.ref1 <- sqrt(var.eta.ref1);
  }
  tmat <- cbind(eta.ref, eta.ref-qnorm(qvalue)*se.eta.ref1, eta.ref+qnorm(qvalue)*se.eta.ref1);
  line <- rep(0, n);
  jj <- match(sort(unique(a[,k])), a[,k]);
  pp <- length(jj);
  mat2 <- c(a[jj[1],k], tmat[jj,1][1], tmat[jj,2][1], tmat[jj,3][1]);
  for (b in 2:pp) {mat2 <- rbind( mat2, c(a[jj[b],k], tmat[jj,1][b], tmat[jj,2][b], tmat[jj,3][b]) );}
  matriz <- matrix( 0, ncol=4, nrow=length(prediction.values) );
  matriz[,1] <- prediction.values;
  for ( qq in 1:length(prediction.values) ) {
    matriz[qq,2] <- approx(x=mat2[,1], y=mat2[,2], xout=matriz[qq,1])[[2]];
    matriz[qq,3] <- approx(x=mat2[,1], y=mat2[,3], xout=matriz[qq,1])[[2]];
    matriz[qq,4] <- approx(x=mat2[,1], y=mat2[,4], xout=matriz[qq,1])[[2]];
  }
  c1 <- conf.level;
  if ( missing(ref.label) ) {mat.name <- c( names(a)[k], "LnOR", paste("lower .", conf.level*100, sep=""), paste("upper .", conf.level*100, sep="") );}
  else {mat.name <- c( ref.label, "LnOR", paste("lower .", conf.level*100, sep=""), paste("upper .", conf.level*100, sep="") );}
  colnames(matriz) <- mat.name;
  rownames(matriz) <- rep( "", length(prediction.values) );
  return(matriz);
}