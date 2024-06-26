% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AICc.R
\name{AICc}
\alias{AICc}
\title{AICc: Calculation for GAM Models}
\usage{
AICc(object)
}
\arguments{
\item{object}{An object of class "Gam" or "gam" representing a fitted GAM model.}
}
\value{
A numeric value representing the AICc for the GAM model.
}
\description{
Calculates AICc (Akaike Information Criterion corrected for small sample sizes)
for Generalized Additive Models (GAM).
}
\details{
This function calculates the AICc value (Akaike Information Criterion corrected
for small sample sizes) for a given GAM model. AICc is a measure of model fit
that penalizes the number of parameters in the model to avoid overfitting.
}
\examples{
# Load dataset
data(PimaIndiansDiabetes2, package="mlbench");

# Fit GAM model
fit <- mgcv::gam(
  diabetes ~ s(age) + s(mass) + s(pedigree) + pressure + glucose,
  data=PimaIndiansDiabetes2,
  family=binomial
);

# Calculate AICc
AICc(fit);

}
\references{
Azevedo, M., Meira-Machado, L., Gude, F., and Araújo, A. (2024).
 Pointwise Nonparametric Estimation of Odds Ratio Curves with R:
 Introducing the flexOR Package. \emph{Applied Sciences}, \bold{14}(9), 1-17.
 \doi{10.3390/app14093897}
}
\seealso{
\code{\link{gam}}, \code{\link{logLik}}, \code{\link{AIC}}
}
\keyword{models}
\keyword{nonlinear}
\keyword{regression}
\keyword{smooth}
