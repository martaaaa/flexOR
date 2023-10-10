# flexOR
The OR curve for a continuous predictor $Z$ in an additive logistic regression model can be written as $OR(Z,z_{ref})=exp(f(z)-f(z_{ref}))$ where $z_{ref}$ is a specific value of the predictor taken as the reference. After taking logarithms for simplicity, the asymptotic variance of $Ln\widehat{OR}(Z,z_{ref})$ can be expressed in terms of the covariance matrix of the spline estimate. Finally, assuming normality, $(1 - \alpha)100$\% pointwise confidence limits can be constructed around the $OR(Z,z_{ref})$ curve.

One disadvantage of using splines for modelling a continuous covariate’s effect is the difficulty in choosing the number and location of the knots between which the smooth line is drawn. An arbitrary choice of number of knots and/or arbitrary knot location can mask important features in the data. While too many knots can lead to oversmoothing, too few can lead to undersmoothing. The package flexor provides an R function that provides the optimal number of degrees of freedom in the multivariable additive logistic model. The optimal degree of smoothing is obtained by minimizing Akaike’s Information Criterion (AIC), a corrected version of this (AICc) or based on the Bayesian Information Criterion (BIC). The AIC and BIC scores have similar forms, differing only in the penalty coefficient. In both scores, the first term rewards goodness of fit whereas the second is a penalty that is an increasing function of the number of estimated parameters (df). The penalties in the expressions of the AIC and BIC ($2 × df$ and $\log(n) \times df$, resp.) discourage overfitting.
