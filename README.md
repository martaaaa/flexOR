# flexOR
Flexible Estimation of Odds Ratio Curves: Introducing the flexOR Package

## Description
Explore the relationship between continuous predictors and binary outcomes with flexOR, an R package designed for robust nonparametric estimation of odds ratio curves. 
Overcome limitations of traditional regression methods by leveraging smoothing techniques, particularly spline-based methods, providing adaptability to complex datasets. 
The package includes options for automatic selection of degrees of freedom in multivariable models, enhancing adaptability to diverse datasets and intuitive visualization functions facilitate the interpretation and presentation of estimated odds ratio curves.

## Installation
If you want to use the release version of the **flexOR** package, you can install the package from CRAN as follows:
```r
install.packages(pkgs="flexOR");
```
If you want to use the development version of the **flexOR** package, you can install the package from GitHub via the [**remotes**](https://remotes.r-lib.org) package:
```r
remotes::install_github(
  repo="martaaaa/flexOR",
  build=TRUE,
  build_manual=TRUE
);
```

## Authors
Marta Azevedo <marta.vasconcelos4@gmail.com> and Luís Meira-Machado <lmachado@math.uminho.pt> \
Maintainer: Marta Azevedo <marta.vasconcelos4@gmail.com>

## Funding
This research was financed by FCT - Fundação para a Ciência e a Tecnologia, under Projects UIDB/00013/2020, UIDP/00013/2020, and EXPL/MAT-STA/0956/2021.
## References

Hosmer, D. W.; Lemeshow, S.; Sturdivant, R. X. Applied Logistic Regression, 3rd ed.;
Wiley, 2013.

Royston, P.; Altman, D.G.; Sauerbrei, W. Dichotomizing continuous predictors in multiple
regression: A bad idea. Statistics in Medicine 2006, 25, 127–141.

Hastie, T. J. and Tibshirani, R. J. Generalized Additive Models; Chapman & Hall/CRC:
New York, USA, 1990.

Wood, S. Generalized Additive Models: An Introduction with R; Chapman & Hall/CRC:
London, UK, 2017.

Akaike, H. A new look at the statistical model identification. IEEE Transactions on
Automatic Control 1974, 19, 716–723.

Hurvich, C. M.; Simonoff, J. S.; Tsai, Ch. L. Smoothing parameter selection in nonpara-
metric regression using an improved akaike information criterion. Journal of the Royal
Statistical Society, Series B 1998, 60, 271–293.

Schwarz, G. E. Estimating the dimension of a model. Annals of Statistics 1978, 6(2),
461–464.

Cadarso-Suárez, C.; Meira-Machado, L.; Kneib, T.; Gude, F. Flexible hazard ratio curves
for continuous predictors in multi-state models. Statistical Modelling 2010, 10(3),
291–314.

Meira-Machado, L.; Cadarso-Su ́arez, C.; Ara ́ujo, A.; Gude, F. smoothHR: An R Pac-
kage for Pointwise Nonparametric Estimation of Hazard Ratio Curves of Continuous
Predictors. Comput. Math. Methods Medicine 2013.
de Boor, C. A Practical Guide to Splines (Rev. Edn); Springer, New York, 2001.

Wood, S.; Pya, N.; , A.; S ̈afken, B. Smoothing parameter and model selection for gene-
ral smooth models (with discussion). Journal of the American Statistical Association
2016, 111, 1548-1575