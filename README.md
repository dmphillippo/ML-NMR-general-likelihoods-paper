# Multilevel network meta-regression for general likelihoods

This repository contains analysis code and data to accompany the manuscript:

> Phillippo, D. M., Dias, S., Ades, A. E., Welton, N. J. (2024). "Multilevel
> network meta-regression for general likelihoods: synthesis of individual and 
> aggregate data with applications to survival analysis".
> *arXiv*:[2401.12640](https://arxiv.org/abs/2401.12640).


Two sets of analysis code are provided, one which implements the analyses using
the [multinma R package](https://dmphillippo.github.io/multinma/), and another 
which calls Stan models directly. We recommend following the multinma scripts 
for most analysis purposes; the latter may be useful for those wishing to 
further modify the models beyond the level of customisation offered by the
multinma package.


## Analyses using multinma

The folder `multinma` contains R code to recreate the analyses in the
paper using the multinma R package:

 * `multinma/simulated_parametric_survival.R` recreates the simulated example,
   fitting parametric proportional hazards models to data simulated from a 
   Weibull distribution.

 * `multinma/ndmm_survival.R` recreates the newly diagnosed multiple myeloma 
   example, fitting models with a flexible M-spline baseline hazard. Synthetic 
   data are included in the multinma R package, along with a walkthrough 
   [vignette](https://dmphillippo.github.io/multinma/articles/example_ndmm.html)
   of this analysis.


## Analyses using Stan directly

The folder `rstan` contains R code to recreate the analyses in the
paper, setting up and calling Stan models directly via the rstan package:

 * `rstan/simulated_parametric_survival.R` recreates the simulated example,
   fitting parametric proportional hazards models to data simulated from a 
   Weibull distribution. The corresponding Stan code is found in
   `rstan/survival_param.stan`.

 * `rstan/ndmm_survival.R` recreates the newly diagnosed multiple myeloma 
   example, fitting models with a flexible M-spline baseline hazard. Synthetic
   data are constructed from the summary information in the folder `ndmm_data`.
   The corresponding Stan code is found in `rstan/survival_mspline.stan`.
