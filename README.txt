###############################################################################
# Analysis code and data to accompany the manuscript:
#
#   Multilevel network meta-regression for general likelihoods: synthesis of 
#   individual and aggregate data with applications to survival analysis
#
###############################################################################

This archive contains two sets of analysis code, one which implements the 
analyses using the multinma R package, and another which calls Stan models 
directly. We recommend following the multinma scripts for most analysis
purposes; the latter may be useful for those wishing to further modify the 
models beyond the level of customisation offered by the multinma package.


# Analyses using multinma

The folder `multinma` contains R code to recreate the analyses in the
paper using the multinma R package:

 * `multinma/simulated_parametric_survival.R` recreates the simulated example,
   fitting parametric proportional hazards models to data simulated from a 
   Weibull distribution.

 * `multinma/ndmm_survival.R` recreates the newly diagnosed multiple myeloma 
   example, fitting models with a flexible M-spline baseline hazard. Synthetic 
   data are included in the multinma R package.


# Analyses using Stan directly

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
