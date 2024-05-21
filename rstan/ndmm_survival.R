################################################################################
# Analysis of Newly Diagnosed Multiple Myeloma PFS using flexible M-splines on 
# the baseline hazard with Multilevel Network Meta-Regression
################################################################################

library(dplyr)     # For data manipulation and plotting
library(readr)
library(tidyr)
library(purrr)
library(ggplot2)

library(simsurv)   # For simulating survival data
library(estmeansd)

library(survival)  # For Kaplan-Meier plots

library(multinma)  # Used here for setting up integration points and RW1 prior

library(rstan)
library(splines2)
library(loo)
library(broom.mixed)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Colour palettes
trt_pal <- c(Pbo = "#113259", Len = "#865800", Thal = "#55A480")

# Set seed
set.seed(54321)


# Set up output directories -----------------------------------------------

if (!dir.exists("./figures")) dir.create("./figures")
if (!dir.exists("./tables")) dir.create("./tables")


# Construct synthetic data ------------------------------------------------

# These synthetic data are available directly as part of the multinma package.
# Here we show how the data are constructed from Generalised-F distributions fit
# to digitized Kaplan-Meier curves and published summary data on covariates and 
# regression coefficients.

# Read in covariate summaries, survival distribution estimates
study_covs <- read_csv("../ndmm_data/covariate_summaries.csv")
surv_pars <- read_csv("../ndmm_data/genf_pars.csv")
surv_maxt_pcens <- read_csv("../ndmm_data/maxt_pcens.csv")

# Reconstruct distribution for age
study_covs <- study_covs %>% 
  rowwise() %>% 
  mutate(age_fit = list(bc.mean.sd(age_min, age_iqr_l, age_median, age_iqr_u, age_max, sample_size)))

purrr::map(study_covs$age_fit, ~summary(.$bc.norm.rvs))
withr::with_par(list(mfrow = c(2, 5)),
                purrr::walk(study_covs$age_fit, ~hist(.$bc.norm.rvs)))

# Simulate covariates
ipd_all <-
  study_covs %>% 
  group_by(study, treatment) %>% 
  transmute(age = list(sample(age_fit[[1]]$bc.norm.rvs, size = sample_size, replace = FALSE)),
            iss_stage3 = list(rbinom(sample_size, 1, iss_stage3)),
            response_cr_vgpr = list(rbinom(sample_size, 1, response_cr_vgpr)),
            male = list(rbinom(sample_size, 1, male))) %>% 
  unnest(cols = c(age, iss_stage3, response_cr_vgpr, 
                  male)) %>% 
  # Add treatment and study factor variables
  mutate(studyf = factor(study, levels = c("Attal2012", "McCarthy2012", "Palumbo2014", "Jackson2019", "Morgan2012")),
         trtf = factor(treatment, levels = c("Pbo", "Len", "Thal")),
         .after = treatment)

by(ipd_all, ipd_all$studyf, summary)


# Simulate survival outcomes using a generalised F model
# The gen-F parameters were obtained from unadjusted fits to the digitized KM curves
Hgenf_reg <- function(t, x, betas, intercept, trteff, sigma, Q, P) {
  xtemp <- as.matrix(x[, names(betas)])
  flexsurv::Hgenf(t, mu = intercept + trteff + xtemp %*% betas, sigma = sigma, Q = Q, P = P)
}

# Prognostic and EM coefficients, informed by Leahy and Walsh (2019)
betas <- -c(
  age = 0.065, 
  iss_stage3 = 0.27,
  response_cr_vgpr = -0.2, 
  male = 0.05,
  `age:trt` = -0.024, 
  `iss_stage3:trt` = 0.26, 
  `response_cr_vgpr:trt` = 0.31, 
  `male:trt` = 0
)

surv_all <- ipd_all %>%
  group_by(study, treatment) %>% 
  mutate(across(age:male, ~.x - mean(.x))) %>% 
  mutate(`age:trt` = if_else(treatment != "Pbo", age, 0),
         `iss_stage3:trt` = if_else(treatment != "Pbo", iss_stage3, 0),
         `response_cr_vgpr:trt` = if_else(treatment != "Pbo", response_cr_vgpr, 0),
         `male:trt` = if_else(treatment != "Pbo", male, 0)) %>% 
  nest(x = age:`male:trt`) %>% 
  left_join(surv_pars) %>% 
  left_join(surv_maxt_pcens) %>% 
  mutate(trteff = recode(treatment,
                         Pbo = 0,
                         Len = treatmentLen,
                         Thal = treatmentThal),
         surv = list(simsurv(cumhazard = Hgenf_reg,
                             x = x[[1]],
                             betas = betas,
                             intercept = mu,
                             sigma = sigma,
                             P = P,
                             Q = Q,
                             trteff = trteff,
                             interval = c(1e-08, 10000),
                             maxt = max_time,
         ))
  ) %>% 
  unnest(cols = c(x, surv)) %>% 
  mutate(# Apply uniform censoring
    status = if_else(runif(n()) <= 0.2, 0L, status)
  )

# Quick check of regression
coxph(Surv(eventtime, status) ~ study + male + (age + iss_stage3 + response_cr_vgpr)*trtf -1, data = surv_all)

# Produce KM plots
kmdat <- surv_all %>% 
  group_by(studyf, trtf) %>% 
  do(broom::tidy(survfit(Surv(eventtime, event = status) ~ 1, data = .))) %>% 
  # Add S(0) = 1
  bind_rows(
    distinct(., studyf, trtf, .keep_all = TRUE) %>%
      mutate(time = 0, n.event = 0, n.censor = 0, estimate = 1, std.error = 0, conf.high = 1, conf.low = 1)
  ) %>%
  arrange(studyf, trtf, time)

ggplot(kmdat, aes(x = time, y = estimate, colour = trtf)) +
  geom_step() +
  geom_point(shape = 3, data = function(x) filter(x, n.censor >= 1)) +
  facet_wrap(~studyf) +
  labs(y = "Survival probability", x = "Time (months)") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual("Treatment", 
                      breaks = c("Pbo", "Len", "Thal"), 
                      values = trt_pal) +
  theme_multinma() +
  theme(legend.text = element_text(size = 9), legend.title = element_text(size = 9),
        legend.position = "top", 
        legend.box.spacing = unit(0, "lines"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
        panel.border = element_rect(colour = "grey70", fill = NA),
        strip.background = element_rect(colour = "grey70", fill = "grey90"))

ggsave("./figures/ndmm_survival_sim_curves.pdf", width = 6, height = 5, scale = 1.2)


# Create aggregate data for Morgan2012 and Jackson2016
agd_covs <- ipd_all %>% filter(studyf %in% c("Morgan2012", "Jackson2019")) %>%
  group_by(study, studyf, treatment, trtf) %>%
  summarise(sample_size = n(),
            as_tibble(as.list(setNames(fivenum(age), c("age_min", "age_iqr_l", "age_median", "age_iqr_h", "age_max")))),
            # Use estmeansd to calculate mean and sd
            setNames(as_tibble(unclass(estmeansd::bc.mean.sd(age_min, age_iqr_l, age_median, age_iqr_h, age_max, n = sample_size))[c("est.mean", "est.sd")]),
                     c("age_mean", "age_sd")),
            iss_stage3 = mean(iss_stage3),
            response_cr_vgpr = mean(response_cr_vgpr),
            male = mean(male))

agd_surv <- surv_all %>% filter(studyf %in% c("Morgan2012", "Jackson2019")) %>%
  select(study, studyf, treatment, trtf, eventtime, status)


# Data for for IPD studies
ipd_surv <- ipd_all %>% filter(!studyf %in% c("Morgan2012", "Jackson2019")) %>% 
  bind_cols(
    surv_all %>% 
      ungroup() %>% 
      filter(!studyf %in% c("Morgan2012", "Jackson2019")) %>% 
      select(eventtime, status)
  )


# ML-NMR ------------------------------------------------------------------


## M-spline 7 knots -------------------------------------------------------

# Ensure that study data are contiguous
ipd_surv <- arrange(ipd_surv, studyf, trtf)
agd_surv <- arrange(agd_surv, studyf, trtf)
agd_covs <- arrange(agd_covs, studyf, trtf)

n_int <- 64  # number of sample points for numerical integration

# Derive weighted z-transformed average covariate correlation matrix from IPD studies
ipd_cors <- dplyr::group_by(ipd_surv, studyf) %>%
  dplyr::group_modify(~tibble::tibble(
    w = nrow(.) - 3,
    r = list(cor(dplyr::select(., age, iss_stage3, response_cr_vgpr, male),
                 method = "pearson",
                 use = "complete.obs"))
  )) %>%
  dplyr::mutate(z = purrr::map2(w, r, ~.x * log((1 + .y) / (1 - .y)) / 2))

w_z <- Reduce(`+`, ipd_cors$z) / sum(ipd_cors$w)
ipd_cor <- (exp(2 * w_z) - 1) / (exp(2 * w_z) + 1)

diag(ipd_cor) <- 1

# Add integration points
agd_covs <- add_integration(agd_covs,
                            age = distr(qgamma, mean = age_mean, sd = age_sd),
                            iss_stage3 = distr(qbern, iss_stage3),
                            response_cr_vgpr = distr(qbern, response_cr_vgpr),
                            male = distr(qbern, male),
                            cor = ipd_cor,
                            n_int = n_int)

agd_xpoints <- unnest_integration(agd_covs)

# Derive global means for centering 
age_gmean <- weighted.mean(c(ipd_surv$age, agd_covs$age_mean), 
                           c(rep_len(1, nrow(ipd_surv)), agd_covs$sample_size))

# Now put all the IPD and AgD integration points into one data frame
stan_xdat <- ipd_surv %>%
  select(-eventtime, -status) %>%
  bind_rows(agd_xpoints) %>%
  mutate(
    # Center the continuous covariates
    age = age - age_gmean,
    # Add in treatment class
    trtclass = factor(if_else(trtf == "Pbo", "Pbo", "Active"), levels = c("Pbo", "Active"))
  )

# Create the model matrix
MLNMR_formula <- ~ -1 + studyf + trtf + 
  age + iss_stage3 + response_cr_vgpr + male +
  (age + iss_stage3 + response_cr_vgpr + male):trtclass

X_all <- model.matrix(MLNMR_formula, data = stan_xdat)

# Ensure columns are in the order: studies, treatments, regression terms
colnames(X_all)

# Thin QR decomposition
X_all_qr <- qr(X_all)
X_all_Q <- qr.Q(X_all_qr) * sqrt(nrow(X_all) - 1)
X_all_R <- qr.R(X_all_qr)[, sort.list(X_all_qr$pivot)] / sqrt(nrow(X_all) - 1)
X_all_R_inv <- solve(X_all_R)

# M-splines for baseline hazard
degree <- 3    # Degree for the M-splines (0 is piecewise constant)
n_iknots <- 7  # Number of internal knots

mspline_dat <- bind_rows(
  ipd = select(ipd_surv, studyf, trtf, eventtime, status),
  agd = select(agd_surv, studyf, trtf, eventtime, status),
  .id = "dtype"
) %>% 
  group_by(studyf, dtype) %>% 
  summarise(
    times = list(eventtime),
    # Boundary knots at 0 and last time in each study
    bknots = list(c(0, max(eventtime))),
    # Internal knots at quantiles of observed event times in each study
    iknots = list(quantile(eventtime[status == 1], probs = seq(1, n_iknots) / (n_iknots + 1))),
    # Construct basis and integrated basis
    basis = list(mSpline(eventtime, knots = unlist(iknots), degree = degree, 
                         intercept = TRUE, Boundary.knots = unlist(bknots))),
    ibasis = list(mSpline(eventtime, knots = unlist(iknots), degree = degree, 
                          intercept = TRUE, Boundary.knots = unlist(bknots), integral = TRUE))
  )

# Softmax and inverse softmax transforms
softmax <- function(x) {
  x0 <- c(0, x)
  exp(x0 - logsumexp(x0))
}

logsumexp <- function(x) {
  maxx <- max(x)
  max(x) + log(sum(exp(x - maxx)))
}

inv_softmax <- function(p) {
  log(p[-1]) - log(p[1])
}

# Derive prior mean and weights on logit spline coefficients
mspline_constant_hazard <- function(basis) {
  df <- ncol(basis)
  ord <- attr(basis, "degree") + 1
  iknots <- attr(basis, "knots")
  bknots <- attr(basis, "Boundary.knots")
  
  # Using approach of Jackson arXiv:2306.03957
  knots <- c(rep(bknots[1], ord), iknots, rep(bknots[2], ord))
  coefs <- (knots[(1:df) + ord] - knots[1:df]) / (ord * (diff(bknots)))
  
  # inverse softmax transform
  inv_softmax(coefs)
}

rw1_prior_weights <- function(basis) {
  nscoef <- ncol(basis)
  ord <- attr(basis, "degree") + 1
  iknots <- attr(basis, "knots")
  bknots <- attr(basis, "Boundary.knots")
  knots <- c(rep(bknots[1], ord), iknots, rep(bknots[2], ord))
  wts <- 1/(ord - 1) * (knots[(ord + 1):(nscoef + ord - 1)] - knots[2:nscoef])
  return(sqrt(wts / sum(wts)))
}

prior_lscoef_location <- map(mspline_dat$basis, mspline_constant_hazard)
prior_lscoef_weight <- map(mspline_dat$basis, rw1_prior_weights)

# Construct stan data
standat_MLNMR <- list(
  
  # Constants
  ns_ipd = n_distinct(ipd_surv$studyf),
  ns_agd = n_distinct(agd_surv$studyf),
  ni_ipd = nrow(ipd_surv),
  ni_agd = nrow(agd_surv),
  nt = n_distinct(stan_xdat$trtf),
  nX = ncol(X_all),
  nint = n_int,
  
  # Survival data
  ipd_status = ipd_surv$status,
  agd_status = agd_surv$status,
  
  # M-splines,
  n_scoef = n_iknots + degree + 1,
  ipd_time = do.call(rbind, filter(mspline_dat, dtype == "ipd")$basis),
  ipd_itime = do.call(rbind, filter(mspline_dat, dtype == "ipd")$ibasis),
  agd_time = do.call(rbind, filter(mspline_dat, dtype == "agd")$basis),
  agd_itime = do.call(rbind, filter(mspline_dat, dtype == "agd")$ibasis),
  
  # Auxiliary parameters (i.e. spline coefficients)
  # These are study ids (i.e. stratified by study)
  aux_id = c(as.numeric(ipd_surv$studyf), as.numeric(agd_surv$studyf)),    # assumes factors have same levels
  
  # Integration point info
  int_id =
    left_join(
      agd_surv,
      agd_covs %>% ungroup() %>% 
        mutate(ag_id = 1:nrow(.)) %>% 
        select(studyf, trtf, ag_id)
    ) %>% 
    pull(ag_id),
  
  # Design matrix
  X = X_all_Q,
  
  # QR decomposition
  QR = TRUE,
  R_inv = X_all_R_inv,
  
  # Priors
  prior_intercept_sd = 100,  # intercept
  prior_reg_sd = 100,        # regression terms
  prior_trt_sd = 100,        # treatment effects
  prior_sigma_sd = 1,        # RW1 smoothing SD
  prior_lscoef_location = prior_lscoef_location,
  prior_lscoef_weight = prior_lscoef_weight)


# Fit using rstan
mspline_MLNMR_stan <- stan("survival_mspline.stan", 
                           data = standat_MLNMR,
                           pars = c("mu", "beta", "gamma", "scoef", "log_lik"),
                           iter = 2000,
                           chains = 4,
                           # Reduce max_treedepth to speed up warm-up (default 10)
                           # Increase if you get max_treedepth warnings
                           control = list(max_treedepth = 7))


# Save stan object
saveRDS(mspline_MLNMR_stan, file = "./ndmm_mspline_MLNMR_stan.rds")
# mspline_MLNMR_stan <- readRDS("./ndmm_mspline_MLNMR_stan.rds")
# check_hmc_diagnostics(mspline_MLNMR_stan)

print(mspline_MLNMR_stan, pars = c("beta", "gamma"))


## M-spline 10 knots ------------------------------------------------------


# M-splines for baseline hazard
degree <- 3    # Degree for the M-splines (0 is piecewise constant)
n_iknots <- 10  # Number of internal knots

mspline10_dat <- bind_rows(
  ipd = select(ipd_surv, studyf, trtf, eventtime, status),
  agd = select(agd_surv, studyf, trtf, eventtime, status),
  .id = "dtype"
) %>% 
  group_by(studyf, dtype) %>% 
  summarise(
    times = list(eventtime),
    # Boundary knots at 0 and last time in each study
    bknots = list(c(0, max(eventtime))),
    # Internal knots at quantiles of observed event times in each study
    iknots = list(quantile(eventtime[status == 1], probs = seq(1, n_iknots) / (n_iknots + 1))),
    # Construct basis and integrated basis
    basis = list(mSpline(eventtime, knots = unlist(iknots), degree = degree, 
                         intercept = TRUE, Boundary.knots = unlist(bknots))),
    ibasis = list(mSpline(eventtime, knots = unlist(iknots), degree = degree, 
                          intercept = TRUE, Boundary.knots = unlist(bknots), integral = TRUE))
  )

prior_lscoef_location <- map(mspline10_dat$basis, mspline_constant_hazard)
prior_lscoef_weight <- map(mspline10_dat$basis, rw1_prior_weights)

# Construct stan data
standat_MLNMR_10kt <- list_modify(standat_MLNMR,
  # M-splines
  n_scoef = n_iknots + degree + 1,
  ipd_time = do.call(rbind, filter(mspline10_dat, dtype == "ipd")$basis),
  ipd_itime = do.call(rbind, filter(mspline10_dat, dtype == "ipd")$ibasis),
  agd_time = do.call(rbind, filter(mspline10_dat, dtype == "agd")$basis),
  agd_itime = do.call(rbind, filter(mspline10_dat, dtype == "agd")$ibasis),
  prior_lscoef_location = prior_lscoef_location,
  prior_lscoef_weight = prior_lscoef_weight)


# Fit using rstan
mspline_MLNMR_10kt_stan <- stan("survival_mspline.stan", 
                           data = standat_MLNMR_10kt,
                           pars = c("mu", "beta", "gamma", "scoef", "log_lik"),
                           iter = 10,#2000,
                           chains = 4,
                           # Reduce max_treedepth to speed up warm-up (default 10)
                           # Increase if you get max_treedepth warnings
                           control = list(max_treedepth = 7))


# Save stan object
saveRDS(mspline_MLNMR_10kt_stan, file = "./ndmm_mspline_MLNMR_10kt_stan.rds")
# mspline_MLNMR_10kt_stan <- readRDS("./ndmm_mspline_MLNMR_10kt_stan.rds")
# check_hmc_diagnostics(mspline_MLNMR_10kt_stan)

print(mspline_MLNMR_10kt_stan, pars = c("beta", "gamma"))



## Compare model fit ------------------------------------------------------

# Get model fit using LOO
(mspline_MLNMR_loo <- loo(mspline_MLNMR_stan))
(mspline_MLNMR_10kt_loo <- loo(mspline_MLNMR_10kt_stan))

loo_compare(list(`7 knots` = mspline_MLNMR_loo,
                 `10 knots` = mspline_MLNMR_10kt_loo))


## Non-PH -----------------------------------------------------------------


# M-splines for baseline hazard
degree <- 3    # Degree for the M-splines (0 is piecewise constant)
n_iknots <- 7  # Number of internal knots

mspline_nph_dat <- bind_rows(
  ipd = select(ipd_surv, studyf, trtf, eventtime, status),
  agd = select(agd_surv, studyf, trtf, eventtime, status),
  .id = "dtype"
) %>% 
  # Generate splines by treatment arm for non-PH model
  group_by(studyf, trtf, dtype) %>% 
  summarise(
    times = list(eventtime),
    # Boundary knots at 0 and last time in each study
    bknots = list(c(0, max(eventtime))),
    # Internal knots at quantiles of observed event times in each study
    iknots = list(quantile(eventtime[status == 1], probs = seq(1, n_iknots) / (n_iknots + 1))),
    # Construct basis and integrated basis
    basis = list(mSpline(eventtime, knots = unlist(iknots), degree = degree, 
                         intercept = TRUE, Boundary.knots = unlist(bknots))),
    ibasis = list(mSpline(eventtime, knots = unlist(iknots), degree = degree, 
                          intercept = TRUE, Boundary.knots = unlist(bknots), integral = TRUE))
  )

prior_lscoef_location <- map(mspline_nph_dat$basis, mspline_constant_hazard)
prior_lscoef_weight <- map(mspline_nph_dat$basis, rw1_prior_weights)

# Construct stan data
standat_MLNMR_nph <- list_modify(standat_MLNMR,
                                 # M-splines
                                 n_scoef = n_iknots + degree + 1,
                                 ipd_time = do.call(rbind, filter(mspline_nph_dat, dtype == "ipd")$basis),
                                 ipd_itime = do.call(rbind, filter(mspline_nph_dat, dtype == "ipd")$ibasis),
                                 agd_time = do.call(rbind, filter(mspline_nph_dat, dtype == "agd")$basis),
                                 agd_itime = do.call(rbind, filter(mspline_nph_dat, dtype == "agd")$ibasis),
                                 prior_lscoef_location = prior_lscoef_location,
                                 prior_lscoef_weight = prior_lscoef_weight,
                                 # aux IDs - now by treatment arm as well as study
                                 aux_id = bind_rows(ipd_surv, agd_surv) %>% group_by(studyf, trtf) %>% group_indices())


# Fit using rstan
mspline_MLNMR_nph_stan <- stan("survival_mspline.stan", 
                                data = standat_MLNMR_nph,
                                pars = c("mu", "beta", "gamma", "scoef", "log_lik"),
                                iter = 10,#2000,
                                chains = 4,
                                # Reduce max_treedepth to speed up warm-up (default 10)
                                # Increase if you get max_treedepth warnings
                                control = list(max_treedepth = 7))


# Save stan object
saveRDS(mspline_MLNMR_nph_stan, file = "./ndmm_mspline_MLNMR_nph_stan.rds")
# mspline_MLNMR_nph_stan <- readRDS("./ndmm_mspline_MLNMR_nph_stan.rds")
# check_hmc_diagnostics(mspline_MLNMR_nph_stan)

print(mspline_MLNMR_nph_stan, pars = c("beta", "gamma"))

# Spline coefficients are now stratified by treatment arm as well as study
print(mspline_MLNMR_nph_stan, pars = "scoef")
print(mspline_MLNMR_stan, pars = "scoef")


## Compare model fit ------------------------------------------------------

# Get model fit using LOO
(mspline_MLNMR_nph_loo <- loo(mspline_MLNMR_nph_stan))

loo_compare(list(PH = mspline_MLNMR_loo,
                 `non-PH` = mspline_MLNMR_nph_loo))


# Results -----------------------------------------------------------------

# Using PH M-spline model with 7 internal knots


## Baseline hazard --------------------------------------------------------

# Plot splines on baseline hazard
as.matrix(mspline_MLNMR_stan, pars = "scoef") %>% 
  as_tibble() %>% 
  mutate(.sample = 1:n(), .before = 1) %>% 
  pivot_longer(starts_with("scoef"), 
               names_pattern = "scoef\\[([0-9]+),([0-9]+)\\]", 
               names_to = c("studyf", ".value"),
               names_transform = list(studyf = function(x) factor(as.integer(x), labels = levels(mspline_dat$studyf)))) %>% 
  left_join(mspline_dat, by = "studyf") %>% 
  rename_with(~paste0("scoef_", .), .cols = matches("[0-9]+")) %>% 
  mutate(scoef = array_tree(select(., starts_with("scoef_")), 1)) %>% 
  rowwise() %>% 
  mutate(newx = list(seq(0, max(unlist(times)), length.out = 100)), 
         ms = list(predict(unlist(basis), unlist(newx)) %*% unlist(scoef))) %>%
  unnest(cols = c(newx, ms)) %>% 
  group_by(studyf, newx) %>% 
  summarise(ms_med = median(ms), ms_lo = quantile(ms, prob = 0.025), ms_hi = quantile(ms, prob = 0.975)) %>% 
  
  ggplot(aes(x = newx, y = ms_med, ymin = ms_lo, ymax = ms_hi)) +
  geom_vline(aes(xintercept = iknots), data = unnest(mspline_dat, cols = iknots), linetype = 2, colour = "grey60") +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  facet_wrap("studyf") +
  xlab("Time (months)") + ylab("Baseline hazard") +
  theme_multinma()

ggsave("./figures/ndmm_mspline_baseline_hazard.pdf", width = 6, height = 5, scale = 1.2)


## Survival curves --------------------------------------------------------

# Create population-average survival curves
get_survival_curve <- function(trtf,      # Treatment to get survival curve for (factor variable)
                               classf,    # Treatment class for EMs (factor variable)
                               studyf,    # Study to get survival curve for (factor variable)
                               regression = ~-1 + studyf + trtf, # Regression model fitted
                               xdat,     # Covariate data (either IPD or integration points)
                               stanfit,  # Fitted stan model object
                               times,    # Vector of time points to evaluate Sbar(t) at
                               basis,    # M-spline basis created with splines2::mSpline()
                               scoef_id = studyf,  # Spline coefficient vector parameter ID
                               ...) {    # Other arguments passed to tidyMCMC
  
  
  # Define survival function S(t)

  S <- function(t, eta, scoef, basis, ...) {
    multinma::pmspline(t, basis = basis, scoef = scoef, rate = exp(eta), lower.tail = FALSE)
  }
  
  # Set up prediction data
  xdat <- mutate(xdat, 
                 trtf = .env$trtf, classf = .env$classf,
                 # Dummy study variable
                 studyf = 1)
  
  # Construct design matrix
  regression <- update(regression, ~.+1)
  X <- model.matrix(regression, data = xdat)
  X <- X[,-1]
  
  post_beta <- as.matrix(stanfit,
                         pars = c(sprintf("mu[%d]", as.numeric(studyf)),
                                  sprintf("gamma[%d]", 1:(nlevels(trtf) - 1)),
                                  sprintf("beta[%d]", 1:(ncol(X) - nlevels(trtf)))))

  # Pull posterior samples of spline coefficient parameters
  post_scoef <- as.matrix(stanfit, pars = sprintf("scoef[%d, %d]", as.numeric(scoef_id), 1:ncol(basis)))
  
  # Evaluate the linear predictor eta at each posterior sample and each
  # point in xdat: eta is a matrix [iter x n_int]
  eta <- tcrossprod(post_beta, X)
  
  # For each individual point in postdat (row of eta):
  #   1. Evaluate S(t) at every point in xdat (column of eta)
  #   2. Take the mean to get Sbar(t)
  
  names(times) <- sprintf("t_%d", seq_along(times))
  
  Sbar <- pmap_dfr(list(.eta = array_tree(eta, 1),  # Rows (iters) of eta
                        .scoef = array_tree(post_scoef, 1)),  # Rows (iters) of scoef
                   function(.eta, .scoef) map(times, function(.t) mean(S(.t, eta = .eta, scoef = .scoef, basis = basis))))
  
  # Produce posterior statistics
  tidyMCMC(Sbar, conf.int = TRUE, ...) %>% 
    mutate(time = times) %>% 
    select(-term)
}

# Apply this function for each treatment in each study
mspline_surv_pred <- crossing(studyf = unique(surv_all$studyf), trtf = unique(surv_all$trtf)) %>% 
  # Add treatment classes
  mutate(classf = factor(if_else(trtf == "Pbo", "Pbo", "Active"), levels = c("Pbo", "Active"))) %>% 
  group_by(studyf, trtf) %>% 
  group_modify(
    ~get_survival_curve(trtf = .$trtf, 
                        classf = .$classf,
                        studyf = .$studyf, 
                        scoef_id = .$studyf,
                        basis = filter(mspline_dat, studyf == .$studyf) %>% unnest(basis) %>% pull(basis),
                        xdat = filter(stan_xdat, studyf == .$studyf),
                        regression = ~-1 + studyf + trtf + 
                                       age + iss_stage3 + response_cr_vgpr + male + 
                                       (age + iss_stage3 + response_cr_vgpr + male):classf,
                        stanfit = mspline_MLNMR_stan,
                        times = do.call(seq, c(filter(mspline_dat, studyf == .$studyf) %>% unnest(bknots) %>% pull(bknots), list(length.out = 50)))),
    .keep = TRUE)

# Plot estimated survival curves
ggplot(mspline_surv_pred, 
       aes(x = time, y = estimate, colour = trtf, fill = trtf)) +
  # KM curves
  geom_step(linewidth = 0.15, data = kmdat) +
  geom_point(stroke = 0.15, shape = 3, data = filter(kmdat, n.censor >= 1)) +
  # ML-NMR curves
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, linetype = 0) +
  geom_line(linewidth = 0.4) +
  facet_wrap("studyf") +
  labs(y = "Survival probability", x = "Time (months)") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual("Treatment", 
                      breaks = c("Pbo", "Len", "Thal"), 
                      values = trt_pal, aesthetics = c("colour", "fill")) +
  theme_multinma() +
  theme(legend.position = "top", legend.box.spacing = unit(0, "lines"))


ggsave("./figures/ndmm_mspline_survival.pdf", width = 6, height = 5, scale = 1.2)


# Hazard ratios -----------------------------------------------------------

# Calculate population-average conditional hazard ratios
get_hazard_ratio <- function(trtf,      # Treatment to get hazard ratio for
                             classf,    # Treatment class for shared EMs
                             regression = ~-1 + studyf + trtf, # Regression model fitted
                             xdat,     # Covariate data (either IPD or integration points)
                             stanfit,  # Fitted stan model object
                             log = FALSE,  # Get logHR rather than HR
                             ...) {    # Other arguments passed to tidyMCMC
  
  if (as.numeric(trtf) == 1) stop("Calculate HR for trt > 1 compared to 1.")
  
  # Set up prediction data
  xdat <- mutate(xdat, 
                 trtf = .env$trtf, classf = .env$classf,
                 # Dummy study variable
                 studyf = 1)
  
  # Construct design matrix
  regression <- update(regression, ~.+1)
  X <- model.matrix(regression, data = xdat)
  X <- X[,-(1:2)]
  
  # Pull posterior samples of model parameters
  post_beta <- as.matrix(stanfit, pars = c("gamma", "beta"))
  
  # Only need to use treatment effect and EM interaction columns
  col_select <- c(1:(nlevels(trtf)-1), 
                  grep("(\\:(trtf|classf))|((trtf|classf).+\\:)", colnames(X)))
  
  # Evaluate the log HR
  logHR <- tcrossprod(post_beta[, col_select], X[, col_select])
  HR <- exp(logHR)
  
  # Integrate/average over the covariate distribution
  logHR <- apply(logHR, 1, mean)
  HR <- apply(HR, 1, mean)
  
  {if (tolower(log) == "both") tibble(HR = HR, logHR = logHR)
    else if (log) tibble(logHR = logHR)
    else tibble(HR = HR)} %>% 
    
    # Produce posterior statistics
    tidyMCMC(conf.int = TRUE, ...)
}

# Get population average hazard ratios in every study
mspline_hrs <- crossing(studyf = unique(surv_all$studyf), trtf = unique(surv_all$trtf)) %>% 
  filter(trtf != "Pbo") %>% 
  # Add treatment classes
  mutate(classf = factor(if_else(trtf == "Pbo", "Pbo", "Active"), levels = c("Pbo", "Active"))) %>% 
  group_by(studyf, trtf) %>% 
  group_modify(
    ~get_hazard_ratio(trtf = .$trtf, classf = .$classf,
                      xdat = filter(stan_xdat, studyf == .$studyf),
                      regression = ~-1 + studyf + trtf + 
                        age + iss_stage3 + response_cr_vgpr + male + 
                        (age + iss_stage3 + response_cr_vgpr + male):classf,
                      stanfit = mspline_MLNMR_stan,
                      log = TRUE),
    .keep = TRUE
  )

mspline_hrs

