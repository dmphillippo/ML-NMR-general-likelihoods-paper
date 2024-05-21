################################################################################
# Simulated parametric survival example calling Stan directly
#   - David Phillippo, University of Bristol
################################################################################

library(dplyr)    # For data manipulation and plotting
library(tidyr)
library(purrr)
library(ggplot2)

library(simsurv)  # For simulating survival data
library(survival) # For Kaplan-Meier plots
library(flexsurv) # For creating model predictions

library(multinma)  # Used here for setting up integration points

library(rstan)
library(broom)
library(broom.mixed)
library(loo)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Colour palette
trt_pal <- c(A = "#113259", B = "#865800", C = "#55A480")



# Simulated example data --------------------------------------------------

# Set seed
set.seed(54321)

# Simulate individuals for two trials
nAB <- 500
nAC <- 400

covdat_AB <- tibble(study = 1,
                    id = 1:nAB,
                    trtB = rep(c(0, 1), each = nAB/2),
                    trtC = 0,
                    x1 = rnorm(nAB, mean = 0, sd = 0.5),
                    x2 = rgamma(nAB, shape = 4, rate = 2),
                    x3 = if_else(runif(nAB) > 0.8, 1, 0)) %>% 
  # Add EM interactions
  bind_cols(model.matrix(~(trtB + trtC):(x1 + x2 + x3), data = .)[,-1] %>% as_tibble())

covdat_AC <- tibble(study = 2,
                    id = 1:nAC,
                    trtB = 0,
                    trtC = rep(c(0, 1), each = nAC/2),
                    x1 = rnorm(nAC, mean = 1, sd = 0.4),
                    x2 = rgamma(nAC, shape = 6, rate = 2),
                    x3 = if_else(runif(nAC) > 0.3, 1, 0)) %>% 
  # Add EM interactions
  bind_cols(model.matrix(~(trtB + trtC):(x1 + x2 + x3), data = .)[,-1] %>% as_tibble())

covdat <- bind_rows(covdat_AB, covdat_AC) %>%
  mutate(trtn = case_when(trtB == 1 ~ 2L, trtC == 1 ~ 3L, TRUE ~ 1L),
         studyf = factor(study, labels = c("AB", "AC")),
         trtf = factor(trtn, labels = c("A", "B", "C")),
         trtclass = recode_factor(trtf, A = 1, B = 2, C = 2))

# True covariate mean values
true_covmeans <- tribble(
  ~study, ~x1, ~x2, ~x3,
  1L, 0, 4/2, 1-0.8,
  2L, 1, 6/2, 1-0.3 
)

true_xbar <- t(crossprod(as.matrix(true_covmeans[, -1]), c(nAB, nAC)) / sum(nAB, nAC))

# Set coefficients
betas <- c(trtB = -1.2, trtC = -0.5, 
           x1 = 0.1, x2 = 0.05, x3 = -0.25,
           `trtB:x1` = -0.2, `trtB:x2` = -0.2, `trtB:x3` = -0.1,
           `trtC:x1` = -0.2, `trtC:x2` = -0.2, `trtC:x3` = -0.1
)

# Set censoring rate
cens_rate <- 0.1

# Set shapes and scales
weib_sim_pars <- c(scaleAB = 6.2, shapeAB = 0.8, scaleAC = 5.8, shapeAC = 1.2)

# Simulate a Weibull outcome
# lambda is scale, gamma is shape
survdat_AB <- simsurv("weibull", 
                      lambdas = weib_sim_pars["scaleAB"], gammas = weib_sim_pars["shapeAB"],
                      x = covdat_AB, betas = betas, maxt = 1) %>%
  mutate(status = if_else(runif(nAB) <= cens_rate, 0L, status))

survdat_AC <- simsurv("weibull", 
                      lambdas = weib_sim_pars["scaleAC"], gammas = weib_sim_pars["shapeAC"],
                      x = covdat_AC, betas = betas, maxt = 1) %>%
  mutate(status = if_else(runif(nAC) <= cens_rate, 0L, status))

# Datasets for ML-NMR and IPD NMA
survdat <- bind_rows(full_join(covdat_AB, survdat_AB),
                     full_join(covdat_AC, survdat_AC)) %>%
  mutate(studyf = factor(study, labels = c("AB", "AC")),
         trtn = case_when(trtB == 1 ~ 2L, trtC == 1 ~ 3L, TRUE ~ 1L),
         trtf = factor(trtn, labels = c("A", "B", "C")),
         trtclass = recode_factor(trtf, A = 1, B = 2, C = 2),
         .after = 1) %>% 
  select(-trtB, -trtC)

survdat_AB <- filter(survdat, studyf == "AB")
survdat_AC <- filter(survdat, studyf == "AC")

survdat_AC_agd <- select(survdat_AC, -x1, -x2, -x3)

covdat_AC_agd <- filter(covdat, studyf == "AC") %>% 
  group_by(studyf, trtf, trtclass) %>% 
  summarise(x1_mean = mean(x1),
            x1_sd = sd(x1),
            x2_mean = mean(x2),
            x2_sd = sd(x2),
            x3 = mean(x3),
            sample_size = n())


# Produce KM plots
kmdat <- survdat %>% 
  group_by(studyf, trtf) %>% 
  group_modify(~with(survfit(Surv(eventtime, event = status) ~ 1, data = .),
                     tibble(time, n.censor, estimate = surv, std.err, upper, lower))) %>% 
  # Add S(0) = 1
  group_modify(~add_row(., time = 0, n.censor = 0, estimate = 1, std.err = 0, upper = 1, lower = 1, .before = 0))

ggplot(kmdat, aes(x = time, y = estimate, colour = trtf)) +
  geom_step() +
  geom_point(shape = 3, data = function(x) filter(x, n.censor >= 1)) +
  facet_grid(~studyf) +
  labs(y = "Survival probability", x = "Time") +
  scale_colour_manual("Treatment", 
                      breaks = c("A", "B", "C"), 
                      values = trt_pal) +
  theme_multinma(base_family = "Source Sans Pro") +
  theme(legend.position = "top", legend.box.spacing = unit(0, "lines"))


# ML-NMR analysis ---------------------------------------------------------


# Create integration points

n_int <- 64 # Number of integration points

covdat_AC_agd <- add_integration(covdat_AC_agd,
                   # Specify marginal distributions for each covariate
                   x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                   x2 = distr(qgamma, mean = x2_mean, sd = x2_sd),
                   x3 = distr(qbern, x3),
                   # Take correlation matrix from IPD study
                   cor = cor(survdat_AB[c("x1", "x2", "x3")], method = "spearman"),
                   n_int = n_int)
  
agd_xpoints <- unnest_integration(covdat_AC_agd)

# Derive global means for centering (for sampling efficiency)
x1_gmean <- weighted.mean(c(survdat_AB$x1, covdat_AC_agd$x1_mean), 
                          c(rep_len(1, nrow(survdat_AB)), covdat_AC_agd$sample_size))
x2_gmean <- weighted.mean(c(survdat_AB$x2, covdat_AC_agd$x2_mean), 
                          c(rep_len(1, nrow(survdat_AB)), covdat_AC_agd$sample_size))


# Now put all the IPD covariates and AgD integration points into one data frame
# Note: Stan code assumes these are IPD studies first, AgD studies second
stan_xdat <- bind_rows(survdat_AB, agd_xpoints) %>% 
  select(study, studyf, trtn, trtf, trtclass, x1, x2, x3) %>% 
  # Center the continuous covariates
  mutate(x1 = x1 - x1_gmean,
         x2 = x2 - x2_gmean)

# Create the model matrix
# Note: Treatments B, C, and D are in the same class here, so we make the shared EM
# assumption and give these the same EM interactions
X_all <- model.matrix(~ -1 + studyf + trtf + x1 + x2 + x3 + (x1 + x2 + x3):trtclass,
                      data = stan_xdat)

# Note: Stan code assumes that columns are ordered as 
#  - study intercepts
#  - treatments
#  - regression terms (prognostic factors, effect modifier interactions)
colnames(X_all)


# Thin QR decomposition (for sampling efficiency)
X_all_qr <- qr(X_all)
X_all_Q <- qr.Q(X_all_qr) * sqrt(nrow(X_all) - 1)
X_all_R <- qr.R(X_all_qr)[, sort.list(X_all_qr$pivot)] / sqrt(nrow(X_all) - 1)
X_all_R_inv <- solve(X_all_R)

# Construct stan data
# Note: inputs also described in Stan model code
standat_MLNMR <- list(
  # Constants
  ns_ipd = 1,               # Number of IPD studies
  ns_agd = 1,               # Number of AgD studies
  ni_ipd = nrow(survdat_AB),  # Number of IPD individuals
  ni_agd = nrow(survdat_AC_agd),  # Number of AgD event times
  nt = 3,                   # Number of treatments
  nX = ncol(X_all),         # Number of columns of design matrix
  nint = n_int,             # Number of integration points
  
  # Survival data
  time = c(survdat_AB$eventtime, survdat_AC_agd$eventtime),  # Event times
  status = c(survdat_AB$status, survdat_AC_agd$status),   # Status (event/censor)
  
  # Auxiliary parameters (i.e. shapes)
  # For a PH/AFT model these are study ids (i.e. stratified by study)
  aux_id = c(as.numeric(survdat_AB$studyf), as.numeric(survdat_AC_agd$studyf)),    # assumes factors have same levels

  # Integration point info
  int_id =
    left_join(
      survdat_AC_agd,
      covdat_AC_agd %>% ungroup() %>% 
        mutate(ag_id = 1:nrow(.)) %>% 
        select(studyf, trtf, ag_id)
    ) %>% 
    pull(ag_id),                             # Integration ID for each AgD arm
  
  # Distribution flag:
  # 1 = Exponential PH
  # 2 = Weibull PH
  # 3 = Gompertz PH
  # 4 = Exponential AFT
  # 5 = Weibull AFT
  # 6 = log Normal AFT
  # 7 = log logistic AFT
  # 8 = Gamma AFT
  # 9 = Generalised Gamma AFT
  dist = 2,
  
  # Design matrix
  X = X_all_Q,
  
  # QR decomposition
  QR = TRUE,
  R_inv = X_all_R_inv,
  
  # Prior SDs
  prior_intercept_sd = 100,
  prior_reg_sd = 100,
  prior_trt_sd = 100,
  prior_aux_sd = 5,
  prior_aux2_sd = 5  # Ignored except for Generalised Gamma model
  )

# Run Stan model
weib_MLNMR_stan <- stan("survival_param.stan", 
                        data = standat_MLNMR,
                        pars = c("mu",          # Study intercepts
                                 "beta",        # Regression terms
                                 "gamma",       # Individual-level treatment effect
                                 "aux", "aux2", # Parameters of survival distribution
                                 "log_lik"),    # log likelihood
                        iter = 2000,
                        chains = 4)

print(weib_MLNMR_stan, pars = c("log_lik", "beta_tilde"), include = FALSE)

# Get model fit using LOO
(weib_MLNMR_loo <- loo(weib_MLNMR_stan))


# Fit Exponential model by changing the `dist` argument in the Stan input data
exp_MLNMR_stan <- stan("survival_param.stan", 
                       data = list_modify(standat_MLNMR, dist = 1),
                       pars = c("mu", "beta", "gamma", "aux", "aux2", "log_lik"),
                       iter = 2000,
                       chains = 4)

print(exp_MLNMR_stan, pars = c("log_lik", "beta_tilde"), include = FALSE)

# Get model fit using LOO
(exp_MLNMR_loo <- loo(exp_MLNMR_stan))


# Fit Gompertz model by changing the `dist` argument in the Stan input data
gomp_MLNMR_stan <- stan("survival_param.stan", 
                       data = list_modify(standat_MLNMR, dist = 3),
                       pars = c("mu", "beta", "gamma", "aux", "aux2", "log_lik"),
                       iter = 2000,
                       chains = 4)

print(gomp_MLNMR_stan, pars = c("log_lik", "beta_tilde"), include = FALSE)

# Get model fit using LOO
(gomp_MLNMR_loo <- loo(gomp_MLNMR_stan))


# Compare model fits
loo_compare(list(Weibull = weib_MLNMR_loo, 
                 Exponential = exp_MLNMR_loo,
                 Gompertz = gomp_MLNMR_loo))



# Results -----------------------------------------------------------------

# Create population-average survival curves
get_survival_curve <- function(trtf,      # Treatment to get survival curve for (factor variable)
                               classf,    # Treatment class for EMs (factor variable)
                               studyf,    # Study to get survival curve for (factor variable)
                               regression = ~-1 + studyf + trtf, # Regression model fitted
                               xdat,     # Covariate data (either IPD or integration points)
                               stanfit,  # Fitted stan model object
                               times,    # Vector of time points to evaluate Sbar(t) at
                               dist = c("exponential", "weibull", "gompertz",
                                        "exponential-aft", "weibull-aft", "lognormal", 
                                        "loglogistic", "gamma", "gengamma"),  # Survival distribution used
                               aux_id = studyf,  # Auxiliary parameter ID
                               ...) {    # Other arguments passed to tidyMCMC
  
  
  # Define survival function S(t)
  dist <- match.arg(dist)
  S <- switch(dist,
    exponential = function(t, eta, ...) {pexp(t, rate = exp(eta), lower.tail = FALSE)},
    weibull = function(t, eta, aux, ...) {flexsurv::pweibullPH(t, shape = aux, scale = exp(eta), lower.tail = FALSE)},
    gompertz = function(t, eta, aux, ...) {flexsurv::pgompertz(t, shape = aux, rate = exp(eta), lower.tail = FALSE)},
    `exponential-aft` = function(t, eta, ...) {pexp(t, rate = exp(-eta), lower.tail = FALSE)},
    `weibull-aft` = function(t, eta, aux, ...) {flexsurv::pweibullPH(t, shape = aux, scale = exp(-aux*eta), lower.tail = FALSE)},
    lognormal = function(t, eta, aux, ...) {plnorm(t, sdlog = aux, meanlog = eta, lower.tail = FALSE)},
    loglogistic = function(t, eta, aux, ...) {flexsurv::pllogis(t, shape = aux, scale = exp(eta), lower.tail = FALSE)},
    gamma = function(t, eta, aux, ...) {pgamma(t, shape = aux, rate = exp(-eta), lower.tail = FALSE)},
    gengamma = function(t, eta, aux, aux2, ...) {flexsurv::pgengamma(t, mu = eta, sigma = aux, Q = 1/sqrt(aux2), lower.tail = FALSE)})

  # Set up prediction data
  xdat <- mutate(xdat, 
                 trtf = trtf, classf = classf,
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
  
  # Pull posterior samples of auxiliary parameter(s)
  # post_aux <- drop(postdat[, sprintf("aux[%d]", studyn)])
  post_aux <- if (!dist %in% c("exponential", "exponential-aft")) drop(as.matrix(stanfit, pars = sprintf("aux[%d]", as.numeric(aux_id)))) else NA
  post_aux2 <- if (dist == "gengamma") drop(as.matrix(stanfit, pars = sprintf("aux2[%d]", as.numeric(aux_id)))) else NA
  
  # Evaluate the linear predictor eta at each posterior sample and each
  # point in xdat: eta is a matrix [iter x n_int]
  eta <- tcrossprod(post_beta, X)
  
  # For each individual point in postdat (row of eta):
  #   1. Evaluate S(t) at every point in xdat (column of eta)
  #   2. Take the mean to get Sbar(t)
  
  names(times) <- sprintf("t_%d", seq_along(times))
  
  Sbar <- pmap_dfr(list(.eta = array_tree(eta, 1),  # Rows of eta
                        .aux = post_aux, 
                        .aux2 = post_aux2),
                  function(.eta, .aux, .aux2) map(times, function(.t) mean(S(.t, eta = .eta, aux = .aux, aux2 = .aux2))))
  
  # Produce posterior statistics
  tidyMCMC(Sbar, conf.int = TRUE, ...) %>% 
    mutate(time = times) %>% 
    select(-term)
}

# Apply this function for each treatment in each study
weib_surv_pred <- crossing(studyf = factor(c("AB", "AC")), trtf = factor(c("A", "B", "C"))) %>% 
  # Add treatment classes
  mutate(classf = recode_factor(trtf, A = 1, B = 2, C = 2)) %>% 
  group_by(studyf, trtf) %>% 
  group_modify(
    ~get_survival_curve(trtf = .$trtf, 
                       classf = .$classf,
                       studyf = .$studyf, 
                       xdat = filter(stan_xdat, studyf == .$studyf) %>% select(x1:x3),
                       regression = ~ -1 + studyf + trtf + x1 + x2 + x3 + (x1 + x2 + x3):classf,
                       stanfit = weib_MLNMR_stan,
                       times = seq(0, 1, by = 0.01),
                       dist = "weibull"),
    .keep = TRUE)

# Plot estimated survival curves
ggplot(weib_surv_pred, 
       aes(x = time, y = estimate, colour = trtf, fill = trtf)) +
  # KM curves
  geom_step(size = 0.15, data = kmdat) +
  geom_point(stroke = 0.15, shape = 3, data = filter(kmdat, n.censor >= 1)) +
  # ML-NMR curves
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, linetype = 0) +
  geom_line(size = 0.4) +
  facet_grid(~studyf) +
  labs(y = "Survival probability", x = "Time") +
  scale_colour_manual("Treatment",
                      breaks = c("A", "B", "C"),
                      values = trt_pal,
                      aesthetics = c("colour", "fill")) +
  theme_multinma()



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
weib_hrs <- crossing(studyf = unique(stan_xdat$studyf), 
                     trtf = unique(stan_xdat$trtf)) %>%
  # Add treatment classes, study and treatment labels
  mutate(classf = recode_factor(trtf, A = 1, B = 2, C = 2)) %>% 
  filter(trtf != "A") %>% 
  group_by(studyf, trtf) %>% 
  group_modify(
    ~get_hazard_ratio(trtf = .$trtf, classf = .$classf,
                      xdat = filter(stan_xdat, studyf == .$studyf),
                      regression = ~ -1 + studyf + trtf + x1 + x2 + x3 + (x1 + x2 + x3):classf,
                      stanfit = weib_MLNMR_stan,
                      log = TRUE),
    .keep = TRUE
  )

weib_hrs


# Full IPD meta-regression ------------------------------------------------


# Create the model matrix
# Note: Treatments B, C, and D are in the same class here, so we make the shared EM
# assumption and give these the same EM interactions
X_all_ipd <- model.matrix(~ -1 + studyf + trtf + x1 + x2 + x3 + (x1 + x2 + x3):trtclass,
                          # Centre the covariates at the same values as the ML-NMR
                          data = mutate(survdat,
                                        x1 = x1 - x1_gmean,
                                        x2 = x2 - x2_gmean))

# Note: Stan code assumes that columns are ordered as 
#  - study intercepts
#  - treatments
#  - regression terms (prognostic factors, effect modifier interactions)
colnames(X_all_ipd)


# Thin QR decomposition (for sampling efficiency)
X_all_ipd_qr <- qr(X_all_ipd)
X_all_ipd_Q <- qr.Q(X_all_ipd_qr) * sqrt(nrow(X_all_ipd) - 1)
X_all_ipd_R <- qr.R(X_all_ipd_qr)[, sort.list(X_all_ipd_qr$pivot)] / sqrt(nrow(X_all_ipd) - 1)
X_all_ipd_R_inv <- solve(X_all_ipd_R)

# Construct stan data
# Note: inputs also described in Stan model code
standat_IPD <- list(
  # Constants
  ns_ipd = 2,               # Number of IPD studies
  ns_agd = 0,               # Number of AgD studies
  ni_ipd = nrow(survdat),  # Number of IPD individuals
  ni_agd = 0,  # Number of AgD event times
  nt = 3,                   # Number of treatments
  nX = ncol(X_all),         # Number of columns of design matrix
  nint = 1,             # Number of integration points
  
  # Survival data
  time = survdat$eventtime,  # Event times
  status = survdat$status,   # Status (event/censor)
  
  # Auxiliary parameters (i.e. shapes)
  # For a PH/AFT model these are study ids (i.e. stratified by study)
  aux_id = as.numeric(survdat$studyf),
  
  # Integration point info
  int_id = integer(),
  
  # Distribution flag:
  # 1 = Exponential PH
  # 2 = Weibull PH
  # 3 = Gompertz PH
  # 4 = Exponential AFT
  # 5 = Weibull AFT
  # 6 = log Normal AFT
  # 7 = log logistic AFT
  # 8 = Gamma AFT
  # 9 = Generalised Gamma AFT
  dist = 2,
  
  # Design matrix
  X = X_all_ipd_Q,
  
  # QR decomposition
  QR = TRUE,
  R_inv = X_all_ipd_R_inv)

# Run Stan model
weib_IPD_stan <- stan("survival_param.stan", 
                        data = standat_IPD,
                        pars = c("mu",          # Study intercepts
                                 "beta",        # Regression terms
                                 "gamma",       # Individual-level treatment effect
                                 "aux", "aux2", # Parameters of survival distribution
                                 "log_lik"),    # log likelihood
                        iter = 2000,
                        chains = 4)

print(weib_IPD_stan, pars = c("log_lik", "beta_tilde"), include = FALSE)

# Get model fit using LOO
(weib_IPD_loo <- loo(weib_IPD_stan))


# Fit Exponential model by changing the `dist` argument in the Stan input data
exp_IPD_stan <- stan("survival_param.stan", 
                       data = list_modify(standat_IPD, dist = 1),
                       pars = c("mu", "beta", "gamma", "aux", "aux2", "log_lik"),
                       iter = 2000,
                       chains = 4)

print(exp_IPD_stan, pars = c("log_lik", "beta_tilde"), include = FALSE)

# Get model fit using LOO
(exp_IPD_loo <- loo(exp_IPD_stan))


# Fit Gompertz model by changing the `dist` argument in the Stan input data
gomp_IPD_stan <- stan("survival_param.stan", 
                        data = list_modify(standat_IPD, dist = 3),
                        pars = c("mu", "beta", "gamma", "aux", "aux2", "log_lik"),
                        iter = 2000,
                        chains = 4)

print(gomp_IPD_stan, pars = c("log_lik", "beta_tilde"), include = FALSE)

# Get model fit using LOO
(gomp_IPD_loo <- loo(gomp_IPD_stan))


# Compare model fits
loo_compare(list(Weibull = weib_IPD_loo, 
                 Exponential = exp_IPD_loo,
                 Gompertz = gomp_IPD_loo))

