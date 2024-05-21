################################################################################
# Simulated parametric survival example using the multinma package
################################################################################

library(dplyr)    # For data manipulation and plotting
library(tidyr)
library(ggplot2)

library(simsurv)  # For simulating survival data
library(survival) # For Kaplan-Meier plots

library(multinma)
options(mc.cores = parallel::detectCores())  # Set up parallel processing

library(loo)  # For model comparison using LOOIC

library(xtable)

# Fonts for plots
library(showtext)
font_add_google("Source Sans Pro")
showtext_auto()

# Colour palette
trt_pal <- c(A = "#113259", B = "#865800", C = "#55A480")


# Set up output directories -----------------------------------------------

if (!dir.exists("./figures")) dir.create("./figures")
if (!dir.exists("./tables")) dir.create("./tables")


# Simulated data ----------------------------------------------------------

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
  group_by(studyf, trtf) %>% 
  summarise(x1_mean = mean(x1),
            x1_sd = sd(x1),
            x2_mean = mean(x2),
            x2_sd = sd(x2),
            x3 = mean(x3))


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

ggsave("./figures/survival_sim_curves.pdf", width = 6, height = 3, scale = 1.2)


# ML-NMR analyses ---------------------------------------------------------

# In the ML-NMR scenario, we have IPD from the AB study
head(survdat_AB)

# And AgD from the AC study which consists of i) event/censoring times
head(survdat_AC_agd)

# and ii) aggregate covariate information
covdat_AC_agd

## Set up the network ------------------------------------------------------

sim_net <- combine_network(
  set_ipd(survdat_AB, 
          trt = trtf,
          study = studyf,
          trt_class = trtclass,
          Surv = Surv(eventtime, status)),
  set_agd_surv(survdat_AC,
               trt = trtf,
               study = studyf,
               trt_class = trtclass,
               Surv = Surv(eventtime, status),
               covariates = covdat_AC_agd)
)

# Add numerical integration points over AgD covariate distribution
sim_net <- add_integration(sim_net,
                           x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                           x2 = distr(qgamma, mean = x2_mean, sd = x2_sd),
                           x3 = distr(qbern, x3))

sim_net


## Fit Weibull model -------------------------------------------------------

weib_MLNMR <- nma(sim_net,
                  regression = ~(x1 + x2 + x3)*.trt,
                  likelihood = "weibull",
                  prior_intercept = normal(0, 100),
                  prior_trt = normal(0, 100),
                  prior_reg = normal(0, 100),
                  prior_aux = half_normal(10),
                  QR = TRUE)

weib_MLNMR

# Get model fit using LOO
(weib_MLNMR_loo <- loo(weib_MLNMR))


## Fit Exponential model ---------------------------------------------------

exp_MLNMR <- nma(sim_net,
                 regression = ~(x1 + x2 + x3)*.trt,
                 likelihood = "exponential",
                 prior_intercept = normal(0, 100),
                 prior_trt = normal(0, 100),
                 prior_reg = normal(0, 100),
                 QR = TRUE) 

exp_MLNMR

# Get model fit using LOO
(exp_MLNMR_loo <- loo(exp_MLNMR))


## Fit Gompertz model ------------------------------------------------------

gomp_MLNMR <- nma(sim_net,
                  regression = ~(x1 + x2 + x3)*.trt,
                  likelihood = "gompertz",
                  prior_intercept = normal(0, 100),
                  prior_trt = normal(0, 100),
                  prior_reg = normal(0, 100),
                  prior_aux = half_normal(10),
                  QR = TRUE) 

gomp_MLNMR

# Get model fit using LOO
(gomp_MLNMR_loo <- loo(gomp_MLNMR))


## Compare model fits ------------------------------------------------------

loo_compare(list(
  Weibull = weib_MLNMR_loo, 
  Exponential = exp_MLNMR_loo, 
  Gompertz = gomp_MLNMR_loo))


# Full IPD NMA ------------------------------------------------------------

# For comparison, we fit a full IPD NMA as the "gold standard" analysis, using
# IPD from both studies


## Set up the network ------------------------------------------------------

sim_net_IPD <- set_ipd(survdat,
                       study = studyf,
                       trt = trtf,
                       trt_class = trtclass,
                       Surv = Surv(eventtime, status))


## Fit IPD NMA -------------------------------------------------------------

weib_IPD <- nma(sim_net_IPD,
                regression = ~(x1 + x2 + x3)*.trt,
                likelihood = "weibull",
                prior_intercept = normal(0, 100),
                prior_reg = normal(0, 100),
                prior_trt = normal(0, 100),
                prior_aux = half_normal(10),
                QR = TRUE)

weib_IPD

exp_IPD <- nma(sim_net_IPD,
                regression = ~(x1 + x2 + x3)*.trt,
                likelihood = "exponential",
                prior_intercept = normal(0, 100),
                prior_reg = normal(0, 100),
                prior_trt = normal(0, 100),
                QR = TRUE)

exp_IPD

gomp_IPD <- nma(sim_net_IPD,
                regression = ~(x1 + x2 + x3)*.trt,
                likelihood = "gompertz",
                prior_intercept = normal(0, 100),
                prior_reg = normal(0, 100),
                prior_trt = normal(0, 100),
                prior_aux = half_normal(10),
                QR = TRUE)

gomp_IPD

## Compare model fits ------------------------------------------------------

(weib_IPD_loo <- loo(weib_IPD))
(exp_IPD_loo <- loo(exp_IPD))
(gomp_IPD_loo <- loo(gomp_IPD))

loo_compare(list(
  Weibull = weib_IPD_loo,
  Exponential = exp_IPD_loo,
  Gompertz = gomp_IPD_loo
))

# Compare model fit of IPD model to ML-NMR
loo_compare(weib_IPD_loo, weib_MLNMR_loo)

# Plot contributions to LOOIC from each model
tibble(IPD = weib_IPD_loo$pointwise[,4], MLNMR = weib_MLNMR_loo$pointwise[,4]) %>% 
  bind_cols(survdat) %>% 
  mutate(statusf = factor(status, labels = c("Censored", "Event"))) %>% 
  ggplot(aes(x = IPD, y = MLNMR, shape = statusf, colour = statusf)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(slope = 1, colour = "grey60") +
  # Plot points for censored events on top
  geom_point(data = function(x) filter(x, status == 1)) +
  geom_point(data = function(x) filter(x, status == 0)) +
  facet_wrap(~studyf + trtf, ncol = 2) +
  scale_shape_manual("Status", values = c(Censored = 3, Event = 1)) +
  scale_colour_manual("Status", values = c(Censored = "#7B3294", Event = "#ABDDA4")) +
  coord_equal() +
  labs(x = "IPD LOOIC contributions", y = "ML-NMR LOOIC contributions") +
  #theme_minimal(base_family = "Source Sans Pro") +
  theme_multinma(base_family = "Source Sans Pro") +
  theme(legend.position = "top", legend.box.spacing = unit(0, "lines"))

ggsave("./figures/looic_contributions.pdf", width = 5, height = 6)


# Standard unadjusted IC --------------------------------------------------

# For a fair comparison with ML-NMR, perform an IPD analysis in each trial,
# adjusting for prognostic factors, and then compare reported log HRs

# For convenience, we fit these in a single run

weib_unadj_IC <- nma(sim_net_IPD,
                     regression = ~(x1 + x2 + x3):.study,
                     likelihood = "weibull",
                     prior_intercept = normal(0, 100),
                     prior_reg = normal(0, 100),
                     prior_trt = normal(0, 100),
                     prior_aux = half_normal(10),
                     QR = TRUE)

weib_unadj_IC


# Compare estimates with true values --------------------------------------

# Lookup table for parameter names
par_lookup <- tribble(
  ~stan_term, ~model_term, ~term_type, ~eqlabel,
  "beta[x1]", "x1", "Prognostic Effect", "beta[1*';'~1]",
  "beta[x2]", "x2", "Prognostic Effect", "beta[1*';'~2]",
  "beta[x3]", "x3", "Prognostic Effect", "beta[1*';'~3]",
  "beta[x1:.trtclass2]", "trtclass:x1", "EM Interaction", "beta[2*';'~1]",
  "beta[x2:.trtclass2]", "trtclass:x2", "EM Interaction", "beta[2*';'~2]",
  "beta[x3:.trtclass2]", "trtclass:x3", "EM Interaction", "beta[2*';'~3]",
  "d[B]", "trtB", "Treatment Effect", "gamma[B]",
  "d[C]", "trtC", "Treatment Effect", "gamma[C]",
  "shape[AB]", "shapeAB", "Shape", "nu[AB]",
  "shape[AC]", "shapeAC", "Shape", "nu[AC]"
)

# True values
truthdat <- 
  tibble(!!! c(betas, weib_sim_pars)) %>% 
  mutate(
    `trtclass:x1` = `trtB:x1`,
    `trtclass:x2` = `trtB:x2`,
    `trtclass:x3` = `trtB:x3`
  ) %>% 
  pivot_longer(cols = everything(), names_to = "model_term", values_to = "truth")

# Uncenter treatment effect estimates for comparison to true values
xb_MLNMR <- apply(as.array(weib_MLNMR, pars = c("beta[x1:.trtclass2]",
                                                "beta[x2:.trtclass2]", 
                                                "beta[x3:.trtclass2]")),
                  MARGIN = 1:2,
                  FUN = function(x) x %*% weib_MLNMR$xbar)
weib_MLNMR_d_uncen <- summary(sweep(as.array(weib_MLNMR, pars = "d"), 
                                    MARGIN = 1:2, STATS = xb_MLNMR, FUN = "-"))


xb_IPD <- apply(as.array(weib_IPD, pars = c("beta[x1:.trtclass2]",
                                            "beta[x2:.trtclass2]", 
                                            "beta[x3:.trtclass2]")),
                MARGIN = 1:2,
                FUN = function(x) x %*% weib_IPD$xbar)
weib_IPD_d_uncen <- summary(sweep(as.array(weib_IPD, pars = "d"), 
                                  MARGIN = 1:2, STATS = xb_IPD, FUN = "-"))


# Plot estimates of IPD NMA and ML-NMR together
weib_est_dat <- 
  bind_rows(
    summary(weib_IPD, pars = c("shape", "beta")) %>%
      as_tibble() %>% 
      add_row(as_tibble(weib_IPD_d_uncen)) %>% 
      mutate(model = "IPD NMA"),
    summary(weib_MLNMR, pars = c("shape", "beta")) %>%
      as_tibble() %>% 
      add_row(as_tibble(weib_MLNMR_d_uncen)) %>% 
      mutate(model = "ML-NMR")
  ) %>% 
  left_join(par_lookup, by = c(parameter = "stan_term")) %>% 
  left_join(truthdat)

ggplot(weib_est_dat,
       aes(y = eqlabel, x = mean, xmin = `2.5%`, xmax = `97.5%`, 
           colour = model, shape = model, linetype = model, fill = model)) +
  geom_pointrange(position = position_dodge(0.6)) +
  geom_point(aes(x = truth), shape = 23, colour = "black", fill = "grey80") +
  scale_colour_manual("Model", 
                      limits = c("ML-NMR", "IPD NMA", "Truth"),
                      values = c(`ML-NMR` = "#113259", `IPD NMA` = "#55A480", Truth = "black"), 
                      drop = FALSE) +
  scale_shape_manual("Model", 
                      limits = c("ML-NMR", "IPD NMA", "Truth"),
                      values = c(`ML-NMR` = 16, `IPD NMA` = 15, Truth = 23), 
                      drop = FALSE) +
  scale_linetype_manual("Model",
                        limits = c("ML-NMR", "IPD NMA", "Truth"),
                        values = c(`ML-NMR` = 1, `IPD NMA` = 1, Truth = 0), 
                        drop = FALSE) +
  scale_fill_manual("Model", 
                    limits = c("ML-NMR", "IPD NMA", "Truth"),
                    values = c(`ML-NMR` = "#113259", `IPD NMA` = "#55A480", Truth = "grey80"), 
                    drop = FALSE) +
  scale_y_discrete(labels = scales::label_parse()) +
  facet_grid(rows = vars(term_type), scales = "free", space = "free") +
  labs(x = "Estimate", y = "") +
  theme_multinma(base_family = "Source Sans Pro") +
  theme(legend.position = "top", legend.box.spacing = unit(0, "lines"))
  
ggsave("./figures/parameters_vs_truth.pdf", width = 5, height = 6)

# Estimates show good agreement, recover the true values well.
# Uncertainty does not seem to differ greatly between IPD NMA and ML-NMR, there 
# is a slight increase in uncertainty for ML-NMR.


# Function for easy formatting of estimates
tabfmt <- function(x, digits = 2, format = "f",
                   trunc.lo = 0, #10^-digits,
                   na.char = "-",
                   prefix = "$", suffix = "$", ...) {
  out <- formatC(x, digits = digits, format = format, ...)
  out[x > 0 & x < trunc.lo] <- paste0("<", trunc.lo)
  out[x < 0 & x > trunc.lo] <- paste0(">-", trunc.lo)
  out[is.na(x) | is.nan(x)] <- na.char
  if (!is.null(prefix) || !is.null(suffix)) out <- paste0(prefix, out, suffix)
  return(out)
}

# Function for only printing (sub) grouping variables on the first row
subfmt <- function(subgroup, by = rep_len(1, length(subgroup))) {
  #if (mode(subgroup) != "character") stop("Only character subgroup supported.")
  
  out <- as.character(subgroup)
  
  for (i in 2:length(subgroup)) {
    if (subgroup[i] == subgroup[i - 1] && by[i] == by[i - 1]) {
      out[i] <- ""
    }
  }
  
  return(out)
}

# Create table of estimates
weib_est_dat %>% 
  mutate(term_typef = ordered(term_type, 
                              levels = c("Treatment Effect", "Prognostic Effect",
                                         "EM Interaction", "Shape")),
         model_termf = ordered(model_term,
                               levels = c("trtB", "trtC",
                                          "x1", "x2", "x3",
                                          "trtclass:x1", "trtclass:x2", "trtclass:x3",
                                          "shapeAB", "shapeAC"),
                               labels = c("$\\gamma_B$", "$\\gamma_C$",
                                          "$\\beta_{1;1}$", "$\\beta_{1;2}$", "$\\beta_{1;3}$",
                                          "$\\beta_{2;1}$", "$\\beta_{2;2}$", "$\\beta_{2;3}$",
                                          "$\\nu_{AB}$", "$\\nu_{AC}$")),
         est_fmt = tabfmt(mean),
         ci_fmt = paste0("(", tabfmt(`2.5%`), ", ", tabfmt(`97.5%`), ")"),
         truth_fmt = tabfmt(truth)) %>% 
  pivot_longer(names_to = "stat", values_to = "stat_fmt", cols = c(est_fmt, ci_fmt)) %>% 
  select(model:stat_fmt) %>% 
  group_by(stat) %>% 
  pivot_wider(names_from = model, values_from = stat_fmt) %>% 
  arrange(term_typef, model_termf, desc(stat)) %>% 
  ungroup() %>% 
  transmute(term_typec = subfmt(as.character(term_typef)), 
            model_termc = subfmt(as.character(model_termf), by = term_typef), 
            Truth = subfmt(truth_fmt, by = model_termf),
            `IPD NMA`,
            `ML-NMR`) %>%
  
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/weib_estimates.tex",
        sanitize.text.function = function(x) {x},
        include.colnames = FALSE,
        include.rownames = FALSE,
        booktabs = TRUE,
        only.contents = TRUE,
        add.to.row = lst(pos = as.list(seq(2, nrow(.) - 1, by = 2)),
                         command = rep_len("[0.5ex]", length(pos))),
        hline.after = NULL,
        comment = FALSE)}
  


# Population-average survival curves --------------------------------------

# Survival estimates can be produced using the predict() function, and plotted
# using the plot() method

plot(predict(weib_MLNMR, type = "survival", times = seq(0, 1, by = 0.02))) +
  # overlay KM data
  geom_km(sim_net) +
  scale_colour_manual("Treatment",
                      breaks = c("A", "B", "C"),
                      values = trt_pal,
                      aesthetics = c("colour", "fill")) +
  theme_multinma(base_family = "Source Sans Pro") +
  theme(legend.position = "top", legend.box.spacing = unit(0, "lines"))
  

ggsave("./figures/survival_est_curves.pdf", width = 6, height = 3, scale = 1.2)


# Population-average hazard ratios ----------------------------------------

# Population-average (conditional) hazard ratios can be calculated using the
# relative_effects() function

(weib_MLNMR_hrs <- relative_effects(weib_MLNMR, all_contrasts= TRUE))
(weib_IPD_hrs <- relative_effects(weib_IPD, all_contrasts= TRUE))

# Create table
weib_hrs <- 
  bind_rows(
    as_tibble(weib_MLNMR_hrs) %>% mutate(model = "ML-NMR"),
    as_tibble(weib_IPD_hrs) %>% mutate(model = "IPD NMA"),
  ) %>% 
  transmute(studyf = .study,
            label = paste0(.trtb, " vs. ", .trta),
            estimate = mean, `2.5%`, `97.5%`,
            model) %>% 
  
  # Add in true values
  bind_rows(
    true_covmeans %>% mutate(studyf = factor(study, labels = c("AB", "AC"))) %>% group_by(studyf) %>% 
      do(
        tribble(
          ~label, ~estimate,
          "B vs. A", betas["trtB"] + c(.$x1, .$x2, .$x3) %*% betas[c("trtB:x1", "trtB:x2", "trtB:x3")],
          "C vs. A", betas["trtC"] + c(.$x1, .$x2, .$x3) %*% betas[c("trtC:x1", "trtC:x2", "trtC:x3")],
          "C vs. B", betas["trtC"] - betas["trtB"]
        )
      ) %>% 
      mutate(model = "Truth")
  ) 

  # Create table
  weib_hrs %>% 
  mutate(modelf = ordered(model, levels = c("Truth", "ML-NMR", "IPD NMA")),# "Standard IC")),
         est_fmt = tabfmt(estimate),
         ci_fmt = paste0("(", tabfmt(`2.5%`), ", ", tabfmt(`97.5%`), ")")) %>% 
  pivot_longer(names_to = "stat", values_to = "stat_fmt", cols = c(est_fmt, ci_fmt)) %>% 
  select(studyf, modelf, label, stat, stat_fmt) %>% 
  group_by(stat) %>% 
  pivot_wider(names_from = label, values_from = stat_fmt) %>% 
  arrange(studyf, modelf, desc(stat)) %>% 
  ungroup() %>% 
  filter(!(modelf == "Truth" & stat == "ci_fmt")) %>% 
  transmute(studyc = subfmt(studyf),
            modelc = subfmt(modelf, by = studyf), 
            `B vs. A`, `C vs. A`, `C vs. B`) %>% 
  
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/weib_log_hazard_ratios.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = FALSE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         add.to.row = lst(pos = list(1, 3, 5, 6, 8), #list(1, 3, 5, 7, 8, 10, 12),
                          command = rep_len("[0.5ex]", length(pos))),
         hline.after = NULL,
         comment = FALSE)}


# Population-average RMST -------------------------------------------------

# For a fair comparison between the non-population-adjusted IC and the
# ML-NMR/IPD NMA, we need to compare compatible quantities. Log HRs are not
# comparable, since they are conditioned in different ways between the
# approaches. Instead, let's compare the RMST on each treatment under each
# approach.
  
(weib_MLNMR_rmst <- predict(weib_MLNMR, type = "rmst", time = 1))
(weib_IPD_rmst <- predict(weib_IPD, type = "rmst", time = 1))
(weib_unadj_IC_rmst <- predict(weib_unadj_IC, type = "rmst", time = 1))
  

  
# Create table
weib_rmst <- 
  bind_rows(
    as_tibble(weib_MLNMR_rmst) %>% mutate(model = "ML-NMR"),
    as_tibble(weib_IPD_rmst) %>% mutate(model = "IPD NMA"),
    as_tibble(weib_unadj_IC_rmst) %>% mutate(model = "Standard IC")
  ) %>% 
  transmute(studyf = .study,
            trtf = .trt,
            estimate = mean, `2.5%`, `97.5%`,
            model)

weib_rmst %>% 
  mutate(modelf = ordered(model, levels = c("ML-NMR", "IPD NMA", "Standard IC")),
         est_fmt = tabfmt(estimate),
         ci_fmt = paste0("(", tabfmt(`2.5%`), ", ", tabfmt(`97.5%`), ")")) %>% 
  pivot_longer(names_to = "stat", values_to = "stat_fmt", cols = c(est_fmt, ci_fmt)) %>% 
  select(studyf, modelf, trtf, stat, stat_fmt) %>% 
  group_by(stat) %>% 
  pivot_wider(names_from = trtf, values_from = stat_fmt) %>% 
  arrange(studyf, modelf, desc(stat)) %>% 
  ungroup() %>% 
  transmute(studyc = subfmt(studyf),
            modelf = subfmt(modelf, by = studyf), 
            A, B, C) %>% 
  
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/weib_rmst.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = FALSE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         add.to.row = lst(pos = list(2, 4, 6, 8, 10),
                          command = rep_len("[0.5ex]", length(pos))),
         hline.after = NULL,
         comment = FALSE)}
  

# Fit statistics under all models -----------------------------------------

tibble(
  surv = rep(c("Exponential", "Gompertz", "Weibull"), times = 2),
  mod = rep(c("ML-NMR", "IPD NMA"), each = 3),
  mod_surv = paste(mod, surv),
  loo = list(exp_MLNMR_loo, gomp_MLNMR_loo, weib_MLNMR_loo, 
             exp_IPD_loo, gomp_IPD_loo, weib_IPD_loo)
  ) %>% 
  rowwise() %>% 
  mutate(loo_stats = list(as.data.frame(loo$estimates) %>% 
                            tibble::rownames_to_column(var = "stat") %>% 
                            as_tibble()),
         loo_compare = case_when(mod == "ML-NMR" ~ list(loo_compare(loo, weib_MLNMR_loo)),
                                 mod == "IPD NMA" ~ list(loo_compare(loo, weib_IPD_loo))),
         loo_all = list(add_row(loo_stats, 
                                stat = "elpd_diff",
                                Estimate = loo_compare[2, "elpd_diff"],
                                SE = loo_compare[2, "se_diff"]))) %>% 
  unnest(loo_all) %>% 
  mutate(est_fmt = tabfmt(Estimate, digits = 1), 
         se_fmt = tabfmt(SE, prefix = "($", suffix = "$)", digits = 1),
         cells_fmt = if_else(surv == "Weibull" & stat == "elpd_diff", 
                             " & ", paste(est_fmt, se_fmt, sep = " & "))) %>%
  select(mod_surv, stat, cells_fmt) %>% 
  pivot_wider(names_from = mod_surv, values_from = cells_fmt) %>% 
  mutate(stat = ordered(stat, 
                        levels = c("looic", "elpd_loo", "p_loo", "elpd_diff"),
                        labels = c("LOOIC", "ELPD", "$p_\\mathrm{LOO}$", "ELPD difference"))) %>% 
  arrange(stat) %>% 
  
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/looic_comparison.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = FALSE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         hline.after = NULL,
         comment = FALSE)}

