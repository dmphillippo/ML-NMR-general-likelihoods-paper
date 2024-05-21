################################################################################
# Analysis of Newly Diagnosed Multiple Myeloma PFS using flexible M-splines on 
# the baseline hazard with Multilevel Network Meta-Regression
################################################################################

# For data manipulation and plotting
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# For KM plots
library(survival)

# multinma for fitting ML-NMR models
library(multinma)
options(mc.cores = parallel::detectCores())

# For model comparison using LOOIC
library(loo)

# For tables
library(xtable)

# Fonts for plots
library(showtext)
font_add_google("Source Sans Pro")
showtext_auto()


# Colour palettes
trt_pal <- c(Pbo = "#113259", Len = "#865800", Thal = "#55A480")


# Set up output directories -----------------------------------------------

if (!dir.exists("./figures")) dir.create("./figures")
if (!dir.exists("./tables")) dir.create("./tables")


# Synthetic data ----------------------------------------------------------

# The NDMM data are available in the multinma package. For details of how
# these are constructed, see rstan_survival/rstan_ndmm_survival.R

# IPD - (synthetic) individual event/censoring times and covariates
head(ndmm_ipd)

# AgD - reconstructed event/censoring times from published KM curves
head(ndmm_agd)

# AgD - summary covariate distributions
ndmm_agd_covs


# Produce KM plots
kmdat <- bind_rows(ndmm_ipd, ndmm_agd) %>% 
  group_by(studyf, trtf) %>% 
  group_modify(~with(survfit(Surv(eventtime, event = status) ~ 1, data = .),
                     tibble(time, n.censor, estimate = surv, std.err, upper, lower))) %>% 
  # Add S(0) = 1
  group_modify(~add_row(., time = 0, n.censor = 0, estimate = 1, std.err = 0, upper = 1, lower = 1, .before = 0)) %>% 
  arrange(studyf, trtf, time)

ggplot(kmdat, aes(x = time, y = estimate, colour = trtf)) +
  geom_step() +
  geom_point(shape = 3, data = function(x) filter(x, n.censor >= 1)) +
  facet_wrap(~studyf) +
  labs(y = "Survival probability", x = "Time") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual("Treatment", 
                      breaks = c("Pbo", "Len", "Thal"), 
                      values = trt_pal) +
  theme_multinma(base_family = "Source Sans Pro") +
  theme(legend.position = "top", legend.box.spacing = unit(0, "lines"))

ggsave("./figures/ndmm_survival_sim_curves.pdf", width = 6, height = 5, scale = 1.2)



# Function for easy formatting of estimates in tables
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


# Table of baseline covariates
bind_rows(
  summarise(ndmm_ipd, 
            N = n(),
            age_mean = mean(age), age_sd = sd(age),
            iss_stage3 = mean(iss_stage3),
            response_cr_vgpr = mean(response_cr_vgpr),
            male = mean(male),
            .by = c(studyf, trtf)),
  transmute(ndmm_agd_covs,
            studyf, trtf,
            N = sample_size,
            age_mean, age_sd, iss_stage3, response_cr_vgpr, male)
) %>% 
  arrange(studyf, trtf) %>% 
  transmute(
    studyf,
    `Study / Treatment` = recode_factor(trtf, Pbo = "\\quad Placebo", Len = "\\quad Lenalidomide", Thal = "\\quad Thalidomide"), 
    `Sample size` = tabfmt(N, 0),
    `Age (years)` = paste(tabfmt(age_mean), tabfmt(age_sd, prefix = "$(", suffix = ")$")),
    `ISS Stage III (\\%)` = tabfmt(iss_stage3 * 100),
    `Response CR/VGPR (\\%)` = tabfmt(response_cr_vgpr * 100),
    `Male (\\%)` = tabfmt(male * 100),
  ) %>% 
  group_by(studyf) %>% 
  group_modify(~add_row(.x, `Study / Treatment` = .y$studyf, .before = 1)) %>% 
  ungroup() %>% 
  select(-studyf) %>% 
  
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/ndmm_baseline_covariates.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = TRUE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         add.to.row = lst(pos = list(3, 6, 9, 12, 15),
                          command = rep_len("[0.5ex]", length(pos))),
         hline.after = 0,
         comment = FALSE)}


# Set up network ----------------------------------------------------------

# Create treatment class variable
ndmm_ipd$trtclass <- case_match(ndmm_ipd$trtf, 
                                "Pbo" ~ "Placebo",
                                c("Len", "Thal") ~ "Active")

ndmm_agd$trtclass <- case_match(ndmm_agd$trtf, 
                                "Pbo" ~ "Placebo",
                                c("Len", "Thal") ~ "Active")

# Define network
ndmm_net <- combine_network(
  set_ipd(ndmm_ipd,
          study = studyf,
          trt = trtf,
          trt_class = trtclass,
          Surv = Surv(eventtime, status)),
  set_agd_surv(ndmm_agd,
               study = studyf,
               trt = trtf,
               trt_class = trtclass,
               Surv = Surv(eventtime, status),
               covariates = ndmm_agd_covs)
)

# Add integration points
ndmm_net <- add_integration(ndmm_net,
                            age = distr(qgamma, mean = age_mean, sd = age_sd),
                            iss_stage3 = distr(qbern, iss_stage3),
                            response_cr_vgpr = distr(qbern, response_cr_vgpr),
                            male = distr(qbern, male))

ndmm_net

# Basic network plot
plot(ndmm_net,
     weight_nodes = TRUE,
     weight_edges = TRUE)

# Customise for paper
plot(ndmm_net, 
    weight_nodes = TRUE, 
    weight_edges = TRUE,
    # Nudge treatment labels away from nodes
    nudge = 0.1,
    # Manual layout
    layout = data.frame(x = c(0, -1, 1),
                        y = c(-0.5, 0, 0))) +
  # Edge numbering
  ggraph::geom_edge_fan(aes(label = .nstudy), colour = NA, strength = 2.25) +
  # Legends and margins
  guides(edge_width = "none", size = "none", 
         edge_colour = guide_legend(override.aes = list(edge_width = 2))) +
  theme(legend.position = "top")

ggsave("./figures/ndmm_network.pdf", width = 3, height = 1.5, scale = 2)


# Fit ML-NMR M-spline models ----------------------------------------------

# We fit ML-NMR models with cubic M-splines on the baseline hazard

ndmm_fit <- nma(ndmm_net, 
                regression = ~(age + iss_stage3 + response_cr_vgpr + male)*.trt,
                likelihood = "mspline",
                prior_intercept = normal(0, 100),
                prior_trt = normal(0, 100),
                prior_reg = normal(0, 100),
                prior_aux = half_normal(1),
                QR = TRUE,
                # Reduce max_treedepth to speed up warm-up (default 10)
                # Increase if you get max_treedepth warnings
                control = list(max_treedepth = 7))
ndmm_fit

saveRDS(ndmm_fit, "./ndmm_fit.rds")
# ndmm_fit <- readRDS("./ndmm_fit.rds")


ndmm_fit_10kt <- nma(ndmm_net, 
                     regression = ~(age + iss_stage3 + response_cr_vgpr + male)*.trt,
                     likelihood = "mspline",
                     n_knots = 10,
                     prior_intercept = normal(0, 100),
                     prior_trt = normal(0, 100),
                     prior_reg = normal(0, 100),
                     prior_aux = half_normal(1),
                     QR = TRUE,
                     # Reduce max_treedepth to speed up warm-up (default 10)
                     # Increase if you get max_treedepth warnings
                     control = list(max_treedepth = 7))

ndmm_fit_10kt

saveRDS(ndmm_fit_10kt, "./ndmm_fit_10kt.rds")
# ndmm_fit_10kt <- readRDS("./ndmm_fit_10kt.rds")


# Compare model fit -------------------------------------------------------

# Compare LOOIC
(ndmm_fit_loo <- loo(ndmm_fit))
(ndmm_fit_10kt_loo <- loo(ndmm_fit_10kt))

loo_compare(list("7 knots" = ndmm_fit_loo,
                 "10 knots" = ndmm_fit_10kt_loo))

# Increasing the number of knots does not improve model fit. The model with 7
# internal knots is adequate.


# Check that no single study has a better fit with a larger number of knots
# (could be washed out by the other studies)

studies_all <- c(ndmm_ipd$study, ndmm_agd$study)
looic_by_study <- 
  cbind(
  `7 knots` = by(ndmm_fit_loo$pointwise[, "looic"], studies_all, sum),
  `10 knots` = by(ndmm_fit_10kt_loo$pointwise[, "looic"], studies_all, sum)
)
looic_by_study

# No, 7 knots is sufficient for all studies


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

# Produce latex table
tibble(
  label = c("7 knots", "10 knots"),
  loo = list(ndmm_fit_loo, ndmm_fit_10kt_loo),
) %>% 
  rowwise() %>% 
  mutate(loo_stats = list(as.data.frame(loo$estimates) %>% 
                            tibble::rownames_to_column(var = "stat") %>% 
                            as_tibble()),
         loo_compare = list(loo_compare(loo, ndmm_fit_loo)),
         loo_all = list(add_row(loo_stats, 
                                stat = "elpd_diff",
                                Estimate = loo_compare[2, "elpd_diff"],
                                SE = loo_compare[2, "se_diff"]))) %>% 
  unnest(loo_all) %>% 
  mutate(est_fmt = tabfmt(Estimate, digits = 1), 
         se_fmt = tabfmt(SE, prefix = "($", suffix = "$)", digits = 1),
         cells_fmt = if_else(label == "7 knots" & stat == "elpd_diff", 
                             " & ", paste(est_fmt, se_fmt, sep = " & "))) %>%
  select(label, stat, cells_fmt) %>% 
  pivot_wider(names_from = label, values_from = cells_fmt) %>% 
  mutate(stat = ordered(stat, 
                        levels = c("looic", "elpd_loo", "p_loo", "elpd_diff"),
                        labels = c("LOOIC", "ELPD", "$p_\\mathrm{LOO}$", "ELPD difference"))) %>% 
  arrange(stat) %>% 
  
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/ndmm_looic_comparison.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = FALSE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         hline.after = NULL,
         comment = FALSE)}


colnames(looic_by_study) <- c("7 knots", "10 knots")
looic_by_study %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Study") %>% 
  mutate_at(vars(-Study), tabfmt, digits = 1)%>% 
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/ndmm_looic_by_study.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = TRUE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         hline.after = 0,
         comment = FALSE)}
  

# Relax PH assumption -----------------------------------------------------

# We can relax the proportional hazards assumption by allowing the spline
# coefficients to vary between treatment arms within each study, using the
# `aux_by` argument to nma().

ndmm_fit_nph <- nma(ndmm_net, 
                    regression = ~(age + iss_stage3 + response_cr_vgpr + male)*.trt,
                    likelihood = "mspline",
                    prior_intercept = normal(0, 100),
                    prior_trt = normal(0, 100),
                    prior_reg = normal(0, 100),
                    prior_aux = half_normal(1),
                    aux_by = c(.study, .trt),
                    QR = TRUE,
                    # Reduce max_treedepth to speed up warm-up (default 10)
                    # Increase if you get max_treedepth warnings
                    control = list(max_treedepth = 7))
ndmm_fit_nph

saveRDS(ndmm_fit_nph, "./ndmm_fit_nph.rds")
# ndmm_fit_nph <- readRDS("./ndmm_fit_nph.rds")

# Calculate LOOIC
(ndmm_fit_nph_loo <- loo(ndmm_fit_nph))

# Compare to PH model
loo_compare(ndmm_fit_loo, ndmm_fit_nph_loo)

# Overall fit for the PH model is better (but not substantially)

# Check that no single study has a better fit with the non-PH model
# (could be washed out by the other studies)

looic_by_study_nph <- 
  cbind(
    PH = by(ndmm_fit_loo$pointwise[, "looic"], studies_all, sum),
    `non-PH` = by(ndmm_fit_nph_loo$pointwise[, "looic"], studies_all, sum)
  )
looic_by_study_nph

# LOOIC is similar or lower for the PH model in all studies.

# Overall, little evidence to suggest that the PH assumption is invalid here
# based on LOOIC alone


# Produce latex tables
tibble(
  label = c("PH", "Non-PH"),
  loo = list(ndmm_fit_loo, ndmm_fit_nph_loo)
) %>% 
  rowwise() %>% 
  mutate(loo_stats = list(as.data.frame(loo$estimates) %>% 
                            tibble::rownames_to_column(var = "stat") %>% 
                            as_tibble()),
         loo_compare = list(loo_compare(loo, ndmm_fit_loo)),
         loo_all = list(add_row(loo_stats, 
                                stat = "elpd_diff",
                                Estimate = loo_compare[2, "elpd_diff"],
                                SE = loo_compare[2, "se_diff"]))) %>% 
  unnest(loo_all) %>% 
  mutate(est_fmt = tabfmt(Estimate, digits = 1), 
         se_fmt = tabfmt(SE, prefix = "($", suffix = "$)", digits = 1),
         cells_fmt = if_else(label == "PH" & stat == "elpd_diff", 
                             " & ", paste(est_fmt, se_fmt, sep = " & "))) %>%
  select(label, stat, cells_fmt) %>% 
  pivot_wider(names_from = label, values_from = cells_fmt) %>% 
  mutate(stat = ordered(stat, 
                        levels = c("looic", "elpd_loo", "p_loo", "elpd_diff"),
                        labels = c("LOOIC", "ELPD", "$p_\\mathrm{LOO}$", "ELPD difference"))) %>% 
  arrange(stat) %>% 
  
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/ndmm_looic_comparison_nph.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = FALSE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         hline.after = NULL,
         comment = FALSE)}

colnames(looic_by_study_nph) <- c("PH", "Non-PH")
looic_by_study_nph %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Study") %>% 
  mutate_at(vars(-Study), tabfmt, digits = 1)%>% 
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/ndmm_looic_by_study_nph.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = TRUE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         hline.after = 0,
         comment = FALSE)}


# Compare to unadjusted NMA -----------------------------------------------

ndmm_fit_nma <- nma(ndmm_net, 
                    likelihood = "mspline",
                    prior_intercept = normal(0, 100),
                    prior_trt = normal(0, 100),
                    prior_aux = half_normal(1))

ndmm_fit_nma

ndmm_fit_nma_nph <- nma(ndmm_net, 
                        likelihood = "mspline",
                        prior_intercept = normal(0, 100),
                        prior_trt = normal(0, 100),
                        prior_aux = half_normal(1),
                        aux_by = c(.study, .trt))

ndmm_fit_nma_nph

# Compare overall model fit
(ndmm_fit_nma_loo <- loo(ndmm_fit_nma))
(ndmm_fit_nma_nph_loo <- loo(ndmm_fit_nma_nph))
loo_compare(ndmm_fit_nma_loo, ndmm_fit_nma_nph_loo)

# Compare model fit by study
looic_by_study_nma <- 
  cbind(
    PH = by(ndmm_fit_nma_loo$pointwise[, "looic"], studies_all, sum),
    `non-PH` = by(ndmm_fit_nma_nph_loo$pointwise[, "looic"], studies_all, sum)
  )
looic_by_study_nma

# PH model is reasonable for every study except Jackson2019


# Produce latex tables
tibble(
  label = c("PH", "Non-PH"),
  loo = list(ndmm_fit_nma_loo, ndmm_fit_nma_nph_loo)
) %>% 
  rowwise() %>% 
  mutate(loo_stats = list(as.data.frame(loo$estimates) %>% 
                            tibble::rownames_to_column(var = "stat") %>% 
                            as_tibble()),
         loo_compare = list(loo_compare(loo, ndmm_fit_nma_loo)),
         loo_all = list(add_row(loo_stats, 
                                stat = "elpd_diff",
                                Estimate = loo_compare[2, "elpd_diff"],
                                SE = loo_compare[2, "se_diff"]))) %>% 
  unnest(loo_all) %>% 
  mutate(est_fmt = tabfmt(Estimate, digits = 1), 
         se_fmt = tabfmt(SE, prefix = "($", suffix = "$)", digits = 1),
         cells_fmt = if_else(label == "PH" & stat == "elpd_diff", 
                             " & ", paste(est_fmt, se_fmt, sep = " & "))) %>%
  select(label, stat, cells_fmt) %>% 
  pivot_wider(names_from = label, values_from = cells_fmt) %>% 
  mutate(stat = ordered(stat, 
                        levels = c("looic", "elpd_loo", "p_loo", "elpd_diff"),
                        labels = c("LOOIC", "ELPD", "$p_\\mathrm{LOO}$", "ELPD difference"))) %>% 
  arrange(stat) %>% 
  
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/ndmm_looic_comparison_nma.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = FALSE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         hline.after = NULL,
         comment = FALSE)}

colnames(looic_by_study_nma) <- c("PH", "Non-PH")
looic_by_study_nma %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Study") %>% 
  mutate_at(vars(-Study), tabfmt, digits = 1)%>% 
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/ndmm_looic_by_study_nma.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = TRUE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         hline.after = 0,
         comment = FALSE)}


# Plot hazards ------------------------------------------------------------

# Baseline hazards can be produced using the predict() function and then plotted

# We'll show the baseline hazard for the "reference" individual
refdat <- tibble(study = ndmm_net$studies, 
                 age = ndmm_fit$xbar["age"],
                 iss_stage3 = 0,
                 response_cr_vgpr = 0,
                 male = 0)

# At evenly spaced times between the boundary knots
tdat <- purrr::imap_dfr(ndmm_fit$basis, 
                        ~tibble(study = factor(.y, levels = ndmm_net$studies), 
                                lower = attr(.x, "Boundary.knots")[1], 
                                upper = attr(.x, "Boundary.knots")[2],
                                times = seq(lower, upper, length = 50)))
  
refdat <- left_join(refdat, tdat, by = "study")

studies <- as.list(setNames(nm = levels(ndmm_net$studies)))

# Also show the knot locations
knotdat7 <- purrr::imap_dfr(ndmm_fit$basis,
                            ~data.frame(Study = .y,
                                        knots = c(attr(.x, "knots"), attr(.x, "Boundary.knots"))))

knotdat10 <- purrr::imap_dfr(ndmm_fit_10kt$basis,
                             ~data.frame(Study = .y,
                                         knots = c(attr(.x, "knots"), attr(.x, "Boundary.knots"))))

hp7 <- plot(predict(ndmm_fit, type = "hazard", level = "individual", 
                    newdata = refdat, study = study, times = times,
                    baseline = studies, aux = studies)) +
  geom_vline(aes(xintercept = knots), data = knotdat7,
             linetype = 2, colour = "grey30", alpha = 0.25, linewidth = 0.25) +
  scale_colour_manual("Treatment", 
                      breaks = c("Pbo", "Len", "Thal"), 
                      values = trt_pal,
                      aesthetics = c("colour", "fill")) +
  theme_multinma(base_family = "Source Sans Pro")

hp10 <- plot(predict(ndmm_fit_10kt, type = "hazard", level = "individual", 
                    newdata = refdat, study = study, times = times,
                    baseline = studies, aux = studies)) +
  geom_vline(aes(xintercept = knots), data = knotdat10,
             linetype = 2, colour = "grey30", alpha = 0.25, linewidth = 0.25) +
  scale_colour_manual("Treatment", 
                      breaks = c("Pbo", "Len", "Thal"), 
                      values = trt_pal,
                      aesthetics = c("colour", "fill")) +
  theme_multinma(base_family = "Source Sans Pro")

# Combining these into a single plot using patchwork
hp7 + facet_grid(rows = vars("7 knots"), cols = vars(Study)) +
  hp10 + facet_grid(rows = vars("10 knots"), cols = vars(Study)) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "top", legend.box.spacing = unit(0, "lines"))

ggsave("./figures/ndmm_est_baseline_hazard_curves.pdf", width = 6, height = 5, scale = 1.2)


# Population-average marginal hazards are also produced using predict()

mhp7 <- plot(predict(ndmm_fit, type = "hazard", level = "aggregate")) +
  geom_vline(aes(xintercept = knots), data = knotdat7,
             linetype = 2, colour = "grey30", alpha = 0.25, linewidth = 0.25) +
  scale_colour_manual("Treatment",
                      breaks = c("Pbo", "Len", "Thal"),
                      values = trt_pal,
                      aesthetics = c("colour", "fill")) +
  theme_multinma(base_family = "Source Sans Pro")

mhp10 <- plot(predict(ndmm_fit_10kt, type = "hazard", level = "aggregate")) +
  geom_vline(aes(xintercept = knots), data = knotdat10,
             linetype = 2, colour = "grey30", alpha = 0.25, linewidth = 0.25) +
  scale_colour_manual("Treatment",
                      breaks = c("Pbo", "Len", "Thal"),
                      values = trt_pal,
                      aesthetics = c("colour", "fill")) +
  theme_multinma(base_family = "Source Sans Pro")

# Combining these into a single plot using patchwork
mhp7 + facet_grid(rows = vars("7 knots"), cols = vars(Study)) +
  mhp10 + facet_grid(rows = vars("10 knots"), cols = vars(Study)) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "top", legend.box.spacing = unit(0, "lines"))

ggsave("./figures/ndmm_est_hazard_curves.pdf", width = 6, height = 5, scale = 1.2)


# Plot survival curves ----------------------------------------------------

# Population-average survival curves are also produced using predict() 
plot(predict(ndmm_fit, type = "survival")) +
  # Overlay the KM data
  geom_km(ndmm_net) +
  scale_colour_manual("Treatment", 
                      breaks = c("Pbo", "Len", "Thal"), 
                      values = trt_pal,
                      aesthetics = c("colour", "fill")) +
  theme_multinma(base_family = "Source Sans Pro") +
  theme(legend.position = "top", legend.box.spacing = unit(0, "lines"))

ggsave("./figures/ndmm_est_survival_curves.pdf", width = 6, height = 3, scale = 1.2)

ggsave("./figures/ndmm_est_survival_curves_PRES.pdf", width = 12, height = 5, scale = 1.2)


# Population-average relative effects -------------------------------------

# Produce population-average (conditional) log hazard ratios using the
# relative_effects() function
(loghr <- relative_effects(ndmm_fit, all_contrasts = TRUE))


# LaTeX table
as_tibble(loghr) %>% 
  transmute(.study,
            label = paste0(.trtb, " vs. ", .trta),
            est_fmt = tabfmt(mean),
            ci_fmt = paste0("(", tabfmt(`2.5%`), ", ", tabfmt(`97.5%`), ")")) %>% 
  pivot_longer(names_to = "stat", values_to = "stat_fmt", cols = c(est_fmt, ci_fmt)) %>% 
  group_by(stat) %>% 
  pivot_wider(names_from = label, values_from = stat_fmt) %>% 
  arrange(.study, desc(stat)) %>% 
  ungroup() %>% 
  mutate(.study = subfmt(.study)) %>% 
  select(-stat) %>% 
  
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/ndmm_log_hazard_ratios.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = FALSE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         add.to.row = lst(pos = list(2, 4, 6, 8, 10),
                          command = rep_len("[0.5ex]", length(pos))),
         hline.after = NULL,
         comment = FALSE)}


# Population-average absolute effects -------------------------------------

# Produce population-average absolute effects such as median survival times
# using the predict() function
(medsurv <- predict(ndmm_fit, type = "median"))

# LaTeX table
as_tibble(medsurv) %>% 
  transmute(.study, .trt,
            est_fmt = tabfmt(mean),
            ci_fmt = paste0("(", tabfmt(`2.5%`), ", ", tabfmt(`97.5%`), ")")) %>% 
  pivot_longer(names_to = "stat", values_to = "stat_fmt", cols = c(est_fmt, ci_fmt)) %>% 
  group_by(stat) %>% 
  pivot_wider(names_from = .trt, values_from = stat_fmt) %>% 
  arrange(.study, desc(stat)) %>% 
  ungroup() %>% 
  mutate(.study = subfmt(.study)) %>% 
  select(-stat) %>% 
  
  # Print to LaTeX table
  xtable() %>%
  {print(., file = "./tables/ndmm_median_survival.tex",
         sanitize.text.function = function(x) {x},
         include.colnames = FALSE,
         include.rownames = FALSE,
         booktabs = TRUE,
         only.contents = TRUE,
         add.to.row = lst(pos = list(2, 4, 6, 8, 10),
                          command = rep_len("[0.5ex]", length(pos))),
         hline.after = NULL,
         comment = FALSE)}


