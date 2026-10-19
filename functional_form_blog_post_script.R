# ==========================================================
# Title:   Checking the Functional Form of a Continuous Covariate
#          in a Cox Proportional Hazards Model (synthetic BMI example)
# Authors: Fan Fan & Jane Liu
# Date:    2025-07-16 (revised 2026-10-19)
# Purpose: Simulate a U-shaped BMI–survival relationship, then compare
#          linear, quartile, quadratic and spline specifications of BMI.
#
# Figures produced (saved to `out_dir`):
#   fig1_u_shape_concept.png    - conceptual U-shaped BMI–risk curve
#   fig2_bmi_vs_time.png        - BMI vs observed time, events vs censored
#   fig3a_quartile_loghr.png    - log-HR by BMI quartile
#   fig3b_age_illustration.png  - illustrative (hand-entered) linear pattern
#   fig4_model_comparison.png   - apparent and cross-validated C-index
#   fig5a_quad_spline.png       - linear vs quadratic vs natural spline
#   fig5b_overfit_spline.png    - flexible spline on a small subsample (overfitting)
# ==========================================================

library(dplyr)
library(survival)
library(ggplot2)
library(splines)

# ---------------------------------------------------------
# 0. Parameters (change here, everything below follows)
# ---------------------------------------------------------
seed            <- 123
n               <- 1000      # number of patients
bmi_mean        <- 27
bmi_sd          <- 5
bmi_optimum     <- 25        # BMI with the longest expected survival
max_time        <- 200
noise_mean      <- 50
noise_sd        <- 5
target_censor   <- 0.20      # target proportion censored
censor_mode     <- "independent"  # "independent" (recommended) or "original"
k_folds         <- 5         # folds for cross-validated C-index
spline_df       <- 4         # df for the natural spline in Fig 5A
overfit_n       <- 150       # small subsample for the overfitting demo (Fig 5B)
overfit_df      <- 20        # df for the overfit spline in Fig 5B
out_dir         <- "figures"

dir.create(out_dir, showWarnings = FALSE)
set.seed(seed)

# ---------------------------------------------------------
# Figure 1: conceptual U-shaped BMI–risk curve (illustration, no data)
#   log relative risk is set to 0 at BMI 22 and rises on both sides.
#   A straight line fitted to the same curve shows what a linear
#   term would capture instead.
# ---------------------------------------------------------
fig1_optimum <- 22
fig1_df <- data.frame(bmi = seq(15, 40, by = 0.1)) %>%
  mutate(log_rr = 0.01 * (bmi - fig1_optimum)^2)  # U-shape on the log scale
# Best straight line on the log scale, weighted toward typical BMIs
fig1_wts    <- dnorm(fig1_df$bmi, mean = bmi_mean, sd = bmi_sd)
fig1_linear <- lm(log_rr ~ bmi, data = fig1_df, weights = fig1_wts)
fig1_df$log_rr_linear <- predict(fig1_linear)

p_fig1 <- ggplot(fig1_df, aes(bmi)) +
  annotate("rect", xmin = -Inf, xmax = 18.5, ymin = -Inf, ymax = Inf,
           fill = "gray85", alpha = 0.5) +
  annotate("rect", xmin = 35, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "gray85", alpha = 0.5) +
  annotate("text", x = 16.75, y = max(fig1_df$log_rr) * 0.97,
           label = "BMI < 18.5", size = 3.5, colour = "gray30") +
  annotate("text", x = 37.5, y = min(fig1_df$log_rr_linear) * 0.9,
           label = "BMI > 35", size = 3.5, colour = "gray30") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
  geom_line(aes(y = log_rr_linear, colour = "Linear assumption"),
            linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = log_rr, colour = "U-shaped relationship"), linewidth = 1.3) +
  geom_point(data = data.frame(bmi = fig1_optimum, log_rr = 0), aes(y = log_rr),
             size = 3) +
  annotate("text", x = fig1_optimum, y = 0, vjust = -1.5,
           label = paste("Lowest risk at BMI", fig1_optimum), size = 3.5) +
  scale_colour_manual(values = c("U-shaped relationship" = "black",
                                 "Linear assumption" = "steelblue")) +
  scale_x_continuous(breaks = seq(15, 40, 5)) +
  labs(x = "BMI", y = paste0("log(Relative risk) vs BMI ", fig1_optimum),
       colour = NULL,
       title = "A U-shaped association between BMI and risk") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_fig1)
ggsave(file.path(out_dir, "fig1_u_shape_concept.png"), p_fig1, width = 7, height = 5)

# ---------------------------------------------------------
# 1. Simulate data
# ---------------------------------------------------------
bmi <- rnorm(n, mean = bmi_mean, sd = bmi_sd)

# U-shaped true survival time: longest at bmi_optimum
noise     <- rnorm(n, mean = noise_mean, sd = noise_sd)
true_time <- round(pmax(max_time - (bmi - bmi_optimum)^2 + noise, 1))

if (censor_mode == "independent") {
  # Censoring time drawn independently of the event time (non-informative
  # censoring, as the Cox model assumes). The exponential rate is tuned so
  # that about `target_censor` of patients end up censored.
  u <- runif(n)
  censor_prop <- function(rate) mean(ceiling(-log(u) / rate) < true_time)
  rate <- uniroot(function(r) censor_prop(r) - target_censor,
                  interval = c(1e-6, 1))$root
  censor_time <- ceiling(-log(u) / rate)
} else {
  # Original scheme: 20% of patients get a censoring time drawn between 1 and
  # their own true time. Note this makes censoring depend on the outcome.
  censor_flag <- rbinom(n, 1, target_censor)
  censor_time <- ifelse(censor_flag == 1,
                        round(runif(n, min = 1, max = true_time)),
                        Inf)
}

surv_data <- data.frame(
  id     = seq_len(n),
  bmi    = bmi,
  time   = pmin(true_time, censor_time),
  status = as.integer(true_time <= censor_time)   # 1 = event, 0 = censored
)

cat("\nCensoring check (0 = censored, 1 = event):\n")
print(prop.table(table(surv_data$status)))
n_events <- sum(surv_data$status)
cat("Number of events:", n_events, "\n")

# ---------------------------------------------------------
# 2. Figure 2: BMI vs observed time
# ---------------------------------------------------------
p_fig2 <- ggplot(surv_data,
                 aes(bmi, time,
                     colour = factor(status, 1:0, c("Event", "Censored")),
                     shape  = factor(status, 1:0, c("Event", "Censored")))) +
  geom_point(alpha = 0.5) +
  scale_colour_manual(values = c(Event = "blue", Censored = "red")) +
  scale_shape_manual(values = c(Event = 19, Censored = 1)) +
  labs(x = "BMI", y = "Observed time", colour = NULL, shape = NULL,
       title = "BMI vs. observed survival time") +
  theme_minimal()
print(p_fig2)
ggsave(file.path(out_dir, "fig2_bmi_vs_time.png"), p_fig2, width = 7, height = 5)

# ---------------------------------------------------------
# 3. Quartile model (Figure 3A)
# ---------------------------------------------------------
quartile_breaks <- quantile(surv_data$bmi, probs = seq(0, 1, 0.25))
surv_data$bmi_quartile <- cut(surv_data$bmi, breaks = quartile_breaks,
                              include.lowest = TRUE,
                              labels = c("Q1", "Q2", "Q3", "Q4"))
# cut() returns a factor whose first level (Q1) is the reference

cox_quartile <- coxph(Surv(time, status) ~ bmi_quartile, data = surv_data)
print(summary(cox_quartile))

quartile_medians <- surv_data %>%
  group_by(bmi_quartile) %>%
  summarize(median_bmi = median(bmi), .groups = "drop")

ci <- confint(cox_quartile)   # on the log-HR scale
quartile_df <- data.frame(
  bmi_quartile = c("Q1", sub("bmi_quartile", "", names(coef(cox_quartile)))),
  log_hr       = c(0, coef(cox_quartile)),
  lower        = c(0, ci[, 1]),
  upper        = c(0, ci[, 2])
) %>%
  left_join(quartile_medians %>% mutate(bmi_quartile = as.character(bmi_quartile)),
            by = "bmi_quartile")

p_fig3a <- ggplot(quartile_df, aes(median_bmi, log_hr)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
  geom_point(size = 3) +
  geom_text(aes(label = bmi_quartile), vjust = -1.2) +
  labs(x = "Median BMI per quartile", y = "log(Hazard ratio) vs Q1",
       title = "Log-hazard ratios across BMI quartiles") +
  theme_minimal()
print(p_fig3a)
ggsave(file.path(out_dir, "fig3a_quartile_loghr.png"), p_fig3a, width = 7, height = 5)

# ---------------------------------------------------------
# 4. Figure 3B: illustrative linear pattern (hand-entered, NOT simulated)
# ---------------------------------------------------------
age_illustration <- data.frame(
  age_quartile = c("Q1", "Q2", "Q3", "Q4"),
  median_age   = c(45, 55, 65, 75),      # illustrative median ages
  HR           = c(1.0, 2.3, 3.3, 6.7),  # illustrative hazard ratios
  lower        = c(1.0, 1.8, 2.5, 5.2),
  upper        = c(1.0, 2.8, 4.1, 8.5)
)

p_fig3b <- ggplot(age_illustration, aes(median_age, log(HR))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray") +
  geom_errorbar(aes(ymin = log(lower), ymax = log(upper)), width = 1) +
  geom_point(size = 3, colour = "red") +
  labs(x = "Median age per quartile", y = "log(Hazard ratio) vs Q1",
       title = "Illustration: an approximately linear pattern (age)") +
  theme_minimal()
print(p_fig3b)
ggsave(file.path(out_dir, "fig3b_age_illustration.png"), p_fig3b, width = 7, height = 5)

# ---------------------------------------------------------
# 5. Compare BMI specifications (Figure 4)
#    Splines are written inside the formula so that predict() reuses the
#    same knots on new data.
# ---------------------------------------------------------
spline_label <- paste0("Natural spline (", spline_df, " df)")
model_formulas <- list(
  Surv(time, status) ~ bmi,
  Surv(time, status) ~ bmi_quartile,
  Surv(time, status) ~ bmi + I(bmi^2),
  as.formula(paste0("Surv(time, status) ~ ns(bmi, df = ", spline_df, ")"))
)
names(model_formulas) <- c("Linear (1 df)", "Quartiles (3 df)",
                           "Quadratic (2 df)", spline_label)
# bquote() writes each formula into the stored call, so predict() on new
# data can rebuild the model frame later
fits <- lapply(model_formulas, function(f)
  eval(bquote(coxph(.(f), data = surv_data))))

# Apparent (in-sample) C-index, read from the fitted models
apparent_c <- sapply(fits, function(f) unname(f$concordance["concordance"]))

# K-fold cross-validated C-index: fit on K-1 folds, score the held-out fold.
# Quartile cut points are re-derived inside each training fold.
cv_cindex <- function(formula, data, k) {
  folds <- sample(rep(seq_len(k), length.out = nrow(data)))
  sapply(seq_len(k), function(i) {
    train <- data[folds != i, ]
    test  <- data[folds == i, ]
    brks  <- quantile(train$bmi, seq(0, 1, 0.25))
    brks[c(1, 5)] <- c(-Inf, Inf)
    lab <- c("Q1", "Q2", "Q3", "Q4")
    train$bmi_quartile <- cut(train$bmi, brks, labels = lab)
    test$bmi_quartile  <- cut(test$bmi,  brks, labels = lab)
    fit <- coxph(formula, data = train)
    test$lp <- predict(fit, newdata = test, type = "lp")
    # reverse = TRUE: higher linear predictor means shorter survival
    unname(concordance(Surv(time, status) ~ lp, data = test,
                       reverse = TRUE)$concordance)
  }) %>% mean()
}
set.seed(seed)
cv_c <- sapply(model_formulas, cv_cindex, data = surv_data, k = k_folds)

model_comparison <- data.frame(
  model      = names(fits),
  df         = sapply(fits, function(f) sum(!is.na(coef(f)))),
  AIC        = sapply(fits, AIC),
  apparent_C = apparent_c,
  cv_C       = cv_c,
  row.names  = NULL
)
cat("\nModel comparison (lower AIC is better; higher C is better):\n")
print(model_comparison, digits = 4)

# Likelihood ratio test: does the spline improve on the nested linear model?
cat("\nLikelihood ratio test, linear vs natural spline:\n")
print(anova(fits[["Linear (1 df)"]], fits[[spline_label]]))

p_fig4 <- model_comparison %>%
  tidyr::pivot_longer(c(apparent_C, cv_C), names_to = "type", values_to = "C") %>%
  mutate(type = recode(type, apparent_C = "Apparent (training data)",
                       cv_C = paste0(k_folds, "-fold cross-validated")),
         model = factor(model, levels = names(fits))) %>%
  ggplot(aes(model, C, fill = type)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", C)),
            position = position_dodge(0.8), vjust = -0.4, size = 3.2) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = NULL, y = "C-index", fill = NULL,
       title = "Prognostic performance by BMI specification") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_fig4)
ggsave(file.path(out_dir, "fig4_model_comparison.png"), p_fig4, width = 8, height = 5)

# ---------------------------------------------------------
# 6. Figure 5A: estimated log-HR curves, relative to BMI = bmi_optimum
# ---------------------------------------------------------
bmi_grid <- data.frame(bmi = seq(min(surv_data$bmi), max(surv_data$bmi),
                                 length.out = 200))

curve_for <- function(fit, label, ref = bmi_optimum) {
  pred <- predict(fit, newdata = bmi_grid, type = "lp", se.fit = TRUE)
  ref_lp <- predict(fit, newdata = data.frame(bmi = ref), type = "lp")
  data.frame(bmi = bmi_grid$bmi, model = label,
             log_hr = pred$fit - ref_lp,
             # band shows pointwise uncertainty of the curve's shape
             lower = pred$fit - ref_lp - 1.96 * pred$se.fit,
             upper = pred$fit - ref_lp + 1.96 * pred$se.fit)
}
curves <- bind_rows(
  curve_for(fits[["Linear (1 df)"]],    "Linear"),
  curve_for(fits[["Quadratic (2 df)"]], "Quadratic"),
  curve_for(fits[[spline_label]], paste0("Natural spline (", spline_df, " df)"))
)
bmi_p05_p95 <- quantile(surv_data$bmi, c(0.05, 0.95))

p_fig5a <- ggplot(curves, aes(bmi, log_hr, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = bmi_p05_p95, linetype = "dotted", colour = "gray40") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
  labs(x = "BMI", y = paste0("log(Hazard ratio) vs BMI ", bmi_optimum),
       colour = NULL, fill = NULL,
       title = "Estimated BMI effect under different functional forms",
       caption = "Dotted lines: 5th and 95th percentiles of BMI") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_fig5a)
ggsave(file.path(out_dir, "fig5a_quad_spline.png"), p_fig5a, width = 8, height = 5)

# ---------------------------------------------------------
# 7. Figure 5B: overfitting with too many df in a small sample
#    A flexible spline is fit to a random subsample of `overfit_n`
#    patients and compared with the quadratic fit from the full data.
# ---------------------------------------------------------
set.seed(seed)
small_data   <- surv_data[sample(nrow(surv_data), overfit_n), ]
small_events <- sum(small_data$status)
cat("\nOverfitting demo: n =", overfit_n, ", events =", small_events,
    ", df =", overfit_df, ", events per df =", round(small_events / overfit_df, 1),
    "(rule of thumb: about 10 or more)\n")

overfit_fit <- eval(bquote(coxph(
  .(as.formula(paste0("Surv(time, status) ~ bs(bmi, df = ", overfit_df, ")"))),
  data = small_data)))

# Plot only within the 5th-95th percentile range, where data are dense
grid_mid <- filter(bmi_grid, bmi >= bmi_p05_p95[1], bmi <= bmi_p05_p95[2])
lp_rel <- function(fit, grid)
  predict(fit, newdata = grid, type = "lp") -
  predict(fit, newdata = data.frame(bmi = bmi_optimum), type = "lp")

fig5b_df <- bind_rows(
  data.frame(grid_mid, log_hr = lp_rel(overfit_fit, grid_mid),
             model = paste0("B-spline, ", overfit_df, " df (n = ", overfit_n, ")")),
  data.frame(grid_mid, log_hr = lp_rel(fits[["Quadratic (2 df)"]], grid_mid),
             model = paste0("Quadratic, 2 df (n = ", n, ")"))
)

p_fig5b <- ggplot(fig5b_df, aes(bmi, log_hr, colour = model, linetype = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
  geom_line(linewidth = 1.1) +
  geom_rug(data = small_data, aes(x = bmi), inherit.aes = FALSE,
           alpha = 0.4, sides = "b") +
  scale_colour_manual(values = c("darkorange", "black")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x = "BMI (5th to 95th percentile)",
       y = paste0("log(Hazard ratio) vs BMI ", bmi_optimum),
       colour = NULL, linetype = NULL,
       title = "Too many degrees of freedom for too few events",
       caption = "Rug marks: BMI values in the small subsample") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_fig5b)
ggsave(file.path(out_dir, "fig5b_overfit_spline.png"), p_fig5b, width = 8, height = 5)

cat("\nDone. Figures saved to:", normalizePath(out_dir), "\n")

