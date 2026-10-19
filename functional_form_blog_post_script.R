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
#   fig4_model_comparison.png   - C-index, linear vs quartile BMI
#   fig5a_quad_spline.png       - linear vs linear + quadratic + spline
#   fig5b_overfit_spline.png    - very flexible spline (df = 100)
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
# 5. Compare linear vs quartile BMI (Figure 4)
# ---------------------------------------------------------
model_formulas <- list(
  "Linear (1 df)"    = Surv(time, status) ~ bmi,
  "Quartiles (3 df)" = Surv(time, status) ~ bmi_quartile
)
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

p_fig4 <- model_comparison %>%
  tidyr::pivot_longer(c(apparent_C, cv_C), names_to = "type", values_to = "C") %>%
  mutate(type = recode(type, apparent_C = "Apparent (training data)",
                       cv_C = paste0(k_folds, "-fold cross-validated")),
         model = factor(model, levels = names(fits))) %>%
  ggplot(aes(model, C, fill = type)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", C)),
            position = position_dodge(0.8), vjust = -0.4, size = 3.5) +
  coord_cartesian(ylim = c(0.5, 0.85)) +
  labs(x = NULL, y = "C-index", fill = NULL,
       title = "Prognostic performance: linear vs quartile BMI") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_fig4)
ggsave(file.path(out_dir, "fig4_model_comparison.png"), p_fig4, width = 7, height = 5)

# ---------------------------------------------------------
# 6. Linear vs linear + quadratic + spline (Figure 5A)
# ---------------------------------------------------------
library(splines)

# Set up
bmi <- surv_data$bmi
bmi_quad <- bmi^2
bmi_spline_basis <- bs(bmi, degree = 3, df = 3)

# Create transformed dataset
surv_data_transformed <- cbind(surv_data,
                               bmi_quad = bmi^2,
                               bmi_cub1 = bmi_spline_basis[,1],
                               bmi_cub2 = bmi_spline_basis[,2],
                               bmi_cub3 = bmi_spline_basis[,3]
)
# Fit both models
cox_linear <- coxph(Surv(time, status) ~ bmi, data = surv_data_transformed)
cox_quad_spline <- coxph(Surv(time, status) ~ bmi + bmi_quad + bmi_cub1 + bmi_cub2 + bmi_cub3,
                         data = surv_data_transformed)
# Note: bmi_cub2 and bmi_cub3 have NA coefficients. The cubic spline basis
# already contains linear and quadratic components, so these two columns are
# redundant with bmi and bmi_quad; coxph drops them automatically.
coef(cox_quad_spline)

# Predict log hazard over range of BMI
X <- seq(min(bmi), max(bmi), length.out = 100)
X_quad <- X^2
X_spline <- predict(bmi_spline_basis, newx = X)   # same knots as the fitted basis

# Coefficients (dropped terms contribute nothing, so treat NA as 0)
coef_spline <- coef(cox_quad_spline)
coef_spline[is.na(coef_spline)] <- 0
coef_linear <- coef(cox_linear)

# Predict linear predictors
pred_spline <- X * coef_spline["bmi"] +
  X_quad * coef_spline["bmi_quad"] +
  X_spline[,1] * coef_spline["bmi_cub1"] +
  X_spline[,2] * coef_spline["bmi_cub2"] +
  X_spline[,3] * coef_spline["bmi_cub3"]

pred_linear <- X * coef_linear["bmi"]

# Shift for comparison
pred_linear <- pred_linear + (mean(pred_spline) - mean(pred_linear))

# Plot and save (Figure 5A)
plot_fig5a <- function() {
  par(mfrow=c(1,1), mar = c(5, 4, 2, 1))
  plot(X, pred_spline, type="l", lwd=3, col="black",
       xlab="BMI", ylab="Log Hazard Ratio",
       main="BMI: Linear vs Linear + Quadratic + Spline")
  lines(X, pred_linear, lwd=3, lty=2, col="blue")
  abline(v = quantile(bmi, c(0.05, 0.95)), col = "gray", lty = 3)
  legend("top", legend = c("Linear + Quadratic + Spline", "Linear"),
         col = c("black", "blue"), lty = c(1, 2), lwd = 3, bty = "n")
}
plot_fig5a()                                   # show on screen
png(file.path(out_dir, "fig5a_quad_spline.png"), width = 8, height = 5, units = "in", res = 300)
plot_fig5a()                                   # save to file
dev.off()


# ---------------------------------------------------------
# 7. Very flexible spline, df = 100 (Figure 5B)
# ---------------------------------------------------------

# ----- Set parameters -----
df_spline <- 100  # change this to increase/decrease degrees of freedom

# Set up
bmi <- surv_data$bmi

# Create spline basis
bmi_spline_basis <- bs(bmi, df = df_spline)
colnames(bmi_spline_basis) <- paste0("bmi_spline", 1:df_spline)

# Add spline basis columns to dataset
surv_data_transformed <- cbind(surv_data, bmi_spline_basis)

# Build formula dynamically for spline terms
formula_terms <- paste(colnames(bmi_spline_basis), collapse = " + ")
cox_formula <- as.formula(paste("Surv(time, status) ~", formula_terms))

# Fit the spline model
cox_spline <- coxph(cox_formula, data = surv_data_transformed)

# Predict log hazard over range of BMI
X <- seq(min(bmi), max(bmi), length.out = 100)
X_spline <- predict(bmi_spline_basis, newx = X)   # same knots as the fitted basis

# Coefficients for spline model (treat any dropped terms as 0)
coef_spline <- coef(cox_spline)
coef_spline[is.na(coef_spline)] <- 0

# Predict linear predictors for spline model
pred_spline <- drop(X_spline %*% coef_spline)

# Predict for linear model
coef_linear <- coef(cox_linear)
pred_linear <- X * coef_linear["bmi"]

# Align for comparison
pred_linear <- pred_linear + (mean(pred_spline) - mean(pred_linear))

# Plot and save comparison (Figure 5B)
plot_fig5b <- function() {
  par(mfrow = c(1,1), mar = c(5, 4, 2, 1))
  plot(X, pred_spline, type = "l", lwd = 3, col = "darkorange",
       xlab = "BMI", ylab = "Log Hazard Ratio",
       main = paste("BMI: Linear vs Spline (df =", df_spline, ")"))
  lines(X, pred_linear, lwd = 3, lty = 2, col = "blue")
  abline(v = quantile(bmi, c(0.05, 0.95)), col = "gray", lty = 3)
  legend("topleft", legend = c(paste("Spline df =", df_spline), "Linear"),
         col = c("darkorange", "blue"), lty = c(1, 2), lwd = 3, bty = "n")
}
plot_fig5b()                                   # show on screen
png(file.path(out_dir, "fig5b_overfit_spline.png"), width = 8, height = 5, units = "in", res = 300)
plot_fig5b()                                   # save to file
dev.off()

cat("\nDone. Figures saved to:", normalizePath(out_dir), "\n")
