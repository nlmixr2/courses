################################################################################
## R/Pharma workshop 2026
## Time-to-event (survival) modelling with nlmixr2
##
## What this script covers:
##   1. Simulate a two-arm trial whose true hazard follows a Gompertz model
##   2. Explore the data: Kaplan-Meier curves and the shape of the hazard
##   3. Fit a (misspecified) exponential model in nlmixr2 and check it with a VPC
##   4. Benchmark the treatment effect with a Cox proportional hazards model
##   5. Fit the correct Gompertz model and compare the two parametric models
##   6. Check the Gompertz model against Kaplan-Meier curves and with a VPC
##
## Key idea: nlmixr2 has no built-in "survival" likelihood. We write the
## log-likelihood for each subject ourselves and pass it in with ll().
################################################################################

# 0. Set up --------------------------------------------------------------------

library(nlmixr2)    # model fitting
library(survival)   # Surv(), survfit(), coxph()
library(survminer)  # ggsurvplot() for Kaplan-Meier plots
library(flexsurv)   # reference parametric survival fits
library(muhaz)      # non-parametric (kernel) hazard estimate
library(dplyr)
library(ggplot2)
library(knitr)      # kable() tables
library(vpc)        # vpc_tte() visual predictive checks


# 1. Simulate data -------------------------------------------------------------
#
# The true model is a Gompertz hazard with a proportional treatment effect:
#
#   h(t) = alpha * exp(gamma * t) * HR^trt
#
# alpha is the hazard at t = 0 and gamma is how fast it grows over time
# (gamma > 0 means risk increases). The drug halves the hazard (HR = 0.5).
# Follow-up stops at one year, so anyone without an event by then is censored.

set.seed(260528)

n_arm       <- 150        # patients per arm
true_alpha  <- 0.002      # baseline hazard at t = 0 (/day)
true_gamma  <- 0.005      # hazard growth rate (/day)
true_log_hr <- log(0.5)   # log hazard ratio, drug vs. placebo
maxt        <- 365        # end of follow-up (days)

# Event times by inverse-CDF sampling: if U ~ Uniform(0, 1), solving
# S(t) = U for t gives a draw from the survival distribution S(t).
gompertz_time <- function(u, alpha, gamma, hr) {
  log(1 + gamma * (-log(u)) / (alpha * hr)) / gamma
}
exponential_time <- function(u, h0, hr) {
  -log(u) / (h0 * hr)
}

n_total  <- 2L * n_arm
test_trt <- rep(0L:1L, each = n_arm)   # 0 = placebo, 1 = drug

u       <- runif(n_total)
t_event <- gompertz_time(u, true_alpha, true_gamma, hr = exp(true_log_hr * test_trt))

tte_sim <- data.frame(
  id       = seq_len(n_total),
  time     = pmin(t_event, maxt),              # observed time (event or censoring)
  event    = as.integer(t_event <= maxt),      # 1 = event, 0 = censored at maxt
  test_trt = test_trt,
  dv       = pmin(t_event, maxt),              # nlmixr2 expects a dv column
  evid     = 0L                                # every row is an observation
)
tte_sim$trt_label <- factor(tte_sim$test_trt, levels = 0:1,
                            labels = c("Placebo", "Drug"))

cat(sprintf(
  "%d patients | Placebo: %d events (%.0f%%) | Drug: %d events (%.0f%%)\n",
  nrow(tte_sim),
  sum(tte_sim$event[tte_sim$test_trt == 0]),
  100 * mean(tte_sim$event[tte_sim$test_trt == 0]),
  sum(tte_sim$event[tte_sim$test_trt == 1]),
  100 * mean(tte_sim$event[tte_sim$test_trt == 1])
))


# 2. Explore the data ----------------------------------------------------------

## 2a. Kaplan-Meier curves by arm

km_fit <- survfit(Surv(time, event) ~ trt_label, data = tte_sim)

ggsurvplot(
  km_fit, data = tte_sim,
  palette      = c("steelblue", "firebrick"),
  legend.labs  = c("Placebo", "Drug"),
  legend.title = "Arm",
  xlab = "Days", ylab = "Survival probability",
  title    = "Kaplan-Meier curves by treatment arm (simulated Gompertz data)",
  conf.int = TRUE, ggtheme = theme_bw()
)

## 2b. Cumulative hazard by arm
#
# H(t) = -log S(t). With a constant hazard (exponential model) H(t) is a
# straight line through the origin. If it curves upward, the hazard is
# increasing over time. Under proportional hazards, the Drug curve is a
# constant multiple (the HR) of the Placebo curve.

ggsurvplot(
  km_fit, data = tte_sim, fun = "cumhaz",
  palette      = c("steelblue", "firebrick"),
  legend.labs  = c("Placebo", "Drug"),
  legend.title = "Arm",
  xlab = "Days", ylab = "Cumulative hazard H(t)",
  title    = "Cumulative hazard by treatment arm",
  conf.int = TRUE, ggtheme = theme_bw()
)

## 2c. What shape is the hazard?
#
# Before choosing a parametric model, compare a non-parametric (kernel) hazard
# estimate with the hazard implied by several standard distributions. The
# model whose curve tracks the kernel estimate (and has the lowest AIC) is a
# good candidate. A flat kernel estimate would suggest an exponential model.

#' @param dat           data frame with one row per subject
#' @param stime,evt     names of the time and event-indicator columns
#' @param cutoff        latest time for the kernel estimate
#' @param bwsmooth      smoothing bandwidth passed to muhaz()
#' @param includeCubic3 also show Royston-Parmar spline models (2 and 3 knots)
plot_hazard <- function(dat, stime, evt, cutoff = max(dat[[stime]]),
                        bwsmooth = 3, includeCubic3 = TRUE) {
  # Kernel estimate, capped at 90% of the time range to avoid the
  # right-boundary artefact
  kernel_haz_est <- muhaz(dat[[stime]], dat[[evt]],
                          max.time = cutoff * 0.9, bw.smooth = bwsmooth)
  kernel_haz <- tibble(
    time = kernel_haz_est$est.grid,
    est  = kernel_haz_est$haz.est
  ) %>% filter(!is.nan(est))

  # Parametric fits with flexsurv
  dists <- c(exp = "Exponential", weibull = "Weibull", gompertz = "Gompertz",
             gamma = "Gamma", lognormal = "Lognormal", llogis = "Log-logistic")

  dat$STIME <- dat[[stime]]
  dat$EVENT <- dat[[evt]]

  hazard_curve <- function(fit, method) {
    out <- summary(fit, type = "hazard", ci = FALSE, tidy = TRUE)
    out$Method <- method
    out$AIC    <- AIC(fit)
    out
  }

  haz_data <- do.call(rbind, lapply(names(dists), function(d) {
    fit <- flexsurvreg(Surv(STIME, EVENT) ~ 1, data = dat, dist = d)
    hazard_curve(fit, dists[[d]])
  }))

  if (includeCubic3) {
    spline_haz <- lapply(2:3, function(k) {
      fit <- flexsurvspline(Surv(STIME, EVENT) ~ 1, data = dat, k = k, scale = "hazard")
      hazard_curve(fit, sprintf("Cubic splines (%d knots)", k))
    })
    haz_data <- rbind(haz_data, do.call(rbind, spline_haz))
  }

  # Repeat the kernel curve in every panel. A layer with no Method column
  # would make every panel share one y range, defeating scales = "free_y".
  kernel_all <- do.call(rbind, lapply(unique(haz_data$Method), function(m) {
    data.frame(time = kernel_haz$time, est = kernel_haz$est, Method = m)
  }))

  aic_labels <- haz_data %>%
    group_by(Method) %>% slice_tail(n = 1) %>% ungroup() %>%
    mutate(label = paste0("AIC: ", round(AIC, 1)))

  ggplot() +
    geom_line(data = kernel_all,
              aes(time, est, linetype = "Kernel density"),
              col = "black", linewidth = 0.8) +
    geom_line(data = haz_data,
              aes(time, est, linetype = "Parametric"),
              col = "blue", linewidth = 0.8) +
    scale_linetype_manual(
      values = c("Kernel density" = "dashed", "Parametric" = "solid"), name = NULL) +
    geom_text(data = aic_labels, aes(x = Inf, y = Inf, label = label),
              hjust = 1.05, vjust = 1.5, size = 2.8, col = "blue",
              inherit.aes = FALSE) +
    scale_x_continuous("Time (d)") + scale_y_continuous("Hazard") +
    facet_wrap(~Method, ncol = 2, scales = "free_y") +
    theme_light() +
    theme(legend.position = "bottom", panel.grid = element_line(linetype = 3))
}

plot_hazard(dat = tte_sim, stime = "time", evt = "event")


# Helpers for the visual predictive checks (VPCs) used below -------------------
#
# A time-to-event VPC simulates many replicate trials from the fitted model,
# applies the same censoring as the real study, and checks that the observed
# Kaplan-Meier curve falls inside the band of simulated curves.

n_vpc_sim <- 50   # replicate trials per VPC

vpc_obs <- tte_sim %>%
  transmute(id = id, t = time, dv = event, trt = trt_label)

# Simulate n_vpc_sim trials. event_time() maps uniform draws and each
# subject's hazard ratio to event times for the fitted model.
simulate_vpc <- function(event_time, log_hr) {
  hr <- exp(log_hr * tte_sim$test_trt)
  bind_rows(lapply(seq_len(n_vpc_sim), function(s) {
    t_evt <- event_time(runif(nrow(tte_sim)), hr)
    data.frame(
      sim = s,
      id  = tte_sim$id,
      t   = pmin(t_evt, maxt),             # same censoring as the real trial
      dv  = as.integer(t_evt <= maxt),
      trt = tte_sim$trt_label
    )
  }))
}

# About 10-15 s per VPC at 50 simulations. Runtime grows faster than
# linearly with the number of simulations (about 3 min at 200).
plot_tte_vpc <- function(sim, title) {
  vpc::vpc_tte(
    obs      = vpc_obs,
    sim      = sim,
    obs_cols = list(id = "id", idv = "t", dv = "dv"),
    sim_cols = list(id = "id", idv = "t", dv = "dv", sim = "sim"),
    stratify = "trt",
    labeller = ggplot2::label_both,   # vpc's default labeller turns text strata into "NA"
    show     = list(bin_sep = FALSE), # with bins = FALSE, bin separators form a solid black bar
    as_percentage = FALSE,
    xlab  = "Days",
    ylab  = "Survival probability",
    title = title
  )
}


# 3. Exponential model (constant hazard) ---------------------------------------
#
# The simplest parametric model assumes the hazard is constant over time:
#
#   h(t) = h0 * HR^trt        H(t) = h(t) * t     (H = cumulative hazard)
#
# For right-censored data, each subject's log-likelihood is
#
#   ll = event * log h(t) - H(t)
#
# Subjects with an event contribute the density, f(t) = h(t) * S(t).
# Censored subjects contribute only survival, S(t) = exp(-H(t)).
#
# Here we write that log-likelihood out and hand it to nlmixr2 with ll().
# bobyqa is used because the model has no random effects (no etas).

exp_model <- function() {
  ini({
    log_h0      <- log(0.004)   # log baseline hazard (/day)
    test_log_hr <- -0.1         # log hazard ratio, drug vs. placebo
  })
  model({
    log_h  <- log_h0 + test_log_hr * test_trt
    h      <- exp(log_h)
    H      <- h * time
    tte_ll <- event * log_h - H
    ll(tte) ~ tte_ll
  })
}

fit_exp <- nlmixr(exp_model, tte_sim, est = "bobyqa",
                  control = bobyqaControl(print = 0))
fit_exp$parFixedDf

theta_e <- fit_exp$theta
h0_exp  <- exp(theta_e["log_h0"])

## 3a. VPC for the exponential model
#
# Look for the misfit: the model can't bend to follow a hazard that rises
# over time, so the observed curve drifts outside the simulated band.

set.seed(8675)
vpc_sim_exp <- simulate_vpc(
  event_time = function(u, hr) exponential_time(u, h0_exp, hr),
  log_hr     = theta_e["test_log_hr"]
)
plot_tte_vpc(vpc_sim_exp, "Exponential TTE VPC by treatment arm")


# 4. Cox proportional hazards model --------------------------------------------
#
# The Cox model estimates the hazard ratio without assuming any shape for the
# baseline hazard, which makes it a useful benchmark. Compare coef (log HR)
# with test_log_hr from each nlmixr2 fit. The exponential estimate is pulled
# toward zero because its baseline hazard is wrong. The Gompertz estimate
# (next section) should almost match Cox.

fit_cox <- coxph(Surv(time, event) ~ trt_label, data = tte_sim)
summary(fit_cox)


# 5. Gompertz model (the true model) -------------------------------------------
#
#   h(t) = alpha * exp(gamma * t) * HR^trt
#   H(t) = (alpha / gamma) * (exp(gamma * t) - 1) * HR^trt
#
# The log-likelihood has the same form as before. Only h(t) and H(t) change.

gompertz_model <- function() {
  ini({
    log_alpha   <- log(0.003)   # log hazard at t = 0 (/day)
    log_gamma   <- log(0.003)   # log growth rate; estimating on the log scale keeps gamma > 0
    test_log_hr <- -0.5         # informed by KM: ~40-50% separation visible
  })
  model({
    alpha  <- exp(log_alpha)
    gamma  <- exp(log_gamma)
    h      <- alpha * exp(gamma * time) * exp(test_log_hr * test_trt)
    H      <- (alpha / gamma) * (exp(gamma * time) - 1) * exp(test_log_hr * test_trt)
    tte_ll <- event * log(h) - H
    ll(tte) ~ tte_ll
  })
}

fit_gompertz <- nlmixr(gompertz_model, tte_sim, est = "bobyqa",
                       control = bobyqaControl(print = 0))
fit_gompertz$parFixedDf

theta_g <- fit_gompertz$theta
alpha_g <- exp(theta_g["log_alpha"])
gamma_g <- exp(theta_g["log_gamma"])
hr_drug <- exp(theta_g["test_log_hr"])

## 5a. Did we recover the simulation parameters?

kable(
  data.frame(
    Parameter = c("alpha (/day)", "gamma (/day)", "HR (drug vs. placebo)"),
    True      = c(true_alpha, true_gamma, exp(true_log_hr)),
    Estimated = round(c(alpha_g, gamma_g, hr_drug), 4),
    Interpretation = c(
      "Initial hazard at t = 0, placebo arm",
      "Hazard growth rate; >0 confirms increasing risk",
      "Instantaneous hazard ratio, drug vs. placebo"
    )
  ),
  digits = 4,
  caption = "Gompertz parameter recovery: true vs. estimated values"
)

## 5b. Exponential vs. Gompertz: which fits better? (lower AIC is better)

aic_table <- AIC(fit_exp, fit_gompertz) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Model") %>%
  mutate(
    Model = c("Exponential", "Gompertz"),
    m2LL  = round(c(-2 * logLik(fit_exp), -2 * logLik(fit_gompertz)), 2)
  )

kable(aic_table, digits = 2,
      caption = "AIC and -2 log-likelihood: exponential vs. Gompertz")


# 6. Check the Gompertz model --------------------------------------------------

# Model-predicted survival curves, S(t) = exp(-H(t))
t_grid <- seq(0, maxt, length.out = 400)

exp_surv      <- function(t)     exp(-h0_exp * t)
gompertz_surv <- function(t, hr) exp(-(alpha_g / gamma_g) * (exp(gamma_g * t) - 1) * hr)

# Kaplan-Meier curve for one arm, starting at S(0) = 1
km_df <- function(trt, label) {
  fit <- survfit(Surv(time, event) ~ 1, data = filter(tte_sim, test_trt == trt))
  data.frame(time = c(0, fit$time), survival = c(1, fit$surv), label = label)
}

## 6a. Both models vs. Kaplan-Meier (placebo arm)

model_pred_df <- bind_rows(
  data.frame(time = t_grid, survival = exp_surv(t_grid),         label = "Exponential"),
  data.frame(time = t_grid, survival = gompertz_surv(t_grid, 1), label = "Gompertz")
)

ggplot() +
  geom_step(data = km_df(0, "Kaplan-Meier"), aes(time, survival, color = label),
            linewidth = 0.8) +
  geom_line(data = model_pred_df, aes(time, survival, color = label), linewidth = 1) +
  scale_color_manual(values = c("Kaplan-Meier" = "black",
                                "Exponential"  = "steelblue",
                                "Gompertz"     = "firebrick")) +
  labs(x = "Days", y = "Survival probability",
       title = "Exponential vs. Gompertz model fit (placebo arm)", color = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() + theme(legend.position = "bottom")

## 6b. Gompertz model vs. Kaplan-Meier, both arms

pred_df <- bind_rows(
  data.frame(time = t_grid, survival = gompertz_surv(t_grid, 1),       label = "Placebo (predicted)"),
  data.frame(time = t_grid, survival = gompertz_surv(t_grid, hr_drug), label = "Drug (predicted)")
)
km_both <- bind_rows(km_df(0, "Placebo (KM)"), km_df(1, "Drug (KM)"))

ggplot() +
  geom_step(data = km_both, aes(time, survival, color = label), linewidth = 0.8) +
  geom_line(data = pred_df, aes(time, survival, color = label),
            linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Placebo (KM)" = "steelblue", "Placebo (predicted)" = "steelblue",
                                "Drug (KM)"    = "firebrick", "Drug (predicted)"    = "firebrick")) +
  labs(x = "Days", y = "Survival probability",
       title = "Gompertz model vs. Kaplan-Meier by arm (simulated data)",
       color = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() + theme(legend.position = "bottom")

## 6c. VPC for the Gompertz model: compare with the exponential VPC in 3a

set.seed(8675)
vpc_sim_gompertz <- simulate_vpc(
  event_time = function(u, hr) gompertz_time(u, alpha_g, gamma_g, hr),
  log_hr     = theta_g["test_log_hr"]
)
plot_tte_vpc(vpc_sim_gompertz, "Gompertz TTE VPC by treatment arm")
