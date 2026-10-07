################################################################################
## R/Pharma workshop 2026 - bonus material
## Binary and categorical endpoints with nlmixr2
##
## What this script covers:
##   1. Binary response: logistic regression
##   2. Unordered categories: multinomial logit model
##   3. Ordered categories: proportional-odds model with a random effect
##
## Each example simulates its own dataset, fits the model from the slides,
## then checks the fit with observed-vs-fitted plots and a categorical VPC.
##
## As in survival.R, we supply the likelihood ourselves:
##   - binary data:      y ~ dbinom(1, p)
##   - categorical data: y ~ c(p1, p2, ...)  (the last category is 1 - sum(p))
##     Categorical outcomes must be coded 1, 2, ..., K.
##
## The name on the left of ~ (RESPONSE, CAT, SEV) is only a label for the
## endpoint. nlmixr2 always reads the observed value from the DV column.
##
## Note: the dose column is called DOSE_MG, not DOSE. rxode2 reserves the
## name "dose" (in any case), and a model covariate called DOSE fails with
## "The following parameter(s) are required for solving: DOSE".
################################################################################

# 0. Set up --------------------------------------------------------------------

library(nlmixr2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vpc)        # vpc_cat() visual predictive checks

dose_levels <- c(0, 50, 100, 200)   # mg, used in every example

# Shared helpers ---------------------------------------------------------------

# Draw one category per row from a matrix of category probabilities
# (rows = subjects, columns = categories 1..K; each row sums to 1).
sample_categories <- function(P) {
  cum_p <- t(apply(P, 1, cumsum))
  u     <- runif(nrow(P))
  1L + rowSums(u > cum_p[, -ncol(P), drop = FALSE])
}

# Stacked bars of observed vs. model-fitted category proportions by dose.
# obs:    data frame with DOSE_MG and a factor column `category`
# fitted: data frame with DOSE_MG, `category` and probability `p`
plot_obs_fitted <- function(obs, fitted, fill_label, title) {
  obs_prop <- obs %>%
    count(DOSE_MG, category) %>%
    group_by(DOSE_MG) %>%
    mutate(p = n / sum(n)) %>%
    ungroup()

  bind_rows(Observed = obs_prop, Fitted = fitted, .id = "source") %>%
    mutate(source = factor(source, levels = c("Observed", "Fitted"))) %>%
    ggplot(aes(factor(DOSE_MG), p, fill = category)) +
    geom_col() +
    facet_wrap(~source) +
    labs(x = "Dose", y = "Category proportion", fill = fill_label, title = title) +
    theme_bw()
}

# Categorical VPC. For each of n_vpc_sim replicate trials, simulate_trial()
# returns a vector of simulated categories, one per row of obs_df.
n_vpc_sim <- 200

run_cat_vpc <- function(obs_df, simulate_trial, title) {
  sim_df <- bind_rows(lapply(seq_len(n_vpc_sim), function(s) {
    data.frame(sim = s, id = obs_df$id, DOSE_MG = obs_df$DOSE_MG,
               dv = simulate_trial())
  }))
  vpc::vpc_cat(
    obs      = obs_df,
    sim      = sim_df,
    obs_cols = list(id = "id", idv = "DOSE_MG", dv = "dv"),
    sim_cols = list(id = "id", idv = "DOSE_MG", dv = "dv", sim = "sim"),
    bins     = c(0, 25, 75, 150, 250),   # one bin per dose level
    smooth   = FALSE,                    # show a step per bin, not a smoothed ribbon
    xlab     = "Dose",
    ylab     = "Category proportion",
    title    = title
  )
}


# 1. Binary response: logistic regression --------------------------------------
#
# Each subject either responds (1) or not (0). The probability of response
# depends on dose through the log-odds (logit) scale:
#
#   logit(p) = logit_p0 + beta_dose * DOSE      p = expit(logit(p))
#
# A binary outcome is binomial with a single trial: RESPONSE ~ dbinom(1, p).

## 1a. Simulate data: 80 subjects per dose level

set.seed(20260930)

n_per_dose_bin <- 80
true_logit_p0  <- -1       # log-odds of response at dose 0 (p ~ 0.27)
true_beta_dose <- 0.008    # change in log-odds per mg

binary_data <- data.frame(
  ID      = seq_len(n_per_dose_bin * length(dose_levels)),
  TIME    = 0,
  DOSE_MG = rep(dose_levels, each = n_per_dose_bin)
)
binary_data$RESPONSE <- rbinom(nrow(binary_data), size = 1,
                               prob = plogis(true_logit_p0 + true_beta_dose * binary_data$DOSE_MG))
binary_data$DV <- binary_data$RESPONSE   # nlmixr2 always reads the observation from DV

## 1b. Explore: observed response rate by dose (with 95% binomial CIs)

binary_data %>%
  group_by(DOSE_MG) %>%
  summarise(x = sum(RESPONSE), n = n()) %>%
  rowwise() %>%
  mutate(p     = x / n,
         lower = binom.test(x, n)$conf.int[1],
         upper = binom.test(x, n)$conf.int[2]) %>%
  ggplot(aes(DOSE_MG, p)) +
  geom_line(col = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 8, col = "firebrick") +
  geom_point(col = "steelblue", size = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Dose", y = "Observed response probability",
       title = "Observed binary response by dose") +
  theme_bw()

## 1c. Fit the logistic regression model

logistic_binary_model <- function() {
  ini({
    # Log odds of response at dose 0
    logit_p0 <- -1
    # Change in log odds for each one-unit increase in dose
    beta_dose <- 0.008
  })
  model({
    # Linear predictor on the log-odds scale
    lp <- logit_p0 + beta_dose * DOSE_MG
    # Inverse-logit transformation to obtain a probability
    p <- expit(lp)
    # A binary response is a binomial endpoint with one trial
    one <- 1
    RESPONSE ~ dbinom(one, p)
  })
}

fit_logistic_binary <- nlmixr(
  logistic_binary_model, binary_data, est = "nlm",
  control = nlmControl(print = 0))
fit_logistic_binary

theta_bin <- fit_logistic_binary$theta

# Odds ratio for a 100 mg increase in dose
exp(100 * theta_bin["beta_dose"])

## 1d. VPC: simulate replicate trials from the fitted model

set.seed(8675)
p_fit_bin <- plogis(theta_bin["logit_p0"] + theta_bin["beta_dose"] * binary_data$DOSE_MG)
bin_labels <- c("No response", "Response")

run_cat_vpc(
  obs_df = data.frame(id = binary_data$ID, DOSE_MG = binary_data$DOSE_MG,
                      dv = bin_labels[binary_data$RESPONSE + 1]),
  simulate_trial = function() bin_labels[rbinom(length(p_fit_bin), 1, p_fit_bin) + 1],
  title = "Categorical VPC for binary logistic regression"
)


# 2. Unordered categories: multinomial logit -----------------------------------
#
# Example: the type of the first adverse event, which is Dermatologic (1),
# Gastrointestinal (2) or Neurologic (3). There is no natural ordering.
#
# Category 1 is the reference. Categories 2 and 3 each get their own logit
# relative to category 1, and the probabilities are
#
#   p_k = exp(lp_k) / (1 + exp(lp_2) + exp(lp_3))     (lp_1 = 0)
#
# Dose is centred at 100 mg so the intercepts describe a typical dose.

## 2a. Simulate data: 60 subjects per dose level

ae_labels <- c("Dermatologic", "Gastrointestinal", "Neurologic")

unordered_probs <- function(cat2_intercept, cat2_dose, cat3_intercept, cat3_dose, dose) {
  lp2 <- cat2_intercept + cat2_dose * (dose - 100)
  lp3 <- cat3_intercept + cat3_dose * (dose - 100)
  cbind(1, exp(lp2), exp(lp3)) / (1 + exp(lp2) + exp(lp3))
}

set.seed(20261001)

n_per_dose_cat <- 60

unordered_data <- data.frame(
  ID      = seq_len(n_per_dose_cat * length(dose_levels)),
  TIME    = 0,
  DOSE_MG = rep(dose_levels, each = n_per_dose_cat)
)
# True model: GI events become slightly more likely with dose, neurologic
# events much more likely, both at the expense of dermatologic events.
unordered_data$CAT <- sample_categories(
  unordered_probs(cat2_intercept = -0.4, cat2_dose = 0.004,
                  cat3_intercept = -0.6, cat3_dose = 0.008,
                  dose = unordered_data$DOSE_MG)
)
unordered_data$DV <- unordered_data$CAT

## 2b. Fit the multinomial logit model

unordered_model <- function() {
  ini({
    # Category 2 versus category 1 at DOSE = 100
    cat2_intercept <- -0.2
    cat2_dose <- 0.004
    # Category 3 versus category 1 at DOSE = 100
    cat3_intercept <- 0.3
    cat3_dose <- -0.003
  })
  model({
    # Logits for categories 2 and 3 relative to reference category 1
    lp2 <- cat2_intercept + cat2_dose * (DOSE_MG - 100)
    lp3 <- cat3_intercept + cat3_dose * (DOSE_MG - 100)
    # Category 1 is represented by exp(0) = 1
    denom <- 1 + exp(lp2) + exp(lp3)
    p1 <- 1 / denom
    p2 <- exp(lp2) / denom
    # Three categories: 1, 2, and the remainder category 3
    CAT ~ c(p1, p2)
  })
}

fit_unordered <- nlmixr(unordered_model, unordered_data, est = "nlm",
                        control = nlmControl(print = 0))
fit_unordered

theta_un <- fit_unordered$theta
fitted_probs_un <- function(dose) {
  unordered_probs(theta_un["cat2_intercept"], theta_un["cat2_dose"],
                  theta_un["cat3_intercept"], theta_un["cat3_dose"], dose)
}

## 2c. Observed vs. fitted category proportions

P_levels_un <- fitted_probs_un(dose_levels)
plot_obs_fitted(
  obs    = mutate(unordered_data, category = factor(ae_labels[CAT], levels = ae_labels)),
  fitted = data.frame(DOSE_MG  = rep(dose_levels, times = 3),
                      category = factor(rep(ae_labels, each = length(dose_levels)),
                                        levels = ae_labels),
                      p        = as.vector(P_levels_un)),
  fill_label = "Category",
  title      = "Observed and fitted unordered category distributions"
)

## 2d. VPC

set.seed(8675)
P_fit_un <- fitted_probs_un(unordered_data$DOSE_MG)

run_cat_vpc(
  obs_df = data.frame(id = unordered_data$ID, DOSE_MG = unordered_data$DOSE_MG,
                      dv = ae_labels[unordered_data$CAT]),
  simulate_trial = function() ae_labels[sample_categories(P_fit_un)],
  title = "Categorical VPC for unordered categories"
)


# 3. Ordered categories: proportional-odds model -------------------------------
#
# Example: symptom severity scored 1 (none) to 4 (severe), recorded at four
# visits per subject. Categories are ordered, so we model the cumulative
# probabilities P(Y <= k) with a single latent severity scale:
#
#   P(Y <= k) = expit(theta_k - lp)      theta_1 < theta_2 < theta_3
#
# A higher lp (more drug, or a more sensitive subject) shifts probability
# toward the more severe categories. The same dose effect applies at every
# threshold, which is the "proportional odds" assumption.
#
# Writing theta_2 = theta_1 + exp(log_delta2) (and so on) keeps the
# thresholds in order. eta_sev lets some subjects score consistently higher
# than others, and repeated visits per subject let us estimate it.

## 3a. Simulate data: 25 subjects per dose level, 4 visits each

ordinal_probs <- function(theta1, log_delta2, log_delta3, beta_dose, dose, eta = 0) {
  lp     <- beta_dose * (dose - 100) + eta
  theta2 <- theta1 + exp(log_delta2)
  theta3 <- theta2 + exp(log_delta3)
  cp     <- cbind(plogis(theta1 - lp), plogis(theta2 - lp), plogis(theta3 - lp), 1)
  cp - cbind(0, cp[, -4])   # cumulative -> category probabilities
}

set.seed(20261002)

n_per_dose_ord <- 25
n_visits       <- 4
true_omega_sev <- 0.5       # variance of eta_sev

n_subj_ord <- n_per_dose_ord * length(dose_levels)
subjects <- data.frame(
  ID      = seq_len(n_subj_ord),
  DOSE_MG = rep(dose_levels, each = n_per_dose_ord),
  eta     = rnorm(n_subj_ord, 0, sqrt(true_omega_sev))
)
ordered_data <- subjects[rep(seq_len(n_subj_ord), each = n_visits), ]
ordered_data$TIME <- rep(seq_len(n_visits), times = n_subj_ord)
ordered_data$SEV  <- sample_categories(
  ordinal_probs(theta1 = -0.5, log_delta2 = log(1.4), log_delta3 = log(1.2),
                beta_dose = 0.008, dose = ordered_data$DOSE_MG, eta = ordered_data$eta)
)
ordered_data$DV  <- ordered_data$SEV
ordered_data$eta <- NULL   # the random effects are unobserved
rownames(ordered_data) <- NULL

## 3b. Fit the proportional-odds model
#
# The random effect enters a non-Gaussian likelihood, so we use the Laplace
# approximation.

ordinal_model <- function() {
  ini({
    # First threshold and positive log-increments to the next thresholds
    theta1 <- -1
    log_delta2 <- log(1.4)
    log_delta3 <- log(1.2)
    # Dose effect on the latent severity scale
    beta_dose <- 0.005
    # Subject-level variability in latent severity. Don't start this too
    # small: from 0.05, the Laplace fit stopped at exactly 5x the initial
    # value (0.25) instead of the true optimum (~0.58 for these data).
    eta_sev ~ 0.2
  })
  model({
    # Linear predictor centered at 100 mg for numerical stability
    lp <- beta_dose * (DOSE_MG - 100) + eta_sev
    # Ordered thresholds
    theta2 <- theta1 + exp(log_delta2)
    theta3 <- theta2 + exp(log_delta3)
    # Cumulative probabilities P(Y <= k)
    cp1 <- expit(theta1 - lp)
    cp2 <- expit(theta2 - lp)
    cp3 <- expit(theta3 - lp)
    # Category probabilities
    p1 <- cp1
    p2 <- cp2 - cp1
    p3 <- cp3 - cp2
    # Four ordered categories: 1, 2, 3, and the remainder category 4
    SEV ~ c(p1, p2, p3)
  })
}

fit_ordinal <- nlmixr(ordinal_model, ordered_data, est = "laplace",
                      control = laplaceControl(print = 0))
fit_ordinal

theta_ord <- fit_ordinal$theta
omega_ord <- fit_ordinal$omega["eta_sev", "eta_sev"]
fitted_probs_ord <- function(dose, eta = 0) {
  ordinal_probs(theta_ord["theta1"], theta_ord["log_delta2"], theta_ord["log_delta3"],
                theta_ord["beta_dose"], dose, eta)
}

## 3c. Observed vs. fitted category proportions
#
# The fitted proportions average over subjects (eta), approximated here by
# simulating 2000 random subjects per dose.

set.seed(4242)
eta_draws <- rnorm(2000, 0, sqrt(omega_ord))
P_levels_ord <- t(sapply(dose_levels, function(d) colMeans(fitted_probs_ord(d, eta_draws))))
sev_levels <- as.character(1:4)

plot_obs_fitted(
  obs    = mutate(ordered_data, category = factor(SEV, levels = 1:4, labels = sev_levels)),
  fitted = data.frame(DOSE_MG  = rep(dose_levels, times = 4),
                      category = factor(rep(sev_levels, each = length(dose_levels)),
                                        levels = sev_levels),
                      p        = as.vector(P_levels_ord)),
  fill_label = "Severity",
  title      = "Observed and fitted ordered categorical distributions"
)

## 3d. VPC
#
# Each replicate trial draws a new eta per subject (shared across that
# subject's visits), then a severity score at each visit.

set.seed(8675)
subj_index <- match(ordered_data$ID, unique(ordered_data$ID))

run_cat_vpc(
  obs_df = data.frame(id = ordered_data$ID, DOSE_MG = ordered_data$DOSE_MG,
                      dv = ordered_data$SEV),
  simulate_trial = function() {
    eta <- rnorm(n_subj_ord, 0, sqrt(omega_ord))[subj_index]
    sample_categories(fitted_probs_ord(ordered_data$DOSE_MG, eta))
  },
  title = "Categorical VPC for ordered categories"
)
