#################################################################################
# ---- Exercises ----------------------------------------------
##  How should the dosing regimen be adujusted to achieve the target          ##
# Which dose gives ≥90% of patients a steady-state average concentration 
# within the target window, and does a weight-based dose help?"               ##
#################################################################################

# 200 children: subject id and body weight WT drawn from a uniform distribution between 10 and 40 (kg)
pedCovs <- data.frame(
  id = 1:200,
  WT = runif(200, 10, 40)
)


# Build one event table for 200 subjects
dose <- et(amt = 10, ii = 24, addl = 10) |>   # 10 mg once daily (every 24 h), 10 additional doses
  et(seq(0, 240, by = 1))  |># hourly observations from 0 to 240 h
  et(id = 1:200) # 200 subjects

# load the saved fit run107 (the model with allometric scaling on WT)
load("run107.Rdata")


# simulate the 200 subjects, the individual covariates (WT) are taken from pedCovs
pedSim <- rxSolve(
  run107,
  dose,
  iCov = pedCovs,
  nSub = nrow(pedCovs),
  returnType = "data.frame"
)


# 5th, 50th (median) and 95th percentile of the simulated cp per time
pedCI <- pedSim |>
  group_by(time) |>
  mutate(p05 = quantile(cp, 0.05),
         p50 = quantile(cp, 0.50),
         p95 = quantile(cp, 0.95))

# median (line) and 5th - 95th percentile (ribbon) of the concentration over time
# (note: the title says 100 mg QD, the dose in the code above is 10 mg)
ggplot(pedCI, aes(x = time)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), fill = "coral", alpha = 0.3) +
  geom_line(aes(y = p50), colour = "coral", linewidth = 1) +
  labs(x = "Time (h)", y = "Concentration (mg/L)",
       title = "Paediatric simulation -- 100 mg QD with allometric scaling",
       subtitle = "WT range 10-40 kg") +
  theme_bw()


# Which dose would be appropriate for children given the same targets?

# 200 new children with newly drawn body weights
nSub <- 200 
pedCovs <- data.frame(
  id = 1:  nSub,
  WT = runif(  nSub, 10, 40)
)


doses  <- c(1,2.5, 5, 7.5, 10)   # doses to compare (mg)
# Build a list of event tables -- one per dose
evList <- lapply(doses, function(d) {
  et(amt = d, ii = 24, addl = 10) |>   # once daily (every 24 h), 10 additional doses
    et(seq(0, 240, by = 1))  |># hourly observations from 0 to 240 h
    et(id = 1:nSub) # 200 subjects
  # the event table is returned
})

# Simulate all doses (one rxSolve call per dose, nSub subjects each, WT from pedCovs)
# and store the dose in the column dose
pedList <- lapply(seq_along(doses), function(i) {#i=1 (value of i to run the code inside the function step by step)
  res <- rxSolve(run107, evList[[i]], nSub = nSub,
                 iCov = pedCovs,
                 returnType = "data.frame")
  res$dose <- doses[i]
  res
})

# combine the simulations of all doses into one data frame
pedAll <- rbindlist(pedList)


# 5th, 50th (median) and 95th percentile of the simulated cp per time and dose
pedCI_all <- pedAll |>
  group_by(time, dose) |>
  mutate(p05 = quantile(cp, 0.05),
         p50 = quantile(cp, 0.50),
         p95 = quantile(cp, 0.95))

# median (line) and 5th - 95th percentile (ribbon) of the concentration over time, per dose
ggplot(pedCI_all, aes(x = time, group = factor(dose), colour = factor(dose),
                      fill = factor(dose)))  +
  geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.15, colour = NA) +
  geom_line(aes(y = p50), linewidth = 1) +
  labs(x = "Time (h)", y = "Concentration (mg/L)",
       colour = "Dose (mg)", fill = "Dose (mg)",
       title = "Paediatric simulation -- Dose comparison -- steady-state concentration profiles", 
       subtitle = "WT range 10-40 kg") +
  theme_bw()

# percentage of subjects per dose with a concentration in the target range (0.5 to 2.5 mg/L) at 240 h
pedtrough <-  pedAll |>
  filter(time == 240) %>%
  group_by(dose) %>%
  summarise(pctAbove = mean(cp >= 0.5 & cp < 2.5 ) * 100
  )

# target attainment per dose
# (note: the y-axis label says ">= 1 mg/L", the code above uses the range 0.5 to 2.5 mg/L)
ggplot(pedtrough, aes(x = factor(dose), y = pctAbove)) +
  geom_col(fill = "steelblue", width = 0.5) +
  geom_hline(yintercept = 90, linetype = "dashed", colour = "red") +
  labs(x = "Dose (mg)", y = "% subjects with trough >= 1 mg/L",
       title = "Target attainment at steady-state trough",
       subtitle = "Dashed line: 90% target attainment criterion") +
  theme_bw()


# ---- 6. would BID (twice daily) dosing help? ------------------------------------


# 200 new children with newly drawn body weights
nSub <- 200 
pedCovs <- data.frame(
  id = 1:  nSub
  )


doses  <- c(1,2.5, 5, 7.5, 10)   # doses to compare (mg), now given every 12 h
# Build a list of event tables -- one per dose
evList <- lapply(doses, function(d) {
  et(amt = d, ii = 12, addl =20) |>   # twice daily (every 12 h), 20 additional doses
    et(seq(0, 240, by = 1))  |># hourly observations from 0 to 240 h
    et(id = 1:nSub)
    et() # 200 subjects
  # the event table is returned
})




# Simulate all doses (one rxSolve call per dose, nSub subjects each, WT from pedCovs)
# and store the dose in the column dose
pedList <- lapply(seq_along(doses), function(i) {#i=1 (value of i to run the code inside the function step by step)
  res <- rxSolve(run107, evList[[i]], nSub = nSub,
                 iCov = pedCovs,
                 returnType = "data.frame")
  res$dose <- doses[i]
  res
})

# combine the simulations of all doses into one data frame
pedAll <- rbindlist(pedList)






# 5th, 50th (median) and 95th percentile of the simulated cp per time and dose
pedCI_all <- pedAll |>
  group_by(time, dose) |>
  mutate(p05 = quantile(cp, 0.05),
         p50 = quantile(cp, 0.50),
         p95 = quantile(cp, 0.95))

# median (line) and 5th - 95th percentile (ribbon) of the concentration over time, per dose
ggplot(pedCI_all, aes(x = time, group = factor(dose), colour = factor(dose),
                      fill = factor(dose)))  +
  geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.15, colour = NA) +
  geom_line(aes(y = p50), linewidth = 1) +
  labs(x = "Time (h)", y = "Concentration (mg/L)",
       colour = "Dose (mg)", fill = "Dose (mg)",
       title = "Paediatric simulation -- Dose comparison -- steady-state concentration profiles", 
       subtitle = "WT range 10-40 kg") +
  theme_bw()

# percentage of subjects per dose with a concentration in the target range (0.5 to 2.5 mg/L) at 240 h
pedtrough <-  pedAll |>
  filter(time == 240) %>%
  group_by(dose) %>%
  summarise(pctAbove = mean(cp >= 0.5 & cp < 2.5 ) * 100
  )

# target attainment per dose
# (note: the y-axis label says ">= 1 mg/L", the code above uses the range 0.5 to 2.5 mg/L)
ggplot(pedtrough, aes(x = factor(dose), y = pctAbove)) +
  geom_col(fill = "steelblue", width = 0.5) +
  geom_hline(yintercept = 90, linetype = "dashed", colour = "red") +
  labs(x = "Dose (mg)", y = "% subjects with trough >= 1 mg/L",
       title = "Target attainment at steady-state trough",
       subtitle = "Dashed line: 90% target attainment criterion") +
  theme_bw()