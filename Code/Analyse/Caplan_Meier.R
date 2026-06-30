# Load required libraries
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(survival)
library(broom)

# Set working directory to project root
setwd(here())

# Load the main data set
data_capm <- readRDS(here("Daten", "Daten_verarbeitet", "data_overweight_cens.rds"))


# Filter to only Study Day 1
data_capm <- subset(data_capm, Study_Day == 1)

# Status coding: 0 = censored, 1 = discharged, 2 = deceased
# Round survival time to whole days
data_capm$surv_0to30 <- round(data_capm$surv_0to30)

# Required packages
library(survival)
library(broom)

# 1) KM for death vs. censored
# Subset for death and censored
# Fit Kaplan-Meier estimator for death
# Prepare data frame for plotting
dat_death     <- subset(data_capm, status_0to30 %in% c(0,2))
km_death      <- survfit(Surv(surv_0to30, status_0to30==2) ~ 1, data = dat_death)
df_death      <- tidy(km_death)[, c("time","estimate","conf.high","conf.low")]
df_death      <- rbind(
  data.frame(time=0, estimate=1, conf.high=1, conf.low=1),
  df_death
)

# 2) KM for discharge vs. censored
# Subset for discharge and censored
# Fit Kaplan-Meier estimator for discharge
# Prepare data frame for plotting
dat_disc      <- subset(data_capm, status_0to30 %in% c(0,1))
km_disc       <- survfit(Surv(surv_0to30, status_0to30==1) ~ 1, data = dat_disc)
df_disc       <- tidy(km_disc)[, c("time","estimate","conf.high","conf.low")]
df_disc       <- rbind(
  data.frame(time=0, estimate=1, conf.high=1, conf.low=1),
  df_disc
)

# 3) Combined plot with geom_step(direction="hv")
p_km <- ggplot() +
  # Death
  geom_step(data = df_death, aes(x = time, y = estimate,     color = "Verstorben"),
            size = 1, direction = "hv") +
  geom_step(data = df_death, aes(x = time, y = conf.low,     color = "Verstorben"),
            linetype = "dashed", direction = "hv") +
  geom_step(data = df_death, aes(x = time, y = conf.high,    color = "Verstorben"),
            linetype = "dashed", direction = "hv") +
  # Discharged
  geom_step(data = df_disc,  aes(x = time, y = estimate,     color = "Entlassen"),
            size = 1, direction = "hv") +
  geom_step(data = df_disc,  aes(x = time, y = conf.low,     color = "Entlassen"),
            linetype = "dashed", direction = "hv") +
  geom_step(data = df_disc,  aes(x = time, y = conf.high,    color = "Entlassen"),
            linetype = "dashed", direction = "hv") +
  # Legend and axes
  scale_color_manual(
    name = "Status\n(95 %-KI)",
    values = c("Entlassen" = "darkblue", "Verstorben" = "red")
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Tag",
    y = "Überlebenswahrscheinlichkeit"
  ) +
  theme_bw() +
  theme(
    text          = element_text(size = 14),
    axis.title    = element_text(size = 14),
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 12)
  )

# Show plot
print(p_km)
