library(dplyr)
library(here)

complete_data <- readRDS(here("Daten", "Daten_roh", "mergedAndCleanedData.Rds"))
### 1. Compute avDPI ###

# filter all observations with BMI >= 25 (Übergewicht oder adipös)
data_subset <- complete_data[complete_data$BMI >= 25, ] # 100452 observations
data_overweight <- data_subset

# calculate needed variables and daily protein intake (DPI), daily calorien intake (DCI)
data_subset$h <- 100 * sqrt(data_subset$Weight / data_subset$BMI)
data_subset$iBW <- ifelse(data_subset$Gender == "Male", 
                          0.9 * data_subset$h - 100,  # Formel für Männer
                          0.9 * data_subset$h - 106)  # Formel für Frauen
data_subset$aBW <- (data_subset$Weight - data_subset$iBW) * 0.25 + data_subset$iBW
data_subset$DPI <- (data_subset$proteinGproKG * data_subset$Weight) / data_subset$aBW
data_subset$DCI <- (data_subset$calproKg * data_subset$Weight) / data_subset$aBW

# filter days 1-4 (time frame of interest)
data_subset <- data_subset %>%
  filter(Study_Day %in% c(1, 2, 3, 4))

# calculate average daily protein intake (avDPI), create dataset with one row per patient
# some variables are added for identification
data_patient_avDPI <- data_subset %>%
  group_by(CombinedID) %>%
  summarise(
    avDPI = mean(DPI, na.rm = TRUE),  # Durchschnitt der Proteinaufnahme
  )

# calculate average daily Calorie intake (avDCI), create dataset with one row per patient
# some variables are added for identification
data_patient_avDCI <- data_subset %>%
  group_by(CombinedID) %>%
  summarise(
    avDCI = mean(DCI, na.rm = TRUE),  # Durchschnitt der Kalorienaufnahme
  )

# prepare data for merging
data_overweight <- data_overweight %>%
  select(-Propofol) # remove variables (optional)

data_overweight <- data_overweight %>%
  left_join(data_patient_avDPI, by = c("CombinedID"))

data_overweight <- data_overweight %>%
  left_join(data_patient_avDCI, by = c("CombinedID"))

# Save
saveRDS(data_overweight, file.path(here("Daten", "Daten_verarbeitet"), "data_overweight.rds"))


### 2. CENSOR BEFORE DAY 9 AND AFTER DAY 30 ###

# We are only interested at in-hospital death or discharge up to 30 days after ICU admission.
# If surv_icu0to60 is >30, the observations are censored.
# 5 cases: already censored after 60 days, death within 30 days, discharge within 30 days, death/discharge before 60 and after 30 days (stay on ICU >30, censored)

# standard value for the status is 0, for surv_0to30 it's 31 (covers all observations which are already censored after 60 days)
data_overweight$status_0to30 <- 0
data_overweight$surv_0to30 <- 31

# All patients with surv_icu_status = 2 have died on ICU within 60 days. All patients with surv_icu_status = 1 have been discharged within 60 days.
# If those patients surv_icu0to60 (corresponding days until death/discharge) is < 31, they are considered as dead/discharged in this study. Else they are censored (status: 0, surv_0to30: 31).
data_overweight <- data_overweight %>%
  mutate(
    status_0to30 = ifelse(surv_icu_status %in% c(1, 2) & surv_icu0to60 <= 30, surv_icu_status, 0),
    surv_0to30 = ifelse(surv_icu_status %in% c(1, 2) & surv_icu0to60 <= 30, surv_icu0to60, 31)
  )%>%
  mutate(
    status_0to30 = ifelse(surv_icu0to60 < 9, 0, status_0to30)
  )%>%
  select(-surv_icu0to60, -surv_icu_status) # remove old information about censoring
# note: event != surv_0to60 if DaysInICU is smaller

# check whether everything works
test_data <- data_overweight %>%
  filter(surv_0to30 != event) %>%
  arrange(surv_0to30)%>%
  select(event, DaysInICU, surv_0to30, status_0to30) %>%
  distinct()


### 3. PROCESS CALORIE INTAKE ###
# Info: calorie category is moved one cat up for oral intake
# categorisation of preproc-data is used (after this procedure)
# calCat: 0 means < 30, 1 means 30 <= caloriesPercentage < 70, 2 means 70 < caloriesPercentage
# if calCat is used: when oral intake then patient moved one categorie higher
# if caloriesCat would be used: original categories by daily caloriesPercentage

data_overweight["calorie_group"] = (data_overweight["calCat3"] * 2) + data_overweight["calCat2"]

# remove variables not needed
data_overweight <- data_overweight %>%
  select(-calCat2, -calCat3, -caloriesCat2,-caloriesCat3, -protCat2, -protCat3, -proteinCat2, -proteinCat3, -Surv0To60, -Disc0To60)

### 4. ANALYSE DATASET ###

# number of censored ids --> 4780
censored <- length(unique(data_overweight$CombinedID[data_overweight$status_0to30 == 0]))

# number of discharged ids --> 3448
discharged <- length(unique(data_overweight$CombinedID[data_overweight$status_0to30 == 1]))

# number of died ids --> 840
died <- length(unique(data_overweight$CombinedID[data_overweight$status_0to30 == 2]))

# save data in data folder
# saveRDS(data_overweight, file.path(here("Daten", "Daten_verarbeitet"), "data_overweight_cens.rds"))
