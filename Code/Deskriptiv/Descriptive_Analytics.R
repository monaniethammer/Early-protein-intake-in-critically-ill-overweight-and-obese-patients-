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
data_overweight_cens <- readRDS(here("Daten", "Daten_verarbeitet", "data_overweight_cens.rds"))

###### OUTCOME DISTRIBUTION (Überlebensstatus 0-30 Tage) ######
# This plot shows the distribution of patient outcomes (censored, discharged, deceased)
# Only the first day per patient is used, as the outcome does not change over time

data_overweight_unique <- data_overweight_cens %>% filter(Study_Day == 1)

ggplot(data_overweight_unique, aes(x = factor(status_0to30), fill = factor(status_0to30))) +
  geom_bar() +
  scale_fill_manual(
    values = c("gray", "darkblue", "red"),
    labels = c("0 = Zensiert", "1 = Entlassen", "2 = Verstorben")
  ) +
  labs(
    title = "Überlebensstatus (0-30 Tage)",
    x = "",
    y = "Absolute Verteilung",
    fill = "Status"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
###### ZEITABHÄNGIGE CONFOUNDER (Grouped Lollipop Chart) ######
# This section calculates which variables change within patients over time
# and which are used as confounders in the model. The result is visualized as a grouped lollipop chart.

# Function to detect which variables change within each patient
# Returns a summary data frame with the percentage of patients for which each variable changes
get_variable_change_summary <- function(data, id_column = "CombinedID") {
  patient_ids <- unique(data[[id_column]])
  summary_counts <- data.frame(
    variable = colnames(data),
    patients_with_changes = 0,
    stringsAsFactors = FALSE
  )
  for (pid in patient_ids) {   # This for loop has been created with help of GitHub Copilot
    patient_data <- data[data[[id_column]] == pid, ]
    if (nrow(patient_data) < 2) next
    for (col in colnames(patient_data)) {
      if (length(unique(patient_data[[col]])) > 1) {
        summary_counts[summary_counts$variable == col, "patients_with_changes"] <- 
          summary_counts[summary_counts$variable == col, "patients_with_changes"] + 1
      }
    }
  }
  summary_counts$percentage <- (summary_counts$patients_with_changes / length(patient_ids)) * 100
  summary_counts <- summary_counts[order(-summary_counts$percentage), ]
  return(summary_counts)
}

# Calculate the change summary for all variables
change_summary <- get_variable_change_summary(data_overweight_cens)

# Define which variables are confounders (from the model)
confounder_vars <- c(
  "Age", "BMI", "ApacheIIScore", "OralIntake2_4", "PN2_4", "inMV2_4", "Propofol2_4", "Year",
  "DiagID2", "AdmCatID", "Gender", "CombinedicuID", "avDPI", "avDCI"
)

# Build the variable_info data frame for plotting
variable_info <- change_summary %>%
  filter(variable %in% c(
    # All variables of interest (changing or confounder)
    "Study_Day", "calproKg", "caloriesPercentage", "caloriesIntake", "proteinGproKG", "proteinAdjustedPercentage", "calorie_group", "EN_Protein", "EN", "incomplete_day", "OralIntake", "PN_Protein", "PN",
    "Age", "BMI", "ApacheIIScore", "OralIntake2_4", "PN2_4", "inMV2_4", "Propofol2_4", "Year", "DiagID2", "AdmCatID", "Gender", "CombinedicuID", "avDPI", "avDCI"
  )) %>%
  mutate(
    changing = percentage > 0,
    confounder = variable %in% confounder_vars
  )

# Prepare for grouped lollipop plot
variable_info$id <- 1:nrow(variable_info)
variable_info_grouped <- variable_info %>%
  mutate(group = ifelse(confounder, "Confounder", "Kein Confounder")) %>%
  arrange(group, -percentage)
variable_info_grouped$variable <- factor(variable_info_grouped$variable, 
                                         levels = variable_info_grouped$variable)

# Plot: Grouped lollipop chart showing the percentage of patients with changes per variable
# Colored by confounder status
ggplot(variable_info_grouped, aes(x = variable, y = percentage, color = group)) +
  geom_segment(aes(x = variable, xend = variable, y = 0, yend = percentage),
               size = 1.5, alpha = 0.7) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Confounder" = "#F8766D", "Kein Confounder" = "#619CFF")) +
  coord_flip() +
  labs(title = "Änderungsrate der Variablen",
       subtitle = "Gruppiert nach Confounder-Status",
       x = "Variable",
       y = "Prozent der Patienten mit Änderungen",
       color = "Kategorie") +
  theme_minimal() +
  theme(legend.position = "top")

###### KORRELATIONSPLOT DER CONFOUNDER ######
# This plot shows the correlation matrix of the main confounder variables (first day per patient)
data_overweight_one <- data_overweight_cens[data_overweight_cens$Study_Day == 1, ]
confounder <- data_overweight_one %>%
  mutate(Year = as.numeric(as.character(Year))) %>%
  mutate(Gender = ifelse(Gender == "Male", 0, 1)) %>%
  select(OralIntake2_4, PN2_4, inMV2_4, Propofol2_4, Gender, Year, ApacheIIScore, Age, BMI, avDPI)
cor_matrix <- cor(confounder)
# Save as PNG so it can be viewed in any environment
png("Grafiken/Julian/Correlation_Plot.png", width = 900, height = 900, res = 150)
corrplot::corrplot(
  cor_matrix,
  method = "square",
  type = "lower",
  order = "original",
  col = colorRampPalette(c("blue", "white", "red"))(200),
  tl.col = "black",
  tl.cex = 0.9,
  addCoef.col = "black",
  number.cex = 0.7,
  title = "Korrelationsmatrix der Confounder-Variablen",
  mar = c(0,0,2,0),
  addgrid.col = "lightgray"
)
dev.off()

###### DIAGNOSE VS. AUFNAHMEKATEGORIE (HEATMAP) ######
# This heatmap shows the distribution of patients by diagnosis and admission category (first day per patient)
data_overweight_one <- data_overweight_cens[data_overweight_cens$Study_Day == 1, ]
table_plot <- table(data_overweight_one$DiagID2, data_overweight_one$AdmCatID)
df_melt <- melt(table_plot)
ggplot(df_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Diagnose", y = "Aufnahmekategorie", fill = "Anzahl",
       title = "Verteilung der Patienten nach Diagnose und Aufnahmekategorie") +
  theme_minimal()
