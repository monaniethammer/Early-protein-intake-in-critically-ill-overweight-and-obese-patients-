# Exploratory BMI effect-modification analysis, supplementary to Modeling.Rmd.
# Run from the project: source("Code/Subgruppenanalyse/Modeling_BMI.R")
# Replace s(BMI) with BMI_group and fit group-specific protein smooths.
# This does not adjust for BMI differences within either group.

library(pammtools)
library(dplyr)
library(mgcv)
library(survival)
library(here)

set.seed(123)

input_file <- here("Daten", "Daten_verarbeitet", "data_overweight_cens.rds")
if (!file.exists(input_file)) {
  stop("Processed data not found: ", input_file)
}
data_original <- readRDS(input_file)
required_columns <- c(
  "Study_Day", "CombinedID", "surv_0to30", "status_0to30", "BMI",
  "Age", "ApacheIIScore", "OralIntake2_4", "PN2_4", "inMV2_4",
  "Propofol2_4", "Year", "DiagID2", "AdmCatID", "Gender",
  "CombinedicuID", "avDPI"
)
missing_columns <- setdiff(required_columns, names(data_original))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

data_bmi <- data_original %>%
  filter(Study_Day == 1)
if (anyNA(data_bmi$CombinedID) || anyDuplicated(data_bmi$CombinedID)) {
  stop("Expected one baseline row per patient with a non-missing CombinedID.")
}
n_baseline <- nrow(data_bmi)
data_bmi <- data_bmi %>%
  filter(is.finite(BMI), BMI >= 25)
n_eligible <- nrow(data_bmi)

# Full baseline cohort, before complete-case exclusions and the day-8 restriction.
# Protein is the prescribed target at admission according to preproc-data.R.
# Its g/kg/day unit and weight denominator need source-dictionary confirmation.
# Protein is descriptive only and is not included in model complete-case filtering.
if (!"Protein" %in% names(data_bmi) || !is.numeric(data_bmi$Protein)) {
  stop("Expected numeric Protein column for prescribed-target summaries.")
}
baseline_descriptives <- data_bmi %>%
  mutate(BMI_group = ifelse(BMI < 30, "Overweight", "Obese"))
summarise_baseline <- function(x, group) {
  bind_rows(lapply(c("BMI", "Protein", "avDPI"), function(variable) {
    values <- x[[variable]]
    observed <- values[is.finite(values)]
    quartiles <- if (length(observed)) {
      unname(quantile(observed, probs = c(.25, .5, .75)))
    } else {
      rep(NA_real_, 3)
    }
    data.frame(
      BMI_group = group, variable = variable, patients = nrow(x),
      observed = length(observed), missing = sum(is.na(values)),
      non_finite = sum(!is.finite(values) & !is.na(values)),
      q1 = quartiles[1], median = quartiles[2], q3 = quartiles[3]
    )
  }))
}
baseline_summary_bmi <- bind_rows(
  summarise_baseline(baseline_descriptives, "Overall"),
  summarise_baseline(filter(baseline_descriptives, BMI_group == "Overweight"), "Overweight"),
  summarise_baseline(filter(baseline_descriptives, BMI_group == "Obese"), "Obese")
)
print(baseline_summary_bmi)
output_dir <- here("Grafiken", "Subgruppenanalyse", "BMI")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(baseline_summary_bmi,
          file.path(output_dir, "baseline_bmi_protein_target_summary.csv"),
          row.names = FALSE)

# Publication table; preserve the numeric summary above for further analysis.
baseline_table_bmi <- data.frame(Characteristic = c(
  "Patients, n", "BMI (kg/m^2)", "Prescribed protein target (g/kg/day)*",
  "Average delivered protein intake, days 1-4 (g/kg adjusted body weight/day)"
), check.names = FALSE)
for (group in c("Overall", "Overweight", "Obese")) {
  rows <- baseline_summary_bmi[baseline_summary_bmi$BMI_group == group, ]
  rows <- rows[match(c("BMI", "Protein", "avDPI"), rows$variable), ]
  baseline_table_bmi[[group]] <- c(
    format(rows$patients[1], big.mark = ",", trim = TRUE),
    sprintf("%.2f [%.2f-%.2f]", rows$median, rows$q1, rows$q3)
  )
}
table_notes_bmi <- c(
  "Values are median [25th-75th percentile] unless stated otherwise.",
  "Overweight: 25 <= BMI < 30 kg/m^2; obese: BMI >= 30 kg/m^2.",
  "One baseline row per patient; before model exclusions and the day-8 restriction.",
  "avDPI is each patient's mean recorded protein intake over days 1-4, summarized across patients.",
  "* Prescribed target at admission, not delivered intake or measured physiological requirement.",
  "The g/kg/day unit and weight denominator require source-data-dictionary confirmation."
)
print(baseline_table_bmi, row.names = FALSE)
write.csv(baseline_table_bmi,
          file.path(output_dir, "baseline_bmi_protein_target_table.csv"), row.names = FALSE)
writeLines(c(
  "# BMI, prescribed protein target, and early delivered protein intake", "",
  paste0("| ", paste(names(baseline_table_bmi), collapse = " | "), " |"),
  "| --- | ---: | ---: | ---: |",
  apply(baseline_table_bmi, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |")),
  "", table_notes_bmi
), file.path(output_dir, "baseline_bmi_protein_target_table.md"))

# Use an explicit complete-case cohort for all model terms before PED expansion.
data_bmi <- data_bmi[complete.cases(data_bmi[, required_columns]), ]
data_bmi <- data_bmi %>%
  mutate(
    BMI_group = factor(ifelse(BMI < 30, "Overweight", "Obese"),
                       levels = c("Overweight", "Obese")),
    CombinedicuID = factor(CombinedicuID)
  ) %>%
  droplevels()
message("Baseline patients: ", n_baseline,
        "; excluded for missing/non-finite BMI or BMI < 25: ", n_baseline - n_eligible,
        "; excluded for missing model variables: ", n_eligible - nrow(data_bmi))
if (!is.numeric(data_bmi$status_0to30) ||
    !all(data_bmi$status_0to30 %in% 0:2) ||
    !is.numeric(data_bmi$surv_0to30) ||
    any(!is.finite(data_bmi$surv_0to30) |
        data_bmi$surv_0to30 <= 0 | data_bmi$surv_0to30 > 30)) {
  stop("Expected numeric status 0/1/2 and survival times in (0, 30].")
}

# Same daily intervals and lag restriction as the main analysis.
# The other cause terminates follow-up in each cause-specific model.
ped_bmi <- as_ped(
  data = data_bmi,
  formula = Surv(surv_0to30, status_0to30) ~ BMI_group + Age +
    ApacheIIScore + OralIntake2_4 + PN2_4 + inMV2_4 + Propofol2_4 +
    Year + DiagID2 + AdmCatID + Gender + CombinedicuID + avDPI,
  cut = 0:30, id = "CombinedID", combine = FALSE
)
if (length(ped_bmi) != 2L) {
  stop("Expected both discharge (1) and death (2) in the competing-risk data.")
}
ped_bmi <- lapply(ped_bmi, function(x) x[x$tend >= 9, ])
names(ped_bmi) <- c("discharged", "died")

model_formula_bmi <- ped_status ~ s(tend, bs = "ps") +
  s(Age, bs = "ps") + BMI_group + ApacheIIScore +
  OralIntake2_4 + PN2_4 + inMV2_4 + Propofol2_4 + Year +
  DiagID2 + AdmCatID + Gender + s(CombinedicuID, bs = "re") +
  s(avDPI, bs = "ps", by = BMI_group)
# The unordered factor creates two protein smooths. Its main effect is needed
# because factor-by smooths are centered. Covariates follow Modeling.Rmd;
# avDCI is not an adjustment term in that model either.

counts_bmi <- bind_rows(lapply(names(ped_bmi), function(outcome) {
  x <- ped_bmi[[outcome]]
  counts <- x %>%
    group_by(BMI_group, .drop = FALSE) %>%
    summarise(patients = n_distinct(CombinedID), events = sum(ped_status),
              person_days = sum(exp(offset)), .groups = "drop")
  if (any(counts$patients == 0 | counts$events == 0)) {
    stop("Both BMI groups must contribute patients and events for ", outcome)
  }
  mutate(counts, outcome = outcome)
}))
print(counts_bmi)

models_bmi <- lapply(ped_bmi, function(x) {
  bam(model_formula_bmi, data = x, family = poisson(), offset = offset,
      method = "fREML", discrete = TRUE, na.action = na.fail)
})
model_discharged_bmi <- models_bmi$discharged
model_died_bmi <- models_bmi$died

output_dir <- here("Grafiken", "Subgruppenanalyse", "BMI")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(models_bmi, file.path(output_dir, "models_bmi.rds"))
write.csv(counts_bmi, file.path(output_dir, "analysis_counts_bmi.csv"),
          row.names = FALSE)

diagnostics_bmi <- lapply(models_bmi, function(model) {
  list(converged = model$converged, rank = model$rank,
       coefficients = length(coef(model)),
       concurvity = concurvity(model, full = TRUE),
       basis_check = k.check(model))
})
for (outcome in names(models_bmi)) {
  model <- models_bmi[[outcome]]
  print(summary(model))
  print(diagnostics_bmi[[outcome]])
  if (!isTRUE(model$converged) || model$rank < length(coef(model))) {
    warning("Check convergence/rank for BMI model: ", outcome)
  }
}
capture.output(
  lapply(models_bmi, summary), diagnostics_bmi, sessionInfo(),
  file = file.path(output_dir, "model_summaries_bmi.txt")
)

# Centered partial effects on the log-hazard scale, with pointwise intervals.
# Smooth-term p-values are not tests of differences between the BMI groups.
save_protein_plot <- function(model, outcome, format) {
  filename <- file.path(output_dir,
                        paste0("protein_smooths_", outcome, "_bmi.", format))
  protein_terms <- which(vapply(model$smooth, function(x) {
    "avDPI" %in% x$term
  }, logical(1)))
  curves <- do.call(rbind, lapply(protein_terms, function(term) {
    smooth <- model$smooth[[term]]
    group <- smooth$by.level
    observed <- model$model$avDPI[model$model$BMI_group == group]
    nd <- model$model[rep(1, 100), , drop = FALSE]
    nd$avDPI <- seq(min(observed), max(observed), length.out = 100)
    nd$BMI_group <- factor(group, levels = model$xlevels$BMI_group)
    prediction <- predict(model, newdata = nd, type = "terms",
                          terms = smooth$label, se.fit = TRUE)
    data.frame(avDPI = nd$avDPI, group = group,
               fit = drop(prediction$fit),
               lower = drop(prediction$fit - 1.96 * prediction$se.fit),
               upper = drop(prediction$fit + 1.96 * prediction$se.fit))
  }))
  curves$group <- factor(curves$group, levels = c("Overweight", "Obese"))
  outcome_color <- if (outcome == "discharged") "darkblue" else "red"
  p <- ggplot2::ggplot(curves, ggplot2::aes(avDPI, fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = outcome_color, alpha = .2) +
    ggplot2::geom_line(color = outcome_color) +
    ggplot2::facet_wrap(~ group, nrow = 1) +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 15)) +
    ggplot2::labs(x = "Protein intake (g/kg adjusted body weight/day)",
                  y = "Centered partial effect on log hazard",
                  title = if (outcome == "discharged") "ICU discharge alive" else "In-ICU death")
  ggplot2::ggsave(filename, p, width = 10, height = 5, dpi = 300)
}
for (format in c("pdf", "png")) {
  invisible(Map(function(model, outcome) {
    save_protein_plot(model, outcome, format)
  }, models_bmi, names(models_bmi)))
}
message("BMI analysis saved to: ", output_dir)
