# Run after Modeling_BMI.R; reads saved models without refitting them.
library(mgcv)
library(ggplot2)
library(here)

output_dir <- here("Grafiken", "Subgruppenanalyse", "BMI")
models <- readRDS(file.path(output_dir, "models_bmi.rds"))
reference_protein <- 1
set.seed(123)

# Derive support from one row per patient in the complete-case risk cohort.
d <- readRDS(here("Daten", "Daten_verarbeitet", "data_overweight_cens.rds"))
variables <- c("CombinedID", "surv_0to30", "status_0to30", "BMI", "Age",
               "ApacheIIScore", "OralIntake2_4", "PN2_4", "inMV2_4",
               "Propofol2_4", "Year", "DiagID2", "AdmCatID", "Gender",
               "CombinedicuID", "avDPI")
d <- d[which(d$Study_Day == 1), ]
stopifnot(!anyDuplicated(d$CombinedID))
d <- d[complete.cases(d[, variables]), ]
d <- d[which(is.finite(d$BMI) & d$BMI >= 25 & d$surv_0to30 > 8), ]
d$BMI_group <- ifelse(d$BMI < 30, "Overweight", "Obese")
support <- t(vapply(split(d$avDPI, d$BMI_group), quantile,
                    numeric(2), probs = c(.05, .95)))
lower <- max(support[, 1])
upper <- min(support[, 2])
if (lower >= upper || reference_protein <= lower || reference_protein >= upper) {
  stop("Choose a reference inside the shared 5th-95th percentile range.")
}
grid <- sort(unique(c(seq(lower, upper, length.out = 200), reference_protein)))
write.csv(data.frame(BMI_group = rownames(support), lower = support[, 1],
                     upper = support[, 2]),
          file.path(output_dir, "protein_support_bmi.csv"), row.names = FALSE)

contrast_summary <- function(L, model) {
  estimate <- drop(L %*% coef(model))
  covariance <- vcov(model, unconditional = !is.null(model$Vc))
  se <- sqrt(pmax(0, rowSums((L %*% covariance) * L)))
  data.frame(avDPI = grid, log_estimate = estimate, se = se,
             ratio = exp(estimate), lower = exp(estimate - 1.96 * se),
             upper = exp(estimate + 1.96 * se))
}

save_plot <- function(plot, name) {
  for (extension in c("pdf", "png")) {
    ggsave(file.path(output_dir, paste0(name, ".", extension)), plot,
           width = 8, height = 5, dpi = 300)
  }
}

for (outcome in names(models)) {
  outcome_color <- if (outcome == "discharged") "darkblue" else "red"
  outcome_title <- if (outcome == "discharged") "ICU discharge alive" else "In-ICU death"
  model <- models[[outcome]]
  if (!isTRUE(model$converged)) stop("Model did not converge: ", outcome)
  # Nuisance terms cancel within each patient-profile contrast, including ICU.
  contrast_matrix <- function(group) {
    nd <- model$model[rep(1, length(grid)), , drop = FALSE]
    nd$avDPI <- grid
    nd$BMI_group <- factor(group, levels = model$xlevels$BMI_group)
    ref <- nd
    ref$avDPI <- reference_protein
    predict(model, nd, type = "lpmatrix") -
      predict(model, ref, type = "lpmatrix")
  }
  L_overweight <- contrast_matrix("Overweight")
  L_obese <- contrast_matrix("Obese")
  curves <- rbind(transform(contrast_summary(L_overweight, model), group = "Overweight"),
                  transform(contrast_summary(L_obese, model), group = "Obese"))
  # Ratio of within-group HRs: separates effect modification from group level.
  L_difference <- L_obese - L_overweight
  difference <- contrast_summary(L_difference, model)

  # Approximate simultaneous 95% band from joint Gaussian coefficient draws.
  # Coverage is over this finite grid, separately for each outcome.
  active <- which(colSums(abs(L_difference)) > 1e-12)
  V <- vcov(model, unconditional = !is.null(model$Vc))[active, active, drop = FALSE]
  eig <- eigen(V, symmetric = TRUE)
  root <- sweep(eig$vectors, 2, sqrt(pmax(eig$values, 0)), "*")
  draws <- root %*% matrix(rnorm(length(active) * 5000), nrow = length(active))
  random_difference <- L_difference[, active, drop = FALSE] %*% draws
  nonzero <- difference$se > 1e-10
  max_z <- apply(abs(sweep(random_difference[nonzero, , drop = FALSE],
                           1, difference$se[nonzero], "/")), 2, max)
  critical <- unname(quantile(max_z, .95))
  difference$sim_lower <- exp(difference$log_estimate - critical * difference$se)
  difference$sim_upper <- exp(difference$log_estimate + critical * difference$se)
  difference$reference <- reference_protein
  difference$sim_critical <- critical
  stopifnot(all(abs(difference$ratio[grid == reference_protein] - 1) < 1e-10),
            all(difference$se[grid == reference_protein] < 1e-10))

  write.csv(curves, file.path(output_dir, paste0("protein_HR_", outcome, ".csv")), row.names = FALSE)
  write.csv(difference, file.path(output_dir, paste0("protein_contrast_", outcome, ".csv")), row.names = FALSE)
  x_label <- "Protein intake (g/kg adjusted body weight/day)"
  p <- ggplot(curves, aes(avDPI, ratio, group = group)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2,
                fill = outcome_color, color = NA) +
    geom_line(aes(linetype = group), color = outcome_color) +
    scale_linetype_manual(values = c(Overweight = "solid", Obese = "dashed"),
                          breaks = c("Overweight", "Obese")) +
    theme_bw() + theme(text = element_text(size = 15), legend.position = "bottom") +
    labs(x = x_label, y = "Hazard ratio relative to 1 g/kg/day",
         title = outcome_title, linetype = "BMI group",
         caption = "Pointwise 95% intervals; shared central exposure range")
  save_plot(p, paste0("protein_HR_", outcome))
  p <- ggplot(difference, aes(avDPI, ratio)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    geom_ribbon(aes(ymin = sim_lower, ymax = sim_upper), fill = outcome_color, alpha = .1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = outcome_color, alpha = .2) +
    geom_line(color = outcome_color) + theme_bw() +
    theme(text = element_text(size = 15)) +
    labs(x = x_label, y = "Ratio of protein HRs (obese / overweight)",
         title = outcome_title,
         caption = "Reference: 1 g/kg/day. Dark: pointwise 95%;\nlight: simultaneous 95% band.")
  save_plot(p, paste0("protein_contrast_", outcome))
  message(outcome, ": simultaneous band excludes 1 anywhere: ",
          any(difference$sim_lower > 1 | difference$sim_upper < 1))
}
capture.output(sessionInfo(), file = file.path(output_dir, "contrast_session_info.txt"))
