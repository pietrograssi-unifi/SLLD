# ==============================================================================
# TITLE: Decoding Secondary Lymphomas: Integrating Macro-Epidemiological Models 
#        and Micro-Genomic Signatures to Disentangle Mutagenesis from Misdiagnosis.
# AUTHORS: Pietro Grassi and Elena Bariletti
# METHODOLOGY: Multiple Imputation (MICE), Factor Analysis of Mixed Data (FAMD), 
#              Fine-Gray Competing Risk Modeling, and Restricted Cubic Splines.
# ==============================================================================

rm(list=ls())
gc()
setwd("/Users/pietro/Desktop/Sant'Anna/SLLD/Progetto")

# ==============================================================================
# PHASE 0: DEPENDENCIES & ENVIRONMENT SETUP
# Rationale: Loading robust libraries for missing data imputation, competing risk 
# survival analysis, and unsupervised learning.
# ==============================================================================

library(data.table)  # High-performance data manipulation
library(mice)        # Multiple Imputation by Chained Equations
library(survival)    # Core survival analysis architecture
library(tidycmprsk)  # Fine-Gray models for competing risks
library(ggsurvfit)   # Cumulative Incidence Function (CIF) visualization
library(cluster)     # Core clustering algorithms
library(factoextra)  # Silhouette statistics and visualization
library(FactoMineR)  # Factor Analysis of Mixed Data (FAMD)
library(lmtest)      # Likelihood Ratio Tests
library(gtsummary)   # Publication-ready clinical tables
library(broom)       # Tidying model outputs
library(ggplot2)     # High-fidelity visualizations
library(dplyr)       # Data manipulation workflows
library(splines)     # Natural (Restricted) Cubic Splines for non-linear modeling

# Ensuring strict computational reproducibility
set.seed(2026)       

# ==============================================================================
# PHASE 1: DATA INGESTION & COMPETING RISK COHORT DEFINITION
# Rationale: Parsing the SEER registry to map patient longitudinal trajectories. 
# "Vital Status" is extracted to define Death as a competing risk to secondary NHL.
# ==============================================================================
cat("\n[PHASE 1] Initializing SEER registry data ingestion...\n")
dt <- fread("TableExport.txt", na.strings = c("NA", "Blank(s)", "Unknown"))

dt[, hist_code := `Histologic Type ICD-O-3`]

# STEP 1A: Define Baseline Cohort (Primary HL strictly as first primary malignancy)
hl_baseline <- dt[(hist_code >= 9650 & hist_code <= 9667) & 
                    (`Sequence number` == "One primary only" | `Sequence number` %like% "1st")]
id_hl_patients <- unique(hl_baseline$`Patient ID`)

# STEP 1B: Extract Longitudinal Trajectories
dt_cohort <- dt[`Patient ID` %in% id_hl_patients]

# STEP 1C: Define Target Outcome and Competing Risk
id_subsequent_nhl <- dt_cohort[
  `Patient ID` %in% id_hl_patients & 
    (`Sequence number` != "One primary only" & !(`Sequence number` %like% "1st")) & 
    ((hist_code >= 9590 & hist_code <= 9642) | (hist_code >= 9670 & hist_code <= 9738)), 
  unique(`Patient ID`)
]

df_model <- hl_baseline
df_model[, Develops_NHL := ifelse(`Patient ID` %in% id_subsequent_nhl, 1, 0)]
df_model[, Vital_Status := `Vital status recode (study cutoff used)`]

# Event Status: 0 = Censored, 1 = Secondary NHL (Event of Interest), 2 = Death (Competing Risk)
df_model[, Event_Status := dplyr::case_when(
  Develops_NHL == 1 ~ 1,
  Develops_NHL == 0 & Vital_Status == "Dead" ~ 2,
  TRUE ~ 0
)]

cat("Total patients with primary HL identified:", nrow(df_model), "\n")
cat("Secondary NHL Events (Target):", sum(df_model$Event_Status == 1), "\n")
cat("Death Events (Competing Risk):", sum(df_model$Event_Status == 2), "\n")

# ==============================================================================
# PHASE 2: FEATURE ENGINEERING & ROBUST MULTIPLE IMPUTATION (MICE)
# Rationale: Using the Nelson-Aalen estimator (H0) to prevent survival bias 
# and treating Event_Status as a factor to correctly guide PMM.
# ==============================================================================
cat("\n[PHASE 2] Executing Feature Engineering and MICE...\n")

df_model[, Age_continuous := as.numeric(gsub("[^0-9]", "", `Age recode with single ages and 85+`))]
df_model[, FollowUp_months := suppressWarnings(as.numeric(as.character(`Survival months`)))]
df_model <- df_model[!is.na(FollowUp_months)]
df_model[, Chemo_flag := factor(ifelse(`Chemotherapy recode (yes, no/unk)` == "Yes", "Treated", "Untreated"), levels = c("Untreated", "Treated"))]
df_model[, Radio_flag := factor(ifelse(`Radiation recode` %like% "Beam" | `Radiation recode` %like% "Radio", "Treated", "Untreated"), levels = c("Untreated", "Treated"))]
df_model[, Stage_advanced := as.factor(ifelse(`Combined Summary Stage with Expanded Regional Codes (2004+)` %like% "Distant", 1, 0))]
df_model[, Sex := as.factor(Sex)]
df_model[, Race := as.factor(`Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)`)]

# Preventing Survival Bias (Nelson-Aalen Estimator)
df_model[, Any_Event := ifelse(Event_Status != 0, 1, 0)]
df_model[, H0 := mice::nelsonaalen(df_model, timevar = "FollowUp_months", statusvar = "Any_Event")]
df_model[, Event_Status_Factor := as.factor(Event_Status)]
model_vars <- c("Event_Status_Factor", "H0", "Age_continuous", "Chemo_flag", "Radio_flag", "Stage_advanced", "Sex", "Race")
df_subset <- df_model[, ..model_vars]

# Multiple Imputation (m=5)
imputed_data <- mice(df_subset, m = 5, method = "pmm", maxit = 5, seed = 2026, printFlag = FALSE)

# Extraction for geometry and clustering
df_imputed <- complete(imputed_data, 1)
setDT(df_imputed)

df_imputed[, Event_Status := as.numeric(as.character(Event_Status_Factor))]
df_imputed[, FollowUp_months := df_model$FollowUp_months]

# ==============================================================================
# PHASE 3: UNSUPERVISED LATENT PHENOTYPING
# Rationale: Discovering unbiased clinical phenotypes in a valid orthogonal space 
# via FAMD. The optimal number of clusters is mathematically derived using the 
# Silhouette statistic to avoid heuristic bias.
# ==============================================================================
cat("\n[PHASE 3] Computing FAMD and Silhouette Optimization...\n")

cluster_features <- c("Age_continuous", "Chemo_flag", "Radio_flag", "Stage_advanced", "Sex", "Race")
df_cluster <- df_imputed[, ..cluster_features]

res.famd <- FAMD(df_cluster, ncp = 3, graph = FALSE)

# Silhouette computation for dynamic K determination
set.seed(2026)
samp_idx <- sample(1:nrow(res.famd$ind$coord), min(2000, nrow(res.famd$ind$coord)))
sil_plot <- fviz_nbclust(res.famd$ind$coord[samp_idx, ], kmeans, method = "silhouette", k.max = 20)

optimal_k <- as.numeric(sil_plot$data$clusters[which.max(sil_plot$data$y)])
cat("Optimal number of latent clinical phenotypes (K):", optimal_k, "\n")

km_res <- cluster::clara(res.famd$ind$coord, k = optimal_k, samples = 50, sampsize = 1000, pamLike = TRUE)
df_imputed[, Phenotype_Cluster := as.factor(km_res$clustering)]

medoid_coords <- km_res$medoids

# ------------------------------------------------------------------------------
# METHODOLOGICAL NOTE: MULTIPLE IMPUTATION PROJECTION (REFERENCE FRAME)
# To preserve Rubin's between-imputation variance (B) without succumbing to the 
# Label Switching Problem typical of independent clustering across m datasets, 
# we establish the m=1 space as a fixed topological Reference Frame. 
# Patients from m=2 to m=5 are subsequently projected into this space. Their 
# cluster assignments will naturally fluctuate depending on their imputed covariates, 
# accurately translating missing-data uncertainty into inflated Standard Errors 
# in the downstream survival models.
# ------------------------------------------------------------------------------

long_data <- mice::complete(imputed_data, action = "long", include = TRUE)
long_data$Phenotype_Cluster <- NA

for (m_idx in 1:5) {
  tmp_df <- mice::complete(imputed_data, m_idx)
  tmp_features <- as.data.frame(tmp_df)[, cluster_features]
  
  tmp_coord <- predict(res.famd, newdata = tmp_features)$coord
  
  closest_cluster <- apply(tmp_coord, 1, function(row) {
    which.min(colSums((t(medoid_coords) - row)^2))
  })
  
  long_data$Phenotype_Cluster[long_data$.imp == m_idx] <- closest_cluster
}

long_data$Phenotype_Cluster[long_data$.imp == 0] <- long_data$Phenotype_Cluster[long_data$.imp == 1]
long_data$Phenotype_Cluster <- as.factor(long_data$Phenotype_Cluster)
long_data$FollowUp_months <- rep(df_model$FollowUp_months, times = 6) 
long_data$Event_Status <- as.numeric(as.character(long_data$Event_Status_Factor))

imputed_data <- as.mids(long_data)

# ==============================================================================
# PHASE 4: MULTIPLE IMPUTATION POOLING (RUBIN'S RULES)
# Rationale: Fitting both Cause-Specific Cox and Fine-Gray across m=5 datasets.
# Note: Phenotype_Cluster is excluded to prevent perfect multicollinearity 
# (infinite coefficient warning) since it is a linear combination of the other variables.
# ==============================================================================
cat("\n[PHASE 4] Fitting Pooled Survival Models (CS-Cox and Fine-Gray)...\n")

# Cause-Specific Cox
cs_models <- with(imputed_data, coxph(Surv(FollowUp_months, Event_Status == 1) ~ 
                                        Chemo_flag * ns(Age_continuous, df = 3) + 
                                        Radio_flag + Sex + Race + Stage_advanced))
pooled_cs <- pool(cs_models)

# Fine-Gray Subdistribution Hazards 
fg_models_list <- lapply(1:5, function(i) {
  tmp_data <- mice::complete(imputed_data, i)
  tidycmprsk::crr(Surv(FollowUp_months, Event_Status_Factor) ~ 
                    Chemo_flag * ns(Age_continuous, df = 3) + 
                    Radio_flag + Sex + Race + Stage_advanced, 
                  data = tmp_data)
})
pooled_fg <- pool(fg_models_list)

# ==============================================================================
# PHASE 4.1: ROBUST SCHOENFELD DIAGNOSTICS (ACROSS ALL IMPUTATIONS)
# ==============================================================================
cat("\n[PHASE 4.1] Testing Proportional Hazards Assumption (Robust CS-Cox)...\n")

zph_pvals <- list()

for(i in 1:5) {
  tmp_data <- mice::complete(imputed_data, i)
  tmp_cox <- coxph(Surv(FollowUp_months, Event_Status == 1) ~ 
                     Chemo_flag * ns(Age_continuous, df = 3) + 
                     Radio_flag + Sex + Race + Stage_advanced, 
                   data = tmp_data)
  
  zph_test_tmp <- cox.zph(tmp_cox)
  zph_pvals[[i]] <- zph_test_tmp$table[, "p"]
}

median_zph_pvals <- apply(do.call(rbind, zph_pvals), 2, median)
cat("\n--- Median Schoenfeld p-values across m=5 imputations ---\n")
print(round(median_zph_pvals, 4))

# ==============================================================================
# PHASE 4.2: MODEL CALIBRATION (10-YEAR COMPETING RISK CALIBRATION PLOT)
# ==============================================================================
cat("\n[PHASE 4.2] Executing Robust 10-Year Calibration Analysis...\n")

calib_results <- list()

for(m in 1:5) {
  tmp_data <- mice::complete(imputed_data, m)
  
  tmp_fg <- tidycmprsk::crr(Surv(FollowUp_months, Event_Status_Factor) ~ 
                              Chemo_flag * ns(Age_continuous, df = 3) + 
                              Radio_flag + Sex + Race + Stage_advanced, 
                            data = tmp_data)
  
  tmp_data$pred_risk_120m <- predict(tmp_fg, times = 120, newdata = tmp_data)[[1]]
  tmp_data$Risk_Decile <- ntile(tmp_data$pred_risk_120m, 10)
  
  tmp_calib <- data.frame(Decile = 1:10, Predicted_Mean = NA, Observed_CIF = NA, Imputation = m)
  
  for(i in 1:10) {
    subset_df <- tmp_data[tmp_data$Risk_Decile == i, ]
    tmp_calib$Predicted_Mean[i] <- mean(subset_df$pred_risk_120m, na.rm = TRUE)
    
    fit_aj <- survfit(Surv(FollowUp_months, Event_Status_Factor) ~ 1, data = subset_df)
    idx <- which(fit_aj$time <= 120)
    if(length(idx) > 0) {
      # Accertiamoci di pescare la probabilità dell'evento "1" (Secondario NHL)
      tmp_calib$Observed_CIF[i] <- fit_aj$pstate[tail(idx, 1), "1"] 
    } else {
      tmp_calib$Observed_CIF[i] <- 0
    }
  }
  calib_results[[m]] <- tmp_calib
}

all_calib <- do.call(rbind, calib_results)
robust_calibration <- all_calib %>%
  group_by(Decile) %>%
  summarize(
    Predicted_Mean_Robust = mean(Predicted_Mean),
    Observed_CIF_Robust = mean(Observed_CIF),
    .groups = 'drop'
  )

p_calib <- ggplot(robust_calibration, aes(x = Predicted_Mean_Robust, y = Observed_CIF_Robust)) +
  geom_abline(slope = 1, intercept = 0, color = "darkgray", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "#2c3e50", linetype = "solid") +
  geom_point(size = 4, color = "#e74c3c") +
  theme_minimal(base_size = 14) +
  labs(title = "Robust 10-Year Model Calibration (Pooled Fine-Gray)",
       subtitle = "Averaged across m=5 MICE datasets",
       x = "Robust Predicted 10-Year Probability",
       y = "Robust Observed 10-Year Cumulative Incidence")

# ==============================================================================
# APPENDIX: EPIDEMIOLOGICAL SENSITIVITY ANALYSES (E-VALUE)
# Rationale: Quantifying the minimum strength of unmeasured confounding required 
# to nullify the observed neutral effect of chemotherapy.
# ==============================================================================
cat("\n[SENSITIVITY] Computing E-values for Unmeasured Confounding...\n")

sum_cs <- summary(pooled_cs)
est_coef <- sum_cs$estimate[sum_cs$term == "Chemo_flagTreated"]
est_se   <- sum_cs$std.error[sum_cs$term == "Chemo_flagTreated"]

HR_chemo <- as.numeric(exp(est_coef))
HR_lower <- as.numeric(exp(est_coef - 1.96 * est_se))
HR_upper <- as.numeric(exp(est_coef + 1.96 * est_se))

cat("Hazard Ratio (Chemotherapy Main Effect):", round(HR_lower,3), round(HR_chemo, 3), round(HR_upper,3), "\n")

calculate_evalue <- function(estimate) {
  if(is.na(estimate)) return(NA)
  if(estimate < 1) estimate <- 1 / estimate 
  if(estimate == 1) return(1)
  return(estimate + sqrt(estimate * (estimate - 1)))
}

cat("E-Value (Point Estimate):", round(calculate_evalue(HR_chemo), 2), "\n")

# ==============================================================================
# MANUSCRIPT VISUALIZATIONS & EXPORT
# ==============================================================================
cat("\n[EXPORT] Generating Publication-Ready Figures & Tables...\n")

out_dir <- "Manuscript_Outputs"
dir.create(out_dir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# TABLES EXPORT (Baseline, Fine-Gray, CS-Cox)
# ------------------------------------------------------------------------------
# Table 1: Baseline Characteristics
df_table1 <- df_imputed[, .(Event_Status_Factor, Age_continuous, FollowUp_months, Chemo_flag, Radio_flag, Stage_advanced, Sex)]
table1 <- tbl_summary(
  df_table1, by = Event_Status_Factor, 
  statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)")
) %>% add_p() %>% add_overall() %>% modify_header(label ~ "**Characteristic**")
as_gt(table1) %>% gt::gtsave(filename = paste0(out_dir, "/Table1_Baseline_Characteristics.docx"))

# Table 2: Fine-Gray Pooled Results (Cumulative Incidence)
sum_fg <- summary(pooled_fg)
df_fg_results <- data.frame(
  Predictor = sum_fg$term,
  HR = exp(sum_fg$estimate),
  LCI = exp(sum_fg$estimate - 1.96 * sum_fg$std.error),
  UCI = exp(sum_fg$estimate + 1.96 * sum_fg$std.error),
  p_value = sum_fg$p.value
)
write.csv(df_fg_results, paste0(out_dir, "/Table2_FineGray_Pooled.csv"), row.names = FALSE)

# Table 3: Cause-Specific Cox Pooled Results (Etiology)
sum_cs <- summary(pooled_cs)
df_cs_results <- data.frame(
  Predictor = sum_cs$term,
  HR = exp(sum_cs$estimate),
  LCI = exp(sum_cs$estimate - 1.96 * sum_cs$std.error),
  UCI = exp(sum_cs$estimate + 1.96 * sum_cs$std.error),
  p_value = sum_cs$p.value
)
write.csv(df_cs_results, paste0(out_dir, "/Table3_Etiology_CS_Cox.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# FIGURE 1: Cumulative Incidence Function (CIF)
# ------------------------------------------------------------------------------
cif_fit <- tidycmprsk::cuminc(Surv(FollowUp_months, Event_Status_Factor) ~ Chemo_flag, data = df_imputed)

fig1 <- cif_fit %>%
  ggcuminc(outcome = "1", size = 1.2) + # <-- FIX: L'esito ora è codificato come "1"
  add_confidence_interval() +
  scale_color_manual(values = c("Untreated" = "#2c3e50", "Treated" = "#e74c3c")) +
  scale_fill_manual(values = c("Untreated" = "#2c3e50", "Treated" = "#e74c3c")) +
  theme_minimal(base_size = 14) +
  labs(title = "Cumulative Incidence of Secondary NHL",
       subtitle = "Fine-Gray Competing Risk Framework (Accounting for Death)",
       x = "Follow-up Time (Months)", y = "Cumulative Probability") +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(paste0(out_dir, "/Figure1_CIF_Trajectories.pdf"), plot = fig1, width = 9, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# FIGURE 2: Etiology Forest Plot
# ------------------------------------------------------------------------------
df_plot <- df_cs_results[!grepl("Intercept", df_cs_results$Predictor) & 
                           df_cs_results$HR > 0 & 
                           df_cs_results$LCI > 0 & 
                           !is.infinite(df_cs_results$HR), ]

fig2 <- ggplot(df_plot, aes(x = HR, y = reorder(Predictor, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_point(size = 3, color = "#2c3e50") +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2, color = "#2c3e50") +
  scale_x_log10() + theme_minimal() +
  labs(title = "Predictors of Mutagenic Clonal Evolution (Cause-Specific Cox)",
       x = "Hazard Ratio (Log Scale) -> Favors Late-Onset Evolution", y = "")
ggsave(paste0(out_dir, "/Figure2_Etiology_ForestPlot.pdf"), plot = fig2, width = 8, height = 6)

# ------------------------------------------------------------------------------
# FIGURE 3: Optimal K Silhouette & FAMD Clusters
# ------------------------------------------------------------------------------
fig3_A <- sil_plot + theme_minimal() + labs(title = "Optimal K Extraction (Silhouette Statistic)")
ggsave(paste0(out_dir, "/Figure3A_Silhouette_Optimal_K.pdf"), plot = fig3_A, width = 7, height = 5)

fig3_B <- fviz_cluster(km_res, data = res.famd$ind$coord, geom = "point", ellipse.type = "convex", 
                       ggtheme = theme_minimal(), main = paste("Latent Clinical Phenotypes (K =", optimal_k, ")"))
ggsave(paste0(out_dir, "/Figure3B_Latent_Clusters.pdf"), plot = fig3_B, width = 8, height = 6)

# ------------------------------------------------------------------------------
# FIGURE 4: Proportional Hazards Check (Schoenfeld Residuals for Stage)
# ------------------------------------------------------------------------------
pdf(paste0(out_dir, "/Figure4_Schoenfeld_Residuals.pdf"), width = 8, height = 6)

zph_plot_model <- coxph(Surv(FollowUp_months, Event_Status == 1) ~ 
                          Chemo_flag + ns(Age_continuous, df = 3) + 
                          Radio_flag + Sex + Race + Stage_advanced, 
                        data = df_imputed)

zph_plot_obj <- cox.zph(zph_plot_model)

stage_var_name <- rownames(zph_plot_obj$table)[grep("Stage", rownames(zph_plot_obj$table))]

if(length(stage_var_name) > 0) {
  plot(zph_plot_obj, var = stage_var_name[1], main = "Schoenfeld Residuals: Advanced Stage") 
  abline(h = 0, col = "red", lty = 2)
} else {
  plot(zph_plot_obj, var = 1, main = "Schoenfeld Residuals (Fallback variable)")
  abline(h = 0, col = "red", lty = 2)
}

dev.off()

# ------------------------------------------------------------------------------
# FIGURE 5: 10-Year Calibration Plot
# ------------------------------------------------------------------------------
ggsave(paste0(out_dir, "/Figure5_Calibration_Plot.pdf"), plot = p_calib, width = 7, height = 7, dpi = 300)

cat("\n[SUCCESS] All manuscript outputs saved to directory.\n")
# ==============================================================================