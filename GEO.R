# ==============================================================================
# TITLE: Decoding Secondary Lymphomas: Integrating Macro-Epidemiological Models 
#        and Micro-Genomic Signatures to Disentangle Mutagenesis from Misdiagnosis.
# AUTHORS: Pietro Grassi and Elena Bariletti
# METHODOLOGY: Nested Sure Independence Screening (SIS), Stability Selection, 
#              Elastic Net, Gene Ontology (GO) Enrichment, and PLS-DA.
# ==============================================================================

rm(list=ls())
gc()
setwd("/Users/pietro/Desktop/Sant'Anna/SLLD/Progetto")

# ==============================================================================
# PHASE 0: DEPENDENCIES & ENVIRONMENT SETUP
# Rationale: Loading specialized statistical learning architectures for high-
# dimensional (n << p) regimes and bioinformatics APIs for omics annotation.
# ==============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

library(GEOquery)        # API for NCBI Gene Expression Omnibus
library(data.table)      # High-performance data manipulation
library(glmnet)          # Penalized Generalized Linear Models (Elastic Net)
library(SIS)             # Sure Independence Screening for ultra-high dimensions
library(caret)           # Cross-validation and data partitioning algorithms
library(pls)             # Partial Least Squares Discriminant Analysis (PLS-DA)
library(pROC)            # Receiver Operating Characteristic (ROC) analytics
library(ggplot2)         # Vector graphics and visualizations
library(pheatmap)        # Unsupervised hierarchical clustering heatmaps
library(gridExtra)       # Multi-panel figure arrangements
library(clusterProfiler) # Over-Representation Analysis (ORA)
library(org.Hs.eg.db)    # Human genome annotation database
library(survival)        # Core survival analysis architecture
library(ggsurvfit)       # High-fidelity Kaplan-Meier visualizations
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("sva")
if (!requireNamespace("impute", quietly = TRUE)) BiocManager::install("impute")
library(impute)

set.seed(2026)       

# ==============================================================================
# PHASE 1: OMICS DATA INGESTION & PROBE ANNOTATION
# Rationale: Retrieving microarray matrices and strictly mapping technical probes 
# to universal Gene Symbols to enable downstream biological pathway enrichment.
# ==============================================================================
cat("\n[PHASE 1] Fetching GSE31312 matrix and executing Gene Symbol mapping...\n")

options(timeout = 600, download.file.method = "libcurl")
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)

gse <- suppressMessages(getGEO("GSE31312", GSEMatrix = TRUE, AnnotGPL = TRUE))
eset <- gse[[1]]
clin_data <- pData(eset)

# Molecular Subtype Identification (Target Extraction)
patient_metadata <- apply(clin_data, 1, paste, collapse = " ")
is_abc <- grepl("ABC", patient_metadata, ignore.case = FALSE)
is_gcb <- grepl("GCB", patient_metadata, ignore.case = FALSE)
valid_patients <- is_abc | is_gcb

# Probe Annotation to Official Gene Symbols
feat_data <- fData(eset)
sym_col <- grep("(?i)^gene symbol$|^symbol$", colnames(feat_data), value = TRUE)

if(length(sym_col) > 0) {
  gene_symbols <- sapply(strsplit(as.character(feat_data[[sym_col[1]]]), "///"), `[`, 1)
} else {
  gene_symbols <- rownames(exprs(eset)) # Fallback to raw IDs if annotation fails
}

Y_Etiology <- ifelse(is_abc[valid_patients], 1, 0) 
Y_Factor <- as.factor(ifelse(Y_Etiology == 1, "ABC_Aggressive", "GCB_Indolent"))

set.seed(2026)
train_idx_global <- createDataPartition(Y_Etiology, p = 0.8, list = FALSE)

exprs_raw <- exprs(eset)[, valid_patients]
exprs_train <- exprs_raw[, train_idx_global]
exprs_test  <- exprs_raw[, -train_idx_global]
imp_train <- impute::impute.knn(exprs_train, k = 10)$data
imp_test  <- impute::impute.knn(exprs_test, k = 10)$data
X_train_global <- t(imp_train)
X_test_internal <- t(imp_test)
valid_genes <- which(!is.na(gene_symbols) & gene_symbols != "")
X_train_global <- X_train_global[, valid_genes]
X_test_internal <- X_test_internal[, valid_genes]

clean_symbols <- make.unique(gene_symbols[valid_genes])
colnames(X_train_global) <- clean_symbols
colnames(X_test_internal) <- clean_symbols
Y_train_global <- Y_Etiology[train_idx_global]
Y_test_internal <- Y_Etiology[-train_idx_global]
Y_Factor_test <- Y_Factor[-train_idx_global]
X_genes <- rbind(X_train_global, X_test_internal)
X_raw <- X_genes

Y_Etiology_reordered <- c(Y_train_global, Y_test_internal)
Y_Factor_reordered <- as.factor(ifelse(Y_Etiology_reordered == 1, "ABC_Aggressive", "GCB_Indolent"))

cat("Data Ingestion Complete | Samples (n):", nrow(X_genes), "| Mapped Features (p):", ncol(X_genes), "\n")

# ==============================================================================
# PHASE 2: NESTED SIS & STABILITY SELECTION
# Rationale: To rigorously prevent high-dimensional Data Leakage, marginal screening 
# is nested within 100 bootstrap iterations. Only highly stable genes (selection 
# frequency > 50%) are retained for the final penalized regression framework.
# ==============================================================================
cat("\n[PHASE 2] Executing Nested SIS & Stability Selection (B = 100)...\n")

B <- 100
gene_selection_tally <- table(character())

pb <- txtProgressBar(min = 0, max = B, style = 3)
for (i in 1:B) {
  set.seed(2026 + i) 
  
  # Stratified out-of-bag sampling
  sub_idx <- createDataPartition(Y_train_global, p = 0.8, list = FALSE)
  x_train <- X_train_global[sub_idx, ]
  y_train <- Y_train_global[sub_idx]
  
  n_theoretical <- floor(nrow(x_train) * 2)
  
  # Nested Marginal Screening
  marginal_cor <- suppressWarnings(abs(cor(x_train, y_train)))
  marginal_cor[is.na(marginal_cor)] <- 0
  
  screened_idx <- order(marginal_cor, decreasing = TRUE)[1:n_theoretical]
  x_screened_train <- x_train[, screened_idx]
  
  # Nested Penalized Regression (Elastic Net)
  cv_fit <- suppressWarnings(cv.glmnet(x_screened_train, y_train, family = "binomial", 
                                       alpha = 0.5, type.measure = "auc", nfolds = 3))
  coefs <- coef(cv_fit, s = "lambda.min")
  selected <- rownames(coefs)[which(coefs != 0)]
  selected <- selected[selected != "(Intercept)"]
  
  if(length(selected) > 0) {
    gene_selection_tally <- c(gene_selection_tally, selected)
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# Core Signature Extraction
gene_freq <- table(gene_selection_tally) / B
core_genes <- names(gene_freq)[gene_freq > 0.5]

# Final Model Construction strictly on the Global Training Set (Zero Data Leakage)
X_train_global_sig <- X_train_global[, core_genes, drop = FALSE]
final_elastic_net <- cv.glmnet(X_train_global_sig, Y_train_global, family = "binomial", 
                               alpha = 0.5, nfolds = 5, keep = TRUE)
optimal_coefs <- coef(final_elastic_net, s = "lambda.min")

cat("\nStability Selection isolated a robust Core Signature of", length(core_genes), "genes.\n")

# ==============================================================================
# PHASE 3: GENE SET ENRICHMENT ANALYSIS (GSEA)
# Rationale: Moving beyond weak ORA by ranking the entire transcriptome based 
# on marginal correlation with etiology to identify systemic biological shifts.
# This prevents extreme FDR penalties on truncated signatures.
# ==============================================================================
cat("\n[PHASE 3] Executing Gene Set Enrichment Analysis (GSEA)...\n")

global_cor <- apply(X_train_global, 2, function(x) cor(x, Y_train_global, use="pairwise.complete.obs"))
global_cor[is.na(global_cor)] <- 0

gene_list <- sort(global_cor, decreasing = TRUE)

names(gene_list) <- gsub("\\.[0-9]+$", "", names(gene_list))

df_genes <- data.frame(gene = names(gene_list), cor = gene_list)
df_genes <- df_genes[!duplicated(df_genes$gene), ]

gene_list_clean <- df_genes$cor
names(gene_list_clean) <- df_genes$gene

set.seed(2026)
gsea_res <- tryCatch({
  gseGO(geneList      = gene_list_clean,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        minGSSize     = 15,
        maxGSSize     = 500,
        pvalueCutoff  = 0.05,
        pAdjustMethod = "BH",
        eps           = 0,
        nPermSimple   = 10000)
}, error = function(e) { NULL })

if(!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
  cat("Successfully executed GSEA: Identified active systemic pathways.\n")
  p_pathway <- dotplot(gsea_res, showCategory = 8, split = ".sign") + 
    facet_grid(.~.sign) +
    theme_minimal(base_size = 12) + 
    theme(plot.title = element_text(face = "bold"))
  
} else {
  cat("[LOG] No significant GSEA pathways detected.\n")
  p_pathway <- ggplot() + annotate("text", x=0, y=0, label="No Significant Pathways Found") + theme_void()
}

# ==============================================================================
# PHASE 4: UNSUPERVISED DIMENSION REDUCTION (PCA)
# Rationale: Visually establishing the natural biological variance captured 
# by the core signature without forcing supervised separation.
# ==============================================================================
cat("\n[PHASE 4] Generating Dimension Reduction Projections...\n")

# Projecting the entire cohort to visualize global phenotypic separation
X_signature_global <- X_genes[, core_genes, drop = FALSE]
pca_res <- prcomp(X_signature_global, scale. = TRUE)

pca_scores <- data.frame(Etiology = Y_Factor_reordered, PC1 = pca_res$x[,1], PC2 = pca_res$x[,2])

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Etiology)) +
  geom_point(size = 3, alpha = 0.8) + stat_ellipse(level = 0.95) +
  scale_color_manual(
    name = "Transcriptomic Phenotype", 
    values = c("GCB_Indolent" = "#0073C2FF", "ABC_Aggressive" = "#EFC000FF"),
    labels = c("Aggressive (ABC-like)", "Indolent (GCB-like)")
  ) +
  theme_minimal(base_size = 14) + 
  labs(title = NULL, subtitle = NULL)

# ==============================================================================
# PHASE 5: STRICT INTERNAL VALIDATION (ON HOLD-OUT SET) & YOUDEN THRESHOLD
# ==============================================================================
cat("\n[PHASE 5] Executing Strict Internal Validation on Hold-Out Set...\n")

# Estraiamo la firma solo sui pazienti vergini
X_test_internal_sig <- X_test_internal[, core_genes, drop = FALSE]

# Previsioni e calcolo AUC puro
cv_preds_internal <- predict(final_elastic_net, newx = X_test_internal_sig, s = "lambda.min", type = "response")
df_results <- data.frame(True_Clinical_Label = Y_Factor_test, Probability_ABC = as.numeric(cv_preds_internal))

roc_obj <- roc(Y_test_internal, df_results$Probability_ABC, quiet = TRUE)
best_thresh <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")$threshold[1]

cat("Strict Internal Validation AUC:", round(auc(roc_obj), 3), "\n")

# ==============================================================================
# PHASE 6: EXTERNAL VALIDATION WITH EMPIRICAL BAYES HARMONIZATION (ComBat)
# Rationale: Testing on GSE10846. To prevent massive batch effects from masking 
# the biological signal, we apply ComBat on the FULL transcriptome before signature extraction.
# ==============================================================================
cat("\n[PHASE 6] External Validation & ComBat Harmonization (GSE10846)...\n")

gse_ext <- suppressMessages(getGEO("GSE10846", GSEMatrix = TRUE, AnnotGPL = TRUE))
eset_ext <- gse_ext[[1]]
clin_data_ext <- pData(eset_ext)

meta_ext <- apply(clin_data_ext, 1, paste, collapse = " ")
valid_ext <- grepl("ABC|GCB", meta_ext, ignore.case = FALSE)
Y_Etiology_ext <- ifelse(grepl("ABC", meta_ext[valid_ext], ignore.case = FALSE), 1, 0)

feat_data_ext <- fData(eset_ext)

sym_col_ext <- grep("(?i)^gene symbol$|^symbol$", colnames(feat_data_ext), value = TRUE)
if(length(sym_col_ext) > 0) {
  gene_symbols_ext <- sapply(strsplit(as.character(feat_data_ext[[sym_col_ext[1]]]), "///"), `[`, 1)
} else {
  gene_symbols_ext <- rownames(exprs(eset_ext))
}

exprs_ext <- exprs(eset_ext)[, valid_ext]
imputed_exprs_ext <- impute::impute.knn(exprs_ext, k = 10)$data
X_raw_ext <- t(imputed_exprs_ext)

colnames(X_raw_ext) <- make.unique(gene_symbols_ext)

common_genes_all <- intersect(colnames(X_train_global), colnames(X_raw_ext))
cat("Found", length(common_genes_all), "total overlapping transcripts for ComBat modeling.\n")

X_train_full <- X_train_global[, common_genes_all]
X_test_full <- X_raw_ext[, common_genes_all]

combined_X_full <- t(rbind(X_train_full, X_test_full))
batch_labels <- c(rep(1, nrow(X_train_full)), rep(2, nrow(X_test_full)))

combat_full <- suppressMessages(sva::ComBat(dat = combined_X_full, batch = batch_labels, 
                                            mod = NULL, par.prior = TRUE, 
                                            prior.plots = FALSE, ref.batch = 1))
combat_full <- t(combat_full)

X_train_harmonized <- combat_full[1:nrow(X_train_full), ]
X_test_harmonized <- combat_full[(nrow(X_train_full)+1):nrow(combat_full), ]

shared_core <- intersect(core_genes, common_genes_all)
X_train_sig <- X_train_harmonized[, shared_core]
X_test_sig <- X_test_harmonized[, shared_core]

preds_ext <- predict(final_elastic_net, newx = X_test_sig, s = "lambda.min", type = "response")
roc_ext <- pROC::roc(Y_Etiology_ext, as.numeric(preds_ext), quiet = TRUE)
auc_ext <- as.numeric(pROC::auc(roc_ext))

cat("\nUniversal Biomarker Validation AUC (Post-ComBat Harmonization):", round(auc_ext, 3), "\n")
if(auc_ext > 0.80) cat("-> SUCCESS: True cross-cohort translational robustness achieved!\n")

# ==============================================================================
# PUBLICATION-READY EXPORT
# ==============================================================================
cat("\n[EXPORT] Generating High-Fidelity Manuscript Visualizations...\n")

out_dir <- "Genomic_Outputs"
dir.create(out_dir, showWarnings = FALSE)

# Figure 1: Unsupervised PCA
ggsave(paste0(out_dir, "/Figure1_Unsupervised_PCA.pdf"), plot = p_pca, width = 8, height = 6, dpi = 300)

# Figure 2: Hierarchical Clustering Heatmap

Clean_Phenotype <- as.character(Y_Factor_reordered)
Clean_Phenotype[Clean_Phenotype == "ABC_Aggressive"] <- "ABC (Aggressive)"
Clean_Phenotype[Clean_Phenotype == "GCB_Indolent"] <- "GCB (Indolent)"

annotation_col <- data.frame(Phenotype = factor(Clean_Phenotype, levels = c("ABC (Aggressive)", "GCB (Indolent)")))
rownames(annotation_col) <- rownames(X_signature_global)

ann_colors <- list(
  Phenotype = c("ABC (Aggressive)" = "#EFC000FF", "GCB (Indolent)" = "#0073C2FF")
)

pheatmap(t(X_signature_global), 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         scale = "row", 
         show_colnames = FALSE, 
         clustering_method = "ward.D2",
         main = "",                 
         fontsize_row = 6,          
         filename = paste0(out_dir, "/Figure2_Genomic_Heatmap.pdf"), 
         width = 8, 
         height = 12)

# Figure 3: ROC Analytics
true_auc <- as.numeric(auc(roc_obj))
fig3 <- ggroc(roc_obj, colour = "#EFC000FF", size = 1.2) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "darkgray") +
  theme_minimal(base_size = 14) +
  labs(x = "Specificity", y = "Sensitivity")
ggsave(paste0(out_dir, "/Figure3_ROC_Curve.pdf"), plot = fig3, width = 7, height = 7, dpi = 300)

# Figure 4: Pathway Analysis
ggsave(paste0(out_dir, "/Figure4_Pathway_Enrichment.pdf"), plot = p_pathway, width = 9, height = 7, dpi = 300)

# Figure 5: External Validation ROC Curve
fig5 <- ggroc(roc_ext, colour = "#27ae60", size = 1.2) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "darkgray") +
  theme_minimal(base_size = 14) +
  labs(x = "Specificity", y = "Sensitivity")
ggsave(paste0(out_dir, "/Figure5_External_Validation_ROC.pdf"), plot = fig5, width = 7, height = 7, dpi = 300)

# Table 1: Stable Signature Extraction
df_coefs <- data.frame(Gene_Symbol = core_genes, Selection_Frequency_Pct = as.numeric(gene_freq[core_genes]) * 100)
df_coefs <- df_coefs[order(df_coefs$Selection_Frequency_Pct, decreasing = TRUE), ]
write.csv(df_coefs, paste0(out_dir, "/Table1_Stable_Genomic_Signature.csv"), row.names = FALSE)

# Table 2: Clinical Reclassification (Occult Aggressiveness)
df_reclass <- table(True_Clinical_Label = df_results$True_Clinical_Label, 
                    Predicted_Genomic_Label = ifelse(df_results$Probability_ABC > best_thresh, "Aggressive", "Indolent"))
write.csv(as.data.frame(df_reclass), paste0(out_dir, "/Table2_Clinical_Reclassification.csv"), row.names = FALSE)

cat("\n[SUCCESS] All analytical outputs rigorously exported.\n")
# ==============================================================================

cat("GSE31312 Raw Samples:", length(patient_metadata), "\n")
cat("GSE31312 Valid ABC/GCB:", sum(valid_patients), "\n")
cat("GSE31312 Train (80%):", length(train_idx_global), "\n")
cat("GSE31312 Hold-out (20%):", sum(valid_patients) - length(train_idx_global), "\n")

cat("GSE10846 Raw Samples:", length(meta_ext), "\n")
cat("GSE10846 Valid ABC/GCB:", sum(valid_ext), "\n")
