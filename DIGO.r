#-----------------------------------------------------------------------------#
#            R Script for Manuscript: GWAS, Fine-mapping, and Haplotype Analysis
#-----------------------------------------------------------------------------#

# This script performs a comprehensive analysis of genetic data, including
# a Genome-Wide Association Study (GWAS) visualization, gene-based analysis,
# fine-mapping of candidate genes, and haplotype-based predictive modeling.

# All required data objects are assumed to be pre-loaded from "MS.Rdata".
# These objects include:
# - combined_gwas: Formatted GWAS results
# - outMM: MAGMA gene-based analysis results
# - top_hits: Top variants from the GWAS
# - genebed: Gene annotations
# - sumstats: Summary statistics for fine-mapping
# - idlist: SNP identifiers
# - ld_matriX: Linkage disequilibrium matrix
# - outHap: Haplotype data frame
# - C1QL2, RSPH4A, RWDD1, BBS1, TTC7B, TOM1L1: SuSiE fine-mapping results for specific genes
# - outHapHap: Haplotype association results
# - model1, model2, model3: Logistic regression models

# Load necessary R packages
# Ensure these packages are installed in your R environment:
# install.packages(c("data.table", "ggplot2", "dplyr", "ggrepel", "susieR",
#                    "Matrix", "cowplot", "haplo.stats", "pROC", "pheatmap", "stringr"))

library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(susieR)
library(Matrix)
library(cowplot)
library(haplo.stats)
library(pROC)
library(pheatmap)
library(stringr) # For str_count in cntHap function

# Load the pre-processed data from MS.Rdata
# This file should contain all data objects required for the script to run.
load("MS.Rdata")

#-----------------------------------------------------------------------------#
# SECTION 1: Genome-Wide Association Study (GWAS) and Manhattan Plot
#-----------------------------------------------------------------------------#

# Prepare data for the Manhattan plot axis and colors
axis_df <- combined_gwas %>%
  group_by(CHROM) %>%
  summarize(center = mean(bp_cum))

# Assign colors to chromosomes for alternating appearance in the plot
combined_gwas$Dataset <- c("steelblue", "skyblue")[combined_gwas$X.CHROM %% 2 + 1]

# Set plot dimensions for the Manhattan plot
options(repr.plot.width = 16, repr.plot.height = 6)

# Generate the Manhattan plot
gwas_manhattan_plot <- ggplot(combined_gwas, aes(x = bp_cum, y = -log10(P))) +
  geom_point(aes(color = Dataset), alpha = 0.6, size = 1.2) + # Plot points, colored by chromosome
  scale_x_continuous(labels = axis_df$CHROM, breaks = axis_df$center) + # Set x-axis labels to chromosome numbers
  geom_text_repel(data = top_hits, aes(label = SYMBOL), # Add gene symbols for top hits
                  size = 3, max.overlaps = 80, nudge_y = 0.3, color = "black") +
  theme_bw() + # Use a black and white theme
  labs(title = "Manhattan Plot of GWAS Results", # Plot title
       x = "Chromosome", y = "-log10(P-value)") + # Axis labels
  theme(legend.position = "none") # Hide the legend for chromosome colors

print(gwas_manhattan_plot)

#-----------------------------------------------------------------------------#
# SECTION 2: Gene-Based Analysis (MAGMA) Plot
#-----------------------------------------------------------------------------#

# Calculate -log10(Adjusted P-value) if not already present in outMM
if (!"neg_log10_AdjP" %in% names(outMM)) {
  outMM$neg_log10_AdjP <- -log10(outMM$AdjP)
}

# Order genes by chromosomal position and convert symbol to factor for plotting order
outMM <- outMM[order(outMM$CHR, outMM$START), ]
outMM$symbol <- factor(outMM$symbol, levels = unique(outMM$symbol))

# Generate the MAGMA gene-based plot
magma_plot <- ggplot(outMM, aes(x = symbol, y = neg_log10_AdjP, size = abs(ZSTAT))) +
  geom_point(aes(color = factor(CHR))) + # Plot points, colored by chromosome
  scale_size_continuous(range = c(2, 8), name = "|ZSTAT|") + # Scale point size by Z-statistic
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) + # Significance threshold
  geom_text_repel(data = outMM[outMM$AdjP < 0.05, ], # Label significant genes
                  aes(label = symbol),
                  size = 3, max.overlaps = 20, nudge_y = 0.2, color = "black") +
  labs(title = "MAGMA Gene-Based Analysis",
       x = "Gene (Ordered by Chromosomal Position)",
       y = "-log10(Adjusted P-value)",
       color = "Chromosome") +
  guides(color = guide_legend(ncol = 2)) + # Legend for chromosome colors
  theme_minimal() + # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Rotate x-axis labels
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "right",
        panel.grid.major.x = element_blank()) # Remove major x-grid lines

print(magma_plot)

#-----------------------------------------------------------------------------#
# SECTION 3: SuSiE Fine-Mapping and Plotting
#-----------------------------------------------------------------------------#

# Function to perform SuSiE fine-mapping for a given gene
getSusie <- function(cfPIP, gene, genebed, sumstats, idlist, ld_matriX) {
  # Subset genebed for the target gene
  tmpgene <- genebed[grep(gene, genebed$V11), ]
  # Subset sumstats for variants within the gene region
  tmpss <- sumstats[sumstats$V1 == gsub("chr", "", tmpgene[, 1]) &
                      sumstats$V2 >= tmpgene[, 2] &
                      sumstats$V2 <= tmpgene[, 3] & sumstats$V8 == "ADD", ]
  # Subset idlist for SNPs matching the sumstats
  tmpID <- idlist[idlist$V2 %in% gsub(":[A-Z].*.", "", tmpss$V3), ]
  # Subset LD matrix for relevant SNPs
  tmpLD <- ld_matriX[idlist$V2 %in% gsub(":[A-Z].*.", "", tmpss$V3),
                     idlist$V2 %in% gsub(":[A-Z].*.", "", tmpss$V3)]
  # Filter sumstats again to ensure consistency with tmpID
  tmpss <- tmpss[grep(paste(tmpID$V2, collapse = "|"), tmpss$V3), ]
  
  # Ensure LD matrix is symmetric
  R_sym <- as.matrix(forceSymmetric((tmpLD + t(tmpLD)) / 2, uplo = "U"))
  
  # Extract Z-scores
  z_scores <- tmpss$V12
  
  # Perform SuSiE fine-mapping
  susie_res <- susie_rss(z = z_scores, R = R_sym, L = 10, n = 74, coverage = 0.95)
  
  # Prepare data for plotting
  plot_data <- data.table(ID = tmpss$V3, POS = tmpss$V2,
                          P = tmpss$V13, PIP = susie_res$pip,
                          loca = paste0(tmpss$V3, "\n     ", tmpss$V2),
                          gene = gene, CF = susie_res$pip > cfPIP)
  return(plot_data)
}

# Function to plot SuSiE fine-mapping results
plotSusie <- function(plot_data, cfPIP) {
  # Top panel: -log10(P) vs. position
  p_pvalue <- ggplot(plot_data, aes(x = POS, y = -log10(P), size = PIP, color = PIP > cfPIP)) +
    geom_point(alpha = 0.8) +
    scale_size_continuous(range = c(2, 8), name = "PIP") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = paste0("PIP > ", cfPIP)) +
    geom_text_repel(data = plot_data[PIP > cfPIP], aes(label = loca),
                    size = 3, max.overlaps = 5, nudge_y = 0, color = "black") +
    labs(title = paste0("SuSiE Fine-Mapping of ", plot_data$gene[1], " (PIP cutoff=", cfPIP, ")"),
         x = NULL, y = "-log10(P-value)") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), # Hide x-axis for alignment
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), # Add rectangle border
          legend.position = "right")

  # Bottom panel: PIP vs. position
  p_pip <- ggplot(plot_data, aes(x = POS, y = PIP, size = 2, color = PIP > cfPIP)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "black"), name = paste0("PIP > ", cfPIP)) +
    labs(x = paste0("Position (CHR ", gsub(":.*.", "", plot_data$ID[1]), ")"), y = "PIP") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), # Add rectangle border
          legend.position = "none") # Hide legend to avoid duplication

  # Combine plots
  return(plot_grid(p_pvalue, p_pip, ncol = 1, align = "v", axis = "lr",
                   rel_heights = c(2, 1)))
}

# Plot fine-mapping results for selected genes
options(repr.plot.width = 14, repr.plot.height = 16)

# Assuming C1QL2, RSPH4A, RWDD1, BBS1, TTC7B, TOM1L1 are generated and loaded from MS.Rdata
# If not, you would need to run the getSusie function for each here:
# C1QL2 <- getSusie(0.2, "C1QL2", genebed, sumstats, idlist, ld_matriX)
# RSPH4A <- getSusie(0.08, "RSPH4A", genebed, sumstats, idlist, ld_matriX)
# RWDD1 <- getSusie(0.15, "RWDD1", genebed, sumstats, idlist, ld_matriX)
# BBS1 <- getSusie(0.15, "BBS1", genebed, sumstats, idlist, ld_matriX)
# TTC7B <- getSusie(0.045, "TTC7B", genebed, sumstats, idlist, ld_matriX)
# TOM1L1 <- getSusie(0.045, "TOM1L1", genebed, sumstats, idlist, ld_matriX)

plot_grid(plotSusie(C1QL2, 0.2), plotSusie(RSPH4A, 0.08), plotSusie(RWDD1, 0.15),
          ncol = 1, align = "v", axis = "lr")
plot_grid(plotSusie(BBS1, 0.15), plotSusie(TTC7B, 0.045), plotSusie(TOM1L1, 0.045),
          ncol = 1, align = "v", axis = "lr")

#-----------------------------------------------------------------------------#
# SECTION 4: Haplotype-Based Predictive Modeling and Forest Plot
#-----------------------------------------------------------------------------#

# Function to count haplotype occurrences
cntHap <- function(x, ii) {
  x[, 4] <- gsub("ref\\.", "", x[, 4]) # Remove "ref." prefix
  x[, 4] <- str_count(x[, 4], pattern = ii) # Count occurrences of haplotype 'ii'
  return(x)
}

# Prepare data for logistic regression models based on significant haplotypes
# Assuming outHapHap is loaded from MS.Rdata and contains p-values for haplotypes
df_p_values <- data.frame(terms = paste0("chr", outHapHap$chr, ":", gsub("hap", "", outHapHap$haplotype)),
                          pvalue = outHapHap$Pr...z..)
significant_terms <- as.character(df_p_values[df_p_values$pvalue < 0.1, ]$terms)

inpHap <- NULL
x <- 0
for (i in sort(significant_terms)) {
  ia <- unlist(strsplit(i, ":"))
  tmphap <- outHap[, c("id", "y", "covar", gsub("chr", "", ia[1]))]
  tmphap <- cntHap(tmphap, ia[2])
  colnames(tmphap)[4] <- i
  if (x < 1) {
    inpHap <- tmphap
  } else {
    inpHap <- merge(inpHap, tmphap[, c(1, 4)], by.x = "id", by.y = "id")
  }
  x <- x + 1
}

colnames(inpHap)[3] <- "AML"

# Build the logistic regression models
# Model 1: Predicts phenotype using AML covariate and significant haplotypes
model1 <- glm(y ~ ., data = inpHap[, -1], family = "binomial")
# Model 2: Predicts phenotype using only the AML covariate
model2 <- glm(y ~ AML, data = inpHap[, -1], family = "binomial")
# Model 3: Predicts phenotype using only the significant haplotypes
model3 <- glm(y ~ ., data = inpHap[, -c(1, 3)], family = "binomial")

# Extract and format coefficients for Model 1 for the forest plot
finalDF <- summary(model1)$coefficients
finalDF <- data.frame(Var = gsub("`|2$", "", row.names(finalDF)), finalDF)
colnames(finalDF) <- c("term", "estimate", "se", "zscore", "pvalue")

# Prepare data for the forest plot
df_plot <- finalDF[!is.na(finalDF$pvalue) & finalDF$term != "(Intercept)", ]
df_plot[, 2:5] <- apply(df_plot[, 2:5], 2, as.numeric)

# Calculate confidence intervals and format labels for the forest plot
df_plot <- df_plot %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    label = sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper)
  ) %>%
  arrange(estimate) %>% # Sort by estimate for better visualization
  mutate(term = factor(term, levels = term)) # Convert term to factor with ordered levels

# Set plot dimensions for the forest plot
options(repr.plot.width = 8, repr.plot.height = 6)

# Generate the Forest Plot
p1 <- ggplot(df_plot, aes(x = estimate, y = term)) +
  geom_point(size = 3) + # Plot point estimates
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) + # Add horizontal error bars for CI
  geom_vline(xintercept = 0, linetype = "dashed") + # Add a vertical dashed line at 0 (no effect)
  geom_text(aes(label = label), hjust = ifelse(df_plot$estimate >= 0, -0.1, 1.1),
            vjust = -0.4, size = 3) + # Add labels for estimates and CIs
  labs(x = "Log Odds Ratio (95% CI)", y = NULL, title = "") + # Axis labels and title
  theme_minimal() + # Minimal theme
  theme(axis.text.y = element_text(size = 10), plot.title = element_text(face = "bold")) +
  xlim(min(df_plot$lower) - 0.5, max(df_plot$upper) + 0.5) # Set x-axis limits

#-----------------------------------------------------------------------------#
# SECTION 5: ROC Curve Comparison
#-----------------------------------------------------------------------------#

# Generate ROC curves and AUC for the three models
roc1 <- roc(inpHap$y, predict(model1, inpHap[, -c(1:2)], type = "response"))
roc2 <- roc(inpHap$y, predict(model2, inpHap[, -1], type = "response"))
roc3 <- roc(inpHap$y, predict(model3, inpHap[, -c(1:3)], type = "response"))

# Define a common FPR grid for smooth plotting of ROC curves
fpr_grid <- seq(0, 1, length.out = 200)

# Function to interpolate True Positive Rate (TPR) at common False Positive Rate (FPR)
get_interp <- function(roc_obj, label, color) {
  fpr <- 1 - roc_obj$specificities
  tpr <- roc_obj$sensitivities
  interp_tpr <- approx(x = fpr, y = tpr, xout = fpr_grid, ties = mean)$y
  data.frame(FPR = fpr_grid, TPR = interp_tpr, Model = label, Color = color)
}

# Combine interpolated ROC data for plotting
df_all <- bind_rows(
  get_interp(roc1, paste0("AML + Haplotypes; AUC = ", round(auc(roc1), 3)), "brown"),
  get_interp(roc2, paste0("AML only; AUC = ", round(auc(roc2), 3)), "darkgreen"),
  get_interp(roc3, paste0("Haplotypes only; AUC = ", round(auc(roc3), 3)), "steelblue")
)

# Generate the ROC plot
p2 <- ggplot(df_all, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1.4) + # Plot ROC curves
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + # Add diagonal reference line
  scale_color_manual(values = c("brown", "darkgreen", "steelblue")) + # Manually set colors
  labs(title = "", x = "False Positive Rate", y = "True Positive Rate", color = "Model") +
  theme_minimal(base_size = 14) + # Minimal theme with larger base font size
  theme(legend.position = c(0.55, 0.2), # Position legend inside the plot
        legend.background = element_rect(fill = "white", color = "black"), # Legend background
        legend.title = element_blank(), # Hide legend title
        legend.key.size = unit(0.6, "lines")) + # Adjust legend key size
  guides(color = guide_legend(ncol = 1)) # Force vertical legend

# Combine the forest plot and ROC plot side-by-side
options(repr.plot.width = 10, repr.plot.height = 4)
plot_grid(p1, p2, ncol = 2, align = "h", axis = "lr", rel_widths = c(3, 2))

# End of script
