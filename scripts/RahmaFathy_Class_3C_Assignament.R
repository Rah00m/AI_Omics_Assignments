#### Load  Data ####
load("Results/GSE183795_processed.RData")

cat("Data loaded - Dimensions:", dim(processed_data), "\n")
cat("Groups:", table(groups), "\n")

#### 1. Probe-to-Gene Mapping ####
cat("Starting probe-to-gene mapping...\n")

# Get probe IDs
probe_ids <- rownames(processed_data)

# Map probe IDs to gene symbols
gene_symbols <- mapIds(
  hugene10sttranscriptcluster.db,
  keys = probe_ids,
  keytype = "PROBEID", 
  column = "SYMBOL",
  multiVals = "first"
)

# Create mapping dataframe
gene_map_df <- data.frame(
  PROBEID = names(gene_symbols),
  SYMBOL = as.character(gene_symbols),
  stringsAsFactors = FALSE
)

# Remove probes without gene symbols
gene_map_df <- gene_map_df[!is.na(gene_map_df$SYMBOL), ]
cat("Probes with gene symbols:", nrow(gene_map_df), "\n")

# Check duplicates
duplicate_summary <- gene_map_df %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(probes_per_gene = n()) %>%
  dplyr::arrange(desc(probes_per_gene))

cat("Probes per gene distribution:\n")
print(table(duplicate_summary$probes_per_gene))

# Merge with expression data - FIXED VERSION
processed_data_df <- processed_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::left_join(gene_map_df, by = "PROBEID") %>%
  dplyr::filter(!is.na(SYMBOL))

# Remove PROBEID column for averaging
expr_only <- processed_data_df %>% 
  dplyr::select(-PROBEID, -SYMBOL)

# Average multiple probes per gene
averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)

cat("Final genes after averaging:", nrow(averaged_data), "\n")

# Convert to matrix
expression_matrix <- as.matrix(averaged_data)

#### 2. Differential Expression Analysis ####
cat("Starting differential expression analysis...\n")

# Create design matrix
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# Fit linear model
fit <- lmFit(expression_matrix, design)

# Define contrast
contrast_matrix <- makeContrasts(
  Tumor_vs_Normal = Tumor - Normal,
  levels = design
)

# Apply contrasts
fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrast)

# Extract results
deg_results <- topTable(fit_ebayes,
                        coef = "Tumor_vs_Normal",
                        number = Inf,
                        adjust.method = "BH")

# Classify genes
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated",
         "Not Significant")
))

# Get significant genes
upregulated <- deg_results %>% dplyr::filter(threshold == "Upregulated")
downregulated <- deg_results %>% dplyr::filter(threshold == "Downregulated")
significant_genes <- rbind(upregulated, downregulated)

cat("DEG Results:\n")
cat("Upregulated genes:", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n")
cat("Total significant DEGs:", nrow(significant_genes), "\n")

# Save results
write.csv(deg_results, "Results/DEGs_Complete_Results.csv")
write.csv(upregulated, "Results/Upregulated_Genes.csv")
write.csv(downregulated, "Results/Downregulated_Genes.csv")
write.csv(significant_genes, "Results/Significant_DEGs.csv")

#### 3. Data Visualization ####
cat("Creating visualizations...\n")

dir.create("Results_Plots", showWarnings = FALSE)

# Volcano Plot
png("Results_Plots/Volcano_Plot_GSE183795.png", width = 2000, height = 1500, res = 300)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue",
                                "Not Significant" = "grey80")) +
  theme_minimal() +
  labs(title = "Volcano Plot - Pancreatic Cancer (Tumor vs Normal)",
       subtitle = paste("Up:", nrow(upregulated), "Down:", nrow(downregulated)),
       x = "log2 Fold Change", 
       y = "-log10(Adjusted P-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue")

dev.off()

# Heatmap of Top DEGs
top_genes <- rownames(significant_genes)[1:min(25, nrow(significant_genes))]
heatmap_data <- expression_matrix[top_genes, ]

png("Results_Plots/Heatmap_Top_DEGs_GSE183795.png", width = 2000, height = 1500, res = 300)

pheatmap(heatmap_data,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = paste("Top", length(top_genes), "Differentially Expressed Genes"),
         annotation_col = data.frame(Group = groups, row.names = colnames(heatmap_data)))

dev.off()

#### 4. Results Summary ####
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS SUMMARY - GSE183795\n")
cat(rep("=", 60), "\n", sep = "")

cat("1. PROBE-TO-GENE MAPPING:\n")
cat("   - Initial probes:", nrow(processed_data), "\n")
cat("   - Probes with gene symbols:", nrow(gene_map_df), "\n")
cat("   - Final unique genes after averaging:", nrow(averaged_data), "\n")
cat("   - Multiple probes mapping to same gene were averaged\n\n")

cat("2. DIFFERENTIAL EXPRESSION ANALYSIS:\n")
cat("   - Comparison: Tumor vs Normal (Pancreatic Cancer)\n")
cat("   - Upregulated genes:", nrow(upregulated), "\n")
cat("   - Downregulated genes:", nrow(downregulated), "\n")
cat("   - Total significant DEGs:", nrow(significant_genes), "\n")
cat("   - Criteria: |logFC| > 1 and adj.P.Val < 0.05\n\n")

cat("3. OUTPUT FILES:\n")
cat("   - DEGs_Complete_Results.csv: All genes with statistics\n")
cat("   - Upregulated_Genes.csv: Significantly upregulated genes\n")
cat("   - Downregulated_Genes.csv: Significantly downregulated genes\n")
cat("   - Volcano_Plot_GSE183795.png: Volcano plot visualization\n")
cat("   - Heatmap_Top_DEGs_GSE183795.png: Heatmap of top DEGs\n")

cat(rep("=", 60), "\n", sep = "")

# Save summary
summary_text <- c(
  "ANALYSIS SUMMARY - GSE183795",
  "============================================================",
  "1. PROBE-TO-GENE MAPPING:",
  paste("   - Initial probes:", nrow(processed_data)),
  paste("   - Probes with gene symbols:", nrow(gene_map_df)),
  paste("   - Final unique genes after averaging:", nrow(averaged_data)),
  "   - Multiple probes mapping to same gene were handled by averaging expression values",
  "",
  "2. DIFFERENTIAL EXPRESSION ANALYSIS:",
  "   - Comparison: Tumor vs Normal (Pancreatic Cancer)",
  paste("   - Upregulated genes:", nrow(upregulated)),
  paste("   - Downregulated genes:", nrow(downregulated)),
  paste("   - Total significant DEGs:", nrow(significant_genes)),
  "   - Criteria: |logFC| > 1 and adj.P.Val < 0.05",
  "",
  "3. NOTE:",
  "   - Significant class imbalance: 241 Tumor vs 3 Normal samples",
  "   - This may affect statistical power and interpretation"
)

writeLines(summary_text, "Results/Analysis_Summary.txt")

cat("\n🎉 DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED! 🎉\n")