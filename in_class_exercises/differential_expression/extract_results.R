setwd("/work/in_class_exercises/differential_expression")

obj1 <- readRDS("data/d1_qc_objects.rds")
obj2 <- readRDS("data/d2_qc_objects.rds")
lim1 <- readRDS("data/d1_limma_results.rds")
edg1 <- readRDS("data/d1_edger_results.rds")
deq1 <- readRDS("data/d1_deseq2_results.rds")
lim2 <- readRDS("data/d2_limma_results.rds")
edg2 <- readRDS("data/d2_edger_results.rds")
deq2 <- readRDS("data/d2_deseq2_results.rds")
raw1 <- readRDS("data/d1_counts_and_meta.rds")
raw2 <- readRDS("data/d2_counts_and_meta.rds")

cat("=== D1 Dataset (GSE233947) ===\n")
cat("Genes before filter:", nrow(raw1$counts), "\n")
cat("Genes after filter:", nrow(obj1$dge_f$counts), "\n")
cat("Samples:", ncol(obj1$dge_f$counts), "\n")
cat("Treatment groups:\n"); print(table(obj1$meta$treatment))

cat("\n=== D1 limma-voom ===\n")
cat("3CTG_vs_NT sig (adj.P<0.05):", sum(lim1$res_3_vs_nt$adj.P.Val < 0.05, na.rm=TRUE), "\n")
cat("  up:", sum(lim1$res_3_vs_nt$adj.P.Val < 0.05 & lim1$res_3_vs_nt$logFC > 0, na.rm=TRUE), "\n")
cat("  down:", sum(lim1$res_3_vs_nt$adj.P.Val < 0.05 & lim1$res_3_vs_nt$logFC < 0, na.rm=TRUE), "\n")
cat("20CTG_vs_NT sig (adj.P<0.05):", sum(lim1$res_20_vs_nt$adj.P.Val < 0.05, na.rm=TRUE), "\n")
cat("  up:", sum(lim1$res_20_vs_nt$adj.P.Val < 0.05 & lim1$res_20_vs_nt$logFC > 0, na.rm=TRUE), "\n")
cat("  down:", sum(lim1$res_20_vs_nt$adj.P.Val < 0.05 & lim1$res_20_vs_nt$logFC < 0, na.rm=TRUE), "\n")
cat("20CTG_vs_3CTG sig (adj.P<0.05):", sum(lim1$res_20_vs_3$adj.P.Val < 0.05, na.rm=TRUE), "\n")
cat("  up:", sum(lim1$res_20_vs_3$adj.P.Val < 0.05 & lim1$res_20_vs_3$logFC > 0, na.rm=TRUE), "\n")
cat("  down:", sum(lim1$res_20_vs_3$adj.P.Val < 0.05 & lim1$res_20_vs_3$logFC < 0, na.rm=TRUE), "\n")

cat("\n=== D1 edgeR QLF ===\n")
cat("3CTG_vs_NT sig (FDR<0.05):", sum(edg1$tab_3_vs_nt$FDR < 0.05, na.rm=TRUE), "\n")
cat("  up:", sum(edg1$tab_3_vs_nt$FDR < 0.05 & edg1$tab_3_vs_nt$logFC > 0, na.rm=TRUE), "\n")
cat("  down:", sum(edg1$tab_3_vs_nt$FDR < 0.05 & edg1$tab_3_vs_nt$logFC < 0, na.rm=TRUE), "\n")
cat("20CTG_vs_NT sig (FDR<0.05):", sum(edg1$tab_20_vs_nt$FDR < 0.05, na.rm=TRUE), "\n")
cat("  up:", sum(edg1$tab_20_vs_nt$FDR < 0.05 & edg1$tab_20_vs_nt$logFC > 0, na.rm=TRUE), "\n")
cat("  down:", sum(edg1$tab_20_vs_nt$FDR < 0.05 & edg1$tab_20_vs_nt$logFC < 0, na.rm=TRUE), "\n")
cat("20CTG_vs_3CTG sig (FDR<0.05):", sum(edg1$tab_20_vs_3$FDR < 0.05, na.rm=TRUE), "\n")
cat("  up:", sum(edg1$tab_20_vs_3$FDR < 0.05 & edg1$tab_20_vs_3$logFC > 0, na.rm=TRUE), "\n")
cat("  down:", sum(edg1$tab_20_vs_3$FDR < 0.05 & edg1$tab_20_vs_3$logFC < 0, na.rm=TRUE), "\n")

cat("\n=== D1 DESeq2 ===\n")
cat("20CTG_vs_NT sig (padj<0.05):", sum(deq1$res_20_vs_nt$padj < 0.05, na.rm=TRUE), "\n")
cat("  up:", sum(deq1$res_20_vs_nt$padj < 0.05 & deq1$res_20_vs_nt$log2FoldChange > 0, na.rm=TRUE), "\n")
cat("  down:", sum(deq1$res_20_vs_nt$padj < 0.05 & deq1$res_20_vs_nt$log2FoldChange < 0, na.rm=TRUE), "\n")
cat("3CTG_vs_NT sig (padj<0.05):", sum(deq1$res_3_vs_nt$padj < 0.05, na.rm=TRUE), "\n")
cat("  up:", sum(deq1$res_3_vs_nt$padj < 0.05 & deq1$res_3_vs_nt$log2FoldChange > 0, na.rm=TRUE), "\n")
cat("  down:", sum(deq1$res_3_vs_nt$padj < 0.05 & deq1$res_3_vs_nt$log2FoldChange < 0, na.rm=TRUE), "\n")
cat("20CTG_vs_3CTG sig (padj<0.05):", sum(deq1$res_20_vs_3$padj < 0.05, na.rm=TRUE), "\n")
cat("  up:", sum(deq1$res_20_vs_3$padj < 0.05 & deq1$res_20_vs_3$log2FoldChange > 0, na.rm=TRUE), "\n")
cat("  down:", sum(deq1$res_20_vs_3$padj < 0.05 & deq1$res_20_vs_3$log2FoldChange < 0, na.rm=TRUE), "\n")

cat("\n=== D2 Dataset (GSE240829) ===\n")
cat("Genes before filter:", nrow(raw2$counts), "\n")
cat("Genes after filter:", nrow(obj2$dge_f$counts), "\n")
cat("Samples:", ncol(obj2$dge_f$counts), "\n")
cat("Status groups:\n"); print(table(obj2$meta$status))

cat("\n=== D2 limma-voom ===\n")
cat("status0_vs_1 sig (adj.P<0.05):", sum(lim2$res_0_vs_1$adj.P.Val < 0.05, na.rm=TRUE), "\n")
cat("  up:", sum(lim2$res_0_vs_1$adj.P.Val < 0.05 & lim2$res_0_vs_1$logFC > 0, na.rm=TRUE), "\n")
cat("  down:", sum(lim2$res_0_vs_1$adj.P.Val < 0.05 & lim2$res_0_vs_1$logFC < 0, na.rm=TRUE), "\n")

cat("\n=== D2 edgeR ===\n")
cat("edg2 names:", paste(names(edg2), collapse=", "), "\n")
# Try QLF first
if ("qlf_0_vs_1" %in% names(edg2)) {
  cat("qlf_0_vs_1 sig (FDR<0.05):", sum(edg2$qlf_0_vs_1$table$FDR < 0.05, na.rm=TRUE), "\n")
} else if ("tab_0_vs_1" %in% names(edg2)) {
  cat("tab_0_vs_1 sig (FDR<0.05):", sum(edg2$tab_0_vs_1$FDR < 0.05, na.rm=TRUE), "\n")
}

cat("\n=== D2 DESeq2 ===\n")
cat("deq2 names:", paste(names(deq2), collapse=", "), "\n")
if ("res_0_vs_1" %in% names(deq2)) {
  cat("status0_vs_1 sig (padj<0.05):", sum(deq2$res_0_vs_1$padj < 0.05, na.rm=TRUE), "\n")
  cat("  up:", sum(deq2$res_0_vs_1$padj < 0.05 & deq2$res_0_vs_1$log2FoldChange > 0, na.rm=TRUE), "\n")
  cat("  down:", sum(deq2$res_0_vs_1$padj < 0.05 & deq2$res_0_vs_1$log2FoldChange < 0, na.rm=TRUE), "\n")
}

# Top DE genes D1 20CTG_vs_NT
cat("\n=== Top 10 D1 DEGs (20CTG_vs_NT, edgeR FDR<0.05, sorted by FDR) ===\n")
top_d1 <- edg1$tab_20_vs_nt[order(edg1$tab_20_vs_nt$FDR),]
print(head(top_d1[c("logFC","FDR")], 10))

# Top DE genes D2 0_vs_1
cat("\n=== Top 10 D2 DEGs (status0_vs_1, limma sorted by adj.P.Val) ===\n")
top_d2 <- lim2$res_0_vs_1[order(lim2$res_0_vs_1$adj.P.Val),]
print(head(top_d2[c("logFC","adj.P.Val")], 10))
