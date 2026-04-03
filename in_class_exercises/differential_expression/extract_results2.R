setwd("/work/in_class_exercises/differential_expression")

edg1 <- readRDS("data/d1_edger_results.rds")
lim1 <- readRDS("data/d1_limma_results.rds")
deq1 <- readRDS("data/d1_deseq2_results.rds")
raw1 <- readRDS("data/d1_counts_and_meta.rds")
raw2 <- readRDS("data/d2_counts_and_meta.rds")

cat("=== D1 mappings structure ===\n")
cat("mappings columns:", paste(names(edg1$counts_mapped_filtered), collapse=", "), "\n")

# Top 10 D1 20CTG_vs_NT by edgeR with gene names
cat("\n=== Top 10 D1 DEGs (20CTG_vs_NT, edgeR) with rownames ===\n")
top_d1 <- edg1$tab_20_vs_nt[order(edg1$tab_20_vs_nt$FDR),]
print(head(top_d1[c("logFC","FDR")], 10))

# Top 10 D1 3CTG_vs_NT
cat("\n=== Top 10 D1 DEGs (3CTG_vs_NT, edgeR) with rownames ===\n")
top_d1_3 <- edg1$tab_3_vs_nt[order(edg1$tab_3_vs_nt$FDR),]
print(head(top_d1_3[c("logFC","FDR")], 10))

# Check D1 mapping stats
cat("\n=== D1 mapping stats ===\n")
cat("raw1 names:", paste(names(raw1), collapse=", "), "\n")
cat("raw1 mappings head:\n"); print(head(raw1$mappings))
cat("total mapped:", sum(!is.na(raw1$mappings$hgnc_symbol) & raw1$mappings$hgnc_symbol != ""), "\n")
cat("total rows:", nrow(raw1$mappings), "\n")

cat("\n=== D2 mapping stats ===\n")
cat("raw2 names:", paste(names(raw2), collapse=", "), "\n")
cat("raw2 mappings head:\n"); print(head(raw2$mappings))
cat("total mapped:", sum(!is.na(raw2$mappings$hgnc_symbol) & raw2$mappings$hgnc_symbol != ""), "\n")
cat("total rows:", nrow(raw2$mappings), "\n")
