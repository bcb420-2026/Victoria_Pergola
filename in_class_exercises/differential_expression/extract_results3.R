setwd("/work/in_class_exercises/differential_expression")

edg1 <- readRDS("data/d1_edger_results.rds")
raw1 <- readRDS("data/d1_counts_and_meta.rds")

# Check structure of the tab objects
cat("tab_20_vs_nt columns:", paste(names(edg1$tab_20_vs_nt), collapse=", "), "\n")
cat("First few rows:\n")
print(head(edg1$tab_20_vs_nt))

# Are rownames the gene IDs?
cat("Rownames (first 5):", paste(head(rownames(edg1$tab_20_vs_nt)), collapse=", "), "\n")

# Map ENSEMBL IDs to HGNC using the mappings
mappings <- raw1$mappings
top10 <- head(edg1$tab_20_vs_nt[order(edg1$tab_20_vs_nt$FDR),], 10)
top10$ensembl_gene_id <- rownames(top10)
top10_named <- merge(top10, mappings[,c("ensembl_gene_id","hgnc_symbol")], by="ensembl_gene_id", all.x=TRUE)
top10_named <- top10_named[order(top10_named$FDR),]
cat("\n=== Top 10 D1 DEGs (20CTG_vs_NT, edgeR) with HGNC ===\n")
print(top10_named[c("hgnc_symbol","logFC","FDR")])

# Also for 3CTG_vs_NT
top10_3 <- head(edg1$tab_3_vs_nt[order(edg1$tab_3_vs_nt$FDR),], 10)
top10_3$ensembl_gene_id <- rownames(top10_3)
top10_3named <- merge(top10_3, mappings[,c("ensembl_gene_id","hgnc_symbol")], by="ensembl_gene_id", all.x=TRUE)
top10_3named <- top10_3named[order(top10_3named$FDR),]
cat("\n=== Top 10 D1 DEGs (3CTG_vs_NT, edgeR) with HGNC ===\n")
print(top10_3named[c("hgnc_symbol","logFC","FDR")])

# D1 mapping proportion
cat("\n=== D1 mapping summary ===\n")
total <- nrow(mappings)
have_hgnc <- sum(!is.na(mappings$hgnc_symbol) & mappings$hgnc_symbol != "")
cat("Total genes:", total, "\n")
cat("With HGNC symbol:", have_hgnc, "\n")
cat("Proportion:", round(have_hgnc/total, 4), "\n")

# D2: check whether HGNC mapping completely failed
raw2 <- readRDS("data/d2_counts_and_meta.rds")
cat("\n=== D2 HGNC mapping ===\n")
cat("Total genes:", nrow(raw2$mappings), "\n")
cat("With HGNC:", sum(!is.na(raw2$mappings$hgnc_symbol) & raw2$mappings$hgnc_symbol != ""), "\n")
cat("Sample IDs from raw2:", head(raw2$mappings$ensembl_gene_id, 5), "\n")
# Check D2 counts matrix IDs
cat("D2 count rownames (first 5):", paste(head(rownames(raw2$counts)), collapse=", "), "\n")
