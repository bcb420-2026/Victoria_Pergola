setwd("/work/in_class_exercises/differential_expression")

edg1 <- readRDS("data/d1_edger_results.rds")
raw1 <- readRDS("data/d1_counts_and_meta.rds")
mappings <- raw1$mappings

# D1 has an "ensembl" column in tab objects
top10 <- head(edg1$tab_20_vs_nt[order(edg1$tab_20_vs_nt$FDR),], 10)
top10_named <- merge(top10, mappings[,c("ensembl_gene_id","hgnc_symbol")],
                     by.x="ensembl", by.y="ensembl_gene_id", all.x=TRUE)
top10_named <- top10_named[order(top10_named$FDR),]
cat("=== Top 10 D1 DEGs (20CTG_vs_NT, edgeR) with HGNC ===\n")
print(top10_named[c("ensembl","hgnc_symbol","logFC","FDR")])

top10_3 <- head(edg1$tab_3_vs_nt[order(edg1$tab_3_vs_nt$FDR),], 10)
top10_3named <- merge(top10_3, mappings[,c("ensembl_gene_id","hgnc_symbol")],
                      by.x="ensembl", by.y="ensembl_gene_id", all.x=TRUE)
top10_3named <- top10_3named[order(top10_3named$FDR),]
cat("\n=== Top 10 D1 DEGs (3CTG_vs_NT, edgeR) with HGNC ===\n")
print(top10_3named[c("ensembl","hgnc_symbol","logFC","FDR")])

# DESeq2 top DEGs
deq1 <- readRDS("data/d1_deseq2_results.rds")
deq1_20 <- as.data.frame(deq1$res_20_vs_nt)
deq1_20$ensembl <- rownames(deq1_20)
deq1_20_sig <- deq1_20[!is.na(deq1_20$padj) & deq1_20$padj < 0.05,]
deq1_20_top <- head(deq1_20_sig[order(deq1_20_sig$padj),], 10)
deq1_20_named <- merge(deq1_20_top, mappings[,c("ensembl_gene_id","hgnc_symbol")],
                        by.x="ensembl", by.y="ensembl_gene_id", all.x=TRUE)
cat("\n=== Top 10 D1 DEGs (20CTG_vs_NT, DESeq2) with HGNC ===\n")
print(deq1_20_named[order(deq1_20_named$padj), c("ensembl","hgnc_symbol","log2FoldChange","padj")])

# D2 HGNC: try stripping version suffix from count IDs
raw2 <- readRDS("data/d2_counts_and_meta.rds")
cat("\n=== D2 IDs with version suffix in counts: ===\n")
cat("First 5 rownames of counts:", paste(head(rownames(raw2$counts)), collapse=", "), "\n")
cat("Mapping ensembl_gene_ids (first 5):", paste(head(raw2$mappings$ensembl_gene_id), collapse=", "), "\n")
# The mappings have stripped IDs but BioMart returned NAs
cat("\n=== D2 HGNC mapping check ===\n")
cat("Total genes with NA:", sum(is.na(raw2$mappings$hgnc_symbol)), "\n")
cat("Total genes:", nrow(raw2$mappings), "\n")
# This means BioMart failed during D2 render and no cache existed
