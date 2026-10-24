# rna_protein_correlation.R
library(boot)

rna <- read.csv("/data/wenjiying/bean/gene_tpm_matrix.csv", row.names=1)
prot <- read.csv("/data/wenjiying/bean/protein_diff_FINAL.csv")
mapping <- read.csv("/data/wenjiying/bean/chip/prot_to_glyma_mapping.csv")

rna_genes <- gsub("\\.Wm82\\.a6\\.v1$","",rownames(rna))
prot$gene <- mapping$gene_id[match(prot$ProteinID, mapping$ProteinID)]

common <- intersect(rna_genes, prot$gene)
cat("Common genes:", length(common), "\n")

df <- data.frame(
  gene = common,
  RNA = rowMeans(rna[match(common, rna_genes),], na.rm=TRUE),
  protein = prot$log2FC[match(common, prot$gene)]
)

rho_all <- cor(df$RNA, df$protein, method="spearman")
cat("Genome-wide rho:", rho_all, "n:", nrow(df), "\n")

rho_fun <- function(data, idx) cor(data$RNA[idx], data$protein[idx], method="spearman")
set.seed(123)
boot_all <- boot(df, rho_fun, R=10000)
ci_all <- boot.ci(boot_all, type="perc")
print(ci_all)

core <- read.csv("/data/wenjiying/bean/core245_protein_level.csv")
core_genes <- unique(core$gene_id)
df_core <- df[df$gene %in% core_genes,]
cat("Core n:", nrow(df_core), "\n")

rho_core <- cor(df_core$RNA, df_core$protein, method="spearman")
cat("Core rho:", rho_core, "\n")

boot_core <- boot(df_core, rho_fun, R=10000)
ci_core <- boot.ci(boot_core, type="perc")
print(ci_core)

tab <- data.frame(
  Group = c("Genome-wide","Core"),
  n = c(nrow(df), nrow(df_core)),
  rho = c(round(rho_all,3), round(rho_core,3)),
  CI_lower = c(round(ci_all$percent[4],3), round(ci_core$percent[4],3)),
  CI_upper = c(round(ci_all$percent[5],3), round(ci_core$percent[5],3))
)
write.csv(tab, "/data/wenjiying/bean/Table_S13_RNA_protein_correlation.csv", row.names=FALSE)
print(tab)
