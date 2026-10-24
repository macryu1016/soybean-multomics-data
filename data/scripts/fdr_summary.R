# fdr_summary.R
add_fdr <- function(file, pcol, out) {
  df <- read.csv(file)
  if (!pcol %in% colnames(df)) {
    cat("Column", pcol, "not found in", file, "\n")
    cat("Available columns:", colnames(df), "\n")
    return(NULL)
  }
  df$FDR <- p.adjust(df[[pcol]], method="BH")
  write.csv(df, out, row.names=FALSE)
  cat(out, ":", sum(df$FDR<0.05, na.rm=TRUE), "significant out of", nrow(df), "\n")
}
add_fdr("/data/wenjiying/bean/rna_diff_results.csv", "P.Value", "/data/wenjiying/bean/Table_S6a_RNA_FDR.csv")
add_fdr("/data/wenjiying/bean/protein_diff_FINAL.csv", "Pvalue", "/data/wenjiying/bean/Table_S6b_protein_FDR.csv")
add_fdr("/data/wenjiying/bean/m6a_diff_results.csv", "P.Value", "/data/wenjiying/bean/Table_S6c_m6A_FDR.csv")
add_fdr("/data/wenjiying/bean/metabolite_sig.csv", "p.value", "/data/wenjiying/bean/Table_S6d_metabolite_FDR.csv")
