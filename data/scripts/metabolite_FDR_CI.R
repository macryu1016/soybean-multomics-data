# metabolite_FDR_CI.R
met <- read.csv("/data/wenjiying/bean/metabolite_sig.csv")
met$FDR <- p.adjust(met$p.value, method="BH")

cat("Total metabolites:", nrow(met), "\n")
cat("p<0.05:", sum(met$p.value<0.05, na.rm=TRUE), "\n")
cat("FDR<0.05:", sum(met$FDR<0.05, na.rm=TRUE), "\n")
cat("FDR<0.01:", sum(met$FDR<0.01, na.rm=TRUE), "\n")
cat("Min FDR:", min(met$FDR, na.rm=TRUE), "\n")

tab <- data.frame(
  Metric = c("Total","p<0.05","FDR<0.05","FDR<0.01","Min FDR","log2FC mean"),
  Value = c(nrow(met),
            sum(met$p.value<0.05, na.rm=TRUE),
            sum(met$FDR<0.05, na.rm=TRUE),
            sum(met$FDR<0.01, na.rm=TRUE),
            format(min(met$FDR, na.rm=TRUE), scientific=TRUE, digits=2),
            round(mean(met$log2FC, na.rm=TRUE),3))
)
write.csv(tab, "/data/wenjiying/bean/Table_S17b_metabolite_summary.csv", row.names=FALSE)
print(tab)

# 完整 FDR + log2FC
out <- met[, c("Name","Description","Fold.change","p.value","FDR","log2FC")]
write.csv(out, "/data/wenjiying/bean/Table_S17_metabolite_FDR.csv", row.names=FALSE)
