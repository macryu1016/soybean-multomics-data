# m6a_FDR_QC.R
df <- read.csv("/data/wenjiying/bean/m6A_site_level_diff_full.csv")
df$FDR <- p.adjust(df$PValue, method="BH")

tab <- data.frame(
  Metric = c("Total sites","Nominal P<0.05","FDR<0.05","FDR<0.1","Minimum FDR"),
  Value = c(nrow(df),
            sum(df$PValue<0.05,na.rm=TRUE),
            sum(df$FDR<0.05,na.rm=TRUE),
            sum(df$FDR<0.1,na.rm=TRUE),
            round(min(df$FDR,na.rm=TRUE),3))
)
write.csv(tab, "/data/wenjiying/bean/Table_S8_m6A_diff_FDR.csv", row.names=FALSE)
print(tab)

edit <- read.csv("/data/wenjiying/bean/m6A_edit_count_matrix.csv", row.names=1)
total <- read.csv("/data/wenjiying/bean/m6A_total_count_matrix.csv", row.names=1)
rate <- edit/total

qc <- data.frame(
  sample = colnames(edit),
  n_sites = colSums(total>=5, na.rm=TRUE),
  mean_rate = colMeans(rate, na.rm=TRUE),
  median_rate = apply(rate,2,median,na.rm=TRUE),
  mean_coverage = colMeans(total, na.rm=TRUE)
)
write.csv(qc, "/data/wenjiying/bean/Table_S9_DARTseq_QC.csv", row.names=FALSE)
print(qc)
