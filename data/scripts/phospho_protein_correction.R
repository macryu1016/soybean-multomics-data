# phospho_protein_correction.R
# Correct phosphosite log2FC for protein abundance changes
# Control = Channel_1 (113) + Channel_5 (117)
# Salt    = Channel_3 (115) + Channel_7 (119)

phos <- read.csv("/data/wenjiying/bean/P14064_phospho_site_matrix_final.csv")
prot <- read.csv("/data/wenjiying/bean/P14064_protein_quant_annotated.csv")

ctrl_cols <- c("Channel_1","Channel_5")
salt_cols <- c("Channel_3","Channel_7")

phos_ctrl <- rowMeans(phos[, ctrl_cols], na.rm=TRUE)
phos_salt <- rowMeans(phos[, salt_cols], na.rm=TRUE)
phos$log2FC_phos <- log2(phos_salt / phos_ctrl)

prot_ctrl <- rowMeans(prot[, ctrl_cols], na.rm=TRUE)
prot_salt <- rowMeans(prot[, salt_cols], na.rm=TRUE)
prot$log2FC_prot <- log2(prot_salt / prot_ctrl)

merged <- merge(phos[, c("ProteinID","Position","log2FC_phos")],
                prot[, c("ProteinID","log2FC_prot")],
                by="ProteinID")
merged$log2FC_corrected <- merged$log2FC_phos - merged$log2FC_prot

write.csv(merged, "/data/wenjiying/bean/Table_S12_phospho_protein_corrected.csv", row.names=FALSE)

tab <- data.frame(
  Metric = c("Merged sites","Mean phos log2FC","Mean prot log2FC",
             "Mean corrected","SD phos","SD prot","SD corrected"),
  Value = c(nrow(merged),
            round(mean(merged$log2FC_phos, na.rm=TRUE),3),
            round(mean(merged$log2FC_prot, na.rm=TRUE),3),
            round(mean(merged$log2FC_corrected, na.rm=TRUE),3),
            round(sd(merged$log2FC_phos, na.rm=TRUE),3),
            round(sd(merged$log2FC_prot, na.rm=TRUE),3),
            round(sd(merged$log2FC_corrected, na.rm=TRUE),3))
)
write.csv(tab, "/data/wenjiying/bean/Table_S12_summary.csv", row.names=FALSE)
print(tab)# phospho_protein_correction.R
phos <- read.csv("/data/wenjiying/bean/P14064_phospho_site_matrix_final.csv")
prot <- read.csv("/data/wenjiying/bean/P14064_protein_quant_annotated.csv")
cat("phos columns:", colnames(phos), "\n")
cat("prot columns:", colnames(prot), "\n")

# 如果有 logFC 列，直接合并
# 否则从 channel 算
