# tf_activity_CI_FDR.R
tf <- read.csv("/data/wenjiying/bean/TF_activity_matrix.csv", row.names=1)
ctrl_cols <- colnames(tf)[c(2,4,6,8,10,12)]
salt_cols <- colnames(tf)[c(1,3,5,7,9,11)]

delta <- rowMeans(tf[,salt_cols],na.rm=TRUE) - rowMeans(tf[,ctrl_cols],na.rm=TRUE)

tab <- data.frame(
  Metric = c("Total TFs","|delta|>3","|delta|>4","90th percentile","95th percentile"),
  Value = c(length(delta), sum(abs(delta)>3,na.rm=TRUE),
            sum(abs(delta)>4,na.rm=TRUE),
            round(quantile(abs(delta),0.90,na.rm=TRUE),3),
            round(quantile(abs(delta),0.95,na.rm=TRUE),3))
)
write.csv(tab, "/data/wenjiying/bean/Table_S10a_TF_activity_threshold.csv", row.names=FALSE)
print(tab)

tfe <- read.csv("/data/wenjiying/bean/TF_enrich_with_RNA.csv")
tab2 <- data.frame(
  Metric = c("Total TFs","q<0.05","q<0.01","p<0.05"),
  Value = c(nrow(tfe), sum(tfe$q_value<0.05,na.rm=TRUE),
            sum(tfe$q_value<0.01,na.rm=TRUE), sum(tfe$p_value<0.05,na.rm=TRUE))
)
write.csv(tab2, "/data/wenjiying/bean/Table_S10b_TF_FDR.csv", row.names=FALSE)
print(tab2)

phos_tf <- read.csv("/data/wenjiying/bean/all_TF_phospho_stats.csv")
tab3 <- data.frame(
  Metric = c("TFs with phospho data","TFs with >=1 site","Background rate"),
  Value = c(nrow(phos_tf), sum(phos_tf$n_sites>0,na.rm=TRUE),
            paste0(round(mean(phos_tf$n_sites>0,na.rm=TRUE)*100,1),"%"))
)
write.csv(tab3, "/data/wenjiying/bean/Table_S10c_TF_phospho_background.csv", row.names=FALSE)
print(tab3)
