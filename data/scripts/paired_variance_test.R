# paired_variance_test.R
paired <- read.csv("/data/wenjiying/bean/protein_diff_paired.csv")
tab14 <- data.frame(n_pairs = nrow(paired), df = nrow(paired) - 1)
write.csv(tab14, "/data/wenjiying/bean/Table_S14_paired_ttest_stats.csv", row.names=FALSE)
print(tab14)

core <- read.csv("/data/wenjiying/bean/core245_protein_level.csv")$logFC
all <- read.csv("/data/wenjiying/bean/protein_diff_FINAL.csv")$logFC
library(car)
lt <- leveneTest(c(core, all),
                 factor(c(rep("core", length(core)), rep("all", length(all)))))
ft <- var.test(core, all)
tab15 <- data.frame(levene_P = lt$`Pr(>F)`[1], Ftest_P = ft$p.value)
write.csv(tab15, "/data/wenjiying/bean/Table_S15_variance_test.csv", row.names=FALSE)
print(tab15)
