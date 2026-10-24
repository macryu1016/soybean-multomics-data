# pathway_joint_Fisher.R
gene_kegg <- read.csv("/data/wenjiying/bean/Table_S16_gene_KEGG_FDR.csv")
met_kegg <- read.csv("/data/wenjiying/bean/fig5C_metabolite_KEGG_top15.csv")

gene_kegg$Pathway_clean <- gsub(" \\[PATH:.*\\]","",gene_kegg$Pathway.Name)
met_kegg$Pathway_clean <- met_kegg$Pathway

common <- intersect(gene_kegg$Pathway_clean, met_kegg$Pathway_clean)
cat("Common pathways:", length(common), "\n")
print(common)

combined <- data.frame()
for (p in common) {
  g_p <- gene_kegg$pvalue[gene_kegg$Pathway_clean==p][1]
  m_p <- met_kegg$Raw.p[met_kegg$Pathway_clean==p][1]
  if (is.na(g_p) || is.na(m_p)) next
  stat <- -2*(log(g_p)+log(m_p))
  pval <- pchisq(stat, df=4, lower.tail=FALSE)
  combined <- rbind(combined, data.frame(
    Pathway=p, Gene_p=g_p, Metabolite_p=m_p,
    Fisher_stat=stat, Fisher_p=pval
  ))
}
combined$Fisher_FDR <- p.adjust(combined$Fisher_p, method="BH")
combined <- combined[order(combined$Fisher_p),]
write.csv(combined, "/data/wenjiying/bean/Table_S16_joint_Fisher.csv", row.names=FALSE)
print(combined)
