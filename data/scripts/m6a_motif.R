# m6a_motif.R
library(Biostrings)

seqs <- readDNAStringSet("/data/wenjiying/bean/core245_m6A_sites_40bp_clean.fa")
cat("Total sequences:", length(seqs), "\n")
cat("Sequence width:", unique(width(seqs)), "\n")

seqs_char <- as.character(seqs)

center <- substr(seqs_char, 19, 23)
cat("Center 5bp head:\n")
print(head(center, 10))

drach_center <- grepl("[AGT][AG]AC[ACT]", center)
cat("DRACH in center:", sum(drach_center), "/", length(center),
    "=", round(mean(drach_center)*100,1), "%\n")

drach_full <- grepl("[AGT][AG]AC[ACT]", seqs_char)
cat("DRACH in full:", sum(drach_full), "/", length(seqs_char),
    "=", round(mean(drach_full)*100,1), "%\n")

tab <- data.frame(
  Metric = c("Total sequences","DRACH in center (5bp)","DRACH in full (40bp)"),
  Value = c(length(seqs), sum(drach_center), sum(drach_full))
)
write.csv(tab, "/data/wenjiying/bean/Table_S7_m6A_motif.csv", row.names=FALSE)
print(tab)
