# chip_peak_CI_FDR.R
ctrl <- c(8241, 17757, 11971)
salt <- c(3166, 2490, 2016)

ci_t <- function(x) {
  m <- mean(x); se <- sd(x)/sqrt(length(x))
  tcrit <- qt(0.975, df = length(x)-1)
  c(mean=m, lower=m-tcrit*se, upper=m+tcrit*se)
}
ci_boot <- function(x, n=10000) {
  set.seed(123)
  m <- replicate(n, mean(sample(x, replace=TRUE)))
  quantile(m, c(0.025, 0.975))
}
diff_boot <- function(a, b, n=10000) {
  set.seed(123)
  d <- replicate(n, mean(sample(a, replace=TRUE)) - mean(sample(b, replace=TRUE)))
  quantile(d, c(0.025, 0.975))
}

ctrl_t <- ci_t(ctrl); salt_t <- ci_t(salt)
ctrl_b <- ci_boot(ctrl); salt_b <- ci_boot(salt)
diff_b <- diff_boot(ctrl, salt)

df <- read.csv("/data/wenjiying/bean/chip/fig1B_diffbind_peaks_final.csv")
fdr05 <- sum(df$FDR < 0.05, na.rm=TRUE)
fdr10 <- sum(df$FDR < 0.1, na.rm=TRUE)
minfdr <- min(df$FDR, na.rm=TRUE)

tab <- data.frame(
  Group = c("Control", "Salt", "Difference"),
  Mean = c(round(ctrl_t["mean"],1), round(salt_t["mean"],1),
           round(ctrl_t["mean"]-salt_t["mean"],1)),
  CI_t = c(paste0(round(ctrl_t["lower"],1)," - ",round(ctrl_t["upper"],1)),
           paste0(round(salt_t["lower"],1)," - ",round(salt_t["upper"],1)), "—"),
  CI_boot = c(paste0(round(ctrl_b[1],1)," - ",round(ctrl_b[2],1)),
              paste0(round(salt_b[1],1)," - ",round(salt_b[2],1)),
              paste0(round(diff_b[1],1)," - ",round(diff_b[2],1)))
)
write.csv(tab, "/data/wenjiying/bean/Table_S1b_H3K27me3_peak_CI.csv", row.names=FALSE)
cat("FDR<0.05:", fdr05, "FDR<0.1:", fdr10, "min FDR:", minfdr, "\n")
print(tab)
