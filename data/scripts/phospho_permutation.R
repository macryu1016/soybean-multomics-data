# phospho_permutation.R
df <- read.csv("/data/wenjiying/bean/fig6E_phospho_heterogeneity.csv", row.names=1)
df <- df[df$n_sites >= 2,]
cat("Multi-site proteins:", nrow(df), "\n")

obs <- mean(df$heterogeneous)
cat("Observed:", obs, "\n")

set.seed(123)
perm_fun <- function(n) {
  v <- sample(c(1,-1), n, replace=TRUE)
  any(v>0) & any(v<0)
}
null <- replicate(10000, mean(sapply(df$n_sites, perm_fun)))

cat("Null mean:", mean(null), "\n")
cat("Null 95% CI:", quantile(null, c(0.025, 0.975)), "\n")
cat("Empirical P:", mean(null >= obs), "\n")

tab <- data.frame(
  observed = obs,
  null_mean = mean(null),
  null_lower = quantile(null, 0.025),
  null_upper = quantile(null, 0.975),
  empirical_P = mean(null >= obs)
)
write.csv(tab, "/data/wenjiying/bean/Table_S11_phospho_permutation.csv", row.names=FALSE)
print(tab)
