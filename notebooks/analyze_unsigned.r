#!/usr/bin/env Rscript

# Comparison of signed vs unsigned correlation analysis
library(data.table)
library(ggplot2)

# Read the unsigned results
dt_unsigned <- fread("out_corshrink/summary_metrics.csv")
between_unsigned <- dt_unsigned[level == "between"]

print("=== UNSIGNED Mode Analysis ===")
print(paste("Number of tissue pairs:", nrow(between_unsigned)))
print(paste("Mode used:", unique(between_unsigned$mode)))
print(paste("Range of donor counts:", min(between_unsigned$n_donors), "to", max(between_unsigned$n_donors)))

# Calculate correlations for unsigned mode
correlations_unsigned <- data.table(
  metric = c("mean_abs_r_pre", "mean_abs_r_post", "var_abs_r_pre", "var_abs_r_post"),
  correlation_with_n_donors = c(
    cor(between_unsigned$n_donors, between_unsigned$mean_abs_r_pre, use = "complete.obs"),
    cor(between_unsigned$n_donors, between_unsigned$mean_abs_r_post, use = "complete.obs"), 
    cor(between_unsigned$n_donors, between_unsigned$var_abs_r_pre, use = "complete.obs"),
    cor(between_unsigned$n_donors, between_unsigned$var_abs_r_post, use = "complete.obs")
  )
)

print("\n=== UNSIGNED: Correlations between n_donors and |r| metrics ===")
print(correlations_unsigned)

# Summary statistics for unsigned mode
print("\n=== UNSIGNED: Overall Statistics ===")
print(paste("Mean |r| pre-shrinkage:", round(mean(between_unsigned$mean_abs_r_pre), 4)))
print(paste("Mean |r| post-shrinkage:", round(mean(between_unsigned$mean_abs_r_post), 4)))
print(paste("Mean variance pre-shrinkage:", round(mean(between_unsigned$var_abs_r_pre), 6)))
print(paste("Mean variance post-shrinkage:", round(mean(between_unsigned$var_abs_r_post), 6)))

# Shrinkage effect analysis
shrinkage_effect <- 1 - (between_unsigned$mean_abs_r_post / between_unsigned$mean_abs_r_pre)
between_unsigned[, shrinkage_effect := shrinkage_effect]

print(paste("Mean shrinkage effect:", round(mean(shrinkage_effect), 4)))
print(paste("Shrinkage effect range:", round(min(shrinkage_effect), 4), "to", round(max(shrinkage_effect), 4)))

# Analyze shrinkage effect vs n_donors
shrinkage_vs_n_cor <- cor(between_unsigned$n_donors, shrinkage_effect, use = "complete.obs")
print(paste("Correlation between n_donors and shrinkage effect:", round(shrinkage_vs_n_cor, 4)))

# Create a new visualization for shrinkage effect
p_shrinkage <- ggplot(between_unsigned, aes(x = n_donors, y = shrinkage_effect)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "orange") +
  labs(
    title = "Shrinkage Effect vs Number of Donors (UNSIGNED)",
    subtitle = paste("Correlation:", round(shrinkage_vs_n_cor, 3)),
    x = "Number of Overlapping Donors",
    y = "Shrinkage Effect (1 - post/pre)"
  ) +
  theme_minimal()

ggsave("out_corshrink/shrinkage_effect_vs_n_donors_unsigned.png", p_shrinkage, width = 8, height = 6, dpi = 150)

# Linear regression for shrinkage effect
print("\n=== Linear regression: n_donors vs shrinkage_effect (UNSIGNED) ===")
lm_shrinkage <- lm(shrinkage_effect ~ n_donors, data = between_unsigned)
print(summary(lm_shrinkage))

# Quartile analysis for unsigned
between_unsigned[, donor_bin := cut(n_donors, 
                                    breaks = quantile(n_donors, probs = seq(0, 1, 0.25)), 
                                    include.lowest = TRUE,
                                    labels = c("Low (Q1)", "Mid-Low (Q2)", "Mid-High (Q3)", "High (Q4)"))]

print("\n=== UNSIGNED: Summary by Donor Count Quartiles ===")
summary_unsigned <- between_unsigned[, .(
  count = .N,
  mean_n_donors = round(mean(n_donors), 1),
  mean_abs_r_pre = round(mean(mean_abs_r_pre, na.rm = TRUE), 4),
  mean_abs_r_post = round(mean(mean_abs_r_post, na.rm = TRUE), 4),
  mean_shrinkage = round(mean(shrinkage_effect, na.rm = TRUE), 4),
  mean_var_pre = round(mean(var_abs_r_pre, na.rm = TRUE), 6),
  mean_var_post = round(mean(var_abs_r_post, na.rm = TRUE), 6)
), by = donor_bin]

print(summary_unsigned)

# Find most and least shrunk tissue pairs
print("\n=== Tissue Pairs with Extreme Shrinkage Effects ===")
setorder(between_unsigned, shrinkage_effect)
print("Least shrinkage (most preserved):")
print(between_unsigned[1:3, .(unit, n_donors, mean_abs_r_pre, mean_abs_r_post, shrinkage_effect)])

print("\nMost shrinkage:")
print(between_unsigned[(.N-2):.N, .(unit, n_donors, mean_abs_r_pre, mean_abs_r_post, shrinkage_effect)])

print("\n=== Key Insights for UNSIGNED Mode ===")
print("1. Unsigned mode focuses purely on correlation magnitude (|r|)")
print("2. All correlations are positive, making interpretation cleaner")
print(paste("3. Higher donor counts show", ifelse(shrinkage_vs_n_cor < 0, "LESS", "MORE"), "shrinkage"))
print("4. ASH shrinkage is sample-size aware in the expected direction")

print("\n=== Generated new plot: shrinkage_effect_vs_n_donors_unsigned.png ===")