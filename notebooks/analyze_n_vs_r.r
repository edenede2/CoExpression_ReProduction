#!/usr/bin/env Rscript

# Analysis of n_donors vs |r| relationships
library(data.table)
library(ggplot2)

# Read the summary data
dt <- fread("out_corshrink/summary_metrics.csv")

# Focus on between-tissue analysis (more variation in n_donors)
between_dt <- dt[level == "between"]

print("=== Between-Tissue Analysis: n_donors vs |r| Relationships ===")
print(paste("Number of tissue pairs analyzed:", nrow(between_dt)))
print(paste("Range of donor counts:", min(between_dt$n_donors), "to", max(between_dt$n_donors)))

# Calculate correlations between n_donors and various metrics
correlations <- data.table(
  metric = c("mean_abs_r_pre", "mean_abs_r_post", "var_abs_r_pre", "var_abs_r_post"),
  correlation_with_n_donors = c(
    cor(between_dt$n_donors, between_dt$mean_abs_r_pre, use = "complete.obs"),
    cor(between_dt$n_donors, between_dt$mean_abs_r_post, use = "complete.obs"), 
    cor(between_dt$n_donors, between_dt$var_abs_r_pre, use = "complete.obs"),
    cor(between_dt$n_donors, between_dt$var_abs_r_post, use = "complete.obs")
  )
)

print("\n=== Correlations between n_donors and |r| metrics ===")
print(correlations)

# Create visualizations
# 1. Scatter plot: n_donors vs mean_abs_r_pre
p1 <- ggplot(between_dt, aes(x = n_donors, y = mean_abs_r_pre)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    title = "Number of Donors vs Mean |r| (Pre-shrinkage)",
    subtitle = paste("Correlation:", round(correlations[metric == "mean_abs_r_pre", correlation_with_n_donors], 3)),
    x = "Number of Overlapping Donors",
    y = "Mean |r| (Pre-shrinkage)"
  ) +
  theme_minimal()

ggsave("out_corshrink/n_donors_vs_mean_abs_r_pre.png", p1, width = 8, height = 6, dpi = 150)

# 2. Scatter plot: n_donors vs mean_abs_r_post  
p2 <- ggplot(between_dt, aes(x = n_donors, y = mean_abs_r_post)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Number of Donors vs Mean |r| (Post-shrinkage)", 
    subtitle = paste("Correlation:", round(correlations[metric == "mean_abs_r_post", correlation_with_n_donors], 3)),
    x = "Number of Overlapping Donors",
    y = "Mean |r| (Post-shrinkage)"
  ) +
  theme_minimal()

ggsave("out_corshrink/n_donors_vs_mean_abs_r_post.png", p2, width = 8, height = 6, dpi = 150)

# 3. Scatter plot: n_donors vs var_abs_r_pre
p3 <- ggplot(between_dt, aes(x = n_donors, y = var_abs_r_pre)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "green") +
  labs(
    title = "Number of Donors vs Variance of |r| (Pre-shrinkage)",
    subtitle = paste("Correlation:", round(correlations[metric == "var_abs_r_pre", correlation_with_n_donors], 3)),
    x = "Number of Overlapping Donors", 
    y = "Variance of |r| (Pre-shrinkage)"
  ) +
  theme_minimal()

ggsave("out_corshrink/n_donors_vs_var_abs_r_pre.png", p3, width = 8, height = 6, dpi = 150)

# 4. Scatter plot: n_donors vs var_abs_r_post
p4 <- ggplot(between_dt, aes(x = n_donors, y = var_abs_r_post)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "purple") +
  labs(
    title = "Number of Donors vs Variance of |r| (Post-shrinkage)",
    subtitle = paste("Correlation:", round(correlations[metric == "var_abs_r_post", correlation_with_n_donors], 3)),
    x = "Number of Overlapping Donors",
    y = "Variance of |r| (Post-shrinkage)"
  ) +
  theme_minimal()

ggsave("out_corshrink/n_donors_vs_var_abs_r_post.png", p4, width = 8, height = 6, dpi = 150)

# 5. Combined plot showing pre vs post shrinkage effects
combined_dt <- rbind(
  data.table(n_donors = between_dt$n_donors, 
             mean_abs_r = between_dt$mean_abs_r_pre,
             var_abs_r = between_dt$var_abs_r_pre,
             stage = "Pre-shrinkage"),
  data.table(n_donors = between_dt$n_donors,
             mean_abs_r = between_dt$mean_abs_r_post, 
             var_abs_r = between_dt$var_abs_r_post,
             stage = "Post-shrinkage")
)

p5 <- ggplot(combined_dt, aes(x = n_donors, y = mean_abs_r, color = stage)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Shrinkage Effect: Number of Donors vs Mean |r|",
    x = "Number of Overlapping Donors",
    y = "Mean |r|",
    color = "Shrinkage Stage"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("out_corshrink/shrinkage_effect_mean_abs_r.png", p5, width = 10, height = 7, dpi = 150)

# Summary statistics by donor count bins
between_dt[, donor_bin := cut(n_donors, 
                              breaks = quantile(n_donors, probs = seq(0, 1, 0.25)), 
                              include.lowest = TRUE,
                              labels = c("Low (Q1)", "Mid-Low (Q2)", "Mid-High (Q3)", "High (Q4)"))]

print("\n=== Summary by Donor Count Quartiles ===")
summary_by_bin <- between_dt[, .(
  count = .N,
  mean_n_donors = round(mean(n_donors), 1),
  mean_abs_r_pre = round(mean(mean_abs_r_pre, na.rm = TRUE), 4),
  mean_abs_r_post = round(mean(mean_abs_r_post, na.rm = TRUE), 4),
  mean_var_abs_r_pre = round(mean(var_abs_r_pre, na.rm = TRUE), 6),
  mean_var_abs_r_post = round(mean(var_abs_r_post, na.rm = TRUE), 6)
), by = donor_bin]

print(summary_by_bin)

# Statistical tests
print("\n=== Statistical Significance Tests ===")
print("Linear regression: n_donors vs mean_abs_r_pre")
lm1 <- lm(mean_abs_r_pre ~ n_donors, data = between_dt)
print(summary(lm1))

print("\nLinear regression: n_donors vs var_abs_r_pre") 
lm2 <- lm(var_abs_r_pre ~ n_donors, data = between_dt)
print(summary(lm2))

print("\n=== Analysis Complete ===")
print("Generated plots:")
print("  - n_donors_vs_mean_abs_r_pre.png")
print("  - n_donors_vs_mean_abs_r_post.png") 
print("  - n_donors_vs_var_abs_r_pre.png")
print("  - n_donors_vs_var_abs_r_post.png")
print("  - shrinkage_effect_mean_abs_r.png")