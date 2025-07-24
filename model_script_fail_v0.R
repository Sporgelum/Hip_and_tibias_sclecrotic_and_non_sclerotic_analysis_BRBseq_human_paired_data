# Combine metadata into a data frame
df <- data.frame(
  age = age,
  sex = sex,
  conditions = conditions,
  batch = batch
)

# Combine with transposed counts
df_counts <- cbind(df, t(counts))

# Convert categorical variables to factors
df_counts$conditions <- factor(df_counts$conditions)
df_counts$batch <- factor(df_counts$batch)
df_counts$sex <- factor(df_counts$sex)

# Convert factors to dummy variables
dummy_vars <- model.matrix(~ conditions + batch + sex - 1, data = df_counts)
df_counts_numeric <- cbind(df_counts[, sapply(df_counts, is.numeric)], dummy_vars)
# Create a scatterplot matrix


# Save the scatterplot matrix to a PNG file
# Increase margins (bottom, left, top, right)
par(mar = c(1, 1, 1, 1)) 

# Try plotting again
png("scatterplot_matrix_with_margins.png", width = 1200, height = 1200)
pairs(df_counts_numeric)
dev.off()


# Calculate the correlation matrix
corr_matrix <- cor(df_counts_numeric, use = "pairwise.complete.obs")

# Plot heatmap
png("heatmap_corr_matrix.png", width = 1200, height = 1200)
heatmap(corr_matrix, main = "Correlation Heatmap")
dev.off()




