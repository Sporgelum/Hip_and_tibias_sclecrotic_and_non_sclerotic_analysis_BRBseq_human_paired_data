# Argument validation
  if (!is.data.frame(counts) || !is.character(conditions) || is.null(ensembl_genes)) {
    stop("Invalid input data.")
  }

# Ensure 'conditions' is a factor
conditions <- factor(conditions)
#non_control_condition <- "non.sclerotic"  # Specify the reference condition
control_condition <- "non.sclerotic"
conditions <- relevel(conditions, ref = control_condition)

# Create the DGEList object
d <- DGEList(counts = counts,genes = ensembl_genes, group = conditions)

# Filter out lowly expressed genes
keep <- rowSums(cpm(d) >= 1) >= 2 & rowSums(d$counts) >= 30
d <- d[keep, , keep.lib.sizes = FALSE]
message("Normalizing using TMM.\n")
d <- calcNormFactors(d, method = "TMM")


#An MDS plot shows the similarities between samples based on their expression profiles.
plotMDS(d,col=as.numeric(as.factor(conditions)))
legend("topright", legend = levels(as.factor(conditions)), col = 1:length(levels(as.factor(conditions))), pch = 16)




normalized_matrix <- cpm(d, normalized.lib.sizes = TRUE)
rownames(normalized_matrix) <- d$genes$genes
# Convert patient or/and batch to factors if they are provided
if (!is.null(patient)) {
  patient <- factor(patient)
}
if (!is.null(batch)) {
  batch <- factor(batch)
  batch <- make.names(batch)
}

# Check d, specify the group column that we are interested n comparing later as the last of the model.
d
# The last variable ofthe design is the one of importance that we want to ask ourselfs the question.. the rest add to the model before.


# MAKE SURE THAT GROUP IS USED HERE #
# comparing the experimental  to the  (control) which is releveled to be the reference as baseline or intercept.
#ref = control.
#logFC condition exp /condition ref -->  if positve logFC the exp gets more expression.


# Construct the design matrix formula and matrix
if (!is.null(batch) && !is.null(patient)) {
  design_formula <- "~batch + patient + group"
  design <- model.matrix(as.formula(design_formula), data = d$samples)
} else if (is.null(patient) && !is.null(batch)) {
  design_formula <- "~batch + group"
  design <- model.matrix(as.formula(design_formula), data = d$samples)
} else if (!is.null(patient) && is.null(batch)) {
  design_formula <- "~patient + group"
  design <- model.matrix(as.formula(design_formula), data = d$samples)
} else {
  design_formula <- "~group"
  design <- model.matrix(as.formula(design_formula), data = d$samples)
}

#design <- model.matrix(~conditions)
# Output the design matrix formula
message(paste0("Using design matrix formula: ", design_formula, "\n"))
print(head(design))

# Fit the NB GLMs with QL methods
message("Fitting the QLNB function.\n")
d <- estimateDisp(d, design)





#2.
# Fit the model and perform the differential expression analysis
library(statmod)
fit <- glmQLFit(d, design,robust = TRUE)
plotQLDisp(fit)
design
# Define contrasts


#the last column is the contrasts


# contrasts <- makeContrasts(
#   sclerotic_vs_non_sclerotic = sclerotic - non.sclerotic,
#   levels = design
# )

# Test for differential expression
################################################################################
#glmQLFTest because the treatment effect is the last coefficient in the model. #
################################################################################
ql_results <- glmQLFTest(fit)#,coef=32) specify which column to check for.
as.data.frame(topTags(ql_results)) |> dplyr::filter(FDR < 0.05)
top_genes <- topTags(ql_results)[,"genes"][[1]][[1]]

#Visualize the desing and samples for one gene
head(design)

conditions

gene="IL11"

as.data.frame(topTags(ql_results,n = Inf)) |> dplyr::filter(FDR < 0.05,genes==gene)

samples_boxplot <- data.frame("samples" = names(normalized_matrix[gene,]))
samples_boxplot$values <- normalized_matrix[gene,]
#samples_boxplot$condition <- d$samples$group
#samples_boxplot$condition <- relevel(samples_boxplot$condition,"non.sclerotic")
samples_boxplot$condition <- conditions
head(samples_boxplot)
color_mapping <- c("non.sclerotic" = "red", "sclerotic" = "blue")

ggplot(samples_boxplot, aes(x = condition, y = values)) +
  geom_boxplot(aes(group = condition), outlier.shape = NA) + # Remove black outlier points Ensure boxplot is black
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, aes(color = condition)) + # Map color inside aes()
  scale_color_manual(values = color_mapping) + # Map colors to conditions
  theme_minimal() +
  labs(title = paste0(gene, " Normalized Values"), x = "Condition", y = paste0("Gene Value of ", gene))

top_genes <- topTags(ql_results)[,"genes"][[1]][[1]]
as.data.frame(topTags(ql_results, n = Inf))
normalized_matrix[top,]

ql_results$comparison



# List of top genes
top_genes <- topTags(ql_results)[,"genes"][[1]][[1]]
# Prepare the data frame for all top genes
samples_boxplot <- data.frame(
  samples = rep(names(normalized_matrix[top_genes[1],]), length(top_genes)),
  values = as.vector(t(normalized_matrix[top_genes,])),
  condition = rep(conditions, length(top_genes)),
  gene = rep(top_genes, each = ncol(normalized_matrix))
)

# Assigning colors to conditions
color_mapping <- c("non.sclerotic" = "red", "sclerotic" = "blue")

# Generate the plot with subplots for each gene
ggplot(samples_boxplot, aes(x = condition, y = values)) +
  geom_boxplot(aes(group = condition), outlier.shape = NA) + # Remove black outlier points
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, aes(color = condition)) + # Map color inside aes()
  scale_color_manual(values = color_mapping) + # Map colors to conditions
  theme_minimal() +
  labs(title = "Top Genes Normalized Values", x = "Condition", y = "Gene Value") +
  facet_wrap(~ gene, scales = "free_y") # Create subplots for each gene


rm(conditions)
rm(d)
rm(design)
rm(fit)
rm(samples_boxplot)
rm(ql_results)
rm(batch)
rm(normalized_matrix)





# Leanring by visualizing models.
# Design matrix for GLM with intercept
design <- model.matrix(~conditions)
dge <- estimateDisp(d, design)  # Estimate dispersions
# Fit GLM
fit <- glmFit(d, design)
# Visualize GLM coefficients (e.g., using plotMD)
plotMD(fit, column = 2, status = d$genes$Status, main = "MD Plot: GLM ~ conditions")


# Design matrix for GLM with intercept and patient
design_pb <- model.matrix(~patient + conditions)
d_pb <- estimateDisp(d, design_pb)  # Estimate dispersions
# Fit GLM
fit_pb <- glmFit(d_pb, design_pb)
# Visualize GLM coefficients (e.g., using plotMD)
plotMD(fit_pb, column = 2, status = d_pb$genes$Status, main = "MD Plot: GLM ~ patient + conditions")





# Design matrix for GLM without intercept
design_no_intercept <- model.matrix(~ 0 + conditions)
d_no_intercept <- estimateDisp(d, design_no_intercept)  # Estimate dispersions
# Fit GLM without intercept
fit_no_intercept <- glmFit(d_no_intercept, design_no_intercept)
# Visualize GLM coefficients (e.g., using plotMD)
plotMD(fit_no_intercept, column = 1, status = d_no_intercept$genes$Status, main = "MD Plot: GLM ~ 0 + conditions")


# Design matrix for GLM without intercept and patient after 
design_no_intercept_patient_after <- model.matrix(~ 0 + conditions + patient)
d_no_intercept_patient_after <- estimateDisp(d, design_no_intercept_patient_after)  # Estimate dispersions
# Fit GLM without intercept
fit_no_intercept_patient_after <- glmFit(d_no_intercept_patient_after, design_no_intercept_patient_after)
# Visualize GLM coefficients (e.g., using plotMD)
plotMD(fit_no_intercept_patient_after, column = 1, status = d_no_intercept_patient_after$genes$Status, main = "MD Plot: GLM ~ 0 + conditions + patient")


# Design matrix for GLM without intercept and patient before
design_no_intercept_patient_before <- model.matrix(~ 0 + patient + conditions)
d_no_intercept_patient_before <- estimateDisp(d, design_no_intercept_patient_before)  # Estimate dispersions
# Fit GLM without intercept
fit_no_intercept_patient_before <- glmFit(d_no_intercept_patient_before, design_no_intercept_patient_before)
# Visualize GLM coefficients (e.g., using plotMD)
plotMD(fit_no_intercept_patient_before, column = 1, status = d_no_intercept_patient_before$genes$Status, main = "MD Plot: GLM ~ 0 + patient + conditions")




# add plot 
#boxplot
#y is gene expression normalzed_matrix["TCN2",]
# x is sample. biological_rep  normalized_matrix[,d$samples]
#plot result

samples_boxplot <- data.frame("samples" = names(normalized_matrix["TCN2",]))
samples_boxplot$values <- normalized_matrix["TCN2",]
samples_boxplot$condition <- d$samples$group
samples_boxplot

ggplot(samples_boxplot, aes(x = condition, y = values)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  theme_minimal() +
  labs(title = "Boxplot of TCN2 Gene Values", x = "Condition", y = "Gene Value")



# keep this in mind in edger, log scale in dge list and log2 in result toptags eextracted
# check toptags
#to check  the reference is lower or higher.



