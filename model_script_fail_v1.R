test_diff_exp <- function(counts, conditions, ensembl_genes, batch = NULL, non_control_condition, patient) {
  
  # Argument validation
  if (!is.data.frame(counts) || !is.character(conditions) || is.null(ensembl_genes)) {
    stop("Invalid input data.")
  }
  
  # Ensure 'conditions' is a factor with 'non_control_condition' as the reference level
  conditions <- factor(conditions)
  conditions <- relevel(conditions, ref = non_control_condition)
  
  # Create DGEList object and filter the genes
  d <- DGEList(counts = counts, group = conditions, genes = ensembl_genes)
  
  # Filter out lowly expressed genes
  keep <- rowSums(cpm(d) >= 1) >= 2 & rowSums(d$counts) >= 30
  d <- d[keep, , keep.lib.sizes = FALSE]
  
  # Perform TMM normalization
  message("Normalizing using TMM.\n")
  d <- calcNormFactors(d, method = "TMM")
  normalized_matrix <- cpm(d, normalized.lib.sizes = TRUE)
  rownames(normalized_matrix) <- d$genes$genes
  
  # Convert patient and batch to factors if provided
  if (!is.null(patient)) {
    patient <- factor(patient)
  }
  if (!is.null(batch)) {
    batch <- factor(batch)
    batch <- make.names(batch)
  }
  
  # Construct the design matrix formula and matrix
  if (!is.null(batch) && !is.null(patient)) {
    design_formula <- "~ conditions + batch + patient"
  } else if (!is.null(batch)) {
    design_formula <- "~ conditions + batch"
  } else if (!is.null(patient)) {
    design_formula <- "~ conditions + patient"
  } else {
    design_formula <- "~ conditions"
  }
  
  # Create the design matrix
  design <- model.matrix(as.formula(design_formula))
  
  # Output the design matrix formula
  message(paste0("Using design matrix formula: ", design_formula, "\n"))
  
  # Fit the model
  message(paste0("Fitting the QLNB function.\n"))
  d <- estimateDisp(d, design)
  fit <- glmQLFit(d, design)
  
  # Generate contrasts programmatically
  condition_levels <- levels(conditions)
  contrasts <- matrix(0, nrow = length(condition_levels) - 1, ncol = ncol(design))
  rownames(contrasts) <- condition_levels[-1]  # Exclude reference level
  colnames(contrasts) <- colnames(design)
  
  for (i in seq_along(condition_levels[-1])) {
    contrasts[i, which(colnames(design) == paste0("conditions", condition_levels[i + 1]))] <- 1
    contrasts[i, which(colnames(design) == paste0("conditions", condition_levels[1]))] <- -1
  }
  
  # Perform differential expression analysis for each contrast
  results_list <- list()
  for (i in seq_len(nrow(contrasts))) {
    contrast <- contrasts[i, ]
    qlf <- glmQLFTest(fit, contrast = contrast)
    
    # Store the results
    comparison_name <- rownames(contrasts)[i]
    message(paste0("Running contrast on: ", comparison_name))
    
    # Calculate row sums for the comparison
    value1 <- strsplit(comparison_name, "-")[[1]][1]
    value2 <- strsplit(comparison_name, "-")[[1]][2]
    
    subset_rows_v1 <- d$samples[d$samples$group == value1, , drop = FALSE]
    subset_rows_v2 <- d$samples[d$samples$group == value2, , drop = FALSE]
    
    sample_names_v1 <- rownames(subset_rows_v1)
    sample_names_v2 <- rownames(subset_rows_v2)
    
    colSum_sample_v1 <- as.matrix(rowSums(d$counts[, sample_names_v1]))
    colSum_sample_v2 <- as.matrix(rowSums(d$counts[, sample_names_v2]))
    
    colnames(colSum_sample_v1) <- paste0("colSums_", subset_rows_v1$group[1])
    colnames(colSum_sample_v2) <- paste0("colSums_", subset_rows_v2$group[1])
    
    colsums_samples <- as.data.frame(d$counts[, c(sample_names_v1, sample_names_v2)])
    colsums_samples$genes <- d$genes$genes
    colsums_samples$C1 <- colSum_sample_v1
    colsums_samples$C2 <- colSum_sample_v2
    
    names_to_add_after_colsums <- colnames(colsums_samples)
    names_to_add_after_colsums[length(names_to_add_after_colsums) - 1] <- paste0("colSums_", subset_rows_v1$group[1])
    names_to_add_after_colsums[length(names_to_add_after_colsums)] <- paste0("colSums_", subset_rows_v2$group[1])
    colnames(colsums_samples) <- names_to_add_after_colsums
    
    # Store the results in the results_list
    results_list[[i]] <- list(name = comparison_name, results = topTags(qlf, n = Inf), colsums_samples = colsums_samples, normalized_matrix = normalized_matrix)
  }
  
  return(results_list)
}
