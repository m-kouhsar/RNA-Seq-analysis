# Load necessary libraries
library(limma)
library(edgeR)
library(stringr)

args <- commandArgs(T)

counts.file <- trimws(args[1])
pheno.file <- trimws(args[2])
var_of_interest <- trimws(args[3])
covariates <- trimws(str_split_1(trimws(args[4]),pattern = ","))
n_perms <- as.numeric(trimws(args[5]))
outliers <- trimws(str_split_1(trimws(args[6]),pattern = ","))
gFilter.min.count <- as.numeric(trimws(args[7]))
gFilter.min.prop <- as.numeric(trimws(args[8]))
OutPrefix <- trimws(args[9])

# ==============================================================================
# 1. SETUP: Define your inputs
# ==============================================================================
# counts: Your matrix of raw counts (Rows = Genes, Cols = Samples)
# pheno:  Your metadata data frame
# var_of_interest: The column name for your Disease/Condition (e.g., "Group")
# covariates: Vector of covariate names (e.g., c("Age", "Sex", "RIN", "PMI"))
# n_perms: Number of permutations (1000 is standard for validation)

run_permutation_test <- function(counts, pheno, var_of_interest, covariates, n_perms=1000) {
  
  cat(paste0("Starting Permutation Test for ", nrow(counts), " genes...\n"))
  
  # -------------------------------------------------------
  # A. Prepare Data (Normalization)
  # -------------------------------------------------------
  # Create DGEList and normalize (TMM is standard)
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge)
  
  # Create the Formula for the Linear Model
  # We construct the formula dynamically based on your inputs
  cov_formula <- paste(covariates, collapse = " + ")
  design_formula_str <- paste("~", var_of_interest, "+", cov_formula)
  design_formula <- as.formula(design_formula_str)
  
  # -------------------------------------------------------
  # B. Run the "True" Model (Observed Data)
  # -------------------------------------------------------
  design <- model.matrix(design_formula, data = pheno)
  v <- voom(dge, design, plot=FALSE) # Voom transforms counts for linear modeling
  
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # Extract the t-statistic for your Condition (usually the 2nd column)
  # NOTE: Check if your Condition is the 2nd column in 'design'. 
  # Usually column 1 is Intercept, column 2 is ConditionCase.
  obs_t_stats <- fit$t[, 2] 
  
  # Store observed p-values for comparison later
  obs_p_values <- fit$p.value[, 2]

  # -------------------------------------------------------
  # C. Run the Permutations (The Loop)
  # -------------------------------------------------------
  # Matrix to store t-stats from random runs (Rows = Genes, Cols = Permutations)
  perm_t_matrix <- matrix(NA, nrow=nrow(counts), ncol=n_perms)
  
  set.seed(123) # Set seed for reproducibility
  
  for(i in 1:n_perms) {
    if(i %% 100 == 0) cat(paste("Permutation:", i, "/", n_perms, "\n"))
    
    # 1. Shuffle the Condition labels ONLY
    # We keep Age/Sex/RIN linked to the specific sample, just swapping who is "Case" vs "Control"
    shuffled_pheno <- pheno
    shuffled_pheno[[var_of_interest]] <- sample(shuffled_pheno[[var_of_interest]])
    
    # 2. Re-create design matrix with shuffled groups
    design_perm <- model.matrix(design_formula, data = shuffled_pheno)
    
    # 3. Re-run model (Fast version)
    # We don't need to re-run voom() every time, the weights are roughly stable
    fit_perm <- lmFit(v, design_perm)
    fit_perm <- eBayes(fit_perm)
    
    # 4. Save the t-statistics
    perm_t_matrix[, i] <- fit_perm$t[, 2]
  }
  
  # -------------------------------------------------------
  # D. Calculate Empirical P-Values
  # -------------------------------------------------------
  # Formula: (Number of random t-stats >= observed t-stat) + 1  /  (Total Perms + 1)
  # We use absolute values for a two-tailed test
  
  emp_p_values <- numeric(nrow(counts))
  
  for(g in 1:nrow(counts)) {
    # How many random t-stats were more extreme than my real result?
    n_extreme <- sum(abs(perm_t_matrix[g, ]) >= abs(obs_t_stats[g]))
    
    # Add pseudocount of 1 to avoid p-value of 0
    emp_p_values[g] <- (n_extreme + 1) / (n_perms + 1)
  }
  
  # -------------------------------------------------------
  # E. Create Result Table
  # -------------------------------------------------------
  results <- data.frame(
    Gene = rownames(counts),
    Observed_T = obs_t_stats,
    Observed_P = obs_p_values,
    Empirical_P = emp_p_values,
    Significant_Perm = ifelse(emp_p_values < 0.05, "YES", "NO")
  )
  
  return(results)
}

# ==============================================================================
# 2. USAGE
# ==============================================================================
counts <- read.table(counts.file , header = T , row.names = 1 , sep = "\t", stringsAsFactors = F, check.names = F)
pheno <- read.csv(pheno.file , row.names = 1 , stringsAsFactors = F)

message("filtering low count genes...")
message("Genes that don't have minimum count of ",gFilter.min.count, " in at least ",(gFilter.min.prop*100) , "% of the samples will be removed.")
keep <- edgeR::filterByExpr(round(counts),group = pheno[,var.trait],min.count = gFilter.min.count, min.prop = gFilter.min.prop)
message(sum(!keep),"/",nrow(counts)," genes removed. Remaining genes:", sum(keep))
counts <- counts[keep,]

if(all(outliers != "")){
  message("Removing outliers...")
  if(! all (outliers %in% colnames(counts))){
    warning("The followin outlier IDs are not exist in the data: ")
    paste(outliers[!(outliers %in% colnames(counts))] , collapse = ";")
  }
  counts <- counts[,!(colnames(counts) %in% outliers)]
  pheno <- pheno[!(rownames(pheno) %in% outliers),]
  
  #paste("Is count and phenotype data are matched?",ifelse(identical(colnames(counts) , rownames(pheno)),"Yes","NO"))
}

results_perm <- run_permutation_test(
   counts = counts,
   pheno = pheno,
   var_of_interest = var_of_interest,
   covariates = covariates,
   n_perms = n_perms
)

write.csv(results_perm , file = paste0(OutPrefix , ".csv"), row.names = F)
# View top results
# head(results_tRNA[order(results_tRNA$Empirical_P), ])


