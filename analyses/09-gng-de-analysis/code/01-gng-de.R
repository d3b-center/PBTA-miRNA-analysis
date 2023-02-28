# Differential Expression Analysis of GNG C1 vs C2 using DESeq2
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia
# Script adapted from ALSF's DE Analysis Tutorial (https://alexslemonade.github.io/refinebio-examples/03-rnaseq/differential-expression_rnaseq_01.html)

## LOAD LIBRARIES
# Library for DE Analysis
library(DESeq2)
# Library data manipulation
#library(dplyr)
# Library for Plotting
library(ggplot2)
# We will need this so we can use the pipe: %>%
library(magrittr)

# Set seed because jitter plot function involves randomness
set.seed(12345)


## SET DIRECTORIES
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "09-gng-de-analysis")

# Set output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

# Make output directories if they don't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Set input directory
data_dir <- file.path(root_dir, "data")

# Declare input file paths
data_file <- file.path(data_dir, "pbta_raw_counts_mrna.rds")
metadata_file <- file.path(data_dir, "pbta_clusters_gng.csv")

#######
# Import metadata and data 
metadata <- readr::read_csv(metadata_file)
expression_df <- readr::read_rds(data_file)

#######
## PREPROCESS THE DATA
# Select only the GNG Sample_IDs
expression_df <- expression_df %>% 
  dplyr::select(metadata$Sample_ID)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$Sample_ID)

# Make cluster a factor and set the levels appropriately
metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    cluster = factor(cluster, levels = c("1", "2"))
  )

# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

## Create DESeq2Dataset
# round all expression counts
gene_matrix <- round(filtered_expression_df)

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~cluster
)

## Run Differential Expression Analysis
deseq_object <- DESeq(ddset)
# Extract results table
deseq_results <- results(deseq_object)

# Use lfcShrink() function to obtain shrunken log fold change estimates based on 
# negative binomial distribution. This will add the estimates to your results table. 
# Using lfcShrink() can help decrease noise and preserve large differences between 
# groups (it requires that apeglm package be installed) (Zhu et al., Bioinformatics 2018).
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)

# Sort and filter DESeq2 results table and convert to dataframe
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

# View results sorted by adjusted p-value
# deseq_df %>%
#   dplyr::arrange(padj)

# Check results by plotting one gene
plotCounts(ddset, gene = "F8A2", intgroup = "cluster")

# Save results as CSV
readr::write_csv(
  deseq_df,
  file.path(
    results_dir,
    "gng_diff_expr_results.csv" # Replace with a relevant output file name
  )
)


## Create volcano plot
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
volcano_plot
# Save volcano plot
ggsave(
  plot = volcano_plot,
  file.path(plots_dir, "gng_de_mrna_volcano_plot.png")
) # Replace with a plot name relevant to your data









