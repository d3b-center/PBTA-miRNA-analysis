# Gene set variation analysis of GNG C1 vs C2 using GSVA
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia
# Script adapted from ALSF's DE Analysis Tutorial (https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_03_gsva.html)

## LOAD LIBRARIES
# Attach the DESeq2 library
library(DESeq2)
# Attach the `qusage` library
library(qusage)
# Attach the `GSVA` library
library(GSVA)
# Human annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)
# We will need this so we can use the pipe: %>%
library(magrittr)



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
expression_df <- expression_df %>%
  # Only keep rows that have total counts above the cutoff
  dplyr::filter(rowSums(.) >= 50) %>%
  # The next DESeq2 functions need the values to be converted to integers
  round()

# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = expression_df, # Our prepped data frame with counts
  colData = metadata, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)

# Perform DESeq2 normalization and transformation
# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)

# Retrieve the normalized data from the `DESeqDataSet`
vst_df <- assay(dds_norm) %>%
  as.data.frame() %>% # Make into a data frame
  tibble::rownames_to_column("gene_symbol") # Make Gene IDs into their own column

# Import Gene Sets
hallmark_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens", # Can change this to what species you need
  category = "H" # Only hallmark gene sets
)
# Make gene sets into a list, which is required for GSVA
hallmarks_list <- split(
  hallmark_gene_sets$entrez_gene, # The genes we want split into pathways
  hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
)
# Look at the gene set list
head(hallmarks_list, n = 2)

# GENE IDENTIFIER CONVERSION
keytypes(org.Hs.eg.db)

# First let's create a mapped data frame we can join to the gene expression values
mapped_df <- data.frame(
  "entrez_id" = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = vst_df$gene_symbol,
    # Replace with the type of gene identifiers in your data
    keytype = "SYMBOL",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
  # drop that from the data frame
  dplyr::filter(!is.na(entrez_id)) %>%
  # Make an `Ensembl` column to store the row names
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(vst_df, by = c("Ensembl" = "gene_symbol"))

# Count up how many Entrez IDs mapped to multiple Ensembl IDs
sum(duplicated(mapped_df$entrez_id))

# Handling duplicate gene identifiers
# First let's determine the gene means
gene_means <- rowMeans(mapped_df %>% dplyr::select(-Ensembl, -entrez_id))

# Let's add this as a column in our `mapped_df`.
mapped_df <- mapped_df %>%
  # Add gene_means as a column called gene_means
  dplyr::mutate(gene_means) %>%
  # Reorder the columns so `gene_means` column is upfront
  dplyr::select(Ensembl, entrez_id, gene_means, dplyr::everything())

filtered_mapped_df <- mapped_df %>%
  # Sort so that the highest mean expression values are at the top
  dplyr::arrange(dplyr::desc(gene_means)) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(entrez_id, .keep_all = TRUE)

# Check for duplicates again
sum(duplicated(filtered_mapped_df$entrez_id))


## Prepare data matrix for GSVA
filtered_mapped_matrix <- filtered_mapped_df %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-Ensembl, -gene_means) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("entrez_id") %>%
  # Now we can convert our object into a matrix
  as.matrix()

## PERFORM GENE SET VARIATION ANALYSIS
# Perform GSVA
gsva_results <- gsva(
  filtered_mapped_matrix,
  hallmarks_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 15,
  # Maximum gene set size
  max.sz = 500,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)
# Print first 6 rows of GSVA results
head(gsva_results[, 1:10])

# Write GSVA results as csv
gsva_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  readr::write_csv(file.path(
    results_dir,
    "gng_gsva_results.csv"
  ))

## Visualizing results with a heatmap
# The pheatmap::pheatmap() will want the annotation data frame to have matching row names to the 
# data we supply it (which is our gsva_results).
annot_df <- metadata %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("Sample_ID") %>%
  arrange(cluster)

# Check if this is in the same order
all.equal(colnames(gsva_results), rownames(annot_df))

# Reorder GSVA results as well
gsva_results <- gsva_results[,rownames(annot_df)]

pathway_heatmap <- pheatmap::pheatmap(gsva_results,
                                      annotation_col = annot_df, # Add metadata labels!
                                      show_colnames = FALSE, # Don't show sample labels
                                      cluster_cols = FALSE,
                                      #cluster_rows=FALSE,
                                      fontsize_row = 10 # Shrink the pathway labels a tad
)

# Print out heatmap here
pathway_heatmap


# Replace file name with a relevant output plot name
heatmap_png_file <- file.path(plots_dir, "gng_gsva_heatmap.png")

# Open a PNG file - width and height arguments control the size of the output
png(heatmap_png_file, width = 700, height = 800)

# Print your heatmap
pathway_heatmap

# Close the PNG file:
dev.off()












