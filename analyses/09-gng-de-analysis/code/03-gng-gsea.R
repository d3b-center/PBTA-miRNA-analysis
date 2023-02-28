# Gene Set Enrichment Analysis of GNG C1 vs C2 using ClusterProfiler
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia
# Script adapted from ALSF's DE Analysis Tutorial (https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html)

## LOAD LIBRARIES
# Attach the library
library(clusterProfiler)
# Package that contains MSigDB gene sets in tidy format
library(msigdbr)
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

# Declare input file paths
dge_results_file <- file.path(results_dir, "gng_diff_expr_results.csv")

# Read in the contents of the differential expression results file
dge_df <- readr::read_csv(dge_results_file)

#######
# Specifying MSigDB gene sets of interest
hs_hallmark_sets <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "H"
)

# Choose other genesets here:https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
c2_canonical_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP"
)

c3_mirna_targets_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C3",
  subcategory = "MIR:MIR_Legacy" # or "MIR:MIRDB"
)

c3_tf_targets_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C3",
  subcategory = "TFT:GTRD" # or "TFT:TFT_Legacy"
)
## OTHERS:
#> 15 C5     "GO:BP"                   7763
#> 16 C5     "GO:CC"                   1035
#> 17 C5     "GO:MF"                   1763
#> 18 C5     "HPO"                     5142
#> 19 C6     ""                         189
#> 20 C7     "IMMUNESIGDB"             4872

#######
## PERFORM GSEA

# 1. Determine our pre-ranked genes list
# Check if there are any duplicate genes present
any(duplicated(dge_df$Gene)) # None

# Create a named vector ranked based on the log2 fold change values
lfc_vector <- dge_df$log2FoldChange
names(lfc_vector) <- dge_df$Gene

# Sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# 2. Run GSEA using the GSEA() function
# Set the seed so our results are reproducible:
set.seed(2020)
# Run GSEA
gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_hallmark_sets,
    gs_name,
    gene_symbol
  )
)

# We can access the results from our `gsea_results` object using `@result`
head(gsea_results@result)

# Convert GSEA results object to dataframe
gsea_result_df <- data.frame(gsea_results@result)

# 3. Visualize GSEA
# Look at the 3 gene sets with the most positive NES
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(NES, n = 3)

# Make GSEA plot
most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_APOPTOSIS",
  title = "HALLMARK_APOPTOSIS",
  color.line = "#0d76ff"
)
most_positive_nes_plot

# Save GSEA enrichment plot as png
ggplot2::ggsave(file.path(plots_dir, "gng_mrna_gsea_most_positive_plot.png"),
                plot = most_positive_nes_plot
)

# Look at the 3 gene sets with the most negative NES
gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(NES, n = 3)

# Make GSEA plot
most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  color.line = "#0d76ff"
)
most_negative_nes_plot

# Save GSEA enrichment plot as png
ggplot2::ggsave(file.path(plots_dir, "gng_mrna_gsea_most_negative_plot.png"),
                plot = most_negative_nes_plot
)

# Save GSEA results
readr::write_csv(
  gsea_result_df,
  file.path(
    results_dir,
    "gng_mrna_gsea_results.csv"
  )
)


HALLMARK_MTORC1_SIGNALING
other_plot <- enrichplot::gseaplot2(
  gsea_results,
  geneSetID = "HALLMARK_MTORC1_SIGNALING",
  title = "HALLMARK_MTORC1_SIGNALING"
)
other_plot

require(DOSE)
dotplot(gsea_results, showCategory=10, split=".sign") + facet_grid(.~.sign)

ridgeplot(gsea_results) + labs(x = "enrichment distribution")

# PubMED terms
terms <- gsea_results$Description[1:3]
pmcplot(terms, 2010:2023, proportion=FALSE)


################################################################################
c5_bp_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C5",
  subcategory = "GO:BP"
)

c5_cc_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C5",
  subcategory = "GO:CC"
)

c5_mf_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C5",
  subcategory = "GO:MF"
)

# Run GSEA
gsea_c5mf_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    c5_mf_sets,
    gs_name,
    gene_symbol
  )
)

# Convert GSEA results object to dataframe
gsea_c5mf_result_df <- data.frame(gsea_c5mf_results@result)

require(DOSE)
dotplot(gsea_c5mf_results, showCategory=10, split=".sign") + facet_grid(.~.sign)




################################################################################

## KEGG GSEA
# Prepare input
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(lfc_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df <- as.data.frame(lfc_vector)
df$SYMBOL <- rownames(df)
rownames(df) <- NULL
df2 <- df[df$SYMBOL %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

kegg_gene_list <- df2$lfc_vector
# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y
# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)



kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "hsa",
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ridgeplot(kk2) + labs(x = "enrichment distribution")

# Visualize MAPK Pathway
library(pathview)

# Produce the native KEGG plot (PNG)
mapk <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04010", species = "hsa")

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04010", species = "hsa", kegg.native = F)




