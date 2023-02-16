# Survival Analysis
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia

# Load libraries
library(survival)
#library(ggpubr)
library(tidyverse)

# Set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "03-survival-analysis")
# Source this script, which contains a wrapper function that can conduct the survival analyses, from OpenPBTA
source(file.path(analysis_dir, "util", "survival_models.R"))
# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Set output directories
input_dir <- file.path(root_dir, "data")
#tma_dir <- file.path(root_dir, "analyses", "05-oncoplot", "input")
results_dir <- file.path(analysis_dir, "output")
plots_dir <- file.path(analysis_dir, "plots")


# Make output directories if they don't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Declare input file paths
metadata_file <- file.path(input_dir, "pbta_clinical_data.csv")
cluster_file <- file.path(input_dir, "pbta_clusters_mirna.csv")

# Declare output file paths
kap_meier_plot_file <- file.path(plots_dir, "KM_mirna_clusters.pdf")
kap_meier_model_file <- file.path(results_dir, "logrank_mirna_clusters.RDS")
table_s2_file <- file.path(results_dir, "Table-S2.xlsx")

# Import metadata and cluster data
metadata <- read.csv(metadata_file)
clusterdata <- read.csv(cluster_file)

# Preprocess the data
metadata_updated <- readr::read_csv(metadata_file) %>%
  # Drop duplicate rows
  distinct(Sample_ID, .keep_all = TRUE) %>%
  # Add Clusters
  inner_join(clusterdata, by="Sample_ID") %>%
  # Convert clusters to factors
  mutate_at(c("Cluster"), as.factor) %>%
  # Keep only necessary columns
  select("Sample_ID", "short_histology", "OS_status", "OS_days", "PFS_days", "Cluster")

## KAPLAN-MEIER FOR CLUSTERS
kap_fit <- survival_analysis(
  metadata  = metadata_updated,
  ind_var = "Cluster",
  test = "kap.meier",
  metadata_sample_col = "Sample_ID"
)

# median and 0.95CI survival days
kap_fit_data <- kap_fit$original_data %>%
  mutate(OS_months = OS_days/30.417,
         OS_status = ifelse(OS_status == 1, "LIVING", "DECEASED")) %>%
  select(-OS_days)

kap_fit_data_summary <- kap_fit_data %>%
  group_by(Cluster) %>%
  summarise_each(funs(median(., na.rm = T), sd(., na.rm = T)), OS_months)



# Make KM survival plot 
surv_plot <- survminer::ggsurvplot(fit = kap_fit$model,
                                   data = kap_fit$original_data,
                                   pval = TRUE,
                                   risk.table = TRUE,
                                   xlim = c(0, 4000),
                                   break.time.by = 500,
                                   ggtheme = theme_minimal(),
                                   risk.table.y.text.col = TRUE,
                                   risk.table.y.text = FALSE,
)

surv_plot$plot <- surv_plot$plot +
  ggtitle(paste0("Kaplan Meier miRNA Clusters")) +
  theme(legend.position = "right")
# Make this plot a combined plot
surv_plot_all_subtype <-
  cowplot::plot_grid(surv_plot[[1]], surv_plot[[2]], nrow = 2,
                     rel_heights = c(3, 2))
# Print it out here
surv_plot_all_subtype

## Save the plot to a file
# We can save the plot like a normal ggplot
cowplot::save_plot(filename = kap_meier_plot_file, plot = surv_plot_all_subtype, base_height = 5, base_width = 8)
# Save the model itself as well as all other output from `survival_analysis()`
readr::write_rds(kap_fit, kap_meier_model_file)


## Analysis: cox regression based on cluster membership
# Run model
fit_save_model(metadata_updated, 
               "Cluster",
               file.path(results_dir, "cox_cluster_membership.RDS"))

# ERRORRRRR
# The fit_save_model() function in survival_models.R is not universal





