# Kaplan Meier Plots
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia


# Install
install.packages("survminer")

# Load
library("survminer")

# Set working directory
#setwd("~/.../survival_plots")

# Load Data
dat <- read.csv("/Users/shehbeel/Documents/OpenPBTA-miRNA-Analysis/data/pbta_clinical_mirna.csv")
dat$OS_days <- as.numeric(dat$OS_days)
dat$PFS_days <- as.numeric(dat$PFS_days)

# Drop GNT, Schwannoma, and Teratoma samples
dat <- dat[dat$short_histology!='GNT',]
dat <- dat[dat$short_histology!='Schwannoma',]
dat <- dat[dat$short_histology!='Teratoma',]

# Overall Survival curves
require("survival")
fit <- survfit(Surv(OS_days, OS_status_boolean) ~ short_histology, data = dat)

# Drawing curves
ggsurvplot(fit, 
           title="Overall Survival (PBTA Cohort)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)

# Progression-free Survival curves
require("survival")
fit <- survfit(Surv(PFS_days) ~ short_histology, data = pbta.dat)

# Drawing curves
ggsurvplot(fit, 
           title="Progression-free Survival (PBTA Cohort)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)

################################################################################
## GENE EXPRESSION CORRELATED WITH SURVIVAL

# Make temporary dataframe containing gene expression column
temp_df <- dat[,c("OS_days", "PFS_days", "OS_status_boolean", "let.7a.5p")]

# Convert gene expression to "HIGH" or "LOW"
temp_df$let.7a.5p_factor <- factor(ifelse(temp_df$let.7a.5p > quantile(temp_df$let.7a.5p, 0.25)[[1]], "HIGH", "LOW"))

# Overall Survival curves
require("survival")
fit <- survfit(Surv(OS_days, OS_status_boolean) ~ let.7a.5p_factor, data = temp_df)

# Drawing curves
ggsurvplot(fit, 
           title="Overall Survival (let-7a-3p)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)

# Progression-free Survival curves
require("survival")
fit <- survfit(Surv(PFS_days) ~ let.7a.5p_factor, data = temp_df)

# Drawing curves
ggsurvplot(fit, 
           title="Progression-free Survival (let-7a-3p)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)



