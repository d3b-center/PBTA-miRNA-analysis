# Large-scale miRNA sequencing analysis reveals four distinct clustering of pediatric brain tumor subtypes

Authors: ...

## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct miRNA analyses for the manuscript noted above.

1. Clone the repository
```
git clone https://github.com/d3b-center/OpenPBTA-miRNA-analysis.git
```

2. Pull the docker container:
```
docker pull shehbeel/pbta-mirna:latest
```

3. Start the docker container, from the `PBTA-miRNA-analysis` folder, run:
```
docker run --name container_name -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/PBTA-miRNA-analysis pbta-mirna:latest
```

4. To execute shell within the docker image, from the `PBTA-miRNA-analysis` folder, run:
```
docker exec -ti container_name bash
```

5. Navigate to an analysis module and run the shell script:
```
cd /home/rstudio/PBTA-miRNA-analysis/analyses/module_of_interest
```


### Below is the main directory structure listing the analyses and data files used in this repository

```
.
├── Dockerfile
├── LICENSE
├── README.md
├── analyses
│   ├── 01-circos-plot
│   ├── 02-consensus-clustering
│   ├── 03-survival-analysis
│   ├── 04-clinical-correlation
│   ├── 05-mutational-burden-correlation
│   ├── 06-tumor-specific-consensus-clustering
│   ├── 07-sankey-plot
│   ├── 08-gng-de-analysis
│   ├── 09-additional_figures
│   └── 10-tables
├── data
│   ├── pbta_clinical_mirna.csv
│   ├── pbta_meta_mirna.rds
│   ├── pbta_meta_mrna.rds
│   ├── pbta_raw_counts_mirna.rds
│   ├── pbta_raw_counts_mrna.rds
│   ├── pbta_clusters_mirna.csv
│   ├── pbta_clusters_gng.csv
│   └── v1
├── utils
├── palletes
├── data
├── download-data.sh
└── scratch
```

## Code Authors

Shehbeel Arif ([@shehbeel](https://github.com/shehbeel))

## Contact

For questions, please submit an issue or send an email to Shehbeel Arif: arifs2@chop.edu
