library(tidyverse)
library(DT)
library(grid)
library(gridExtra)
library(kableExtra)
library(here)
library(jsonlite)
library(plyr)
library(purrr)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ReactomePA)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(data.table)
library(tictoc)
library(fs)
source("code/utilities_drug.R")
source("code/patient_get_network_details.R")
source("code/01_analyze_patient.R")
source("code/format_target_data.R")
source("code/patient_plots.R")
source("code/patient_map_drugs.R")


## Load Global data files
## Load coohrt survival data
pheno <- read_csv("../data/TCGA_Survival_Gbm.csv") %>%
  arrange(desc(GuanScore))

## Load cohort expression data
cohort_expr <- read.csv("/Volumes/omics4tb2/SYGNAL/GBM-Serdar/MINER_MicroLowessRNATMM.08.24.2020/GbmMicroRNAMergedWithIDsZScored.csv")

## get cohort expression data genes
cohort_genes <- cohort_expr %>%
  dplyr::pull(X)

## Load gene identifier mappings
identifiers <- read.delim("../data/identifier_mappings.txt", sep = "\t") %>%
  filter(Source == "Gene Name")

## load additional evidence file
evidence <- read.csv("../data/Additional_Evidences_Gene.txt", sep = "\t",fileEncoding="latin1")

## load clinical trials info file
clinical_trials <- read.delim("../data/Glioblastoma_ClinicalTrials_01092022.txt", header = T, sep = "\t")

## Load patient specific files
patient_folders_path <- fs::dir_info("/Volumes/omics4tb2/SYGNAL/XCures/patients_network_activities/")
patient_folders <- fs::path_file(patient_folders_path$path)
# Filter samples that will not be anlyzed
patient_folders_current <- patient_folders[grep("TL-", patient_folders)]
patient_folders_current <- patient_folders_current[grep("TL-18-31DCF0", patient_folders_current, invert = T)]

## Loop through patient folders
for(patient_id in patient_folders_current[1]){
  # Patient analysis
  #tmp1 <- analyze_patient(sample.name = patient_sample)

  ## Patient network data file
  filename <- paste("output/",patient_id, "/", patient_id, "_drug_therapy_activity.csv",sep="")

  ## load network file for patient
  network_info <- read_csv(filename) %>%
    filter(`Drug Constrained Regulon Activity` > 0.8)
  ## Load patient's expression data for plotting and then z-score
  patient_exp_filename <- paste0("/Volumes/omics4tb2/SYGNAL/XCures/patients_processed_data/",patient_id, "/RNA/", patient_id, ".genes.results")
  patient_exp <- read.delim(patient_exp_filename, sep="\t") %>%
    separate(gene_id, into = c("Gene_ID", "Gene_Symbol"), sep = "_") %>%
    dplyr::select(Gene_ID, Gene_Symbol, TPM) %>%
    filter(Gene_ID %in% cohort_genes) %>%
    mutate(zscore = (TPM - mean(TPM))/sd(TPM))

  rmarkdown::render(input = "analysis/Patient_Analysis.Rmd",
                    output_format = "html_document",
                    output_file = paste0(gsub(".html","", patient_id), ".html"),
                    output_dir = "docs",
                    output_yaml = "analysis/_site.yml"
                    )

}


