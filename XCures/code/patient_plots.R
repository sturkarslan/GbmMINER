library(tidyverse)
library(here)
library(DT)
library(gridExtra)
library(jsonlite)
library(plyr)
library(purrr)
library(ggplot2)
library(reshape2)
library(data.table)

#Load data files
# get program info
programs <- read_json("../data/MINER_MicroLowessRNATMM.08.24.2020/transcriptional_programs.json", simplifyVector = T)

## Read regulon gene list
regulons <- read_json("../data/MINER_MicroLowessRNATMM.08.24.2020/regulons.json",simplifyVector = T)

# reformat regulons
programs_df <- ldply(programs, data.frame) %>%
  dplyr::rename("Program" = 1, "Regulon" = 2)

# get drug information
drug_info <- read_csv("../data/PreliminaryDrugTableForTherapyPrioritizationGBMMINER.csv")

# modify drug information
drug_info_pr <- tibble()
for(i in 1: length(drug_info$MutationSymbol)){
  my_row = drug_info[i,]
  regulon = as.character(drug_info[i,"Regulon"])

  program = programs_df %>%
    filter(Regulon == regulon) %>%
    pull(Program)

  drug_info_pr <- bind_rows(drug_info_pr, bind_cols(my_row, Program = as.integer(program)))
}


## Load cohort regulon activity data
all_patient_reg_activity <- read_csv("data/drug_constrained_network_activity_allcohort.csv")

## Load coort survival data
pheno <- read_csv("../data/TCGA_Survival_Gbm.csv") %>%
  arrange(desc(GuanScore))

# # load cohort gene expression data
# cohort_expression <- read_csv("../data/MINER_MicroLowessRNATMM.08.24.2020-ST/data/GbmMicroRNAMergedWithIDsZScored.csv")
#
# # load cohort gene expression data
# patient_expression1 <- read_delim("/Volumes/omics4tb2/SYGNAL/XCures/P76156/P76156_3/RNA/results_RSEM/P76156_3.genes.results", delim = "\t")
# patient_expression1 <- patient_expression1 %>%
#   separate(gene_id, into = c("Ensembl_id","Gene_Symbol"), sep = "_") %>%
#   filter(Ensembl_id != "ENSG00000202354") %>%
#   filter(Ensembl_id != "ENSG00000202354") %>%
#   filter(Ensembl_id != "ENSG00000201098") %>%
#   filter(Ensembl_id != "ENSG00000206585") %>%
#   filter(!grepl("RNA",Gene_Symbol)) %>%
#   filter(!grepl("RNY",Gene_Symbol)) %>%
#   filter(!grepl("RNU",Gene_Symbol)) %>%
#   mutate(zscore = (TPM - mean(TPM))/sd(TPM)) %>%
#   select(Ensembl_id, Gene_Symbol,zscore)
#
#
# # load cohort gene expression data
# patient_expression2<- read_delim("/Volumes/omics4tb2/SYGNAL/XCures/P76156/P76156_6/RNA/results_RSEM/P76156_6.genes.results", delim = "\t")
# patient_expression2 <- patient_expression2 %>%
#   separate(gene_id, into = c("Ensembl_id","Gene_Symbol"), sep = "_") %>%
#   filter(Ensembl_id != "ENSG00000202354") %>%
#   filter(Ensembl_id != "ENSG00000202354") %>%
#   filter(Ensembl_id != "ENSG00000201098") %>%
#   filter(Ensembl_id != "ENSG00000206585") %>%
#   filter(!grepl("RNA",Gene_Symbol)) %>%
#   filter(!grepl("RNY",Gene_Symbol)) %>%
#   filter(!grepl("RNU",Gene_Symbol)) %>%
#   mutate(zscore = (TPM - mean(TPM))/sd(TPM)) %>%
#   select(Ensembl_id, Gene_Symbol,zscore)
#
#



#### Function for histogram
plot_activity_histogram <- function(drug){
  ## extract target constrained regulon activity for patient 1st timepoint
  patient_score1 <- read_csv(paste0("output/P76156_3/P76156_3_regulon_drug_therapy_activity.csv"))

  ## extract target constrained regulon activity for patient 2nd timepoint
  patient_score2 <- read_csv(paste0("output/P76156_6/P76156_6_regulon_drug_therapy_activity.csv"))

  patient_score1 %>%
    dplyr::filter(Drug==drug) %>%
    dplyr::select(`Drug Constrained Regulon Activity`) %>%
    unlist() ->
    p_reg_score1

  patient_score2 %>%
    dplyr::filter(Drug==drug) %>%
    dplyr::select(`Drug Constrained Regulon Activity`) %>%
    unlist() ->
    p_reg_score2

  p_reg_score <- tibble(Time = c("2020","2021"), Activity=c(p_reg_score1,p_reg_score2))

  activity <-  all_patient_reg_activity
  drug_act <- activity[activity$Drug==drug, pheno$Patient_ID]
  drug_act2 <- melt(drug_act)

  ppp <- ggplot() +
    geom_histogram(data=drug_act2, aes(x=value)) +
    ggtitle(drug) +
    geom_vline(data=p_reg_score, aes(xintercept=Activity, color=Time),linetype="dotted", size=1) +
    geom_point(data = p_reg_score, aes(x=Activity, y=0, color=Time), size=4) +
    labs(x="Drug Constrained Regulon Activity", y="Count")

  ppp

}


plot_expression <- function(drug){


}










