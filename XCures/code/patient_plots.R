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
drug_info <- read.csv("../data/GBM_Master_Drugs_Mapped_CMFlows_120921.csv")

# modify drug information
# drug_info_pr <- tibble()
# for(i in 1: length(drug_info$MutationSymbol)){
#   my_row = drug_info[i,]
#   regulon = as.character(drug_info[i,"Regulon"])
#
#   program = programs_df %>%
#     filter(Regulon == regulon) %>%
#     pull(Program)
#
#   drug_info_pr <- bind_rows(drug_info_pr, bind_cols(my_row, Program = as.integer(program)))
# }


## Load cohort regulon activity data
all_patient_reg_activity <- read_csv("data/drug_constrained_network_activity_allcohort.csv")

## Load coort survival data
pheno <- read_csv("../data/TCGA_Survival_Gbm.csv") %>%
  arrange(desc(GuanScore))


cohort_expr <- read.csv("/Volumes/omics4tb2/SYGNAL/GBM-Serdar/MINER_MicroLowessRNATMM.08.24.2020/GbmMicroRNAMergedWithIDsZScored.csv")

cohort_genes <- cohort_expr %>%
  dplyr::pull(X)

identifiers <- read.delim("../data/identifier_mappings.txt", sep = "\t") %>%
  filter(Source == "Gene Name")


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
  patient_score1 <- read_csv(paste0("output/P76156_3/P76156_3_drug_therapy_activity.csv"))

  ## extract target constrained regulon activity for patient 2nd timepoint
  patient_score2 <- read_csv(paste0("output/P76156_6/P76156_6_drug_therapy_activity.csv"))

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
    #ggtitle(drug) +
    geom_vline(data=p_reg_score, aes(xintercept=Activity, color=Time),linetype="dotted", size=1) +
    geom_point(data = p_reg_score, aes(x=Activity, y=0, color=Time), size=4) +
    labs(x="Drug Constrained Regulon Activity", y="Count")

  ppp

}


plot_expression <- function(target){
  targets <- strsplit(target, split = ",")[[1]]
  targets_ids <- identifiers %>%
    filter(Name %in% targets) %>%
    pull(Preferred_Name)

exp_data_c <- cohort_expr %>%
  filter(X %in% targets_ids) %>%
  unique()

exp_data_c <- reshape2::melt(exp_data_c)

exp_data_d <- left_join(exp_data_c, identifiers, by=c("X" = "Preferred_Name"))

patient_zscore <- patient_exp %>%
  filter(Gene_ID %in% targets_ids) %>%
  select(Gene_ID,zscore)
patient_zscore_m <- reshape2::melt(patient_zscore)
patient_zscore_n <- left_join(patient_zscore_m, identifiers, by=c("Gene_ID" = "Preferred_Name"))




p <- ggplot(exp_data_d, aes(x=Name, y=value, group=Name, fill=Name))
p <- p + geom_violin(alpha=0.5)
p <- p + geom_hline(data=patient_zscore_n, aes(yintercept=value,color=Name),linetype="dotted", size=1)
p <- p + guides(color=guide_legend(title="Patient Expression"))
p <- p + guides(fill=guide_legend(title="Gene(s)"))
p <- p + labs(x="Gene(s)", y="Gene Expression z-score")
#p <- p + geom_hline(yintercept=patient_zscore,linetype="dotted", size=1, color=patient_zscore)
p

}










