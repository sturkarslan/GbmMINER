library(tidyverse)
library(here)
library(DT)
library(gridExtra)
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
################ Mapping function   ################
map_drugs <- function(drug_info_pr, sample.id = sample.name, type="regulon", disease=TRUE, network_dir=network_dir_input,mutation_files=mutation_files){

  all.samples.activities <- data.frame()

  # Load Regulon and Program activities for all or Disease relevant
  if (disease == TRUE) {
      # Regulon activity
    reg_out <- read_csv(paste0(network_dir,"/",sample.id,"_disease_rel_regulon_activity.csv")) %>%
      dplyr::rename("Regulon" = 1, "activity" = 2)

    # Program activity
    prog_out <- read_csv(paste0(network_dir, "/",sample.id, "_disease_rel_program_activity.csv")) %>%
      dplyr::rename("Program" = 1, "activity" = 2)

    cat("Loaded: Disease relevant: ", dim(reg_out)[1], " regulon activity and ", dim(prog_out)[1], " program activity.\n")

  } else{
    # Regulon activity
    reg_out <- read_csv(paste0(network_dir,"/",sample.id,"_all_regulon_activity.csv")) %>%
      dplyr::rename("Regulon" = 1, "activity" = 2)

    # Program activity
    prog_out <- read_csv(paste0(network_dir, "/",sample.id, "_all_program_activity.csv")) %>%
      dplyr::rename("Program" = 1, "activity" = 2)

    cat("Loaded: All: ", dim(reg_out)[1], " regulon activity and ", dim(prog_out)[1], " program activity.\n")

  }


  # create directory to collect results for each sample
  dir_path <- paste0("output/", sample.id)
  if(!dir.exists(dir_path)) {
    dir.create(dir_path)
  }


  ## generate drug regulon/program activity
  drug_reg <- getDrugTherapyActivity(drug_info_pr, regulon=reg_out, program=prog_out)
  drug_reg$patient = sample.id

  ## Mutations
  # Get patients mutation data
  patient_mutation_files <- grep(sample.id, mutation_files, value = T)

  # Get mutated gene names
  patient_mutations <- patient_mutation_files %>%
    map(read_delim, delim="\t") %>%
    purrr::reduce(rbind) %>%
    #mutate(AF = round(AF*100, digits = 2)) %>%
    separate(col=hgncGene,sep = ",", into=c("GENE"), remove = F) %>%
    filter(Consequence %in% c("missense","frameshift","stop_gained")) %>%
    pull(GENE) %>%
    base::unique()

  drug_reg_mut <- drug_reg %>%
    mutate(MutatedInPatient = if_else(MutationGene %in% patient_mutations, paste("TRUE"), paste("FALSE"))) %>%
    mutate(TargetMutatedInPatient = if_else(strsplit(DrugTargetAll,split = ",",fixed = T) %in% patient_mutations,paste("TRUE"), paste("FALSE")))

  #write_csv(drug_reg_mut, paste0("output/", sample.id,"/",sample.id,"_drug_therapy_activity.csv"))
  all.samples.activities <- bind_rows(all.samples.activities, drug_reg_mut)
  return(all.samples.activities)

}
