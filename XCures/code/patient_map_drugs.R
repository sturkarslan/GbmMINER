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
map_drugs <- function(type="regulon", disease=TRUE, network_dir=network_dir_input){
  all.samples.activities <- data.frame()

  # Regulon activities
  if(type == "regulon"){
    if(disease == TRUE){
      samples <- dir(network_dir, pattern = "*_disease_rel_regulon_activity.csv",full.names = T)
    } else{
      samples <- dir(network_dir, pattern = "*_all_regulon_activity.csv",full.names = T)
    }
  }

  # Program activities
  if(type == "program"){
    if(disease == TRUE){
      samples <- dir(network_dir, pattern = "*_disease_rel_program_activity.csv",full.names = T)
    } else{
      samples <- dir(network_dir, pattern = "*_all_program_activity.csv",full.names = T)
    }
  }

  ### Loop through each sample
  for(sample in samples){
    cat("processing ", sample, "\n")
    # Get sample name
    #sample.name1 <- strsplit(sample, split = "_disease_rel_regulon_activity.csv")[[1]][1]
    if(disease == FALSE){
      sample.name2 <- strsplit(sample, split="/")[[1]][9]
      sample.name <- strsplit(sample.name2, split="_all")[[1]][1]
    } else {
      sample.name2 <- strsplit(sample, split="/")[[1]][9]
      sample.name <- strsplit(sample.name2, split="_disease")[[1]][1]
    }

    # create directory to collect results for each sample
    dir_path <- paste0("output/", sample.name)
    if(!dir.exists(dir_path)) {
      dir.create(dir_path)
    }

    # Regulon and program activity files
    reg_out <- read_csv(sample) %>%
      dplyr::rename("Regulon" = 1, "activity" = 2)

    if(disease == TRUE){
      prog_out <- read_csv(paste0(network_dir, "/",sample.name, "_disease_rel_program_activity.csv")) %>%
        dplyr::rename("Program" = 1, "activity" = 2)
    } else {
      prog_out <- read_csv(paste0(network_dir, "/",sample.name, "_all_program_activity.csv")) %>%
        dplyr::rename("Program" = 1, "activity" = 2)
    }

    ## generate drug regulon/program activity
    drug_reg <- getDrugTherapyActivity(drug_info_pr, regulon=reg_out, program=prog_out)
    drug_reg$patient = sample.name

    ## Mutations
    # Get patients mutation data
    patient_mutation_files <- grep(sample.name, mutation_files, value = T)

    # Get mutated gene names
    patient_mutations <- patient_mutation_files %>%
      map(read_delim, delim="\t") %>%
      purrr::reduce(rbind) %>%
      #mutate(AF = round(AF*100, digits = 2)) %>%
      separate(col=hgncGene,sep = ",", into=c("GENE"), remove = F) %>%
      filter(Consequence %in% c("missense","frameshift","stop_gained")) %>%
      pull(GENE) %>%
      base::unique()
    #print(patient_mutations)

    drug_reg_mut <- drug_reg %>%
      mutate(MutatedInPatient = if_else(MutationSymbol %in% patient_mutations, paste("TRUE"), paste("FALSE"))) %>%
      #mutate(TargetMutatedInPatient = str_match(DrugTargetAll, z)) %>%
      #group_by(Drug,TargetMutatedInPatient) %>%
      #summarise(tt = paste0(TargetMutatedInPatient, collapse=","))
      mutate(TargetMutatedInPatient = if_else(strsplit(DrugTargetAll,split = ",",fixed = T) %in% patient_mutations,paste("TRUE"), paste("FALSE")))

    # get network details
    # zz <- drug_network_details(drug_reg_mut)

    # Join network details with main drug table
    # drug_reg_final <- left_join(drug_reg_mut, zz, by="Drug")

    #write_csv(drug_reg_final, paste0("output/", sample.name,"/",sample.name,"_",type,"_drug_therapy_activity.csv"))
    all.samples.activities <- bind_rows(all.samples.activities, drug_reg_mut)

  }

  return(all.samples.activities)

}
