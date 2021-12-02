library(tidyverse)
library(here)
library(DT)
library(data.table)
library(gridExtra)
library(jsonlite)
library(plyr)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ReactomePA)
library(AnnotationDbi)
library(org.Hs.eg.db)
source("code/utilities_drug.R")
source("code/drug_network_details.R")
source("code/regulon_enrichment.R")
###############################################################
programs <- read_json("../data/MINER_MicroLowessRNATMM.08.24.2020/transcriptional_programs.json", simplifyVector = T)
## Read regulon gene list
regulons <- read_json("../data/MINER_MicroLowessRNATMM.08.24.2020/regulons.json",simplifyVector = T)

programs_df <- ldply(programs, data.frame) %>%
  dplyr::rename("Program" = 1, "Regulon" = 2)

drug_info <- read_csv("../data/PreliminaryDrugTableForTherapyPrioritizationGBMMINER.csv")

drug_info_pr <- tibble()
for(i in 1: length(drug_info$MutationSymbol)){
  my_row = drug_info[i,]
  regulon = as.character(drug_info[i,"Regulon"])

  program = programs_df %>%
    filter(Regulon == regulon) %>%
    pull(Program)

  drug_info_pr <- bind_rows(drug_info_pr, bind_cols(my_row, Program = as.integer(program)))

}


network_dir_input = "/Volumes/omics4tb2/SYGNAL/GBM-Serdar/XCures/network_activities/P76156"
dna_dir = "/Volumes/omics4tb2/SYGNAL/XCures/P76156"

# Get list of all mutation files from XCures folder
mutation_files <- dir(dna_dir, pattern = "*_somatic.tsv",full.names = T,recursive = T)

## Mapping function
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
      reduce(rbind) %>%
      #mutate(AF = round(AF*100, digits = 2)) %>%
      separate(col=hgncGene,sep = ",", into=c("GENE"), remove = F) %>%
      filter(Consequence %in% c("missense","frameshift","stop_gained")) %>%
      pull(GENE) %>%
      unique()

      drug_reg_mut <- drug_reg %>%
      mutate(MutatedInPatient = if_else(MutationSymbol %in% patient_mutations, paste("TRUE"), paste("FALSE"))) %>%
      #mutate(TargetMutatedInPatient = str_match(DrugTargetAll, z)) %>%
      #group_by(Drug,TargetMutatedInPatient) %>%
      #summarise(tt = paste0(TargetMutatedInPatient, collapse=","))
      mutate(TargetMutatedInPatient = if_else(strsplit(DrugTargetAll,split = ",",fixed = T) %in% patient_mutations,paste("TRUE"), paste("FALSE")))

     # get network details
      zz <- drug_network_details(drug_reg_mut)

      # Join network details with main drug table
      drug_reg_final <- left_join(drug_reg_mut, zz, by="Drug")

    write_csv(drug_reg_final, paste0("output/", sample.name,"/",sample.name,"_",type,"_drug_therapy_activity.csv"))
    all.samples.activities <- bind_rows(all.samples.activities, drug_reg)

  }

  return(all.samples.activities)

} # function


drugs_all_regulon <- map_drugs(type = "regulon", disease = FALSE)
drugs_all_program <- map_drugs(type = "program", disease = FALSE)


drugs_all_regulon_pathways <- data.frame()

# patient_3
drugs_all_regulon_3 <- drugs_all_regulon %>%
  filter(patient == "P76156_3")

for(i in drugs_all_regulon_3$Drug){
  # select the drug row
  my.row = drugs_all_regulon_3 %>%
    filter(Drug == i)

  # get all the regulons
  my.regulons1 = drugs_all_regulon_3 %>%
    filter(Drug == i) %>%
    pull(DrugRegulonOverActive) %>%
    unique()
  # split regulons if multiple
  my.regulons2 <- strsplit(my.regulons1, split = ",")

  ## For each regulon get the regulon genes
  enrichment.list <- list()
  all.genes <- c()
  ii = 0
  total = length(my.regulons2[[1]])
  for(regulon in my.regulons2[[1]]){
    ii = ii + 1
    cat("\nProcessing regulon: ",regulon, " ", i, "/",total,"\n")
    regulon.genes <- regulons[[regulon]]
    cat("Total ", length(regulon.genes), "genes in regulon", regulon)

    # combine all regulon genes
    all.genes <- append(all.genes, regulon.genes)
    #

  }
  # remve duplicate genes
  all.genes = unique(all.genes)

  ## Convert gene symbols
  gene.entrez <-   AnnotationDbi::select(org.Hs.eg.db, keys= all.genes, column="ENTREZID", keytype="ENSEMBL", multiVals="first")$ENTREZID
  gene.symbols <-   AnnotationDbi::select(org.Hs.eg.db, keys= all.genes, column="SYMBOL", keytype="ENSEMBL", multiVals="first")$SYMBOL

  enrichment.list <- try(regulon_enrichment(my.genes = gene.symbols))
  #reactome.enrich <- enrichPathway(gene = gene.entrez, pvalueCutoff = 0.05, pAdjustMethod = "BH")

  drugs_all_regulon_pathways <- rbind(drugs_all_regulon_pathways,
                                      cbind(my.row,
                                            Reactome = paste0(enrichment.list$table8$Description, collapse=":"),
                                            GO = paste0(enrichment.list$table6$Description, collapse=":"),
                                            Hallmarks =  paste0(enrichment.list$table1$Description, collapse=":")
                                              ))

}










drugs_all_regulon %>%

  select()
 AnnotationDbi::select(org.Hs.eg.db, keys= genes.input, column="ENTREZID", keytype="SYMBOL", multiVals="first") %>%

my.entrez <- my.entrez$ENTREZID
#my.entrez <- as.character(sort(as.numeric(my.entrez), decreasing = T))
reactome.enrich <- enrichPathway(gene = my.entrez, pvalueCutoff = 0.05, pAdjustMethod = "BH")













#
#
#
#
#
# ## Reorder drugs based on their IC50 values
# drug.order <- all.samples.activities %>%
#   arrange(DrugConstrainedRegulonActivity) %>%
#   pull(Drug)
#
#
# ## Heatmap plot for IC50
# q <- ggplot(all.samples.activities, aes(y=sample, x=factor(Drug), fill=IC50))
# q <- q + geom_tile() +
#   scale_fill_gradientn(limits=c(1e-10,1e-05), colours=c("red", "white", "blue"))
# q <- q + theme(axis.text.x = element_text(angle=90))
# q
#
#
#
# scale_fill_distiller(palette = "YlGn", direction=2)
# "navyblue", "darkmagenta", "darkorange1"
#
# scale_fill_gradient2(low = "blue",mid = "red",high = "white", midpoint = 1e-5,)
# q
#
# scale_fill_gradient2(colours = c("white", "red", "blue"), values = c(1e-01,,0))
# q
#
# scale_fill_distiller(name = "Legend title", palette = "Reds", direction = -1, na.value = "transparent")
# q
#
# ## Heatmap plot for IC50
# s <- ggplot(all.samples.activities, aes(y=sample, x=factor(Drug), fill=DrugConstrainedRegulonActivity))
# s <- s + geom_tile() +
#   scale_fill_gradient2(low = "blue",mid = "white",high = "red", midpoint = 0)
# s <- s + theme(axis.text.x = element_text(angle=90))
#
# scale_fill_gradientn(colours = c("white", "red", "blue"), values = c(0,-1,1))
#
# scale_fill_distiller(name = "Legend title", palette = "Blues", direction = 1, na.value = "transparent")
# s
#
#
# #
# # z <- ggplot()
# # z <- z + geom_tile(data=all.samples.activities, aes(y=sample, x=factor(Drug), fill=DrugConstrainedRegulonActivity, alpha=0.3))
# # z <- z + geom_tile(data=all.samples.activities, aes(y=sample, x=factor(Drug), fill=log(IC50), alpha=0.3))
# # z
#
#
#
# #Arrange them in a grid
# gg1 <- ggplot_gtable(ggplot_build(q))
# gg2 <- ggplot_gtable(ggplot_build(s))
#
# grid.arrange(gg1, gg2, ncol=2)
#
#
#
#
#
# # Load HTP screening IC50 values
# hts.ic50 <- read_csv(file = "data/HTS_IC50_PDGSC_R01.csv")
# drugs <- toupper(unique(hts.ic50$`DRUG CANDIDATE`))
# write.csv(drugs, file="data/HTS_Drug_list.csv")
#
# xx <- drug_reg %>%
#   filter(Drug %in% drugs)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# activity.table <- tibble()
# for(sample in samples){
#   sample.name <- strsplit(sample, split = "/Regulon_Activity_")[[1]][2]
#   sample.name <- strsplit(sample.name, split=".csv")[[1]][1]
#
#   tmp1 <- read_csv(sample) %>%
#     mutate(Sample = sample.name) %>%
#     rename("Activity" = 2, "Regulon" = 1)
#
#   activity.table <- bind_rows(activity.table, tmp1)
# }
#
#
#
# reg_out <- read_csv(here("output/scca/regulon_activity_scca2_py.csv"))
# prog_out <- read_csv(here("output/scca/program_activity_scca2_py.csv"))
#
# subjects <- colnames(reg_out[,-1])
# for(sub in subjects) {
#
#   dir_path <- paste0(here("output/scca/scca_drug_therapy/"), sub)
#   if(!dir.exists(dir_path)) {
#     dir.create(dir_path)
#   }
#   reggie <- reg_out[,c("regulon", sample)]
#   proggie <- prog_out[,c("program", sub)]
#   write_csv(reggie, paste0(dir_path, "/", sub, "_regulon_activity.csv"))
#   write_csv(proggie, paste0(dir_path, "/", sub, "_program_activity.csv"))
#
#   ## generate drug regulon/program activity
#   drug_reg <- getDrugTherapyActivity(drug_info, regulon=reggie, program=proggie)
#
#   write_csv(drug_reg, paste0(dir_path, "/", sub, "_drug_therapy_activity.csv"))
# }
# #####################################################################
