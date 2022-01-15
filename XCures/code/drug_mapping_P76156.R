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
library(tictoc)
source("code/utilities_drug.R")
source("code/patient_get_network_details.R")
source("code/regulon_enrichment.R")
source("code/patient_map_drugs.R")
source("code/patient_drug_pathway_enrichment.R")
source("code/patient_rankorder_drug_table.R")
###############################################################
sample.name = "P76156_3"
type = "regulon"

analyze_patient <- function(sample.name)
{

  # define input directories
  network_dir_input = paste0("/Volumes/omics4tb2/SYGNAL/GBM-Serdar/XCures/network_activities/", sample.name)

  # define mutations directory
  dna_dir = paste0("/Volumes/omics4tb2/SYGNAL/XCures/","P76156")

  # get program info
  programs <- read_json("../data/MINER_MicroLowessRNATMM.08.24.2020/transcriptional_programs.json",simplifyVector = T)
  cat("Reading `Program information` is done\n")

  ## Read regulon gene list
  regulons <- read_json("../data/MINER_MicroLowessRNATMM.08.24.2020/regulons.json",simplifyVector = T)
  cat("Reading `Regulon information` is done\n")

  # reformat regulons
  programs_df <- ldply(programs, data.frame) %>%
    dplyr::rename("Program" = 1, "Regulon" = 2)

  # get disease relevant drug information
  #drug_info <- read_csv("../data/PreliminaryDrugTableForTherapyPrioritizationGBMMINER.csv")

  #drug_info <- read.csv("../data/GBM_Master_Drugs_Mapped_CMFlows_120921.csv")
  drug_info_pr <- read.csv("../data/GBM_Master_Drugs_Mapped_CMFlows_120921_with_regulons.csv") %>%
    mutate(Regulon_ID = as.character(Regulon_ID))

  # # modify drug information
  # drug_info_pr <- tibble()
  # for(i in 1: length(drug_info$Drug)){
  #   my_row = drug_info[i,]
  #   regulon = as.character(drug_info[i,"Regulon_ID"])
  #
  #   program = programs_df %>%
  #     filter(Regulon == regulon) %>%
  #     pull(Program)
  #
  #   drug_info_pr <- bind_rows(drug_info_pr, bind_cols(my_row, Program = as.integer(program)))
  # }
  # cat("Reading and formnatting `Drug information` is done\n")
  #write.csv(drug_info_pr, "../data/GBM_Master_Drugs_Mapped_CMFlows_120921_with_regulons.csv")

  # Get list of all mutation files from XCures folder
  mutation_files <- dir(dna_dir, pattern = "*_somatic.tsv",full.names = T,recursive = T)


  ################### Mapping ##################################
  # drug mapping for regulons
  tic("Mapping drugs...")
  drugs_all <- map_drugs(drug_info_pr,sample.id = sample.name, disease = FALSE, network_dir = network_dir_input, mutation_files = mutation_files)
  toc()

  ################## Get Network Details ##################################
  # get summarized network details
  tic("Getting Network details...")
  drugs_all_details <- drug_network_details(data_mat = drugs_all)
  toc()

  # Filter for non-zero regulon and program activity
  drugs_all_details_filt <- drugs_all_details %>%
    filter(`Drug Constrained Regulon Activity` > 0.8)

  # dd <- drugs_all_details_filt %>%
  #   #group_by(`Target Gene(s)`) %>%
  #   mutate(DD = ddply(drugs_all_details_filt, .(`Target Gene(s)`), summarize, Drug=paste(Drug, collapse=",")))
  # ################### Get Functional Enrichment ##################################
  # get pathway information
  tic("Performing pathway enrichment...")
  drugs_all_pathways <- drug_pathway_enrichment(input.df = drugs_all_details_filt)
  toc()

  drugs_all_pathways_ordered <- drugs_all_pathways %>%
    group_by(`Target Gene(s)`) %>%
    arrange(desc(max_glioblastoma.multiforme_phase,`Drug Constrained Regulon Activity`,`Other FDA Appr.`))

  ################### Write results into output file##################################

  write_csv(drugs_all_pathways_ordered, paste0("output/", sample.name,"/",sample.name,"_drug_therapy_activity.csv"))
}


# Analyze patient P76156_3
P76156_3 <- analyze_patient(sample.name = "P76156_3")
P76156_6 <- analyze_patient(sample.name = "P76156_6")






d1 <- read_csv("../data/opentargets_gbm.csv")

d2 <- d1 %>%
  separate_rows(c("approved_name","target_id"), sep = ":")
write_csv(d2, "../data/opentargets_gbm_longer.csv")






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
