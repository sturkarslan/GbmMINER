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
source("code/utilities_drug.R")
source("code/patient_get_network_details.R")
source("code/regulon_enrichment.R")
source("code/patient_map_drugs.R")
source("code/patient_drug_pathway_enrichment.R")
source("code/patient_rankorder_drug_table.R")
###############################################################
# define input directories
network_dir_input = "data"

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

################### Mapping ##################################
# drug mapping for regulons
drugs_all_regulon <- map_drugs(type = "regulon", disease = FALSE)

# drug mapping for programs
drugs_all_program <- map_drugs(type = "program", disease = FALSE)


################### Patient Filtering ##################################
# patient_3
drugs_all_regulon_3 <- drugs_all_regulon %>%
  filter(patient == "P76156_3")

drugs_all_program_3 <- drugs_all_program %>%
  filter(patient == "P76156_3")

# patient_6
drugs_all_regulon_6 <- drugs_all_regulon %>%
  filter(patient == "P76156_6")

drugs_all_program_6 <- drugs_all_program %>%
  filter(patient == "P76156_6")



################### Get Network Details ##################################
# get summarized network details
#source("code/drug_network_details.R")
drugs_all_regulon_3_details <- drug_network_details(data_mat = drugs_all_regulon_3)
drugs_all_regulon_6_details <- drug_network_details(data_mat = drugs_all_regulon_6)

drugs_all_program_3_details <- drug_network_details(data_mat = drugs_all_program_3)
drugs_all_program_6_details <- drug_network_details(data_mat = drugs_all_program_6)

################### Get Functional Enrichment ##################################
# get pathway information
drugs_all_regulon_3_pathways <- drug_pathway_enrichment(input.df = drugs_all_regulon_3_details)
drugs_all_regulon_6_pathways <- drug_pathway_enrichment(input.df = drugs_all_regulon_6_details)

drugs_all_program_3_pathways <- drug_pathway_enrichment(input.df = drugs_all_program_3_details)
drugs_all_program_6_pathways <- drug_pathway_enrichment(input.df = drugs_all_program_6_details)

################### Rank ordering and Final Files  ##################################
drugs_all_regulon_3_final <- rank_order_drug_table(drugs_all_regulon_3_pathways)
drugs_all_regulon_6_final <- rank_order_drug_table(drugs_all_regulon_6_pathways)

drugs_all_program_3_final <- rank_order_drug_table(drugs_all_program_3_pathways)
drugs_all_program_6_final <- rank_order_drug_table(drugs_all_program_6_pathways)

sample.name = "P76156_3"
type = "regulon"
write_csv(drugs_all_regulon_3_final, paste0("output/", sample.name,"/",sample.name,"_",type,"_drug_therapy_activity.csv"))

sample.name = "P76156_6"
type = "regulon"
write_csv(drugs_all_regulon_6_final, paste0("output/", sample.name,"/",sample.name,"_",type,"_drug_therapy_activity.csv"))










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
