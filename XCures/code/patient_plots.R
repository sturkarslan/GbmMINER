
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
drug_info_pr <- drug_info_pr %>%
  select(Drug,)


# Regulon and program activity files
all_reg_activity <- read_csv("data/regulon_activity_all.csv") %>%
  dplyr::rename("Regulon" = 1)

all_prog_activity <- read_csv("data/program_activity_all.csv") %>%
  dplyr::rename("Program" = 1)

samples <- all_reg_activity %>%
  select(!matches("Regulon")) %>%
  colnames()

all_patient_reg_activity <- tibble(Drug = unique(drug_info_pr$Drug))
all_patient_prog_activity <- tibble()
for(sample in samples){
  cat(sample, "\n")
  reg.df <- all_reg_activity %>%
    dplyr::select(Regulon, ends_with(sample)) %>%
    dplyr::rename("Regulon" = 1, "activity" = 2)

  prog.df <- all_prog_activity %>%
    dplyr::select(Program, ends_with(sample)) %>%
    dplyr::rename("Program" = 1, "activity" = 2)

  drug_reg <- getDrugTherapyActivity(drug_info_pr, regulon=reg.df, program=prog.df)
  drug_reg <- drug_reg %>%
    dplyr::select(Drug,DrugConstrainedRegulonActivity) %>%
    base::unique() %>%
    ungroup() %>%
    dplyr::rename_with(drug_reg, .fn = ~paste0(sample), .cols = 2)


  all_patient_reg_activity <- left_join(all_patient_reg_activity, drug_reg, by="Drug")

}

write_csv(all_patient_reg_activity, "data/drug_constrained_network_activity_allcohort.csv")

## extract target constrained regulon activity
patient_score <- read_csv(paste0("output/P76156_3/P76156_3_regulon_drug_therapy_activity.csv"))

patient_score %>%
  dplyr::filter(Drug=="BEVACIZUMAB") %>%
  dplyr::select(`Drug Constrained Regulon Activity`) %>%
  unlist() ->
  p_reg_score

activity <-  all_patient_reg_activity
drug_act <- activity[activity$Drug=="BEVACIZUMAB", pheno$Patient_ID]
drug_act2 <- melt(drug_act)

ppp <- ggplot(drug_act2, aes(x=value)) + geom_histogram() +
  ggtitle("BEVACIZUMAB") + geom_vline(xintercept=p_reg_score, linetype="dotted", color="red", size=1) +
  geom_point(aes(x=p_reg_score, y=0, color="red"), size=3) +
  xlab("drug constrained regulon activity") + labs(color="patient")







pheno <- read_csv("../data/TCGA_Survival_Gbm.csv") %>%
  arrange(desc(GuanScore))
activity <- network_info %>%
  select(Drug, Drug.Constrained.Regulon.Activity)


  ## target: Name of Drug
## pheno: ia12 phenotype file. Samples are ordered by Guan risk top=high, bottom =low
## activity: data.frame of constrained activity with first column being
##           Drug.  It must match the type that is defined in 'target'
## patient_score: table of drug therapy info from mapping patient regulon
##                and program activity to our master drug therapy info file
## sub_pheno: vector of subtype labels that map to pheno
## sub_anno: vector of subtype labels for printing
plotDrugActivityHistogram <- function(target, pheno, activity,
                                      patient_score,
                                      sub_pheno=subtypes,
                                      sub_anno=sub_labels) {

  ## extract target constrained regulon activity
  patient_score %>%
    dplyr::filter(Drug==target) %>%
    dplyr::select(DrugConstrainedRegulonActivity) %>%
    unlist() ->
    p_reg_score

  colnames(activity)[1] <- "Target"
  drug_act <- activity[activity$Target==target, pheno$sample]

  drug_hist <- NULL
  for(sub in 1:length(sub_pheno)) {
    pheno %>%
      dplyr::filter(!!as.name(sub_pheno[sub])==1) %>%
      dplyr::select(sample) %>%
      unlist() ->
      sub_samps

    drug_hist <- rbind(drug_hist,
                       data.frame(subtype=sub_anno[sub],
                                  activity=unlist(drug_act[1, sub_samps])))
  } ## end for sub

  ppp <- ggplot(drug_hist, aes(x=activity, fill=subtype)) + geom_histogram() +
    ggtitle(target) + geom_vline(xintercept=p_reg_score, linetype="dotted", color="red", size=1) +
    geom_point(aes(x=p_reg_score, y=0, color="red"), size=3) +
    xlab("drug constrained regulon activity") + labs(color="patient")

  print(ppp)
}
