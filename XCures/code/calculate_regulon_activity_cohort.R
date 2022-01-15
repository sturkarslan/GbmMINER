## Collect regulon activity for all patients in the cohort
# Regulon and program activity files
all_reg_activity <- read_csv("data/regulon_activity_all.csv") %>%
  dplyr::rename("Regulon" = 1)

all_prog_activity <- read_csv("data/program_activity_all.csv") %>%
  dplyr::rename("Program" = 1)

samples <- all_reg_activity %>%
  dplyr::select(!matches("Regulon")) %>%
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

