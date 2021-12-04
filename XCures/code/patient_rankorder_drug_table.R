library(tidyverse)
rank_order_drug_table <- function(input.df){
  ordered_table <- input.df %>%
    separate(DrugTrialPhaseGBM, into = c("1","2","GBM_Trial_Phase"), fill = "left", remove = F) %>%
    mutate(GBM_Trial_Phase = if_else(GBM_Trial_Phase == "None", 0, as.numeric(GBM_Trial_Phase))) %>%
    mutate(GBM_Trial_Phase_Score = scales::rescale(GBM_Trial_Phase, to=c(0,1))) %>%
    mutate(FDA_Approved_Score = if_else(`Other FDA Appr.` == TRUE, 1,0)) %>%
    mutate(Target_Mutation_Score = if_else(`TargetMutatedInPatient` == TRUE, -1,0)) %>%
    mutate(Target_Type_Score = if_else(`TargetType` == "RegulatorGene", 0,0)) %>%
    mutate(Regulon_RMST_Score = scales::rescale(`Mean RMST Regulon`, to=c(0,1)) ) %>%
    mutate(Overall_Score = `Drug Constrained Regulon Activity` +
             `Drug Constrained Program Activity` +
             GBM_Trial_Phase_Score +
             FDA_Approved_Score +
             Regulon_RMST_Score +
             Target_Mutation_Score +
             Target_Type_Score) %>%
    mutate(Overall_Score = scales::rescale(Overall_Score, to=c(0,1)) ) %>%
    arrange(desc(`Drug Constrained Regulon Activity`,Overall_Score))
  return(ordered_table)
}
