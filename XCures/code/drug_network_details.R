
drug_network_details <- function(data_mat){
  #cat("Collecting network details info\n")

  t1 <- tibble()

  for(drug in unique(data_mat$Drug)){
    cat("Processing ", drug,"\n")
    ## SOC drugs
    drug.1 <- data_mat %>%
      filter(Drug == drug) %>%
      dplyr::select(-IfCancerDriverMutation,RegulatorType,-MutationRegulatorEdge,-MutationRegulonEdge,-SpearmanCorrelationTFAndEigenGene,-`IfCRISPR-Cas9SignificantTF`,-IfDETFForTCGAGBM,-RNAiTFGeneEssentialityAchilles,-Miner_Target_Type) %>%
      base::unique()
    #print(drug.1)

    target.count <-  drug.1 %>%
      dplyr::select(DrugTargetAll) %>%
      separate_rows(DrugTargetAll,sep=",") %>%
      base::unique() %>%
      dplyr::count()

    regulon.count <- drug.1 %>%
      dplyr::select(DrugRegulonAll) %>%
      separate_rows(DrugRegulonAll,sep=",") %>%
      base::unique() %>%
      dplyr::count()

    overactive.regulon.count <- drug.1 %>%
      dplyr::select(DrugRegulonOverActive) %>%
      filter(DrugRegulonOverActive != "NA") %>%
      separate_rows(DrugRegulonOverActive,sep=",") %>%
      base::unique() %>%
      dplyr::count()

    underactive.regulon.count <- drug.1 %>%
      dplyr::select(DrugRegulonUnderActive) %>%
      filter(DrugRegulonUnderActive != "NA") %>%
      separate_rows(DrugRegulonUnderActive,sep=",") %>%
      base::unique() %>%
      dplyr::count()

    program.count <- drug.1 %>%
      dplyr::select(DrugProgramAll) %>%
      separate_rows(DrugProgramAll,sep=",") %>%
      base::unique() %>%
      dplyr::count()

    overactive.program.count <- drug.1 %>%
      dplyr::select(DrugProgramOverActive) %>%
      filter(DrugProgramOverActive != "NA") %>%
      separate_rows(DrugProgramOverActive,sep=",") %>%
      base::unique() %>%
      dplyr::count()

    underactive.program.count <- drug.1 %>%
      dplyr::select(DrugProgramUnderActive) %>%
      filter(DrugProgramUnderActive != "NA") %>%
      separate_rows(DrugProgramUnderActive,sep=",") %>%
      base::unique() %>%
      dplyr::count()

    mutations <- drug.1 %>%
      dplyr::select(MutationSymbol) %>%
      mutate(Mutations = paste0(unique(MutationSymbol), collapse=":")) %>%
      dplyr::pull() %>%
      base::unique()

    regulators <- drug.1 %>%
      dplyr::select(RegulatorSymbol) %>%
      mutate(Regulators = paste0(unique(RegulatorSymbol), collapse=":")) %>%
      dplyr::pull() %>%
      base::unique()

    regulator.type <- drug.1 %>%
      dplyr::select(RegulatorType) %>%
      mutate(RegulatorType = paste0(unique(RegulatorType), collapse=":")) %>%
      dplyr::pull() %>%
      base::unique()


    # if(drug %in% previous_treatments()){
    #   past.treatment <- paste(tags$span(style = "color: blue;","YES"))#paste(icon("check-circle", class="text-success"))
    # } else {
    #   past.treatment <- paste(tags$span(style = "color: gray;","NO"))#paste(icon("times-circle", class="text-danger"))
    # }

    regulon.summary <- paste("O:",overactive.regulon.count,
                             "U:",underactive.regulon.count,
                             "A:",regulon.count)
    program.summary <- paste("O:",overactive.program.count,
                             "U:",underactive.program.count,
                             "A:",program.count)

    # drug.class <- drug.1 %>%
    #   select(TargetClass)
    # #drug.class <- str_replace_all(string = drug.class$TargetClass, pattern = "__", replacement=",")
    # drug.class <- sub(x=drug.class$TargetClass, pattern = "__", replacement=",")
    #
    # drug.target <- drug.1 %>%
    #   pull(DrugTargetAll)

    # drug.target <- strsplit(drug.target, split = ",")[[1]]


    #print(drug.target)
    GBM.trials <- drug.1 %>%
      dplyr::select(DrugTrialPhaseGBM) %>%
      mutate(DrugTrialPhaseGBM =  dplyr::if_else(DrugTrialPhaseGBM == "#N/A", "None", paste0(DrugTrialPhaseGBM)))


    antiCancer.phaseIV <- drug.1 %>%
      dplyr::select(Cancer_PhaseIV) %>%
      base::unique() %>%
      mutate(`Cancer_PhaseIV` = dplyr::if_else(Cancer_PhaseIV == "AntiCancerPhaseIV", paste("TRUE"), paste("FALSE") )) %>%
      pull()#paste(collapse=",", sep="")


    t1 <- rbind(t1, cbind(
      Drug = drug,
      `Drug Mechanism` = unique(drug.1$DrugMechOfAction),
      `Target Gene(s)` = unique(drug.1$DrugTargetAll),
      `Drug Constrained Regulon Activity` = unique(drug.1$DrugConstrainedRegulonActivity),
      `Drug Constrained Program Activity` = unique(drug.1$DrugConstrainedProgramActivity),
      `All Regulon(s)` = unique(drug.1$DrugRegulonAll),
      `Overactive Regulon(s)` = unique(drug.1$DrugRegulonOverActive),
      `Underactive Regulon(s)` = unique(drug.1$DrugRegulonUnderActive),
      `All Program(s)` = unique(drug.1$DrugProgramAll),
      `Overactive Program(s)` = unique(drug.1$DrugProgramOverActive),
      `Underactive Program(s)` = unique(drug.1$DrugProgramUnderActive),
      `Regulon Activity Summary` = regulon.summary,
      `Program Activity Summary` = program.summary,
      `Mean RMST Regulon` = round(drug.1$DrugMeanRMSTRegulon,digits = 2),
      Mutations = mutations,
      Regulators = regulators,
      `Regulator Type` = regulator.type,
      `Max GBM Trial Phase` = GBM.trials,
      `Other FDA Appr.` = antiCancer.phaseIV,
      `TargetMutatedInPatient` = paste0(unique(drug.1$TargetMutatedInPatient), collapse = ":"),
      `MutatedInPatient` = paste0(unique(drug.1$MutatedInPatient),collapse = ":")
      )
      )


  }
  # soc <- if_else(input$soc_select, "TRUE", "FALSE")
  # mmTrials <- if_else(input$myelomaTrial_select, "TRUE", "FALSE")
  # anticancerApproved <- if_else(input$anticancer_select, "TRUE", "FALSE")

  return(unique(t1))
}

