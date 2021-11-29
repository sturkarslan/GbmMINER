
drug_network_details <- function(data_mat){
  #cat("Collecting network details info\n")

  t1 <- tibble()

  for(drug in unique(data_mat$Drug)){
    cat("Processing ", drug,"\n")
    ## SOC drugs
    drug.1 <- data_mat %>%
      filter(Drug == drug) %>%
      select(-MutationSymbol,-IfCancerDriverMutation,-RegulatorSymbol,RegulatorType,-MutationRegulatorEdge,-MutationRegulonEdge,-SpearmanCorrelationTFAndEigenGene,-`IfCRISPR-Cas9SignificantTF`,-IfDETFForTCGAGBM,-RNAiTFGeneEssentialityAchilles,-MutatedInPatient,-TargetMutatedInPatient,-Miner_Target_Type) %>%
      unique()
    #print(drug.1)

    target.count <-  drug.1 %>%
      select(DrugTargetAll) %>%
      separate_rows(DrugTargetAll,sep=",") %>%
      unique() %>%
      dplyr::count()

    regulon.count <- drug.1 %>%
      select(DrugRegulonAll) %>%
      separate_rows(DrugRegulonAll,sep=",") %>%
      unique() %>%
      dplyr::count()

    overactive.regulon.count <- drug.1 %>%
      select(DrugRegulonOverActive) %>%
      filter(DrugRegulonOverActive != "NA") %>%
      separate_rows(DrugRegulonOverActive,sep=",") %>%
      unique() %>%
      dplyr::count()

    underactive.regulon.count <- drug.1 %>%
      select(DrugRegulonUnderActive) %>%
      filter(DrugRegulonUnderActive != "NA") %>%
      separate_rows(DrugRegulonUnderActive,sep=",") %>%
      unique() %>%
      dplyr::count()

    program.count <- drug.1 %>%
      select(DrugProgramAll) %>%
      separate_rows(DrugProgramAll,sep=",") %>%
      unique() %>%
      dplyr::count()

    overactive.program.count <- drug.1 %>%
      select(DrugProgramOverActive) %>%
      filter(DrugProgramOverActive != "NA") %>%
      separate_rows(DrugProgramOverActive,sep=",") %>%
      unique() %>%
      dplyr::count()

    underactive.program.count <- drug.1 %>%
      select(DrugProgramUnderActive) %>%
      filter(DrugProgramUnderActive != "NA") %>%
      separate_rows(DrugProgramUnderActive,sep=",") %>%
      unique() %>%
      dplyr::count()

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
      select(DrugTrialPhaseGBM) %>%
      mutate(DrugTrialPhaseGBM =  dplyr::if_else(DrugTrialPhaseGBM == "#N/A", "None", paste0(DrugTrialPhaseGBM)))


    antiCancer.phaseIV <- drug.1 %>%
      select(Cancer_PhaseIV) %>%
      unique() %>%
      mutate(`Cancer_PhaseIV` = dplyr::if_else(Cancer_PhaseIV == "AntiCancerPhaseIV", paste("TRUE"), paste("FALSE") )) %>%
      pull()#paste(collapse=",", sep="")

    t1 <- rbind(t1, cbind(
      Drug = drug,
      #`Target Gene(s)` = paste0(drug.target, collapse=", "),
      #`Drug Class` = drug.class,
      `Program Activity Summary` =program.summary,
      `Regulon Activity Summary` = regulon.summary,
      #`Mean RMST Program` = round(drug.1$DrugMeanRMSTProgram,digits = 2),
      #`Mean RMST Regulon` = round(drug.1$DrugMeanRMSTRegulon,digits = 2),
      #`Mean RMST Regulon p-value` = signif(drug.1$DrugRMSTRegulonMeanPvalue,digits = 2),
      #`Network Constrained Regulon Activity` = round(drug.1$DrugConstrainedRegulonActivity,digits = 2),
      `Max GBM Trial Phase` = GBM.trials,
      `Other FDA Appr.` = antiCancer.phaseIV
      #SOC = drug.1$soc,
      #`Add Comments` = "Comments"
      )
      )


  }
  # soc <- if_else(input$soc_select, "TRUE", "FALSE")
  # mmTrials <- if_else(input$myelomaTrial_select, "TRUE", "FALSE")
  # anticancerApproved <- if_else(input$anticancer_select, "TRUE", "FALSE")

  return(unique(t1))
}

