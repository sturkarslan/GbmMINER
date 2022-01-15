
format_target_data <- function(target, count){
  targets <- strsplit(target, split = ",")[[1]]
  # We add target again for multiple targets to be able to match
  targets <- append(targets, target)
  # Load clinical trials data
  ct <- read.delim("../data/Glioblastoma_ClinicalTrials_01092022.txt", header = T, sep = "\t")

  # Load addtional evidence data
   evidence1 <- read.csv("../data/Additional_Evidences_Gene.txt", sep = "\t",fileEncoding="latin1") %>%
     filter(Gene %in% targets) %>%
     mutate(Title = if_else(Type == "Publication",
                            paste0('<a href="https://pubmed.ncbi.nlm.nih.gov/', PMID,'">', Title,'</a>'),
                            if_else(Type == "Preprint",
                                    paste0('<a href="https://www.biorxiv.org/content/',PMID,'">', Title,'</a>'),
                            paste0('<a href="https://patents.google.com/patent/', PMID,'">', Title,'</a>'))))

   if(length(evidence1$Gene != 0)){
     evidence2 <- evidence1
   } else {
     evidence2 <- data.frame(NULL)
   }


   # if(length(evidence1$Gene != 0)){
   #   evidence2 <- evidence1 %>%
   #     mutate(evidence.title = paste0('<a href="https://pubmed.ncbi.nlm.nih.gov/',PMID,'">',Evidence,'</a>', collapse=", ")) %>%
   #     pull(evidence.title)
   # } else{
   #   evidence2 <- paste0("No additional evidence")
   # }

  # Select target specific information
  network_selected <- network_info %>%
    filter(`Target Gene(s)` == target) %>%
    separate_rows(`Drug Action`,sep = ":") %>%
    unique()

  # get target drugs
    target.drugs <- network_selected %>%
      pull(Drug)

  # Regulon activity counts
  activity.counts <- network_selected%>%
    separate(`Regulon Activity Summary`, sep=" ", into=c("i1","Over","i3","Under","i5","All"), remove = F) %>%
    select("Over","Under","All") %>%
    unique()

  # Regulon activity score
  activity.score <- network_selected %>%
    pull(`Drug Constrained Regulon Activity`) %>%
    unique()

  # List of target genes
  target.list <- network_selected %>%
    pull(`Target Gene(s)`) %>%
    unique()

  # List of separated target genes
  target.list.seperated <- network_selected %>%
    separate_rows(`Target Gene(s)`, sep=",") %>%
    pull(`Target Gene(s)`) %>%
    unique()

  # GBM Trials info
  gbmtrials <- network_selected %>%
    pull(max_glioblastoma.multiforme_phase)

  # Target Mutations in patient
  target_mutated <- network_selected %>%
    pull(TargetMutatedInPatient) %>%
    unique()

  # FDA approval info
  anticancerPhaseIVs <- network_selected %>%
    select(`Other FDA Appr.`) %>%
    unique() %>%
    mutate(`Other FDA Appr.` = if_else(`Other FDA Appr.` == TRUE, paste("**Phase IV** trial for some other cancer(s)"), paste(" **NOT** in Phase IV trials for any cancer."))) %>%
    pull()

  # Hallmarks information
  hallmarks.list <- network_selected %>%
    dplyr::select(Hallmarks) %>%
    unique() %>%
    pull()

  # Hallmarks seperated
  hallmarks.sep1 <- str_split(hallmarks.list, pattern = ":")[[1]]
  hallmarks.sep2 <- unlist(lapply(hallmarks.sep1, function(x) x[x != "NA"]))
  hallmarks.sep = if_else(!is.na(hallmarks.sep2), paste0('<a href="https://www.gsea-msigdb.org/gsea/msigdb/cards/',hallmarks.sep2,'">', hallmarks.sep2,'</a>', collapse=" "), paste0("None"))

  # Clinical trials
  clinical_trials_raw <- tibble()
  for(target.drug in target.drugs){
    ct1 <- filter(ct, grepl(target.drug, toupper(Interventions)))
    clinical_trials_raw <- bind_rows(clinical_trials_raw, ct1)
  }

  clinical_trials <- clinical_trials_raw %>%
    mutate(Title = paste0('<a href="',URL,'">', Title, '</a>')) %>%
    select(NCT.Number,Title,Status,Interventions,Phases,First.Posted)




  # Build the drug summary table
  target.table <- network_selected %>%
    select(Drug, max_glioblastoma.multiforme_phase, `Other FDA Appr.`, `Drug Action`, `Toxicity Class`, `MEDDRA Code`) %>%
    mutate(Drug_temp = Drug) %>%
    mutate(Drug = paste0('<a name=#',Drug,' href="https://www.ebi.ac.uk/chembl/g/#search_results/all/query=', Drug,'">', Drug, '</a>')) %>%
    mutate(Action = `Drug Action`) %>%
    mutate(`MEDDRA Code` = sub(".0", "", `MEDDRA Code`)) %>%
    mutate(`Toxicity Class` = mapply(
      function(i, j)
        paste0('<a class="btn btn-warning btn-sm" role="button" href="http://purl.bioontology.org/ontology/MEDDRA/', i, '" target="_blank">', j, '</a>', collapse = " "),
      str_split(`MEDDRA Code`, pattern = ":"),
      str_split(`Toxicity Class`, pattern = ":")
      )) %>%
    mutate(`Toxicity Class` = dplyr::if_else(`MEDDRA Code` != "NA", `Toxicity Class`,  paste0('<a class="btn btn-warning btn-sm disabled" role="button" href=" " target="_blank">', 'NA', '</a>', collapse = " ")      )) %>%
    mutate(`Max GBM Phase` = ifelse( is.na(max_glioblastoma.multiforme_phase), 0 ,
      max_glioblastoma.multiforme_phase
      )) %>%
    mutate(`Max GBM Phase` = as.numeric(`Max GBM Phase`)) %>%
    mutate(`FDA Approved` = if_else(`Other FDA Appr.` == "TRUE", "TRUE", "FALSE")) %>%
    mutate(`Clinical Trials` = if_else(
      `Max GBM Phase` != 0,
      paste0(
        '<a href="https://beta.clinicaltrials.gov/search?patient=',
        Drug_temp,
        ' AND  Glioblastoma&locStr=&distance=0"> <button type="button" class="btn btn-primary btn-sm">Clinical Trials</button></a>'
      ),
      paste0(
        '<a href="https://beta.clinicaltrials.gov/search?patient=',
        Drug_temp,
        ' AND Glioblastoma&locStr=&distance=0"> <button type="button" class="btn btn-primary btn-sm disabled">Clinical Trials</button></a>'
      )
    )) %>%
    select(-max_glioblastoma.multiforme_phase,
           -`Other FDA Appr.`,
           -Drug_temp,
           -`MEDDRA Code`,
           -Action)


  ## Build the UI

  # Target
  cat('\n')
  #cat(sep="", paste0('<a name=#', target, '</a>')
  cat(sep="","# ", paste0(count), ". Target: ", paste0('<a name=', target, ' href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=',target,'">',target,'</a>', collapse=", "), "\n")

  # Drug
  cat('\n')
  cat(paste0('<button type="button" class="btn btn-info">  Drugs <span class="badge badge-light"> ',length(target.drugs),' </span> </button> \n'))

  # Summary buttons
  # Active regulons
  cat('\n')
  cat(paste0('<button type="button" class="btn btn-primary"> Active Regulons <span class="badge badge-light"> ',activity.counts$Over,' </span> </button>\n'))

  # Regulon Activity
  cat('\n')
  cat(paste0('<button type="button" class="btn btn-primary"> Regulon Activity <span class="badge badge-light"> ',round(activity.score, digits=3),' </span> </button>\n'))

  # Hallmarks
  cat('\n')
  cat(unique(paste0('<button type="button" class="btn btn-primary"> Hallmarks <span class="badge badge-light"> ',if_else(!is.na(hallmarks.sep2),length(hallmarks.sep), as.integer(0)),' </span> </button>\n')))

  # Target Mutated in Patient
  cat('\n')
  if("TRUE" %in% target_mutated){
    cat(unique(paste0('<button type="button" class="btn btn-warning"> Mutated in patient <span class="badge badge-light"> ','Yes',' </span> </button>\n')))
  } else{
    cat(unique(paste0('<button type="button" class="btn btn-primary disabled"> Mutated in patient <span class="badge badge-light"> ','No',' </span> </button>\n')))
  }


  # Summary Table
  cat("### Drugs \n")
  cat("\n\n")
  cat(knitr::knit_print(DT::datatable(target.table,
                                     #caption = "Drugs matching to targets",
                                     escape=F,
                                     width = "100%",
                                     height = "100%",
                                     options = list(dom='')
                                     ) %>%
                             formatStyle(
                               "FDA Approved",
                               target = "cell",
                               color = styleEqual(
                                 c("TRUE","FALSE"), c('green', 'gray')
                               )
                             )

  ))
  cat("\n\n")



  # Evidence summary
  cat("### Model Evidence \n")
  cat(sep="", "**",activity.counts$Over,"** regulons containing the target(s), **", paste0(target.list, collapse=","), "** are overactive in this patient. \n")
  cat("\n\n")
  #cat(sep="", "* ##### **Evidence:** ", "**",activity.counts$Over,  " ** regulons containing the target(s), **", paste0(target.list, collapse=","), "** are overactive in this patient.\n")

  # Hallmarks summary
  cat("### Hallmarks Enrichment \n")
  if(!is.na(hallmarks.sep[1])){
    cat(hallmarks.sep[1], "\n")
  } else{
    cat("None \n")
  }
  cat("\n\n")

  # Additional Evidence
  cat("### Literature Evidence \n")
  if(length(evidence2$Gene) != 0){
    cat(knitr::knit_print(DT::datatable(evidence2,
                                        escape=F,
                                        width = "100%",
                                        height = "100%",
                                        options = list(dom='')
                                        )
                          )
        )
  } else {
    cat("No additional evidence was found \n")
    }
  cat("\n\n")


  # Clinical Trials
  cat("### Glioblastoma Clinical Trials \n")
  if(length(clinical_trials$NCT.Number) != 0){
    cat(knitr::knit_print(DT::datatable(clinical_trials,
                                        escape=F,
                                        width = "100%",
                                        height = "100%",
                                        options = list(dom='p',
                                                       pageLength = 5)
    )
    )
    )
  } else {
    cat("No glioblastoma clinical trials were found for these drugs. \n")
  }
  cat("\n\n")


  # if(length(evidence2) != 0){
  #   cat(sep="", "* ##### **Literature Evidence: **",  paste0(evidence2, collapse = "\n", sep=""),  "\n")
  # } else {
  #   cat(sep="", "* ##### **Literature Evidence: **",  "No additional evidence was found",  "\n")
  #
  # }

  # ID conversion for plots
  # target_id <- identifiers %>%
  #   filter(Name == target) %>%
  #   pull(Preferred_Name)

    # print(
    #   network_selected %>%
    #     mutate(`MEDDRA Code` = sub(".0", "", `MEDDRA Code`)) %>%
    #     kable() %>%
    #     column_spec(1, color = "blue",
    #           link = paste0('https://www.ebi.ac.uk/chembl/g/#search_results/all/query=', Drug)) %>%
    #     column_spec(1, color = "blue",
    #           link = mutate(`Toxicity Class` = mapply(
    #       function(i, j)
    #         paste0("http://purl.bioontology.org/ontology/MEDDRA/', i", collapse = " "),
    #       str_split(`MEDDRA Code`, pattern = ":"),
    #       str_split(`Toxicity Class`, pattern = ":")
    #     ))) %>%
    #      kable_styling()
    #
    # )








  cat("### Network Plots \n")
  # Plots
  # Activity histogram
  ppp.act <- plot_activity_histogram(target.drugs[1])

  # target gene expression
  ppp.exp <- plot_expression(target=target)

  # combined plot
  ppp <- grid.arrange(ppp.act,ppp.exp,ncol=2)

  cat("\n\n")
  cat("\n")
  cat("***")
  #cat('<div class="tocify-extend-page" data-unique="tocify-extend-page" style="height: 0;"></div>"')


}
