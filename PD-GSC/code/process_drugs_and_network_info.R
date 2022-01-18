library(progress)
library(AnnotationDbi)
options(connectionObserver = NULL)
library(org.Hs.eg.db)

##### Read and process HTP Drug screening info
hts.drugs <- read_csv(here("data/HTS_AUC_data_10-1-2021.csv")) %>%
  rename(Drug = 1) %>%
  mutate(Drug = toupper(Drug)) %>%
  pull(Drug)

##### Load Opentargets Drug mappings based on ell ensembl protein coding genes
ot.drugs <- read_csv(here("data/opentargets_Sep2721.csv"))

#####Kavya's Drug list for GBM
gbm.drugs <- read_csv(here("data/MINER_MicroLowessRNATMM.08.24.2020/AllDrugTargetCombosForSigMinerCausalMappingsOpenTargetsAndXCures.csv")) %>%
  pull(Drug)


hts_drugs_notinOT <- read_csv(here("data/HTP_Drugs_notInOT.csv")) %>%
  separate_rows(Targets, sep = ":") %>%
  filter(!is.na(Targets))

hts_added_drugs <- data.frame()
for(i in 1:length(hts_drugs_notinOT$molecule_name)){
  # get the target gene for the drug
  target <- hts_drugs_notinOT[i,"Targets"] %>%
    pull(Targets)
  print(target)
  # get the emsembl id for the gene
  ensemblid <- try(AnnotationDbi::mapIds(org.Hs.eg.db, keys= target, column="ENSEMBL", keytype="SYMBOL", multiVals="first")[[1]])

  # Create a new row for binding
  new.row <- bind_cols(hts_drugs_notinOT[i,], target_id=ensemblid)
  # bind them together
  hts_added_drugs <- bind_rows(hts_added_drugs, new.row)

}

ot.drugs1 <-ot.drugs %>%
  separate_rows(target_id, sep = ":") %>%
  filter(!is.na(target_id))

pb1 <- progress_bar$new(total = length(ot.drugs1$target_id),
                       format ="[:bar] :current/:total (:percent)")
ot.drugs2 <- data.frame()
for(i in 1:length(ot.drugs1$target_id)){
  pb1$tick()
  Sys.sleep(1 / 100)
  target_id <- ot.drugs1[i,"target_id"] %>%
    pull(target_id)
  my.row <- ot.drugs1[i,]
  target_symbol <- try(AnnotationDbi::mapIds(org.Hs.eg.db, keys= target_id, column="SYMBOL", keytype="ENSEMBL", multiVals="first")[[1]])
  ot.drugs2 <- bind_rows(ot.drugs2, bind_cols(my.row, target_symbol))
}

ot.drugs3 <- ot.drugs2 %>%
  rename('approved_name' = 'target_symbol') %>%
  rename("...18" = 'approved_name')


hts_added_drugs1 <- hts_added_drugs %>%
  dplyr::rename(approved_name = Targets ) %>%
  group_by(molecule_name) %>%
  #mutate(target_id = paste0(unique(target_id), collapse=":")) %>%
  #mutate(approved_name = paste0(unique(approved_name), collapse = ":")) %>%
  ungroup() %>%
  unique()

## Merge large OT drug list with HTP drugs absent in the original list
ot.drugs4 <- full_join(ot.drugs3, hts_added_drugs1)

## Separate multiple targetids into additional rows
#ot.drugs.exp <- ot.drugs1 %>%
#  separate_rows(target_id, sep = ":")


##### HTP drugs that are not in OT Drug list
hts.drugs[!hts.drugs %in% ot.drugs]
##### HTP drugs that are not in GBM Drug list
hts.drugs[hts.drugs %in% gbm.drugs]

##### Load regulons for GBM MINER run
regulons <- read_json(here("data/MINER_MicroLowessRNATMM.08.24.2020/regulons.json"), simplifyVector = T)
##### Load Programsfor GBM MINER run
programs <- read_json(here("data/MINER_MicroLowessRNATMM.08.24.2020/transcriptional_programs.json"), simplifyVector = T)
##### Load Regulator to Regulon relationshipfor GBM MINER run
regulonsDF <- read_csv(here("data/MINER_MicroLowessRNATMM.08.24.2020/regulonDf.csv"))



#### Add regulon and program informatio for eac drug in HTS



###################
### Loop through each row to add regulon and program info
pb <- progress_bar$new(total = length(ot.drugs.exp$target_id),
                       format ="[:bar] :current/:total (:percent)")

ot.drugs.network <- data.frame()
for(i in 1:length(ot.drugs4$target_id)){
  pb$tick()
  Sys.sleep(1 / 100)

  # get the target Gene ID
  my.target <- ot.drugs4[i,"target_id"]
  my.row <- ot.drugs4[i,]

  # get regulons containing target gene id
  my.regulons <- unique(names(regulons[str_detect(regulons, my.target)]))
  my.regulator <- unique(names(regulonsDF[str_detect(regulons, my.target)]))

  zz1 <- bind_cols(my.row, Regulon = my.regulons)

  zz2 <- zz1 %>%
    rowwise() %>%
    mutate(Program = paste(Program = names(programs[str_detect(programs,paste0("\\b",Regulon,"\\b"))]), collapse=":")) %>%
    separate_rows(Program, sep=":")

  ot.drugs.network <- bind_rows(ot.drugs.network, zz2)
}


## Modfy RegulonsDF for merung with larget drug info
regulonsDF2 <- regulonsDF %>%
  dplyr::select(Regulon_ID, Regulator) %>%
  unique()

## Combine egulators with regulons
ot.drugs.network2 <- ot.drugs.network %>%
  mutate(Regulon = as.numeric(Regulon)) %>%
  left_join(regulonsDF2, by = c("Regulon" = "Regulon_ID")) %>%
  rowwise() %>%
  mutate(RegulatorSymbol = AnnotationDbi::mapIds(org.Hs.eg.db, keys= Regulator, column="SYMBOL", keytype="ENSEMBL", multiVals="first")[[1]])

write_csv(ot.drugs.network2, file="data/opentargets_hts_drug_details_network_activity.csv")



#### COncordance
## cutoffs for responder, non-responder
cut_up <- 2/3    ## 0.666...
cut_down <-1/3   ## 0.777

doResponseMap <- function(reggie, network_all,
                          scale_up=cut_up, scale_down=cut_down,
                          net_up=cut_up, net_down=cut_down,
                          labby=" ") {

  res <- reggie[, .(responder=sum(`1minusAUC_scaled`>=scale_up &
                                    network_all>=net_up),
                    nonresponder=sum(`1minusAUC_scaled` <= scale_down &
                                       network_all <= net_down),
                    neutral=sum(`1minusAUC_scaled` > scale_down &
                                  `1minusAUC_scaled` < scale_up &
                                  network_all > net_down &
                                  network_all < net_up),
                    sample=length(`1minusAUC_scaled`)),
                by=.(Drug)]
  res %<>%
    dplyr::mutate(total=responder + nonresponder + neutral)

  res %>%
    dplyr::select(Drug, sample, everything()) %>%
    dplyr::arrange(-responder, nonresponder) ->
    count

  ## create column labelling concordance between sygnal and auc
  reggie$concordance <- NA
  reggie$concordance[reggie$`1minusAUC_scaled`>=scale_up &
                       reggie$network_all>=net_up] <- "responder"
  reggie$concordance[reggie$`1minusAUC_scaled` <= scale_down &
                       reggie$network_all <= net_down] <- "nonresponder"
  reggie$concordance[reggie$`1minusAUC_scaled` > scale_down &
                       reggie$`1minusAUC_scaled` < scale_up &
                       reggie$network_all > net_down &
                       reggie$network_all < net_up] <- "neutral"

  res %<>%
    dplyr::arrange(-responder, -nonresponder)
}



kk <- read_csv(here("data/NCGA_OfIC50Drugs_TargetsAveraged_PDGSC.csv"))
kk_melt <- melt(kk)

kk_mat <- as.matrix(kk %>%
  dplyr::select(Drug, !ends_with("_rank")))
rownames(kk_mat) <- kk_mat[,1]
kk_mat <- kk_mat[,-1]

pheatmap(kk)

p <- ggplot(kk, aes())



























### Loop through each row to add regulon and program info
pb <- progress_bar$new(total = length(ot.drugs.exp$target_id),
                       format ="[:bar] :current/:total (:percent)")
ot.drugs.new <- data.frame()
for(i in 1:length(ot.drugs.exp$target_id)){
  pb$tick()
  Sys.sleep(1 / 100)

  # get the target Gene ID
  my.target <- ot.drugs.exp[i,"target_id"] %>%
    pull(target_id)

  # get regulons regulated by this target
  my.regulated.regulons <- regulonsDF %>%
    filter(Regulator %in% my.target) %>%
    pull(Regulon_ID) %>%
    unique()
  my.regulated.regulons <- paste0(my.regulated.regulons, collapse=":")

  # get regulons containing target gene id
  my.regulons <- unique(names(regulons[str_detect(regulons, my.target)]))

  # get programs containing these regulons
  my.programs <- vector()
  my.regulators <- vector()
  for(my.regulon in my.regulons){
    my.program <- names(programs[str_detect(programs, my.regulon)])
    my.programs <- append(my.programs, my.program) %>%
      unique()

    # get regulators of these regulons
    my.regulator <- regulonsDF %>%
      filter(Regulon_ID %in% my.regulon) %>%
      pull(Regulator) %>%
      unique()
    my.regulators <- append(my.regulators, my.regulator) %>%
      unique()

  }


  my.row <- ot.drugs.exp[i,] %>%
    mutate(Regulon = if_else(
      length(my.regulons) != 0,
      paste0(my.regulons,collapse = ":"),"")) %>%
    mutate(Program = if_else(
      length(my.programs) != 0,
      paste0(my.programs, collapse = ":"),"")) %>%
    mutate(Regulates_Regulons = my.regulated.regulons)

  ot.drugs.new <- bind_rows(ot.drugs.new, my.row)

}

write_csv(ot.drugs.new, file="data/opentargets_hts_drug_details.csv")
t1 <- drug_info %>%
  separate_rows(Regulon, sep = ":") %>%
  separate_rows(Program, sep=":") %>%
  separate_rows(Regulates_Regulons, sep=":")

