## Load patient specific files
patient_folders_path <- fs::dir_info("/Volumes/omics4tb2/SYGNAL/XCures/patients_network_activities/")
patient_folders <- fs::path_file(patient_folders_path$path)
# Filter samples that will not be anlyzed
#patient_folders_current <- patient_folders[grep("TL-", patient_folders)]
patient_folders_current <- patient_folders[grep("TL-18-31DCF0", patient_folders, invert = T)]
patient_folders_current <- patient_folders[grep("abundance", patient_folders, invert = T)]



network_activities <- fs::dir_ls("/Volumes/omics4tb2/SYGNAL/XCures/patients_network_activities_merged//",regexp="_all_regulon_activity.csv", recurse = T)


patient_activities <- network_activities %>%
  map(read_delim, delim = ",") %>%
  purrr::reduce(cbind)

patient_activities$Regulon <- patient_activities$X1

patient_activities_matrix <- tibble(.name_repair = 'unique', patient_activities) %>%
  dplyr::select(!starts_with("X1")) %>%
  dplyr::select(!starts_with("abundance")) %>%
  dplyr::relocate(Regulon)


common_active <- patient_activities_matrix %>%
  filter(if_all(.fns = ~. == 1)) #277


library(jsonlite)
library(gridExtra)
library(ggplot2)
library(tidyverse)

## Data

regulon <- "113"

patient_id= "TL-18-31DCF0" # -1
patient_id2= "TL-19-0B9A1B" # +1
patient_id3= "TL-19-61DB85" # 0


## Read regulon gene list
regulons <- read_json("../data/MINER_MicroLowessRNATMM.08.24.2020/regulons.json",simplifyVector = T)

## Read cohort expression
cohort_expr <- read.csv("/Volumes/omics4tb2/SYGNAL/GBM-Serdar/MINER_MicroLowessRNATMM.08.24.2020/GbmMicroRNAMergedWithIDsZScored.csv")

# Function to plot Single patient single regulon pliot
patient_exp_single < function(patient_id, regulon){
  # get regulon genes
  my.regulon.genes <- regulons[regulon][[1]]

  ## Load patient's expression data for plotting and then z-score
  patient_exp_filename1 <- paste0("/Volumes/omics4tb2/SYGNAL/XCures/patients_processed_data/allpatients_withcohort_zscores.csv")

  ## read patient exp
  patient_exp <- read_csv(patient_exp_filename1)

  # Get single patient value to combine with cohort for plotting
  patients_exp_single <- patient_exp %>%
    select(X1, paste(patient_id))

  # Get TPM values
  patients_exp_data <- patient_exp %>%
    select(X1, paste(patient_id)) %>%
    filter(X1 %in% my.regulon.genes) %>%
    rename("zscore" = paste(patient_id)) %>%
    mutate(patient= patient_id)

  # filter cohort expression data for regulon genes
  exp_data_cohort <- cohort_expr %>%
    filter(X %in% my.regulon.genes) %>%
    unique()

  ## Add patient data to cohort data
  exp_data_with_patient <- inner_join(exp_data_cohort, patients_exp_single, by=c("X" = "X1") )

  ## melt
  exp_data_cohort_melt <- reshape2::melt(exp_data_with_patient)

  ## Add median and color
  exp_data_final <-  exp_data_cohort_melt %>%
    group_by(variable) %>%
    mutate(median=median(value)) %>%
    mutate(Q4=quantile(value)["100%"]) %>%
    mutate(Color=if_else(variable == patient_id, "red", "gray"))


  ### GENE PLOT
  # plot of cohort expression for regulon genes across patients
  q <- ggplot(exp_data_cohort_melt, aes(x=reorder(X,value,na.rm = TRUE), y=value, group=X))
  q <- q + geom_violin(alpha=0.7, show.legend = F, fill="gray")
  q <- q + geom_point(data=patients_exp_data, aes(x=X1, y=zscore), color="red", inherit.aes = FALSE)
  q <- q + geom_line(data=patients_exp_data, aes(x=X1, y=zscore, group=2),color="red", inherit.aes = FALSE)
  q <- q + theme(axis.text.x = element_text(size=7, angle=45))
  q <- q + labs(x="Gene(s)", y="Gene Expression z-score")
  q

  ### PATIENT PLOT
  # plot of regulon gene exopression across patients
  p <- ggplot(exp_data_final,
              aes(x=reorder(variable,value,na.rm = TRUE),
                  y=value, group=variable,fill=Color,color=Color)
              )
  p <- p + geom_violin(show.legend = F)
  p <- p + geom_vline(xintercept = patient_id, color="red", linetype="dotted")
  p <- p + geom_point(data=exp_data_final, aes(x=reorder(variable,value,na.rm = TRUE), y=median, group=variable),color="blue", inherit.aes = FALSE)
  p <- p + scale_fill_identity(aesthetics = c("fill","color"))
  p <- p + theme(axis.text.x = element_blank())
  p <- p + labs(x="Patients", y="Gene Expression z-score")
  p



}








patient_exp < function(patient_id1, patient_id2, patient_id3, regulon){
  my.regulon <- regulons[regulon][[1]]


  ## Load patient's expression data for plotting and then z-score
  patient_exp_filename1 <- paste0("/Volumes/omics4tb2/SYGNAL/XCures/patients_processed_data/",patient_id1, "/RNA/", patient_id1, ".genes.results")
  patient_exp_filename2 <- paste0("/Volumes/omics4tb2/SYGNAL/XCures/patients_processed_data/",patient_id2, "/RNA/", patient_id2, ".genes.results")
  patient_exp_filename3 <- paste0("/Volumes/omics4tb2/SYGNAL/XCures/patients_processed_data/",patient_id3, "/RNA/", patient_id3, ".genes.results")
  patient_exp_filenames <- c(patient_exp_filename1,patient_exp_filename2,patient_exp_filename3)

  patient_exp <- data_frame(filename = patient_exp_filenames) %>%
    mutate(file_contents = map(filename, ~ read_delim(.,delim = "\t")) ) %>%
    unnest(c(file_contents))

  # Get TPM values
  patients_exp_data <- patient_exp %>%
    separate(filename, sep = "/", into = c("1","2","3","4","5","6","Patient")) %>%
    separate(gene_id, sep="_", into=c("Ensembl_Id","Symbol")) %>%
    select(Patient, Ensembl_Id,TPM) %>%
    pivot_wider(id_cols=Ensembl_Id, names_from = Patient, values_from=TPM )


  patient_plot_data <- patients_exp_data %>%
    dplyr::filter(Ensembl_Id %in% my.regulon)

  pp1 <- melt(patients_exp_data, id.vars = "Ensembl_Id", measure.vars = colnames(patient_plot_data)[2:4]) %>%
    group_by(Ensembl_Id, variable) %>%
    mutate(z_score_group = (value - mean(value)) / sd(value)) %>%   # You can also use scale(value) as pointed out by @RuiBarradas
    ungroup %>%
    mutate(z_score_ungrouped = (value - mean(value)) / sd(value))

    rowwise() %>%
    dplyr::mutate(lapply(zscore = (mean() - mean()/sd()))



  p <- ggplot(pp1, aes(x=Ensembl_Id, y=z_score_ungrouped, group=variable, fill=variable))
  p <- p + geom_bar(alpha=0.5, stat="identity")
  p <- p + guides(color=guide_legend(title="Patient"))
  p <- p + guides(fill=guide_legend(title="Genes"))
  p <- p + labs(x="Gene(s)", y="Gene Expression z-score", title = paste0("patient: ", " regulon: ", regulon))
  #p <- p + geom_hline(yintercept=patient_zscore,linetype="dotted", size=1, color=patient_zscore)
  p

}

write_csv(patient_activities_matrix, file="patient_activities_matrix.csv")




dd = combined_rnaseq %>%
  is.na(select(matches("TL-")))

dd <- read_csv("/Volumes/GoogleDrive/My Drive/Manuscripts/MINER Manuscript/MINER/results_minCorrelation_0o2_50_allFiles/regulons_activity_heatmap.csv")

cc <- read_csv(("/Volumes/omics4tb2/SYGNAL/GBM-Serdar/MINER_MicroLowessRNATMM.08.24.2020/GbmMicroRNAMergedWithIDsZScored.csv"))


library(ComplexHeatmap)


data3 <- data_frame(filename = network_activities) %>%
  mutate(file_contents = map(filename,
                             ~ read_csv(file.path()))
  )
data4 <- unnest(data3)

for (patient in patient_folders_current){



}
