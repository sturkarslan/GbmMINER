## Read cohort RNASeq data
rnaseq2 <- read.delim("~/Downloads/gdac.broadinstitute.org_GBM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/GBM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt", sep = "\t")

## Get column names
colnames.rnaseq2 <- colnames(rnaseq2)
myseq <- seq(3, length(colnames.rnaseq2), by = 3)

## Multiply this value by 1e6 to get the TPM
rnaseq5 <- rnaseq2 %>%
  select(1,myseq) %>%
  slice(-1) %>%
  type_convert() %>%
  separate(Hybridization.REF, sep = "[|]", into = c("Symbol", "Entrez_Id")) %>%
  mutate(across(where(is.numeric), ~ .x * 1e+06))

## Load Gene mapping identifiers
identifiers <- read.delim("~/Documents/GitHub/GbmMINER/data/identifier_mappings.txt", sep="\t")

# join ids with cohort data
rnaseq_ids <- left_join(rnaseq5, identifiers, by=c(Entrez_Id = 'Name')) %>%
  filter(!is.na(`Preferred_Name`)) %>%
  select(-Source) %>%
  relocate(`Preferred_Name`) %>%
  rename("Ensembl_Id" = `Preferred_Name` )


### Read all patient data
filenames_df <- data_frame(filename = list.files("/Volumes/omics4tb2/SYGNAL/XCures/patients_processed_data", pattern = ".genes.results",recursive = T,full.names = T)) %>%
  mutate(file_contents = map(filename, ~ read_delim(.,delim = "\t")) ) %>%
  unnest(c(file_contents))

# Get TPM values
patients_data <- filenames_df %>%
  separate(filename, sep = "/", into = c("1","2","3","4","5","6","Patient")) %>%
  separate(gene_id, sep="_", into=c("Ensembl_Id","Symbol")) %>%
  select(Patient, Ensembl_Id,TPM) %>%
  pivot_wider(id_cols=Ensembl_Id, names_from = Patient, values_from=TPM )

# combine patient and cohort data
combined_rnaseq <- left_join(rnaseq_ids, patients_data, by="Ensembl_Id") %>%
  #filter(!is.na(TPM)) %>%
  select(-Symbol, -Entrez_Id) %>%
  relocate(matches("TL-"), .after = Ensembl_Id)

combined_rnaseq[is.na(combined_rnaseq)] <- 0

## Write into a file for network activioty calculations
write_csv(combined_rnaseq, file="/Volumes/omics4tb2/SYGNAL/XCures/patients_processed_data/allpatients_withcohort_tpm.csv")









rnaseq6 <- melt(combined_rnaseq[,c(3:174,176)], id.vars = c("Entrez_Id"), measure.vars = colnames(combined_rnaseq)[c(3:174,176)]) %>%
  filter (variable != "Entrez_Id")


plot_data <- rnaseq6 %>%
  filter(variable %in% unique(rnaseq6[,"variable"])[c(1:5,172)]) %>%
  select(variable, value)
plot_data$value <- as.numeric(plot_data$value)

p <- ggplot(plot_data, aes(x=variable, y=value, group=variable))
p <- p + geom_boxplot()
p <- p + lims(y=c(0,5000))
p






rnaseq6 <- melt(rnaseq5[,2:173], id.vars = c("Entrez_Id"), measure.vars = colnames(rnaseq5[2:173])) %>%
  filter (variable != "Entrez_Id")

plot_data <- rnaseq6 %>%
  filter(variable %in% unique(rnaseq6[,"variable"])[1:10]) %>%
  select(variable, value)
plot_data$value <- as.numeric(plot_data$value)

p <- ggplot(plot_data, aes(x=variable, y=value, group=variable))
p <- p + geom_boxplot()
p <- p + lims(y=c(0,5000))
p



separate(as.tibble(rnaseq2[5,1]), "|", into = c("Symbol", "Entrez_Id") )


rnaseq4 <- rnaseq2
colnames(rnaseq4) <- rnaseq2[1,]
rnaseq4 <- rnaseq4[-1,]


zz <- rnaseq2 %>%
  select(contains(c("TCGA.06.0125","Hybridization")))

yy <- rnaseq %>%
  select(contains(c("TCGA.06.0125", "HUGO")))

ss <- rnaseq3 %>%
  select(contains(c("TCGA.06.0125", "Hybridization")))



identifiers <- read.delim("~/Documents/GitHub/GbmMINER/data/identifier_mappings.txt", sep="\t")

patient_data <- read.delim("/Volumes/omics4tb2/SYGNAL/XCures/patients_processed_data/TL-19-87E3E1/RNA/TL-19-87E3E1.genes.results", sep="\t") %>%
  select(gene_id, TPM) %>%
  separate(gene_id, "_", into=c("Ensembl_Id", "Symbol"))

rnaseq$entrez <- as.character(rnaseq$Entrez_Gene_Id)

rnaseq_ids <- left_join(rnaseq5, identifiers, by=c(Entrez_Id = 'Name')) %>%
  filter(!is.na(`Preferred_Name`)) %>%
  select(-Source) %>%
  relocate(`Preferred_Name`) %>%
  rename("Ensembl_Id" = `Preferred_Name` )

combined_rnaseq <- left_join(rnaseq_ids, patient_data, by="Ensembl_Id") %>%
  filter(!is.na(TPM))

write_csv(combined_rnaseq, file="/Volumes/omics4tb2/SYGNAL/XCures/patients_processed_data/TL-19-87E3E1/RNA/TL-19-87E3E1_withcohort_tpm.csv")



#rnaseq <- read.delim("~/Downloads/gbm_tcga/data_mrna_seq_v2_rsem.txt", sep = "\t")

#rnaseq3 <- read.delim("~/Downloads/gdac.broadinstitute.org_GBM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/GBM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", sep = "\t")

# patient_files <- fs::dir_ls("/Volumes/omics4tb2/SYGNAL/XCures/patients_processed_data", regexp = "*.genes.results", recurse = T) %>%
#   map(read_delim, delim = "\t") %>%
#   purrr::reduce(rbind)
#
# patient_data <- patient_files %>%
#   distinct() %>%
#   select(gene_id,TPM)

