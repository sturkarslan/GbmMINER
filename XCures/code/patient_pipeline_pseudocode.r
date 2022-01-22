 ### Patient Analysis Pipeline
 ## STEP 0 Run wrapper for patient analysis
 # File: 01_analyze_patient.R
 # sample.id: sample.name
 analyze_patient <- function(sample.id, activity_filter=0.8){
    
    ## STEP 1. Maps network activities to drugs
    # File: patient_map_drugs.R
    # drug_info_pr: Master drug list mapped to regulons
    # sample.id : sample name
    # type : base mapping on the regulon or program activity
    # disease: Wehether to use all regulons or disease relevant ones
    # network_dir: Directory containing patients regulon and program activities
    # mutation_files: patients mutation information if available
    # analysis_type: type of the molecular profiling. Used for mutation formatting.
    map_drugs <- function(drug_info_pr, sample.id = sample.id, type="regulon", disease=TRUE, network_dir=network_dir_input,mutation_files=mutation_files, analysis_type="TEMPUS"){
  
        # Given the regulon and program activity, gets drug therapy activity
        # File: utilities.R
        # regulon: regulon activity matrix
        # program: program activity matrix
        getDrugTherapyActivity <- function(drug_info_pr, regulon=reg_out, program=prog_out){

        }
    }

    ## STEP 2. Get network details for each target
    # File: patient_get_network_details.R
    # data_mat: mapped drug and regulon data frame from map_drugs function
    drug_network_details <- function(data_mat = drugs_all){

    }

    ## STEP 3. Pathway enrichment for the active regulons of each drug 
    # File: patient_drug_pathway_enrichment.R
    # input.df: output from filtered network details from STEP 2.
    drug_pathway_enrichment <- function(input.df = drugs_all_details_filt){

    }

}
 
 
 