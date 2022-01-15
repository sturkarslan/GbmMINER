calculate_stats <- function(data_file){
  # get targets
  targets_field <- unique(data_file$`Target Gene(s)`)
  targets <- unique(unlist(lapply(targets_field, function(i) strsplit(i, split = ",")[[1]])))
  # get drugs
  drugs <- unique(data_file$Drug)

  # get overactive regulons
  overactive_regulons <- data_file$`Overactive Regulon(s)`
  overactive_regulons <- unique(unlist(str_split(overactive_regulons,pattern = ",")))

  # get underactive regulons
  underactive_regulons <- data_file$`Underactive Regulon(s)`
  underactive_regulons <- unique(unlist(str_split(underactive_regulons,pattern = ",")))

  # get all regulons
  all_regulons <- data_file$`All Regulon(s)`
  all_regulons <- unique(unlist(str_split(all_regulons,pattern = ",")))

  return(list(targets=unique(targets), drugs=drugs, oa_regs = overactive_regulons, ua_regs = underactive_regulons, all_reg = all_regulons))

}
