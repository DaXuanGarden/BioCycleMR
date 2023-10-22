#' Use a keyword or ID to retrieve GWAS datasets, then analyze each matching outcome_id
#'
#' @param input Can be a keyword, used to match trait in the GWAS dataset; or it can be an outcome_id
#' @return Save the analysis results of each matched outcome_id as a text file
#' @examples
#' \dontrun{
#' get_gwas_id("Myocardial infarction")
#' get_gwas_id("finn-b-N14_ENDOMETRIOSIS")
#' }
#' @export
get_gwas_id <- function(input) {
  # Load the built-in GWAS dataset list
  #load(system.file("data", "available_datasets.rda", package = "BioCycleMR"))

  # Determine whether the input is a keyword or an outcome_id
  if (any(input == available_datasets$id)) { # input is outcome_id
    matched_datasets <- available_datasets[available_datasets$id == input, ]
  } else { # input is a keyword
    matched_datasets <- available_datasets[grep(input, available_datasets$trait, ignore.case = TRUE), ]
  }

  # Get all matched outcome_ids
  outcome_ids <- matched_datasets$id

  # Analyze each outcome_id
  for (outcome_id in outcome_ids) {
    write.table(analyze_outcome_immune(outcome_id = outcome_id, max_retries = 50), file = paste0(outcome_id, ".txt"), sep = "\t", quote = F, row.names = F)
  }
}
