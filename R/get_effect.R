#' get_effect function
#'
#' This function processes .rda files from a specified directory and generates
#' effect size estimates. For each .rda file, a subdirectory is created if it
#' does not already exist. The effect size estimates are then calculated for each
#' subdirectory and written to a .csv file. The function returns the effect size
#' estimates for all .rda files as a data frame.
#'
#' @param directory A character string specifying the directory path that
#' contains the .rda files.
#' @param immune_ref A data frame that is used for merging with the final effect size
#' estimates data frame.
#' @return A data frame containing the effect size estimates for all .rda files.
#'
#' @examples
#' \dontrun{
#' get_effect("~/path_to_directory", immune_ref_data)
#'}
#' @export
get_effect <- function(directory, immune_ref) {
  finn_list <- list.files(directory, full.names = TRUE, recursive = FALSE)
  finn_name_list <- sapply(list.files(directory, full.names = FALSE, recursive = FALSE),FUN = function(x){sub("\\.rda", "",x)})

  for (i in finn_name_list){
    if (!dir.exists(i)){
      dir.create(i)
    }
  }

  process_outcomes(finn_name_list, immune_ref)
}

#' Process Outcomes from Finn Name List and Immune References
#'
#' This function processes the outcomes by using the provided `finn_name_list` and `immune_ref`
#' parameters. It is designed to facilitate the analysis and transformation of outcomes
#' based on specific criteria or references.
#'
#' @param finn_name_list A character vector containing the list of Finn names
#' to be processed. Each name in the list should correspond to a specific outcome
#' that is to be processed.
#'
#' @param immune_ref A data frame or similar structure containing the immune reference
#' data. This reference data is used to map and transform the outcomes from
#' `finn_name_list` during processing.
#'
#' @return A processed list or data frame (you should specify what the function returns
#' based on its logic, e.g., "A data frame containing the processed outcomes.")
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- process_outcomes(sample_finn_names, sample_immune_ref)
#' }
process_outcomes <- function(finn_name_list, immune_ref)
{
  effect_list = list()

  for (i in seq_along(finn_name_list)) {
    outcome_name = finn_name_list[i]
    load(paste0(outcome_name,"/SNP.rda"))
    dataList = split(har_true, har_true$id.exposure)

    effect <- read.table(paste0(outcome_name,"/effect.csv"),sep=",",header = T)
    table1 <- effect %>%
      filter(FDR < 0.05, method %in% c("Wald ratio","Inverse variance weighted")) %>%
      generate_odds_ratios() %>%
      mutate(`OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)",or,or_lci95, or_uci95))
    write.csv(table1, file=paste0(outcome_name,"/fdr.csv"), row.names=F)
    effect_list = list.append(effect_list,table1)
    cat("finished: ",outcome_name,"  ",dim(table1)[1],"\n")
  }

  effect_dat = do.call(rbind,effect_list)
  effect_dat$rank = NULL
  effect_dat = merge(effect_dat, immune_ref, by = "id.exposure")
  write.csv(effect_dat, file="./effect_immune.csv", row.names=F)

  return(effect_dat)
}
