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
#' get_effect("~/path_to_directory", immune_ref_data)
#'
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

process_outcomes <- function(finn_name_list, immune_ref) {
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
