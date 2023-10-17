#' Harmonise and process FinnGen R9 data.
#'
#' This function processes and harmonises the FinnGen R9 data. The function reads each file,
#' formats and harmonises it, and then saves the processed data in a specified directory.
#' It uses parallel computation to process the list of files.
#'
#' @param finn_dir The directory where the FinnGen data files are located.
#' @param save_dir Directory where processed data should be saved. Default is "finn_r".
#' @param cores The number of cores to use for parallel computation. Default is 2.
#'
#' @examples
#' \dontrun{
#' get_finn(finn_dir = "finn", save_dir = "finn_r", cores = 30)
#' }
#' @export
get_finn <- function(finn_dir, save_dir = "finn_r", cores = 30) {

  # Load necessary libraries and data
  library(data.table)
  library(parallel)

  finn_info <- load(system.file("data", "Finn_R9_data.rda", package = "BioCycleMR"))

  # Define the process_finn function
  process_finn <- function(trait) {
    finn <- fread(file = trait, data.table = F)

    finndata <- format_data(finn, type = 'outcome', snp_col = "rsids", beta_col = "beta", se_col = "sebeta",
                            eaf_col = "af_alt", effect_allele_col = "alt", other_allele_col = "ref",
                            pval_col = "pval", gene_col = "nearest_genes", chr_col = "#chrom",
                            pos_col = "pos")
    trait = gsub(finn_dir,"",trait)
    trait_row <- finn_info[grepl(trait, finn_info$path_https),]

    finndata$ncase.outcome <- trait_row$num_cases
    finndata$ncontrol.outcome <- trait_row$num_controls
    finndata$samplesize.outcome <- trait_row$num_cases + trait_row$num_controls

    outcome <- trait_row$name
    finndata$outcome <- outcome
    trait_file_name = sub("\\.gz", "", trait)

    # Save the data
    save(finndata, file = paste0(save_dir, "/", "finndataR9_", trait_file_name, ".rda"))
    cat("Finished:", trait_file_name, "\n")
  }

  finn_list <- list.files(finn_dir, full.names = TRUE, recursive = FALSE)
  parallel::mclapply(finn_list, process_finn, mc.cores = cores)
}

# Usage Instructions:
# 1. Ensure you have the required libraries (`data.table` and `parallel`) installed.
# 2. Replace `"YourPackageName"` with the name of your R package (if you're using one).
# 3. Call the function using the given example or by specifying the `finn_dir`, `save_dir`, and `cores` as needed.
