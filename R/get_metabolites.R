#' Get Metabolites Data Processing
#'
#' This function processes metabolite data using a series of helper functions.
#' It loads necessary libraries, defines helper functions, reads a list of
#' file identifiers, processes the data in parallel, and then saves the
#' results to an `.rda` file.
#'
#' @param index_file The path to the file that contains the list of identifiers.
#' @param out_dir The directory where metabolite data files are stored.
#' @param pvalue The p-value threshold for subset.
#' @param clump_kb The clumping kilobase parameter.
#' @param clump_r2 The clumping r2 parameter.
#' @param count_try_max Maximum number of tries for clumping.
#' @param num_cores Number of cores to use for parallel processing.
#'
#' @return NULL. Results are saved to an `.rda` file.
#' @export
#'
#' @examples
#' \dontrun{
#' get_metabolites(index_file = "metabolitis_id.txt")
#' }
get_metabolites <- function(use_preprocessed = TRUE, out_dir = "metabolitis/",
                            pvalue = 5e-6, clump_kb = 10000, clump_r2 = 0.001,
                            count_try_max = 50, num_cores = 64) {

  if (use_preprocessed) {
    load( "metabolites_raw.rda")
    return(NULL)
  }

  # Load the index from the rda file
  load("metabolitis_id.rda")

  clump_id <- function(exposure_dat, clump_kb, clump_r2, count_try_max) {
    count_try <- 0
    repeat {
      count_try <- count_try + 1
      try({
        cat(count_try, "Clumping for exposure IV", "\n")
        exposure_dat1 <- TwoSampleMR::clump_data(exposure_dat, clump_kb = clump_kb, clump_r2 = clump_r2)
      })
      if ((exists("exposure_dat1") && !is.null(exposure_dat1) && nrow(exposure_dat1) > 0) || count_try == count_try_max) {
        cat("after clump, the number of exposure IV:", dim(exposure_dat1)[1], "\n")
        break
      }
      Sys.sleep(0.5)
    }
    return(exposure_dat1)
  }

  loop_clump <- function(index) {
    meta_list <- list()
    for (i in index) {
      q <- read.table(file = paste0(out_dir, i), header = TRUE)
      w <- tidyverse::subset(q, q$P.value < pvalue)
      id <- substr(i, 1, 8)
      if (length(w[, 1]) > 0) {
        utils::write.csv(w, "w.csv")
        cat("Clumping local data for:", id, "\n")
        exposure_dat <- TwoSampleMR::read_exposure_data("w.csv", sep = ",",
                                                        snp_col = "MarkerName",
                                                        beta_col = "Effect",
                                                        se_col = "StdErr",
                                                        eaf_col = "Freq1",
                                                        effect_allele_col = "Allele1",
                                                        other_allele_col = "Allele2",
                                                        pval_col = "P.value", clump = FALSE)
        cat("The number of exposure IV:", dim(exposure_dat)[1], "\n")
        exposure_dat <- clump_id(exposure_dat, clump_kb, clump_r2, count_try_max)
        meta_list[[id]] <- exposure_dat
        cat(tidyverse::str_c(id, " has been processed!"), "\n")
      }
    }
    return(meta_list)
  }

  loop_clump_parallel <- function(index) {
    process_each_index <- function(i) {
      tryCatch({
        return(loop_clump(i))
      }, error = function(e) {
        message(paste("Error processing index", i, ":", e$message))
        return(list())
      })
    }
    result_list <- c(parallel::mclapply(index, process_each_index, mc.cores = num_cores, mc.preschedule = FALSE))
    return(result_list)
  }

  calculate_list_length <- function(lst) {
    if (is.null(lst) || length(lst) == 0) {
      cat("The input list is empty or NULL.\n")
      return(NULL)
    }
    num_elements <- length(lst)
    stats <- data.frame(TotalElements = num_elements)
    return(stats)
  }

  # Process data
  meta_list_all <- loop_clump_parallel(index)
  SNP_stats_all = calculate_list_length(meta_list_all)

  # Save results
  utils::save(meta_list_all, SNP_stats_all, file = "data/metabolites_raw_0.rda")

  # Return
  return(NULL)
}
