#' Calculate F-values and MAF for given data.
#'
#' This function calculates F-values and MAF for a given dataset. It provides options to
#' specify the F-value threshold, MAF threshold, and sample size.
#'
#' @param dat_object An R object containing the dataset to be processed.
#' @param F_value A numeric threshold for the F-value. Default is 10.
#' @param maf_threshold A numeric threshold for the minor allele frequency (MAF). Default is 0.01.
#' @param samplesize A numeric value specifying the sample size. If NULL, the function attempts
#'                   to derive it from the dat_object. Default is NULL.
#' @param use_preprocessed A logical value indicating whether to use preprocessed data. Default is FALSE.
#'
#' @return A list containing processed datasets and statistics.
#' @examples
#' \dontrun{
#' result <- get_f_maf(dat_object = immu_cell_raw)
#' }
#' @export
get_f_maf <- function(dat_object,
                      F_value = 10,
                      maf_threshold = 0.01,
                      samplesize = NULL,
                      use_preprocessed = FALSE) {

  # If not using preprocessed data, load the necessary data
  if (!use_preprocessed) {
    immu_cell_raw <- dat_object
  }

  # Function definitions
  get_f_maf_function <- function(dat, F_value = 10, samplesize = NULL) {
    log <- is.na(dat$eaf.exposure)
    log <- unique(log)

    if (length(log) == 1) {
      if (log == TRUE) {
        cat("Data does not contain eaf and cannot calculate F-statistic.\n")
        return(dat)
      }
    }

    if (is.null(dat$beta.exposure[1]) || is.na(dat$beta.exposure[1])) {
      cat("Data does not contain beta and cannot calculate F-statistic.\n")
      return(dat)
    }

    if (is.null(dat$se.exposure[1]) || is.na(dat$se.exposure[1])) {
      cat("Data does not contain se and cannot calculate F-statistic.\n")
      return(dat)
    }

    if (is.null(dat$samplesize.exposure[1]) || is.na(dat$samplesize.exposure[1])) {
      if (is.null(samplesize) || is.na(samplesize)) {
        cat("Data does not contain samplesize and cannot calculate F-statistic.\n")
        return(dat)
      }
    }

    if ("FALSE" %in% log && !is.null(dat$beta.exposure[1]) && !is.na(dat$beta.exposure[1]) &&
        !is.null(dat$se.exposure[1]) && !is.na(dat$se.exposure[1])) {

      if (is.null(samplesize) || is.na(samplesize)) {
        samplesize <- dat$samplesize.exposure[1]
      }

      R2 <- (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$beta.exposure^2)) /
        ((2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$beta.exposure^2)) +
           (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$se.exposure^2) * samplesize))

      F <- (samplesize - 2) * R2 / (1 - R2)
      MAF <- ifelse(dat$eaf.exposure < 0.5, dat$eaf.exposure, (1 - dat$eaf.exposure))
      dat$R2 <- R2
      dat$F <- F
      dat$maf <- MAF
      return(dat)
    }
  }

  calculate_df_row_stats <- function(df_list) {
    if (length(df_list) == 0) {
      cat("The input list is empty.\n")
      return(NULL)
    }

    num_rows <- sapply(df_list, function(df) nrow(df))
    max_rows <- max(num_rows)
    min_rows <- min(num_rows)
    median_rows <- median(num_rows)
    num_rows_equal_to_1 <- sum(num_rows == 1)
    stats <- data.frame(Max = max_rows, Min = min_rows, Median = median_rows, NumRowsEqualToOne = num_rows_equal_to_1)
    return(stats)
  }

  select_f_maf_function <- function(dat, F, maf) {
    dat = dat[dat$F > F,]
    dat = dat[dat$maf > maf,]
    return(dat)
  }

  # Setup parallel processing
  cl <- makeCluster(30)
  registerDoParallel(cl)
  clusterExport(cl, list = ls())

  # Perform calculations
  immu_cell_f <- parLapply(cl, immu_cell_raw, function(data) get_f_maf_function(data, F_value, samplesize))
  immu_cell_f_select <- parLapply(cl, immu_cell_f, function(data) select_f_maf_function(data, F_value, maf_threshold))

  # Stop parallel processing
  stopCluster(cl)

  # Calculate statistics
  SNP_stats_f <- calculate_df_row_stats(immu_cell_f_select)

  result_list <- list(
    immu_cell_f = immu_cell_f,
    immu_cell_f_select = immu_cell_f_select,
    SNP_stats_f = SNP_stats_f
  )

  return(result_list)
}
