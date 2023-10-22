#' Extract Instruments for Given Outcomes
#'
#' This function extracts instruments based on given outcomes using data from `ieugwasr` package.
#'
#' @param outcomes A character vector of outcomes.
#' @param p1 A numeric threshold for p-value.
#' @param clump A logical indicating whether to perform clumping.
#' @param p2 A secondary p-value threshold for clumping.
#' @param r2 A numeric threshold for the r-squared value.
#' @param kb A numeric indicating the window size in kilobases for clumping.
#' @param access_token A token to access the `ieugwasr` database. By default, it checks for a valid access token.
#' @param force_server A logical indicating whether to force the server connection.
#'
#' @return A data frame with extracted instrument variables.
#' @export
#'
#' @examples
#' \dontrun{
#'   result <- extract_instruments(outcomes = c("BMI", "LDL"), p1 = 5e-8)
#' }
extract_instruments <- function (outcomes, p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001,
                                 kb = 10000, access_token = ieugwasr::check_access_token(),
                                 force_server = FALSE)
{
  outcomes <- ieugwasr::legacy_ids(unique(outcomes))
  d <- ieugwasr::tophits(outcomes, pval = p1, clump = clump,
                         r2 = r2, kb = kb, force_server = FALSE, access_token = access_token)
  if (nrow(d) == 0)
    return(NULL)
  d$phenotype <- paste0(d$trait, " || id:", d$id)
  d <- format_data(d, type = "exposure", snps = NULL, phenotype_col = "phenotype",
                   snp_col = "rsid", chr_col = "chr", pos_col = "position",
                   beta_col = "beta", se_col = "se", eaf_col = "eaf", effect_allele_col = "ea",
                   other_allele_col = "nea", pval_col = "p", samplesize_col = "n",
                   min_pval = 1e-200, id_col = "id")
  d$data_source.exposure <- "igd"
  return(d)
}
