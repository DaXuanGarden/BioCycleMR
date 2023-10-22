#' Convert local VCF file to Two-Sample MR input format and export it as a CSV file
#'
#' @param vcf_file The local VCF file to be converted
#' @param type The type of GWAS study, "exposure" or "outcome"
#' @return Export the converted data as a CSV file
#' @examples
#' \dontrun{
#' get_local("ieu-a-2.vcf.gz", "exposure")
#' }
#' @export
get_local <- function(vcf_file, type) {
  # Downloaded VCF file should be placed into the working path in advance
  bim_VCF <- VariantAnnotation::readVcf(vcf_file)

  # Convert VCF to TwoSampleMR format
  bmi = gwasglue::gwasvcf_to_TwoSampleMR(vcf=bim_VCF, type=type) # Change to "exposure" or "outcome" based on the data

  # Check if the conversion is successful
  print(head(bmi))

  # Export the data
  fwrite(bmi, "bmi.csv")
}
