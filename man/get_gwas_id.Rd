% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_gwas_id.R
\name{get_gwas_id}
\alias{get_gwas_id}
\title{Use a keyword or ID to retrieve GWAS datasets, then analyze each matching outcome_id}
\usage{
get_gwas_id(input)
}
\arguments{
\item{input}{Can be a keyword, used to match trait in the GWAS dataset; or it can be an outcome_id}
}
\value{
Save the analysis results of each matched outcome_id as a text file
}
\description{
Use a keyword or ID to retrieve GWAS datasets, then analyze each matching outcome_id
}
\examples{
\dontrun{
get_gwas_id("Myocardial infarction")
get_gwas_id("finn-b-N14_ENDOMETRIOSIS")
}
}
