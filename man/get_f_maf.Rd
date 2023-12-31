% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_f_maf.R
\name{get_f_maf}
\alias{get_f_maf}
\title{Calculate F-values and MAF for given data.}
\usage{
get_f_maf(
  dat_object,
  F_value = 10,
  maf_threshold = 0.01,
  samplesize = NULL,
  use_preprocessed = FALSE
)
}
\arguments{
\item{dat_object}{An R object containing the dataset to be processed.}

\item{F_value}{A numeric threshold for the F-value. Default is 10.}

\item{maf_threshold}{A numeric threshold for the minor allele frequency (MAF). Default is 0.01.}

\item{samplesize}{A numeric value specifying the sample size. If NULL, the function attempts
to derive it from the dat_object. Default is NULL.}

\item{use_preprocessed}{A logical value indicating whether to use preprocessed data. Default is FALSE.}
}
\value{
A list containing processed datasets and statistics.
}
\description{
This function calculates F-values and MAF for a given dataset. It provides options to
specify the F-value threshold, MAF threshold, and sample size.
}
\examples{
\dontrun{
result <- get_f_maf(dat_object = immu_cell_raw)
}
}
