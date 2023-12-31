% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_metabolites.R
\name{get_metabolites}
\alias{get_metabolites}
\title{Get Metabolites Data Processing}
\usage{
get_metabolites(
  use_preprocessed = TRUE,
  out_dir = "metabolitis/",
  pvalue = 5e-06,
  clump_kb = 10000,
  clump_r2 = 0.001,
  count_try_max = 50,
  num_cores = 64
)
}
\arguments{
\item{out_dir}{The directory where metabolite data files are stored.}

\item{pvalue}{The p-value threshold for subset.}

\item{clump_kb}{The clumping kilobase parameter.}

\item{clump_r2}{The clumping r2 parameter.}

\item{count_try_max}{Maximum number of tries for clumping.}

\item{num_cores}{Number of cores to use for parallel processing.}

\item{index_file}{The path to the file that contains the list of identifiers.}
}
\value{
NULL. Results are saved to an `.rda` file.
}
\description{
This function processes metabolite data using a series of helper functions.
It loads necessary libraries, defines helper functions, reads a list of
file identifiers, processes the data in parallel, and then saves the
results to an `.rda` file.
}
\examples{
\dontrun{
get_metabolites(index_file = "metabolitis_id.txt")
}
}
