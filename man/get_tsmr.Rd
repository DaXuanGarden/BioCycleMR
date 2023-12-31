% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tsmr.R
\name{get_tsmr}
\alias{get_tsmr}
\title{get_tsmr function}
\usage{
get_tsmr(immu_cell_f_select, finn_r_dir, cores = 64)
}
\arguments{
\item{immu_cell_f_select}{Input exposure data, expected to be a list of data frames.}

\item{finn_r_dir}{Directory containing .rda files for each outcome.}

\item{cores}{Number of cores to use for parallel operations. Default is 64.}
}
\value{
This function does not return anything but creates multiple directories, .rda and .csv files as a part of its operation.
}
\description{
This function is intended to be used as a part of BioCycleMR package.
It performs multiple complex operations on the input data.
}
\examples{
\dontrun{
get_tsmr(immu_cell_f_select, finn_r_dir, cores = 64)
}

}
