% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{run_ELViS}
\alias{run_ELViS}
\title{Run ELViS using input raw depth matrix}
\usage{
run_ELViS(
  X,
  N_cores = min(10L, detectCores()),
  reduced_output = TRUE,
  verbose = FALSE,
  save_intermediate_data = FALSE,
  save_dir = "save_dir",
  overwrite = FALSE
)
}
\arguments{
\item{X}{Raw depth matrix of position x samples}

\item{N_cores}{The number of cores to use (Default : \code{min(10L,detectCores())})}

\item{reduced_output}{logical indicating whether to return only reduced output}

\item{verbose}{logical whether to print out information for debugging}

\item{save_intermediate_data}{logical indicating whether to save intermediate data as rds file. (default : FALSE)}

\item{save_dir}{Name of the directory to save intermediate files in. (default : "save_dir")}

\item{overwrite}{logical indicating whether to overwrite intermediate files. (default : FALSE)}
}
\value{
list containing ELViS run results
}
\description{
Run ELViS using input raw depth matrix
}
\examples{

data(mtrx_samtools_reticulate)
th<-50
filtered_mtrx <- filt_samples(mtrx_samtools_reticulate,th=th,smry_fun=max)

result <- run_ELViS(filtered_mtrx[,seq_len(5)],N_cores=1L)

}
