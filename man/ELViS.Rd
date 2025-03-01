% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ELViS.R
\docType{package}
\name{ELViS}
\alias{ELViS-package}
\alias{ELViS}
\title{ELViS : An R Package for Estimating Copy Number Levels of Viral Genome Segments Using Base-Resolution Read Depth Profile.}
\description{
Base-resolution copy number analysis of viral genome. Utilizes base-resolution read depth data
over viral genome to find copy number segments with two-dimensional segmentation approach.
Provides publish-ready figures, including histograms of read depths, coverage line plots
over viral genome annotated with copy number change events and viral genes, and heatmaps
showing multiple types of data with integrative clustering of samples.
}
\section{Functions}{

\itemize{
\item get_depth_matrix : Generate a read depth matrix of positions x samples from input BAM files list.
\item run_ELViS : Run ELViS using input raw depth matrix.
\item raw depth matrix.
}
}

\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/hyochoi/ELViS}
  \item Report bugs at \url{https://github.com/hyochoi/ELViS/issues}
}

}
\author{
\strong{Maintainer}: Jin-Young Lee \email{jlee307@uthsc.edu} (\href{https://orcid.org/0000-0002-5366-7488}{ORCID}) [copyright holder]

Authors:
\itemize{
  \item Hyo Young Choi \email{hyoyoung.choi@uthsc.edu} (\href{https://orcid.org/0000-0002-7627-8493}{ORCID}) [copyright holder]
  \item Katherine A. Hoadley \email{hoadley@med.unc.edu} (\href{https://orcid.org/0000-0002-1216-477X}{ORCID})
  \item D. Neil Hayes \email{Neil.Hayes@uthsc.edu} (\href{https://orcid.org/0000-0001-6203-7771}{ORCID}) [funder, copyright holder]
}

Other contributors:
\itemize{
  \item Xiaobei Zhao \email{xiaobei_zhao@uthsc.edu} (\href{https://orcid.org/0000-0002-5277-0846}{ORCID}) [contributor]
  \item Jeremiah R. Holt \email{jholt41@uthsc.edu} (\href{https://orcid.org/0000-0002-5201-5015}{ORCID}) [contributor]
}

}
