% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{get_new_baseline}
\alias{get_new_baseline}
\title{Get new baselines according to criteria user designates}
\usage{
get_new_baseline(result, mode = "longest")
}
\arguments{
\item{result}{Run result}

\item{mode}{Indicate how new baseline should be set ("longest","shortest")}
}
\value{
a integer vector indicating new baseline index for each sample
}
\description{
Get new baselines according to criteria user designates
}
\examples{

# its usage example is given in vignette in detail

data(ELViS_toy_run_result)
result <- ELViS_toy_run_result

get_new_baseline(result,mode="longest")

}
