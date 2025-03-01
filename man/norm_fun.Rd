% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{norm_fun}
\alias{norm_fun}
\title{Normalization - scaling by certain quantile}
\usage{
norm_fun(x, probs = 0.75)
}
\arguments{
\item{x}{numeric vector to normalize.}

\item{probs}{a single numeric value of probabilities in \verb{[0,1]} used for normalization.(Default = \code{0.75})}
}
\value{
numeric vector of normalized values
}
\description{
Normalization - scaling by certain quantile
}
\examples{

norm_fun(seq_len(5))
# [1] 0.25 0.50 0.75 1.00 1.25

}
