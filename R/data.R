#' ELViS Toy Example - Metadata
#'
#' Metadata of samples in the toy examples
#'
#' @format ## `toy_example`
#' A data frame with 120 rows and 6 columns:
#' \describe{
#'   \item{VarType}{Variant type. Set to control if there is no variant}
#'   \item{Copies_Altered}{The number of copies duplicated or deleted}
#'   \item{Event_Size}{Variant size}
#'   \item{Mean_Depth}{Mean read depth}
#' }
#' @usage
#' data(toy_example)
"toy_example"

#' ELViS Toy Example - Total Aligned Base
#'
#' Total aligned base both to host and viral genome.
#'
#' @format ## `total_aligned_base__host_and_virus`
#' A vector of length 120
#' @usage
#' data(total_aligned_base__host_and_virus)
"total_aligned_base__host_and_virus"

#' ELViS Toy Example - Base-Resolution Raw Read Depth
#'
#' Base-resolution raw read depth profile over viral genome
#'
#' @format ## `mtrx_samtools_reticulate`
#' A matrix with 7906 rows and 120 columns
#' @usage
#' data(mtrx_samtools_reticulate)
"mtrx_samtools_reticulate"


#' ELViS Toy Example - Run Result
#'
#' List containing ELViS run result
#'
#' @format ## `ELViS_toy_run_result`
#' A list of
#' \describe{
#'   \item{is_reduced_output}{Indicates whether this results is a reduced form}
#'   \item{final_output}{Final base-resolution segmentation output}
#'   \item{final_call}{Indices of samples in which copy number variants were detected}
#'   \item{new_Y_p2}{Normalized read depth}
#' }
#' @usage
#' data(ELViS_toy_run_result)
"ELViS_toy_run_result"
