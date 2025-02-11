#' ELViS : An R Package for Estimating Copy Number Levels of Viral Genome Segments Using Base-Resolution Read Depth Profile.
#'
#' Base-resolution copy number analysis of viral genome. Utilizes base-resolution read depth data
#' over viral genome to find copy number segments with two-dimensional segmentation approach.
#' Provides publish-ready figures, including histograms of read depths, coverage line plots
#' over viral genome annotated with copy number change events and viral genes, and heatmaps
#' showing multiple types of data with integrative clustering of samples.
#'
#' @section Functions:
#' - get_depth_matrix : Generate a read depth matrix of positions x samples from input BAM files list.
#' - run_ELViS : Run ELViS using input raw depth matrix.
#' - raw depth matrix.
#'
#' @name ELViS
#' @rdname ELViS
"_PACKAGE"
