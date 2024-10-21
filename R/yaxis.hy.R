#' from R package SCISSOR  : https://github.com/hyochoi/SCISSOR/blob/377789106d34200baa1d4d826952e2ad5313ba3a/R/yaxis.hy.R
#' @noRd
yaxis.hy <- function(mat){
  #  mat : d by n matrix
  tempmax <- max(mat) ;
  tempmin <- min(mat) ;
  templen <- tempmax-tempmin ;
  return(c(tempmin-0.002*templen, tempmax+0.002*templen)) ;
}
