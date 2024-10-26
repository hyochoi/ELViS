#' from SCISSOR https://github.com/hyochoi/SCISSOR/blob/377789106d34200baa1d4d826952e2ad5313ba3a/R/rhoHuber.R#L4
#' @noRd
rhoHuber <- function(x,c=2.1){
  # x is a univariate sample
  # c is the tuning constant
  # output is rho(x)
  #
  rho <- (x/c)^2
  rho[rho > 1] <- 1
  1.54^2*rho
}
