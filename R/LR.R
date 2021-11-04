#' LR performs BDA with a logistic regression
#'
#' @param dataLR Temporary binary matrix
#' @param covariates Temporary samplesheet
#' @importFrom stats glm
#'
#' @keywords internal
#' @usage NULL
#'
LR <- function(dataLR, covariates){
  if(length(covariates) == 0){
    covs <- ""
  }else{
    covs <- sprintf("+ %s", paste0(covariates, collapse = "+"))
  }
  res <- glm(sprintf("Gene ~ status %s",covs),family = "binomial",data = dataLR)
  sum <- summary(res)
  return(sum$coefficients[2,])
}
