#' chisq performs BDA with a Chi squared test
#'
#' @param dataChisq Temporary binary matrix
#' @importFrom stats chisq.test
#'
#' @keywords internal
#' @usage NULL
#'
chisq <- function(dataChisq){
  contingency_table <- table(dataChisq$Gene,dataChisq$status)
  if(sum(dim(contingency_table))!=4){
    return(c("X-squared" = NA, "P" = 1))
  }else{
    res <- chisq.test(contingency_table)
    return(c("X-squared" = unname(res$statistic), "P" = unname(res$p.value)))
  }
}

