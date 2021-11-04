#' fisher performs BDA with a fisher exact
#'
#' @param dataFisher Temporary binary matrix
#' @importFrom stats fisher.test
#'
#' @keywords internal
#' @usage NULL
#'
fisher <- function(dataFisher){
  contingency_table <- table(dataFisher$Gene,dataFisher$status)
  if(sum(dim(contingency_table))!=4){
    return(c("OR" = NA, "CI_low" = NA,"CI_high" = NA,"P" = NA))
  }else{
    res <- fisher.test(contingency_table)
    return(c("OR" = unname(res$estimate), "CI_low" = unname(res$conf.int[1]),
             "CI_high" = unname(res$conf.int[2]),"P" = unname(res$p.value)))
  }
}


