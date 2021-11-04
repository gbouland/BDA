#' phi performs BDA with a Pearson's correlation
#'
#' @param dataPhi Temporary binary matrix
#' @importFrom stats cor.test
#'
#' @keywords internal
#' @usage NULL
#'
phi <- function(dataPhi){
  res <- cor.test(as.numeric(dataPhi$Gene),as.numeric(dataPhi$status))
  return(c("Phi" = unname(res$estimate),
           "t" = unname(res$statistic),
           "P" = unname(res$p.value)))
}

