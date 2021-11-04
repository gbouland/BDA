#' scBDA single core BDA
#'
#' @param tmp_binaryMatrix Temporary binary matrix
#' @param tmp_samplesheet Temporary samplesheet
#' @param contrast.by Contrast column
#' @param ident.1 First cell population of contrasting cell populations.
#' @param ident.2 Second cell population of contrasting cell populations.
#' @param covariates covariates
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @keywords internal
#' @usage NULL
#'
scBDA <- function(tmp_binaryMatrix ,tmp_samplesheet, contrast.by, ident.1, ident.2, covariates, method){


  popu.results <- matrix(nrow = nrow(tmp_binaryMatrix),
                         ncol = switch(method,"LR" = 4,
                                              "fisher" = 4,
                                              "chisq" = 2,
                                              "phi" = 3))

  pb <-  txtProgressBar(min = 0, max = nrow(tmp_binaryMatrix), initial = 0,
                        style = 3)

  tmp_binaryMatrix <- as.matrix(tmp_binaryMatrix)
  for(i in 1:nrow(popu.results)){
    setTxtProgressBar(pb,i)
    tmpData <- data.frame(status = tmp_samplesheet[,contrast.by],
                          Gene = tmp_binaryMatrix[i,],
                          tmp_samplesheet[,covariates])
    colnames(tmpData) <- c("status","Gene",covariates)
    tmpData$status <- factor(tmpData$status, levels = c(ident.1,ident.2))
    popu.results[i,] <- switch(method,
                               "LR" = LR(tmpData,covariates),
                               "fisher" = fisher(tmpData),
                               "chisq" = chisq(tmpData),
                               "phi" = phi(tmpData))

  }
  return(popu.results)
}
