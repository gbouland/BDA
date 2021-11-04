#' mcBDA multi core BDA
#'
#' @param tmp_binaryMatrix Temporary binary matrix
#' @param tmp_samplesheet Temporary samplesheet
#' @param contrast.by Contrast column
#' @param ident.1 First cell population of contrasting cell populations.
#' @param ident.2 Second cell population of contrasting cell populations.
#' @param covariates covariates
#' @param cores Number of cores
#' @import foreach
#' @import doParallel
#' @importFrom parallel makeCluster stopCluster
#' @keywords internal
#' @usage NULL
#'
mcBDA <- function(tmp_binaryMatrix ,tmp_samplesheet, contrast.by, ident.1,
                  ident.2, covariates, cores, method){

  geneList <- split(sample(rownames(tmp_binaryMatrix)),
                    sort(1:length(rownames(tmp_binaryMatrix))%%cores))
  subCounts <- vector("list", length = cores)
  for(i in 1:cores){
    subCounts[[i]] <- as.matrix(tmp_binaryMatrix[geneList[[i]],])
  }
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  res <- foreach(sub=subCounts) %dopar% {
    scBDA(sub,
          tmp_samplesheet,
          contrast.by,
          ident.1,
          ident.2,
          covariates,
          method)
  }
  stopCluster(cl)
  res <- do.call(what = "rbind",res)
  return(res)

}

