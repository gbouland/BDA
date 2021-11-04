#' Performs Binary Differential analysis on a seurat object
#'
#' @param object  A seurat object
#' @param contrast.by Column in which the contrasting cell populations are
#' defined.
#' @param ident.1 First cell population of contrasting cell populations.
#' @param ident.2 Second cell population of contrasting cell populations.
#' @param group.by Column in which groups are defined. For each
#' group Binary Differential analysis will be performed between the contrasting
#' populations.
#' @param cell.populations Groups of cells for which the differential dropout
#' analysis will be performed separately.
#' @param min.pct Only tested genes that were expressed in at least min.pct% of
#' the cells in either of the respective groups of interest.
#' @param covariates Covariates
#' @param p.adjust.method Multiple testing correction method.
#' @param n.cores Number of cores. When multicore is enabled, progress bar
#' cannot be displayed.
#' @param method which BDA test to use; Logistic regression (LR),
#'  chi squared test (chisq) Fisher exact (fisher) or  Pearson's Correlation
#' (phi). Of Note only LR allows to adjust for covariates.
#' @export
#' @import SeuratObject
#' @import Matrix
#' @importFrom stats na.omit p.adjust
#' @examples \dontrun{
#' library(BDA)
#' results <- BDA(object = data,
#' contrast.by = "Status", ident.1 = "COVID", ident.2 = "Healthy",
#' group.by = "cell.type.coarse", cell.populations = c("RBC","B"),
#' p.adjust.method = "fdr")}
#' @return List of differential dropout analysis results.

BDA <- function(object, contrast.by, ident.1, ident.2, group.by = NULL,
                cell.populations = NULL, min.pct = 0.1, covariates = NULL,
                p.adjust.method, n.cores = 1, method = "LR") {

  ## Param checks ##

  if(class(object)[1] != "Seurat"){
    stop("Data object is not a Seurat object, please provide Seurat object.")
  }

  if(!contrast.by %in% colnames(object@meta.data)){
    stop(sprintf("%s not present in seurat object, please provide a name that
                  is a valid meta.data column name.",contrast.by))
  }

  if(!ident.1 %in% unique(object[[contrast.by]][,1])){
    stop(sprintf("%s is not a defined cell population of %s",ident.1,contrast.by))
  }

  if(!ident.2 %in% unique(object[[contrast.by]][,1])){
    stop(sprintf("%s is not a defined cell population of %s",ident.2,contrast.by))
  }

  if(!is.null(group.by)){
    if(!group.by %in% colnames(object@meta.data)){
      stop(sprintf("%s not present in seurat object, please provide a name that
                   is a valid meta.data column name.",group.by))
    }
  }

  if(!is.null(cell.populations)){
    if(!all(cell.populations %in% unique(object[[group.by]][,1]))){
      notIn <- cell.populations[!cell.populations %in% unique(object[[group.by]][,1])]
      notIn <- paste0(notIn, collapse = ", ")
      stop(sprintf("%s is/are not defined in %s", notIn, group.by))
    }
  }

  if(!is.null(covariates)){
    if(!all(covariates %in% colnames(object@meta.data))){
      covs_not_in <- covariates[!covariates %in% colnames(object@meta.data)]
      covs_not_in <- paste0(covs_not_in,collapse = ", ")
      stop(sprintf("%s not present in seurat object, please provide a name that
                   is a valid meta.data column name.",covs_not_in))
    }
  }
  if(!p.adjust.method %in% c("holm", "hochberg", "hommel",
                              "bonferroni", "BH", "BY","fdr")){
    stop(sprintf("%s is not a valid adjustment method", p.adjust.method))
  }
  if(!method %in% c("LR","chisq","phi","fisher")){
    stop(sprintf("%s is not a valid BDA method", method))
  }
  if(method != "LR" & !is.null(covariates)){
    stop(sprintf("You added covariates while selecting %s as BDA method.
                 Only LR can handle covariates. Change mehtod to LR or
                 remove covariates.",method))
  }

  ## Prepare data ##
  message("\nSubsetting data...")
  samplesheet <- object@meta.data
  if(is.null(group.by)){
    samplesheet$group <- sprintf("%s_vs_%s",ident.1,ident.2)
    group.by <- "group"
    cell.populations <- sprintf("%s_vs_%s",ident.1,ident.2)
  }
  samplesheet <- samplesheet[samplesheet[,contrast.by] %in% c(ident.1,ident.2) &
                               samplesheet[,group.by] %in% cell.populations,]
  countMatrix <- object@assays$RNA[,rownames(samplesheet)]

  message("\nBinarizing count matrix...")
  binaryMatrix <- (countMatrix >= 1)*1

  ## Run tests ##
  res <- lapply(cell.populations,function(population){
    message(sprintf("\n%s",population))
    tmp_samplesheet <- samplesheet[samplesheet[,group.by] == population,]
    tmp_binaryMatrix <- binaryMatrix[,rownames(tmp_samplesheet)]

    tmp_binaryMatrix <- tmp_binaryMatrix[rowSums(tmp_binaryMatrix) > 0, ]
    tmp_binaryMatrix <- tmp_binaryMatrix[rowSums(tmp_binaryMatrix) < ncol(tmp_binaryMatrix), ]


    cols.ident.1 <- rownames(tmp_samplesheet[tmp_samplesheet[,contrast.by] == ident.1,])
    cols.ident.2 <- rownames(tmp_samplesheet[tmp_samplesheet[,contrast.by] == ident.2,])
    pct.expr <- data.frame(
      "A" = rowSums(tmp_binaryMatrix[,cols.ident.1]) / length(cols.ident.1),
      "B" = rowSums(tmp_binaryMatrix[,cols.ident.2]) / length(cols.ident.2))
    pct.expr <- pct.expr[pct.expr$A >= min.pct | pct.expr$B >= min.pct, ]
    tmp_binaryMatrix <- tmp_binaryMatrix[rownames(pct.expr),]

    if(n.cores == 1){
      popu.results <- scBDA(tmp_binaryMatrix,
                                  tmp_samplesheet,
                                  contrast.by,
                                  ident.1,
                                  ident.2,
                                  covariates,
                                  method)
    }else{
      popu.results <- mcBDA(tmp_binaryMatrix,
                                  tmp_samplesheet,
                                  contrast.by,
                                  ident.1,
                                  ident.2,
                                  covariates,
                                  n.cores,
                                  method)

    }
    popu.results <- data.frame(popu.results)
    rownames(popu.results) <- make.unique(rownames(tmp_binaryMatrix))


    colnames(popu.results) <- switch(method,
                               "LR" = c("Estimate","Std.Error","Z","P"),
                               "fisher" = c("OR","CI_Low","CI_high","P"),
                               "chisq" = c("OR","P"),
                               "phi" = c("Phi","t","P"))
    popu.results <- na.omit(popu.results)
    popu.results$padjust <- p.adjust(popu.results$P,method= p.adjust.method)
    popu.results <- popu.results[order(popu.results$padjust),]
    return(popu.results)
  })

  names(res) <- cell.populations
  return(res)
}
