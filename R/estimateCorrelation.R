#' Calculates the correlation coefficient between logOR and logFC of a
#' BDA and DEA respectively
#'
#' @param BDA Results of differential dropout analysis performed with BDA
#' @param DEA Results of differential expression analysis performed with the
#' FindMarkers function of Seurat
#' @param method a character string indicating which correlation coefficient
#' is to be computed. See `cor()`
#' @importFrom stats cor.test
#' @export
#' @return jaccard index
estimateCorrelation <- function(BDA, DEA, method  = "pearson"){
  commonlyTested <- intersect(rownames(BDA), rownames(DEA))
  BDA <- BDA[commonlyTested,]
  DEA <- DEA[commonlyTested,]
  cor.test(BDA[,"Z"], DEA[,"avg_log2FC"],method = method)
}
#' Calculates the correlation coefficient between logOR and logFC of a
#' BDA and DEA respectively
#'
#' @param BDA Results of differential dropout analysis performed with BDA
#' @param DEA Results of differential expression analysis performed with the
#' FindMarkers function of Seurat
#' @param method a character string indicating which correlation coefficient
#' is to be computed. See `cor()`
#' @importFrom stats cor.test
#' @export
#' @return jaccard index
estimateCorrelation <- function(BDA, DEA, method  = "pearson"){
  commonlyTested <- intersect(rownames(BDA), rownames(DEA))
  BDA <- BDA[commonlyTested,]
  DEA <- DEA[commonlyTested,]
  cor.test(BDA[,"Z"], DEA[,"avg_log2FC"],method = method)
}
