#' Calculates the jaccard index between BDA and DEA based on significantly
#' detected genes
#'
#' @param BDA Results of differential dropout analysis performed with BDA
#' @param DEA Results of differential expression analysis performed with the
#' FindMarkers function of Seurat
#' @param signThreshold Threshold at which genes are considered significant
#' @param p.adjust.method Multiple testing corrextion method see `p.adjust()`
#' @importFrom stats p.adjust
#' @export
#' @return jaccard index
compareJaccard <- function(BDA, DEA, signThreshold, p.adjust.method){
  commonlyTested <- intersect(rownames(BDA), rownames(DEA))
  BDA <- BDA[commonlyTested,]
  DEA <- DEA[commonlyTested,]
  BDA$padjust <- p.adjust(BDA$P, method = p.adjust.method)
  DEA$padjust <- p.adjust(DEA$p_val, method = p.adjust.method)
  BDA_sign <- BDA[BDA$padjust <= signThreshold,]
  DEA_sign <- DEA[DEA$padjust <= signThreshold,]
  int <- length(intersect(rownames(BDA_sign), rownames(DEA_sign)))
  union <- length(unique(c(rownames(BDA_sign), rownames(DEA_sign))))
  return(int / union)
}
