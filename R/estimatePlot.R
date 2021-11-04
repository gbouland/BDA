#' Generates a ggplot of the logOR versus the logFC of a BDA and DEA respectively
#'
#' @param BDA Results of differential dropout analysis performed with BDA
#' @param DEA Results of differential expression analysis performed with the
#' FindMarkers function of Seurat
#' @import ggplot2
#' @export
#' @return ggplot
estimatePlot <- function(BDA, DEA){
  commonlyTested <- intersect(rownames(BDA), rownames(DEA))
  BDA <- BDA[commonlyTested,]
  DEA <- DEA[commonlyTested,]
  plotData <- data.frame("logOR" = BDA[,"Estimate"],
                         "logFC" = DEA[,"avg_log2FC"])
  max_est <- max(abs(c(BDA[,"Estimate"],DEA[,"avg_log2FC"])))
  max_est <- ceiling(max_est)
  ggplot(plotData,aes(logOR,logFC)) + geom_point() + geom_smooth(method = "lm") +
    theme_minimal() + xlim(-max_est,max_est) + ylim(-max_est,max_est) +
    geom_hline(yintercept = 0,linetype = "dotdash") +
    geom_vline(xintercept = 0,linetype = "dotdash")
}
