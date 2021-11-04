# Binary Differential Analysis (BDA)
Gerard Bouland, Ahmed Mahfouz, Marcel Reinders

4 November, 2021

## 1.0 Installation

Depends:
  R (>= 4.0.5), SeuratObject(>= 4.0.0)

``` r
install.packages("devtools")
devtools::install_github("gbouland/BDA")
```
## 2.0 Usage
`BDA()` performs a Binary Differential Analysis with a Seurat object.
``` r
library(BDA)
library(SeuratObject)
library(Seurat)
```

```r
data("dataAD")
object <- CreateSeuratObject(dataAD$counts)
object[["celltype"]] <- dataAD$meta.data$celltype
object[["status"]] <- dataAD$meta.data$status
```
```r
results <- BDA(object = object,
               contrast.by = "status", ident.1 = "ct", ident.2 = "AD",
               group.by = "celltype", cell.populations = c("neuron","mg","OPC"),
               p.adjust.method = "fdr", min.pct = 0.1, n.cores = 4, method = "LR")
```
