---
title: "Shiny preparation"
author: ""
date: "2019-11-11"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: yeti
    highlight: tango
    code_folding: show
    keep_md: true
editor_options: 
  chunk_output_type: console
---



# Introduction

This script prepares the data to be used for shiny app. 

# Load packages


```r
suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(limma)
  library(edgeR)
  library(reshape2)
  library(SingleCellExperiment)
  library(S4Vectors)
})
```

# Load data


```r
options(ucscChromosomeNames = FALSE)
sg <- se$sg
st <- se$st
```

# Gene models


```r
create_genemodels <- function(genemodels) {
    idx <- match(c("transcript_id", "gene_id", "exon_id"), 
                 colnames(mcols(genemodels)))
    colnames(mcols(genemodels))[idx] <- c("transcript", "gene", "exon")
    mcols(genemodels)$symbol <- mcols(genemodels)$transcript
    subset(genemodels, type == "exon")
}

if (!is.null(gtffile)) {
    genemodels <- create_genemodels(genemodels)
} else {
    genemodels <- NULL
}
```

# Vector with bigWig file names 


```r
if (!is.null(bigwigdir)) {
    bwfiles <- normalizePath(list.files(bigwigdir, pattern = "\\.bw$", 
                                        full.names = TRUE))
    names(bwfiles) <- gsub("_Aligned.sortedByCoord.out.bw", "", basename(bwfiles))
} else {
    bwfiles <- NA
}
```

# edgeR - gene-level MDS


```r
logcpms <- assay(sg, "logcpm")
mds <- limma::plotMDS(logcpms, top = 500, labels = NULL, pch = NULL,
                      cex = 1, dim.plot = c(1, 2), ndim = min(7, ncol(logcpms) - 1),
                      gene.selection = "common",
                      xlab = NULL, ylab = NULL, plot = FALSE)$cmdscale.out
colnames(mds) <- paste0("MDS", seq_len(ncol(mds)))
mds <- as.data.frame(mds) %>% tibble::rownames_to_column(var = "names") %>%
    dplyr::full_join(data.frame(colData(sg)), by = "names")
```

# SingleCellExperiment on gene level

The `rowData` of `sce_gene` includes the gene information and the result tables
from `edgeR` and `DRIMSeq`. Each result table is stored as a column, and the
column name is composed by `edgeR:` or `DRIMSeq:` and the name of the contrast
used. 

The `colData` of `sce_gene` stores the sample information, the bigWig file names
and condition information

The multidimensional scale data is stored in `reducedDims`.


```r
nam <- colData(sg)$names

## low dimensional representation
reducedData <- mds %>%
    dplyr::arrange(match(names, nam)) %>%
    as.data.frame() %>% 
    dplyr::mutate(namestmp = names) %>%
    tibble::column_to_rownames("namestmp") %>%
    dplyr::select(-one_of(colnames(colData(sg))))
reducedData <- as.matrix(reducedData)

## column data
colData(sg)$bwFiles <- bwfiles[nam]

sce_gene <- SingleCellExperiment(assays = assays(sg), 
                                 rowData = rowData(sg),
                                 colData = colData(sg),
                                 metadata = list(geneModels = genemodels),
                                 reducedDims = SimpleList(MDS = reducedData))
```

# SingleCellExperiment on transcript level

The `rowData` of `sce_tx` includes the information of genes and transcripts,
and the result table on the transcript level from `DRIMSeq`.

The `colData` of `sce_tx` stores the sample information, the bigWig file names
and condition information.


```r
nam <- colData(st)$names

## column data
colData(st)$bwFiles <- bwfiles[nam]

sce_tx <- SingleCellExperiment(assays = assays(st), 
                               rowData = rowData(st),
                               colData = colData(st),
                               metadata = list(geneModels = genemodels))
```

# Output results


```r
saveRDS(list(sce_tx = sce_tx, 
             sce_gene = sce_gene),
        file = "shiny_sce.rds")
```

# Session info

The analyses above were performed with the following package versions:


```r
sessionInfo()
```

```
## R version 3.6.0 (2019-04-26)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 9 (stretch)
## 
## Matrix products: default
## BLAS/LAPACK: /usr/lib/libopenblasp-r0.2.19.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] SingleCellExperiment_1.6.0  SummarizedExperiment_1.14.1
##  [3] DelayedArray_0.10.0         BiocParallel_1.18.1        
##  [5] matrixStats_0.55.0          Biobase_2.44.0             
##  [7] GenomicRanges_1.36.1        GenomeInfoDb_1.20.0        
##  [9] IRanges_2.18.3              S4Vectors_0.22.1           
## [11] BiocGenerics_0.30.0         reshape2_1.4.3             
## [13] edgeR_3.26.8                limma_3.40.6               
## [15] tidyr_1.0.0                 dplyr_0.8.3                
## [17] tibble_2.1.3                rmarkdown_1.16             
## [19] BiocManager_1.30.9         
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_0.2.5         locfit_1.5-9.1          
##  [3] xfun_0.10                purrr_0.3.3             
##  [5] lattice_0.20-38          vctrs_0.2.0             
##  [7] htmltools_0.4.0          rtracklayer_1.44.4      
##  [9] yaml_2.2.0               XML_3.98-1.20           
## [11] rlang_0.4.1              pillar_1.4.2            
## [13] glue_1.3.1               plyr_1.8.4              
## [15] GenomeInfoDbData_1.2.1   lifecycle_0.1.0         
## [17] stringr_1.4.0            zlibbioc_1.30.0         
## [19] Biostrings_2.52.0        evaluate_0.14           
## [21] knitr_1.25               Rcpp_1.0.2              
## [23] backports_1.1.5          XVector_0.24.0          
## [25] Rsamtools_2.0.3          digest_0.6.22           
## [27] stringi_1.4.3            grid_3.6.0              
## [29] tools_3.6.0              bitops_1.0-6            
## [31] magrittr_1.5             RCurl_1.95-4.12         
## [33] crayon_1.3.4             pkgconfig_2.0.3         
## [35] zeallot_0.1.0            Matrix_1.2-17           
## [37] assertthat_0.2.1         R6_2.4.0                
## [39] GenomicAlignments_1.20.1 compiler_3.6.0
```

```r
date()
```

```
## [1] "Mon Nov 11 11:29:15 2019"
```

