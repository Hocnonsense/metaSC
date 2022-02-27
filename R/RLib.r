###
#* @Date: 2022-02-27 13:23:46
#* @LastEditors: Hwrn
#* @LastEditTime: 2022-02-27 13:29:52
#* @FilePath: /metaSC/R/RLib.r
#* @Description:
###


library(data.table)
library(zeallot)

library(reshape2)
library(plyr)
library(dplyr)
library(vegan)
library(ape)
library(stringr)

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(scales)
library(ggpubr)

suppressMessages(library(pathview))
#install.packages("roxygen2", depend = TRUE)
library(roxygen2)

expand.source.abs <- function(x, ...) {
  dir <- "."
  if (sys.nframe() > 0) {
    frame <- sys.frame(1)
    if (!is.null(frame$ofile)) {
      dir <- dirname(frame$ofile)
    }
  }
  return(file.path(dir, x, ...))
}
abs.dir = expand.source.abs("RLib/")
roxygenize(expand.source.abs("RLib/"))
system(paste("R CMD build", abs.dir))
#system(paste("R CMD check", abs.dir))
system(paste("R CMD INSTALL", abs.dir))

library(RLib)