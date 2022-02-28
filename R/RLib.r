###
#* @Date: 2022-02-27 13:23:46
#* @LastEditors: Hwrn
#* @LastEditTime: 2022-02-28 10:44:56
#* @FilePath: /metaSC/R/RLib.r
#* @Description:
#     source me!
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

reload.RLib <- function() {
  dir <- "."
  if (sys.nframe() > 0) {
    frame <- sys.frame(1)
    if (!is.null(frame$ofile)) {
      dir <- dirname(frame$ofile)
    }
  }
  abs.dir = file.path(dir, "RLib/")
  roxygenize(abs.dir)
  system(paste("R CMD build", abs.dir))
  #system(paste("R CMD check", abs.dir))
  system(paste("R CMD INSTALL", abs.dir))
}
reload.RLib()

library(RLib)
