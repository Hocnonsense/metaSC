###
#' @Date: 2022-02-27 13:23:46
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-11-09 17:29:23
#' @FilePath: /metaSC/R/RLib/R/RLib.r
#' @Description:
#'    source me!
###

suppressWarnings(library(tidyverse))
options(dplyr.summarise.inform = FALSE)
suppressWarnings(library(patchwork))
suppressWarnings(library(vegan))

#' library(zeallot)
#'
#' library(stringr)
#' library(plyr)
#' library(dplyr)
#' library(vegan)
#' library(ape)
#' library(multcomp)
#'
#' library(pheatmap)
#' library(ggExtra)
#' library(ggrepel)
#' library(scales)
#' library(ggpubr)
#'
#' library(DESeq2)
#' library(edgeR)
#'
#' suppressMessages(library(pathview))

if (grepl("RLib$", getwd())) {
  #### LOAD general functions                                               ####
  . <- plyr::.
  message.print <- function(...) { # nolint: object_name_linter.
    capture.output(...) %>%
      paste0(collapse = "\n") %>%
      {
        base::message(.)
      }
  }

  #### LOAD general vars                                                    ####
  argv <- commandArgs(trailingOnly = TRUE)
} else {
  make.package <- function(dir) { # nolint: object_name_linter.
    suppressWarnings(roxygen2::roxygenize(dir))
    suppressWarnings(install.packages(dir, quiet = TRUE, repos = NULL))
  }

  make.package("../../RLib")
  library(RLib)
}
