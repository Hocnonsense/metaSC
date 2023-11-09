###
#' @Date: 2023-11-09 14:03:40
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-11-09 14:05:07
#' @FilePath: /meta-snakemake-minimal/src/libs/metaSC/R/RLib/R/vars.r
#' @Description:
###


ggscale_label_format_sci <-
  . %>%
  format(scientific = TRUE) %>%
  str_replace("e\\+0+$", "") %>%
  str_replace("e\\+0", "%*%10^") %>%
  str_replace("e\\+", "%*%10^") %>%
  str_replace("e", "%*%10^") %>%
  parse(text = .)
