###
#* @Date: 2022-02-27 13:20:26
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-11-09 17:21:17
#' @FilePath: /metaSC/R/RLib/R/taxon.r
#* @Description:
###

library(dplyr)
library(stringr)


taxon_levels <- factor(
  c("domain", "phylum", "class", "order", "family", "genus", "species"),
  levels = c("domain", "phylum", "class", "order", "family", "genus", "species")
)
taxon.levels <- taxon_levels # nolint


#' @title change the string of taxon level to number
#'
#' @param taxon
#'        if int: do nothing and return
#'        if char: match it to taxon.levels
#' @return domain: 1, phylum: 2, class: 3, ..., species: 7
#' @export
#'
#' @examples
#' taxon.as.num(2)
#' # [1] 2
#' taxon.as.num("d__")
#' # [1] 1
taxon_as_num <- function(taxon) {
  if (class(taxon) == "character") {
    taxon <- which(tolower(strsplit(taxon, "")[[1]][1]) == c(
      "d", "p", "c", "o",
      "f", "g", "s"
    ))
  }
  return(taxon)
}
taxon.as.num <- taxon_as_num # nolint

#' @title split taxon
#'
#' @param taxon.full the full string of taxon, splited by ";"
#' @param start,end the index to cut taxon,
#'                  if is string, will change to int by taxon.as.num
#'                  if end is not given, will return the specified level
#'                    as end = start
#' @return splited taxon prefix (or string), already changed to factor
#'          order by taxon.full
#' @export
#'
#' @examples
#' taxon.split("Bacteria;Bacteroidota;Bacteroidia;Flavobacteriales;Flavobacteriaceae;Zunongwangia", 3) # nolint
#' # Bacteria;Bacteroidota;Bacteroidia;Flavobacteriales;Flavobacteriaceae;Zunongwangia                 # nolint
#' #                                                                       Bacteroidia                 # nolint
#' # Levels: Bacteroidia
taxon_split <- function(taxon_full, start, end = NULL) {
  if (is.null(end)) {
    end <- start
  }

  start <- taxon_as_num(start)
  end <- taxon_as_num(end)

  taxon_new_factor <- sapply(
    as.character(taxon_full),
    function(x) {
      if (is.na(x)) {
        return(x)
      }
      tax_order <-
        unlist(strsplit(as.character(x), ";")) %>%
        gsub("^.__(.+)$", "\\1", .) %>% # "d__Archaea", "D__Archaea" # nolint
        gsub("^\\(.__(.+)\\)$", "(\\1)", .) # "(d__Archaea)", "(D__Archaea)"
      return(paste(c(tax_order, rep("", end))[start:end], collapse = ";"))
    }
  ) %>%
    factor(levels = unique(.[order(taxon_full)]))

  return(taxon_new_factor)
}
taxon.split <- function(taxon.full, start, end = NULL) taxon_split(taxon.full, start, end) # nolint

.taxon_split_last <- function(x) {
  a <- unlist(strsplit(as.character(x), ";"))
  b <- a[!grepl("^(.__|)$", a)]
  ifelse(length(a) == length(b),
    gsub("^.__", "", a[length(a)]), b[length(b)]
  )
}
#' @title split taxon to the last string
#'
#' @param taxon.full the full string of taxon, splited by ";"
#' @param fill.na string to cover anything unannotated thins
#' @return last annotation taxon name, already changed to factor
#'          order by taxon.full
#' @export
taxon_split_last <- function(taxon_full, fill_na = "others") {
  sort.last.levels <- sapply(
    unique(sort(taxon_full)),
    .taxon_split_last
  )
  name_new_factor <-
    sapply(taxon_full, .taxon_split_last) %>%
    factor(levels = unique(c(sort.last.levels, fill_na))) %>%
    {
      .[is.na(.)] <- fill_na
      .
    }

  return(name_new_factor)
}
taxon.split.last <- function(taxon.full, fill.na = "others") taxon_split_last(taxon.full, fill_na) # nolint

#' @title find common prefix of a vector of string
#'
#' @param names a vector of string
#' @return single string of common prefix of string
common_prefix <- function(names) {
  # end of common part
  end <- 0
  if (length(unique(names)) <= 1) {
    return(names[1])
  }

  while (length(unique(substr(names, start = 1, stop = end + 1))) == 1) {
    end <- end + 1
  }

  # Returns common part of the name
  substr(names[1], 1, end)
}
