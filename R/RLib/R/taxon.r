###
#* @Date: 2022-02-27 13:20:26
#* @LastEditors: Hwrn
#* @LastEditTime: 2022-02-28 14:45:56
#* @FilePath: /metaSC/R/RLib/R/taxon.r
#* @Description:
###


taxon.levels = factor(c("domain", "phylum", "class", "order", "family", "genus", "species"),
                      levels = c("domain", "phylum", "class", "order", "family", "genus", "species"))

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
taxon.as.num <- function(taxon) {
  if (class(taxon) == "character") {
    taxon = which(tolower(strsplit(taxon, "")[[1]][1]) == c("d", "p", "c", "o",
                                                            "f", "g", "s"))
  }
  return(taxon)
}
#' @title split taxon
#'
#' @param taxon.full the full string of taxon, splited by ";"
#' @param start,end the index to cut taxon,
#'                  if is string, will change to int by taxon.as.num
#'                  if end is not given, will return the specified level as end = start
#' @return splited taxon prefix (or string), already changed to factor order by taxon.full
#' @export
#'
#' @examples
#' taxon.split("Bacteria;Bacteroidota;Bacteroidia;Flavobacteriales;Flavobacteriaceae;Zunongwangia", 3)
#' # Bacteria;Bacteroidota;Bacteroidia;Flavobacteriales;Flavobacteriaceae;Zunongwangia
#' #                                                                       Bacteroidia
#' # Levels: Bacteroidia
taxon.split <- function(taxon.full, start, end = NULL) {
  if (is.null(end))
    end = start

  start = taxon.as.num(start)
  end = taxon.as.num(end)

  taxon.new = sapply(as.character(taxon.full), function(x) {
    if (is.na(x)) return(x)
    if (strsplit(x, "^.__") == x) {
      tax_order = unlist(strsplit(as.character(x), ";"))
    } else {
      tax_order = unlist(strsplit(unlist(strsplit(x, "^.__"))[2], ";.__"))
    }
    return(paste(c(tax_order, rep("", end))[start:end], collapse = ";"))
  })

  taxon.new.factor = factor(taxon.new,
                            levels = unique(taxon.new[order(taxon.full)]))
  return(taxon.new.factor)
}


.taxon.split.last <- function(x) {
  a = unlist(strsplit(x, ";"))
  b = a[!grepl("^(.__|)$", a)]
  ifelse(length(a) == length(b),
         gsub("^.__", "", a[length(a)]), b[length(b)])
}
#' @title split taxon to the last string
#'
#' @param taxon.full the full string of taxon, splited by ";"
#' @param fill.na string to cover anything unannotated thins
#' @return last annotation taxon name, already changed to factor order by taxon.full
#' @export
taxon.split.last <- function(taxon.full, fill.na = "others") {
  sort.last.levels = sapply(unique(sort(taxon.full)),
                            .taxon.split.last)
  name.new = sapply(taxon.full, .taxon.split.last)

  name.new.factor = factor(name.new, levels = unique(c(sort.last.levels, fill.na)))
  name.new.factor[is.na(name.new.factor)] = fill.na
  return(name.new.factor)
}
