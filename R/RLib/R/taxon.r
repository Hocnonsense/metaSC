###
#* @Date: 2022-02-27 13:20:26
#* @LastEditors: Hwrn
#* @LastEditTime: 2022-02-27 13:20:26
#* @FilePath: /metaSC/R/taxon.r
#* @Description:
###

taxon.as.num <- function(taxon) {
  if (class(taxon) == "character") {
    taxon = which(tolower(strsplit(taxon, "")[[1]][1]) == c("d", "p", "c", "o",
                                                            "f", "g", "s"))
  }
  return(taxon)
}


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
