###
#* @Date: 2022-02-27 13:19:01
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-11-07 22:55:55
#' @FilePath: /meta-snakemake-minimal/src/libs/metaSC/R/RLib/R/load.r
#* @Description:
###


#' @title check if a function / var is defined
#'
#' @param var a var, should NOT be a string of name
#'
#' @return True if is defined else False
#' @export
#'
#' @examples
#' # please check "load_if_not"
is_defined <- function(var) {
  env <- parent.frame()
  sym <- deparse(substitute(var))
  return(exists(sym, env))
}


#' @title never do the same work again
#'
#' @param var a var, should NOT be a string of name
#'            a function named "load__"var must exist.
#'            If var haven't be loaded to environemt, then return the result
#'            of load__"var(...)
#' @param ... vars for "load__"var
#' @param .force.reolad if True, then force run load__"var(...)
#'
#' @return if not load, load it
#' @export
#'
#' @examples
#' load__a <- function(a) {
#'   a + 1
#' }
#' print(load_if_not(a, 1))
#' print(a)
#' print(load_if_not(a, 2))
#' print(a)
#' print(load_if_not(a, 2, .force.reload = TRUE))
#' print(a)
load_if_not <- function(var, ..., .force_reload = FALSE) {
  env <- parent.frame()
  sym <- deparse(substitute(var))
  if (!exists(sym, env)) {
    message(paste0("var '", sym, "' don't exists, loading"))
  } else if (.force_reload) {
    message(paste0("var '", sym, "' reloading"))
  } else {
    return(var)
  }
  env[[sym]] <- env[[paste0("load__", sym)]](...)
  return(env[[sym]])
}


#' @title source third file when source the other file
#'
#' @param x the other file to source
#' @param ... used in source
#'
#' @return
#' @export
#' @description or just use this:
#'              source(here::here('load.r'))
source_here <- function(x, ...) {
  dir <- "."
  if (sys.nframe() > 0) {
    frame <- sys.frame(1)
    if (!is.null(frame$ofile)) {
      dir <- dirname(frame$ofile)
    }
  }
  source(file.path(dir, x), ...)
}
