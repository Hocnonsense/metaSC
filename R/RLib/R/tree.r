###
#' @Date: 2023-11-09 17:21:07
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-11-09 17:21:11
#' @FilePath: /metaSC/R/RLib/R/tree.r
#' @Description:
###

library(dplyr)

#' @title load a phylogenetic tree and fold at outgroup
#'
#' @param tree_file <str> file to out_group
#' @param out_group <c[str]> an array of outgroup labels
#'
#' @return a tree that reroot at outgroup, and sort as iTol does.
load__tree <- function(tree_file, out_group) {
  tre <-
    tree_file %>%
    ape::read.tree()

  genome_a <-
    tre$tip.label %>% .[. %in% out_group]
  genome_b <-
    tre$tip.label %>% .[!. %in% genome_a]

  ape::is.monophyletic(tre, genome_a) %>%
    assertthat::assert_that(msg = "outgroup is not monophyletic")

  roots <- c(A = ape::getMRCA(tre, genome_a), B = ape::getMRCA(tre, genome_b))
  mid_root_len <-
    ape::root.phylo(tre, genome_b) %>%
    {
      .$edge.length[
        apply(
          .$edge, 2, . %>%
            {
              . %in% roots
            }
        ) %>%
          apply(1, all)
      ]
    }

  tre_reroot <-
    ape::root.phylo(ape::unroot(tre), genome_a, resolve.root = TRUE) %>%
    {
      ape::bind.tree(
        ape::drop.tip(., genome_a) %>%
          {
            .$root.edge <- mid_root_len
            .
          },
        ape::keep.tip(., genome_a) %>%
          {
            .$root.edge <- mid_root_len
            .
          },
        position = mid_root_len
      )
    }

  TreeTools::SortTree(tre_reroot)
}


#' @title sort taxonomy according to label order
#'
#' @param taxon taxonomic annotation of of any label.
#'
#'    we assume that the taxonomic annotation given from general to specific,
#'    and seperated by @{taxon_sep}
#'
#'    example: "Bacteria;Chloroflexota;Dehalococcoidia;SAR202"
#' @param label_taxon the labels ordered with @{taxon}
#' @param label_phylo (a subset of) labels sorted by some methods.
#'    such as a phylogenetic tree.
#' @param taxon_sep <str> seperator of taxon. Only defined as a str of length 1
#'
#'    example: ";"
#'
#' @return <c[str]> ordered taxon.
#'    Any taxon without label given phylogentic order will be added to the list
#'    according to its name.
#'
#'    example:
#'        1.  An order of c("A;C;Dee", "A;C;Hii", "A;B") was given.
#'        2.  An "A;C;Fgg" will be added as:
#'          c("A;C;Dee", "A;C;Hii", "A;C;Fgg", "A;B")
#' @example
#' infer_phylo_order(
#'  c("A;b", "A;C;Dee", "A;C;Fgg", "A;C;Hii"),
#'  c(1, 2, 3, 4),
#'  c(2, 4, 1)
#' )
#' ## [1] "A;C;Dee;" "A;C;Hii;" "A;C;Fgg;" "A;b;"
infer_phylo_order <-
  function(taxon, label_taxon, label_phylo, taxon_sep = ";") {
    taxon_sep_gsub <- paste0(taxon_sep, "[^", taxon_sep, "]*$")
    taxon_order <- as.character(taxon) %>% paste0(taxon_sep)

    sort_taxon_order <-
      rev(unique(taxon_order[match(label_phylo, label_taxon)]))
    for (i in taxon_order %>% .[!. %in% sort_taxon_order]) { # nolint: object_usage_linter, line_length_linter.
      prefix_match <- sapply(
        sort_taxon_order, . %>% c(i) %>% common_prefix() # nolint: object_usage_linter, line_length_linter.
      ) %>%
        gsub(taxon_sep_gsub, "", .) %>%
        stringr::str_length() %>%
        which.max()
      sort_taxon_order <- c(
        sort_taxon_order[1:prefix_match - 1],
        i,
        sort_taxon_order[prefix_match:length(sort_taxon_order)]
      )
    }

    gsub(paste0(taxon_sep, "$"), "", unique(rev(sort_taxon_order)))
  }
