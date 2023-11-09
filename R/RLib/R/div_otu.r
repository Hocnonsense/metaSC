###
#* @Date: 2022-02-27 16:52:29
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-11-08 15:52:03
#' @FilePath: /meta-snakemake-minimal/src/libs/metaSC/R/RLib/R/div.otu.r
#* @Description:
###

suppressMessages(library(data.table))
suppressMessages(library(vegan))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))


#' @title get most abundant otu
#'
#' @param div_otu table of diversity
#'
#'                colnames(div_otu) -> sample
#'
#'                rownames(div_otu) -> otu
#'
#' @param topi if an otu is the topi of one of all samples, return it
#' @param min_value,min_sample
#'      filter before selecting topi.
#'      By default, otus with total existance < 10
#'      or found in single sample will be filtered.
#' @return name of otu (from rownames)
#' @export
otu_div_top <- function(div_otu,
                        topi = 10,
                        min_value = 10,
                        min_sample = 2) {
  # filter zero samples
  div_otu_1 <- div_otu[apply(div_otu, 1, function(x) {
    (sum(x) >= min_value & sum(x > 0) >= min_sample)
  }), ]

  # formatting topi
  if (topi < 1 || topi > length(rownames(div_otu_1))) {
    topi <- length(rownames(div_otu_1))
  }

  div_top <- unique(as.vector(apply(div_otu_1, 2, function(x) {
    names(x[order(x, decreasing = TRUE)][1:topi])
  })))
  return(div_top)
}
otu.div.top <- function(div.otu, topi = 10, min.value = 10, min.sample = 2) { # nolint
  otu_div_top(div.otu, topi, min.value, min.sample)
}


.default_beta_div_color <- c(
  "#3C5488B2", "#00A087B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2",
  "#DC0000B2", "#7E6148B2", "yellow", "darkolivegreen1",
  "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick",
  "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon",
  "darkgoldenrod1", "darkseagreen", "darkorchid"
)
#' @title show diversity of otu table
#'
#' @param div.otu table of diversity
#'
#'                colnames(div.otu) -> sample
#'
#'                rownames(div.otu) -> otu
#'
#' @param pname description of table. or name of div.otu
#' @param method: method of Dimension Reduction
#'                choices: pcoa, nmds
#' @param dist: method of distance between samples
#'              choices: see vegdist. default: jaccard
#' @param binary: if jaccard: set to TRUE, otherwise FALSE
#'                see vegdist
#' @param area: method to note different Location
#'
#'              Location: recognized automatically by the first string
#'                        seperated by "_" of sample (and sample will be cut)
#'                        e.g. the column LOC_SAMPLE1 will colored by Location
#'                             LOC and noted as sample SAMPLE1
#' @param draw_labels: draw the labels by repel
#'
#' @return ggplot
#' @export
plot_beta_div <- function(div_otu,
                          pname = NA,
                          method = c("pcoa", "nmds"),
                          dist = "jaccard", binary = NA,
                          area = c("none", "ellipse", "polygon"),
                          draw_labels = TRUE) {
  # >>->> argParse
  method <- match.arg(method)
  if (is.na(binary)) {
    binary <- tryCatch(
      expr = {
        match.arg(dist, c("jaccard")) == "jaccard"
      },
      error = function(e) FALSE
    )
  }
  area <- match.arg(area)
  if (is.na(pname)) {
    pname <- deparse(substitute(div_otu, ))
  }
  # <<-<<                                                                 <<-<<
  div_otu_t <- t(div_otu)
  title_test_method <- paste0(
    method, " plot of ",
    ifelse(binary, "binary ", ""), dist, " distance"
  )

  # >>->> adonis2
  group <- data.frame(Location = sapply(
    rownames(div_otu_t),
    function(x) {
      unlist(strsplit(x, "\\_"))[1]
    }
  ))
  group_adonis2 <- adonis2(formula("div_otu_t ~ Location"), group,
    method = dist, binary = binary, by = "margin"
  )
  # <<-<<                                                                 <<-<<
  title_adonis_sgnf <- paste0(
    "ADONIS",
    " R^2=", round(group_adonis2$F[1], 4),
    " p(Pr(>F))=", group_adonis2$`Pr(>F)`[1]
  )

  # >>->> Dimensionality reduction
  if (method == "pcoa") {
    pcoa_16s <-
      div_otu_t %>% # nolint: object_usage_linter.
      vegdist(method = dist, binary = binary) %>%
      ape::pcoa()
    div_otu_point_ <- data.frame(pcoa_16s$vectors[, 1:2])
    xylab <- paste0(
      "PCo", 1:2,
      " [", round(pcoa_16s$values$Relative_eig[1:2] * 100, 2), "%]"
    )
  } else {
    nmds_dis <-
      div_otu_t %>% # nolint: object_usage_linter.
      vegdist(method = dist, binary = binary) %>%
      metaMDS(trace = 0)
    if (nmds_dis$stress >= 0.2) {
      warning("应力函数值 >= 0.2, 不合理")
      pname <- paste0(
        pname, ", stress=", as.character(round(nmds_dis$stress, 6))
      )
    }
    div_otu_point_ <- data.frame(nmds_dis$points)
    xylab <- paste0("NMDS ", 1:2)
  }
  # <<-<<                                                                 <<-<<
  div_otu_point <-
    div_otu_point_ %>% # nolint: object_usage_linter.
    `names<-`(c("Axis.1", "Axis.2")) %>%
    mutate( # nolint: object_usage_linter.
      Location = sapply(
        rownames(.), # nolint: object_usage_linter.
        . %>%
          strsplit("\\_") %>%
          unlist() %>%
          .[1]
      ), label = sapply(
        rownames(.),
        . %>%
          strsplit("\\_") %>%
          unlist() %>%
          .[-1] %>%
          paste(collapse = "_")
      )
    )


  p <- ggplot(data = div_otu_point) +
    geom_point(
      aes_string(
        x = "Axis.1", y = "Axis.2",
        color = "Location"
      ),
      size = 2, alpha = 0.65
    ) +
    scale_color_manual(
      values =
        .default_beta_div_color[seq_along(unique(div_otu_point$Location))]
    ) +
    scale_fill_manual(
      values =
        .default_beta_div_color[seq_along(unique(div_otu_point$Location))]
    ) +
    labs(
      title = paste(title_test_method, title_adonis_sgnf, pname,
        sep = "\n"
      ),
      x = xylab[1], y = xylab[2]
    ) +
    theme(
      panel.grid.major = element_line(color = "gray", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent"),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank(),
      legend.key = element_blank()
    )

  # >>->> add plugins
  if (draw_labels) {
    p <- p +
      geom_text_repel(aes_string(x = "Axis.1", y = "Axis.2", label = "label"))
  }

  if (area == "none") {
    p <- p
  } else if (area == "ellipse") {
    p <- p +
      stat_ellipse(
        aes_string(
          x = "Axis.1", y = "Axis.2",
          fill = "Location"
        ),
        type = "norm", geom = "polygon",
        alpha = 0.15, level = 0.95,
        linetype = "dashed", size = 3
      )
  } else if (area == "polygon") {
    p <- p +
      geom_polygon(
        data = . %>% # nolint
          {
            Reduce(
              rbind,
              lapply(
                split(., .$Location), # nolint
                function(x) x[chull(x[c("Axis.1", "Axis.2")]), ]
              )
            )
          },
        aes_string(
          x = "Axis.1", y = "Axis.2",
          fill = "Location", color = "Location"
        ),
        alpha = 0.15, linetype = 3
      )
  }
  # <<-<<                                                                 <<-<<

  return(p)
}

plot.beta.div <- function(div.otu, # nolint
                          pname = NA,
                          method = c("pcoa", "nmds"),
                          dist = "jaccard", binary = NA,
                          area = c("none", "ellipse", "polygon"),
                          draw_labels = TRUE) {
  plot_beta_div(div.otu, pname, method, dist, binary, area, draw_labels)
}
