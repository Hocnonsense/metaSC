###
#* @Date: 2022-02-27 16:52:29
#* @LastEditors: Hwrn
#* @LastEditTime: 2022-02-28 15:32:33
#* @FilePath: /metaSC/R/RLib/R/div.otu.r
#* @Description:
###

library(data.table)


#' @title get most abundant otu
#'
#' @param div.otu table of diversity
#'                colnames(div.otu) -> sample
#'                rownames(div.otu) -> otu
#' @param topi if an otu is the topi of one of all samples, return it
#' @param min.value,min.sample
#'      filter before selecting topi.
#'      By default, otus with total existance < 10
#'      or found in single sample will be filtered.
#' @return name of otu (from rownames)
#' @export
otu.div.top <- function(div.otu,
                          # index: taxon; column: sample
                          topi = 10,
                          min.value = 10,
                          min.sample = 2) {
  # filter zero samples
  div.otu = div.otu[apply(div.otu, 1, function(x)
    (sum(x) >= min.value & sum(x > 0) >= min.sample)), ]

  # formatting topi
  if (topi < 1 || topi > length(rownames(div.otu)))
    topi = length(rownames(div.otu))

  div.top = unique(as.vector(apply(div.otu, 2, function(x)
    names(x[order(x, decreasing = T)][1:topi]))))
  return(div.top)
}


bar.pct.annot <- function(div.otu,
                          taxon.sort,
                          cutoff = 0.5,
                          convert_to_pct = TRUE) {
  taxon.sort.g = sapply(sort(taxon.sort),
                        function(x) {
                          x = gsub("^(\\d)", "X\\1", gsub("[^a-zA-Z0-9_]", ".", x))
                          return(x)
                        })
  if (convert_to_pct) {
    div.otu.pct = div.otu / apply(div.otu, 1, sum) * 100
  } else {
    div.otu.pct = div.otu
  }

  div.grouped = reshape2::melt(
    data.frame(div.otu.pct,
               sample = rownames(div.otu.pct),
               stringsAsFactors = FALSE),
    id.vars = "sample",
    variable.name = "name",
    value.name = "annot.percent"
  )

  div.grouped$location = sapply(
    div.grouped$sample,
    function(x) {unlist(strsplit(x, "\\_"))[1]})
  div.grouped$sample = sapply(
    div.grouped$sample,
    function(x) {paste(unlist(strsplit(x, "\\_"))[-1], collapse = "_")})

  div.grouped$name = sapply(div.grouped$name,
                            function(x)
                              sort(taxon.sort)[which(x == taxon.sort.g)[1]])
  div.grouped$name[is.na(div.grouped$name)] = "others"
  div.grouped = div.grouped[!is.na(div.grouped$annot.percent) &
                              div.grouped$annot.percent > 0,]
  div.grouped$index = sapply(div.grouped$name,
                             function(x)
                               which(x == sort(taxon.sort))[1])

  ce = arrange(div.grouped, location, sample, index, annot.percent)
  ce = ddply(ce,
             "sample",
             transform,
             label.y = cumsum(annot.percent) - annot.percent / 2)
  ce$label.text = ifelse(ce$annot.percent >= cutoff,
                         as.character(ce$name), NA)

  return(ce)
}


plot.beta.div <- function(div.otu,
                            pname = NA,
                            method = c("pcoa", "nmds"),
                            dist = c("bray", "jaccard"),
                            area = c("none", "ellipse", "polygon")) {
  method <- match.arg(method)
  dist <- match.arg(dist)
  area <- match.arg(area)
  if (is.na(pname)) {
    pname = deparse(substitute(div.otu, ))
  }
  if (method == "pcoa") {
    pcoa.16s = pcoa(vegdist(div.otu, method = dist))
    div.otu.point = data.frame(pcoa.16s$vectors[, 1:2])
    colnames(div.otu.point) = c("Axis.1", "Axis.2")
    xylab = paste0(
      "PCo", 1:2, " [", round(pcoa.16s$values$Relative_eig[1:2] * 100, 2), "%]")
  } else {
    nmds.dis = metaMDS(vegdist(div.otu, method = dist))
    if (nmds.dis$stress >= 0.02) {
      warning("应力函数值 >= 0.2, 不合理")
      pname = paste0(
        pname, ", stress=", as.character(round(nmds.dis$stress, 6)))
    }
    #nmds.dis.species = wascores(nmds.dis$points, t(div.otu))
    div.otu.point = data.frame(nmds.dis$points)
    colnames(div.otu.point) = c("Axis.1", "Axis.2")
    xylab = paste0("MDS", 1:2)
  }
  div.otu.point$location = sapply(
    rownames(div.otu.point),
    function(x) {unlist(strsplit(x, "\\_"))[1]})
  rownames(div.otu.point) = sapply(
    rownames(div.otu.point),
    function(x) {paste(unlist(strsplit(x, "\\_"))[-1], collapse = "_")})

  color = c("#3C5488B2", "#00A087B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2",
            "#DC0000B2", "#7E6148B2", "yellow", "darkolivegreen1",
            "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick",
            "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon",
            "darkgoldenrod1", "darkseagreen", "darkorchid")

  p = ggplot(data = div.otu.point) +
    geom_point(aes_string(x = "Axis.1", y = "Axis.2",
                          color = "location"),
               size = 2) +
    geom_text_repel(aes_string(x = "Axis.1", y = "Axis.2",
                               label = "rownames(div.otu.point)")) +
    scale_color_manual(
      values = color[1:length(unique(div.otu.point$location))]) +
    scale_fill_manual(
      values = color[1:length(unique(div.otu.point$location))]) +
    labs(title = paste(paste0(method, " plot of ", dist, " distance"),
                       pname, sep = "\n"),
         x = xylab[1], y = xylab[2]) +
    theme(
      panel.grid.major = element_line(color = 'gray', size = 0.2),
      panel.grid.minor = element_blank(),
      #panel.background = element_blank(),
      panel.background = element_rect(color = 'black', fill = 'transparent'),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    )

  if (area == "none") {
    p = p
  } else if (area == "ellipse") {
    p = p +
      stat_ellipse(aes_string(x = "Axis.1", y = "Axis.2",
                              color = "location"),
                   level = 0.95, show.legend = FALSE, size = 1)
  } else if (area == "polygon") {
    p = p +
      geom_polygon(
        data = Reduce(rbind, lapply(split(div.otu.point, div.otu.point$location),
                                    function(x) x[chull(x[c("Axis.1", "Axis.2")]),])),
        aes_string(x = "Axis.1", y = "Axis.2",
                   fill = "location", color = "location"),
        alpha = 0.1, linetype = 3)
  }

  return(p)
}
