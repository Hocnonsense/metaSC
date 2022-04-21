###
#* @Date: 2022-02-27 16:52:29
#* @LastEditors: Hwrn
#* @LastEditTime: 2022-04-21 11:55:18
#* @FilePath: /metaSC/R/RLib/R/div.otu.r
#* @Description:
###

library(data.table)


#' @title get most abundant otu
#'
#' @param div.otu table of diversity
#'
#'                colnames(div.otu) -> sample
#'
#'                rownames(div.otu) -> otu
#'
#' @param topi if an otu is the topi of one of all samples, return it
#' @param min.value,min.sample
#'      filter before selecting topi.
#'      By default, otus with total existance < 10
#'      or found in single sample will be filtered.
#' @return name of otu (from rownames)
#' @export
otu.div.top <- function(div.otu,
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


#' @title annot an otu table for geom_bar
#'
#' @param div.otu table of diversity
#'
#'                colnames(div.otu) -> sample
#'
#'                rownames(div.otu) -> otu
#'
#' @param taxon.sort a subset of rowname
#'                   usually generate by "otu.div.top(div.otu, topi)"
#' @param cutoff for a bar plot, detect the bar height to only show the name of the majority otus
#' @param convert_to_pct convert the otu table to percent format
#'                       if (1. you provide a table already transfered to percent format)
#'                           or (2. you want to keep the raw count):
#'                         set it to FALSE
#' @param do_cumsum set to FALSE to avoid time-consuming-caluclation
#' @return long format of this table with these columns:
#'
#'         1.  sample: colnames (see location)
#'
#'         2.  name: rownames
#'
#'                   WARNING: rownames shall not different by any word except [a-zA-Z0-9_], which will be changed to "."
#'
#'         3.  annot.percent: percent * 100 of otu in each sample
#'
#'         4.  location: recognized automatically by the first string
#'                       seperated by "_" of sample (and sample will be cut)
#'                       e.g. the column LOC_SAMPLE1 will result in location LOC and sample SAMPLE1
#'
#'         5.  index
#'
#'         6.  label.y: center height of each bar
#'
#'         7.  label.text: name if annot.percent > cutoff else NA
#'
#' @export
#' @examples
#' se = bar.pct.annot(div.otu, otu.div.top(div.otu, topi))
#' ggpolt(data = se) +
#'   geom_bar(
#'     mapping = aes_string(x = "sample", y = "annot.percent", fill = "name"),
#'     position = position_stack(reverse = TRUE),
#'     stat = "identiy",
#'     col = "black"
#'   )
bar.pct.annot <- function(div.otu,
                          taxon.sort,
                          cutoff = 0.5,
                          convert_to_pct = TRUE,
                          do_cumsum = TRUE) {
  taxon.sort.g = sapply(sort(taxon.sort),
                        function(x) {
                          x = gsub("^(\\d)", "X\\1",
                                   gsub("[^a-zA-Z0-9_]", ".", x))
                          return(x)
                        })
  if (convert_to_pct) {
    div.otu.pct = t(div.otu) / apply(div.otu, 2, sum) * 100
  } else {
    div.otu.pct = t(div.otu)
  }

  div.grouped = reshape2::melt(
    data.frame(div.otu.pct,
               sample = rownames(div.otu.pct),
               stringsAsFactors = FALSE),
    id.vars = "sample",
    variable.name = "name",
    value.name = "annot.percent"
  )

  div.location_ = sapply(
    rownames(div.otu.pct),
    function(x) {unlist(strsplit(x, "\\_"))[1]})
  div.location = factor(div.location_, levels = unique(div.location_))

  div.sample_ = sapply(
    rownames(div.otu.pct),
    function(x) {paste(unlist(strsplit(x, "\\_"))[-1], collapse = "_")})
  div.sample = factor(div.sample_, levels = unique(div.sample_))

  div.grouped$location = div.location[div.grouped$sample]
  div.grouped$sample = div.sample[div.grouped$sample]

  div.grouped$name = sapply(div.grouped$name,
                            function(x)
                              sort(taxon.sort)[which(x == taxon.sort.g)[1]])
  if (any(is.na(div.grouped$name))) {
    div.grouped$name[is.na(div.grouped$name)] = "others"
  }
  div.grouped = div.grouped[!is.na(div.grouped$annot.percent) &
                              div.grouped$annot.percent > 0,]
  div.grouped$index = sapply(div.grouped$name,
                             function(x)
                               which(x == sort(taxon.sort))[1])

  ce = arrange(div.grouped, location, sample, index, annot.percent)
  if (do_cumsum) {
      ce = ddply(ce,
             "sample",
             transform,
             label.y = cumsum(annot.percent) - annot.percent / 2)
  }
  ce$label.text = ifelse(ce$annot.percent >= cutoff,
                         as.character(ce$name), NA)

  return(ce)
}


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
#' @param area: method to note different location
#'
#'              location: recognized automatically by the first string
#'                        seperated by "_" of sample (and sample will be cut)
#'                        e.g. the column LOC_SAMPLE1 will colored by location LOC
#'                             and noted as sample SAMPLE1
#' @param draw_labels: draw the labels by repel
#'
#' @return ggplot
#' @export
plot.beta.div <- function(div.otu,
                          pname = NA,
                          method = c("pcoa", "nmds"),
                          dist = "jaccard", binary = NA,
                          area = c("none", "ellipse", "polygon"),
                          draw_labels = TRUE) {
  # >>->> >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> >>->> argParse
  method <- match.arg(method)
  if (is.na(binary)) {
    binary = tryCatch(expr = {match.arg(idist, c("jaccard")) == "jaccard"},
                      error = function(e) FALSE)
  }
  area <- match.arg(area)
  if (is.na(pname)) {
    pname = deparse(substitute(div.otu, ))
  }
  # <<-<<                                                                 <<-<<
  div.otu.t = t(div.otu)
  title.test.method = paste0(method, " plot of ",
                             ifelse(binary, "binary ", ""), dist, " distance")

  # >>->> >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> >>->> adonis2
  group = data.frame(location = sapply(
    rownames(div.otu.t),
    function(x) {unlist(strsplit(x, "\\_"))[1]})
  )
  group.adonis2 = adonis2(formula("div.otu.t~location"), group,
                          method = dist, binary = binary)
  # <<-<<                                                                 <<-<<
  title.adonis.sgnf = paste0("ADONIS",
                             " R^2=", round(group.adonis2$F[1], 4),
                             " p(Pr(>F))=", group.adonis2$`Pr(>F)`[1])

  # >>->> >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> >>->> Dimensionality reduction
  if (method == "pcoa") {
    pcoa.16s = pcoa(vegdist(div.otu.t, method = dist, binary = binary))
    div.otu.point = data.frame(pcoa.16s$vectors[, 1:2])
    xylab = paste0(
      "PCo", 1:2,
      " [", round(pcoa.16s$values$Relative_eig[1:2] * 100, 2), "%]")
  } else {
    nmds.dis = metaMDS(vegdist(div.otu.t, method = dist, binary = binary),
                       trace = 0)
    if (nmds.dis$stress >= 0.02) {
      warning("应力函数值 >= 0.2, 不合理")
      pname = paste0(
        pname, ", stress=", as.character(round(nmds.dis$stress, 6)))
    }
    #nmds.dis.species = wascores(nmds.dis$points, t(div.otu))
    div.otu.point = data.frame(nmds.dis$points)
    xylab = paste0("NMDS ", 1:2)
  }
  # <<-<<                                                                 <<-<<
  colnames(div.otu.point) = c("Axis.1", "Axis.2")
  div.otu.point$location = sapply(
    rownames(div.otu.point),
    function(x) {unlist(strsplit(x, "\\_"))[1]})
  div.otu.point$label = sapply(
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
               size = 2, alpha = 0.65) +
    scale_color_manual(
      values = color[1:length(unique(div.otu.point$location))]) +
    scale_fill_manual(
      values = color[1:length(unique(div.otu.point$location))]) +
    labs(title = paste(title.test.method, title.adonis.sgnf, pname,
                       sep = "\n"),
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

  # >>->> >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> >>->> add plugins
  if (draw_labels) {
    p = p +
      geom_text_repel(aes_string(x = "Axis.1", y = "Axis.2", label = "label"))
  }

  if (area == "none") {
    p = p
  } else if (area == "ellipse") {
    p = p +
      stat_ellipse(aes_string(x = "Axis.1", y = "Axis.2",
                              fill = "location"),
                   type = "norm", geom = "polygon",
                   alpha = 0.15, level = 0.95,  # show.legend = FALSE,
                   linetype = 'dashed', size = 3)
  } else if (area == "polygon") {
    p = p +
      geom_polygon(
        data = Reduce(rbind,
                      lapply(split(div.otu.point, div.otu.point$location),
                             function(x) x[chull(x[c("Axis.1", "Axis.2")]),])),
        aes_string(x = "Axis.1", y = "Axis.2",
                   fill = "location", color = "location"),
        alpha = 0.15, linetype = 3)
  }
  # <<-<<                                                                 <<-<<

  return(p)
}
