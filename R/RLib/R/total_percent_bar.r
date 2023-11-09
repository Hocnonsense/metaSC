###
#' @Date: 2022-06-28 15:37:09
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-11-08 15:08:19
#' @FilePath: /meta-snakemake-minimal/src/libs/metaSC/R/RLib/R/total_percent_bar.r
#' @Description:
###

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


#' @title get total abundance of div.otu
#'
#' @param div_raw a long table of diversity, like:
#'
#'                colnames: <sample.name>, <value.name>, <*keey.cols>
#'                values  : <any>        , <numeric>   , <*any>
#'                example : Sample       , Reads       , Group, Layer, ...
#'
#'                warning: these clumns name will be overwritten:
#'                    ..Sample, ..All, Total
#'
#' @return a smaller row of total diversity, like:
#'
#'                colnames: (<index>), <sample.name>, <value.name>, <*keey.cols>
#'                values  : (<any>  ), <any>        , <numeric>   , <*any>
#'                example : (Sample ), Sample       , Total       , Group, Layer
#'
get_total_ce <-
  function(div_raw, value_name, sample_name = "Sample", keep_cols = c()) {
    keep_div_otu <-
      div_raw %>%
      mutate(..Sample = get(sample_name)) %>%
      .[c("..Sample", sample_name, keep_cols)] %>% # nolint: object_usage_linter
      unique()

    div_raw %>%
      group_by(..Sample = get(sample_name)) %>%
      dplyr::summarise(Total = sum(get(value_name))) %>%
      merge(keep_div_otu) %>%
      column_to_rownames("..Sample") # nolint: object_usage_linter
  }


get_total_plot <- function(
    div_raw, value_name,
    sample_name = "Sample", fill = "Group", keep_cols = c(),
    labs_x = "", labs_y = "",
    font_size_1 = 15, font_size_2 = 13, font_size_3 = 10,
    axis_ticks_length = 0.1) {
  total_ce <- get_total_ce(
    div_raw,
    value_name = value_name, sample_name = sample_name,
    keep_cols = unique(c(fill, keep_cols))
  )

  p_total_mapped_rate <-
    ggplot(
      data = total_ce,
      mapping = aes_string(x = sample_name, y = "Total", fill = fill)
    ) +
    geom_bar(
      stat = "identity",
      col = "black", size = 0.3
    ) +
    scale_x_discrete(limits = rownames(total_ce), labels = NULL) +
    labs(title = "", x = labs_x, y = labs_y) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.text = element_text(
        size = font_size_2, colour = "black", face = "bold"
      ),
      axis.title = element_text(
        size = font_size_1, face = "bold", colour = "black"
      ),
      axis.ticks.length = unit(0, "cm"),
      legend.title = element_text(size = font_size_2, face = "bold"),
      legend.key = element_blank()
    )

  return(p_total_mapped_rate)
}


get_relative_ce <- function(
    div_raw, value_name, sample_name = "Sample",
    total_ce = NULL) {
  if (is.null(total_ce)) {
    total_ce <- get_total_ce(div_raw, value_name, sample_name = sample_name)
  }
  div_raw %>%
    as.data.frame() %>%
    {
      .[, value_name] / total_ce[as.character(.[, sample_name]), "Total"] * 100 # nolint: object_usage_linter, line_length_linter.
    }
}

#' @title draw a picture of relative abundance
#'
#' @param div_raw a long table of diversity, like:
#'
#'                colnames: <sample_name>, <value_name>, <fill_name>, <*keey.cols>      # nolint: line_length_linter.
#'                values  : <any>        , <numeric>   , <str|taxon>, <*any>            # nolint: line_length_linter.
#'                example : Sample       , Reads       , Taxonomy   , Group, Layer, ... # nolint: line_length_linter.
#'
#' @param value_name abundance column name
#' @param fill_name name to group by, and to be colored
#' @param sample_name name of sample to define row names
#' @return plot of relative abundance in each environment
get_percent_plot <- function(
    div_raw, value_name, fill_name, sample_name = "Sample",
    labs_x = "", labs_y = "",
    TOP_N_TAXON_PER_LAYER = NULL, MIN_REPORT_TAXON_PERCENT = 0, # nolint: object_name_linter, line_length_linter.
    font_size_1 = 15, font_size_2 = 13, font_size_3 = 10,
    axis_ticks_length = 0.1,
    geom_bar_color = "#454545") {
  total_ce <- get_total_ce(div_raw, value_name, sample_name)

  ce <-
    div_raw %>%
    as.data.frame() %>%
    mutate(
      Abundance = get_relative_ce(., value_name, sample_name, total_ce),
      name = get(fill_name)
    ) %>%
    {
      top_taxon <-
        if (!is.null(TOP_N_TAXON_PER_LAYER)) {
          reshape2::acast(
            ., formula(paste0("name ~ ", sample_name)),
            value.var = "Abundance", fill = 0, fun.aggregate = sum
          ) %>%
            otu_div_top(TOP_N_TAXON_PER_LAYER, min_value = 0) %>% # nolint: object_usage_linter, line_length_linter.
            sort()
        } else if (MIN_REPORT_TAXON_PERCENT > 0) {
          split(.[c("Abundance", "name")], .[sample_name]) %>%
            lapply(
              . %>%
                group_by(name = get("name")) %>%
                summarise(Abundance = sum(get("Abundance")))
            ) %>%
            bind_rows(.id = sample_name) %>%
            {
              .$name[.$Abundance >= MIN_REPORT_TAXON_PERCENT]
            } %>%
            unique() %>%
            as.character() %>%
            sort()
        } else {
          NULL
        }
      if (!is.null(top_taxon)) {
        .$name <- .$name %>%
          as.character() %>%
          {
            ifelse(. %in% top_taxon, ., "others")
          } %>%
          factor(levels = c(top_taxon, "others"))
      }
      .
    }

  p_relative_abundance <-
    ggplot(
      data = ce,
      mapping = aes_string(
        x = sample_name, y = "Abundance",
        fill = "name"
      )
    ) +
    geom_bar(
      position = position_stack(reverse = TRUE),
      stat = "identity",
      color = geom_bar_color, size = 0.3
    ) +
    paletteer::scale_fill_paletteer_d("ggsci::default_igv") +
    scale_x_discrete(limits = rownames(total_ce)) +
    # guides(fill = "none") +
    guides(fill = guide_legend(title = fill_name, ncol = 1, reverse = TRUE)) +
    # theme(axis.text.x = element_text(colour = axis.x.color)) +
    labs(title = "", x = labs_x, y = labs_y)

  p1 <-
    p_relative_abundance +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
    ) +
    theme(
      axis.line = element_line(colour = "black"),
      axis.text = element_text(
        size = font_size_2, colour = "black", face = "bold"
      ),
      axis.title = element_text(
        size = font_size_1, face = "bold", colour = "black"
      ),
      axis.ticks.length = unit(axis_ticks_length, "cm"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    theme(
      legend.title = element_text(size = font_size_2, face = "bold"),
      legend.text = element_text(size = font_size_3),
      legend.key = element_blank()
    ) +
    theme(text = element_text(
      family = "Arial",
      size = font_size_1,
      hjust = 0.5,
      lineheight = 0.5
    )) +
    theme(plot.title = element_text(hjust = 0.5))

  return(p1)
}
