###
#* @Date: 2022-04-25 00:51:52
#' @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2023-08-10 20:04:58
#' @FilePath: /metaSC/R/RLib/R/ggpoint.signif.r
#* @Description:
###
# FUNCTION FOR SIGNIFICANT TEST >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
suppressMessages(library(ggExtra))
suppressMessages(library(multcomp))
suppressMessages(library(tidyr))
suppressMessages(library(patchwork))
library(ggplot2)


#' @title test difference between `group` on `value`
#'
#' @param data: long format, include 2 cols (value, group)
#' @param value: one of colname in data, numeric value to compare
#' @param group: one of colname in data, group the value
#' @return: data.frame(
#'
#'            index = group,
#'
#'            columns = c(
#'
#'              char,  # marker
#'
#'              locat,  # max locat of the value
#'
#'              locat.mark,  # place marker shall locat
#'
#'              group
#'
#'            )
#'
#'          )
#'
#' @export
group_signif_mark <- function(
    data, value, group, signif_alpha = 0.05,
    max_locat_scale = 0.1, glht_mcp_test = "Tukey") {
  test_b <- cld(
    glht(
      aov(
        formula(paste0(value, " ~ .glht_mcp_factor")),
        cbind(data,
          .glht_mcp_factor = as.factor(data[, group])
        )
      ),
      linfct = mcp(.glht_mcp_factor = glht_mcp_test)
    ),
    alpah = signif_alpha
  )
  stat <-
    data %>%
    {
      split(.[, value], .[, group]) # nolint: object_usage_linter.
    } %>%
    lapply(
      . %>%
        {
          c("min" = min(.), "mean" = mean(.), "max" = max(.)) # nolint.
        }
    ) %>%
    bind_rows(.id = "group") # nolint: object_usage_linter.
  signif_mark <- data.frame(char = test_b$mcletters$Letters, stat)
  signif_mark["locat.mark"] <- signif_mark$max * (1 + max_locat_scale)
  signif_mark[group] <- rownames(signif_mark)
  signif_mark
}

#' @title add signif figure on a given figure
#'
#' @param p a ggplot geom_point plot
#'          g must define `mapping` in `ggplot` but not `geom_point`.
#'
#' @param igroup group and color information
#'               such as `list("env" = sample_meta_col)`
#' @param x_geom,y_geom geom plot type
#' @param static_mark if NULL, show nothing
#'
#'                    if "*", show paired wilcox test between any pairs
#'                            in the group
#'
#'                    if numeric, show character difference given by
#'                                turkey test with adjust P-value lower
#'                                than that level.
#' @param p_width,p_height width and hight of main ggplot p,
#'                         and the total figure is shown as
#'                         (p_width + 1) * (p_height + 1)
#' @return a patchwork picture
#' @export
ggpoint_signif <- function(
    p, igroup,
    x_geom = geom_boxplot, y_geom = geom_boxplot,
    static_mark = "*",
    p_width = 3, p_height = 3, p_theme = NULL,
    scale_x = scale_x_continuous, scale_y = scale_y_continuous) {
  # extract data from plot OBJECT >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
  pbuild <- ggplot_build(p)

  pdata <- pbuild$plot$data
  xname <- quo_name(pbuild$plot$mapping$x)
  yname <- quo_name(pbuild$plot$mapping$y)
  # extract data from plot OBJECT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

  # get groups AND comparisons >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
  fillname <- names(igroup)
  comparisons <- combn(names(igroup[[fillname]]), 2, FUN = list)
  # get groups AND comparisons <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

  # draw plot at side >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
  px <- ggplot(
    data = pdata, mapping = aes_string(y = fillname, fill = fillname, x = xname)
  ) +
    x_geom() +
    scale_fill_manual(values = igroup[[fillname]])
  py <- ggplot(
    data = pdata, mapping = aes_string(x = fillname, fill = fillname, y = yname)
  ) +
    y_geom() +
    scale_fill_manual(values = igroup[[fillname]])
  # add static marks such as '*' or 'a'/'b' ================================== #
  if (is.null(static_mark)) {
  } else if (static_mark == "*") {
    px <- px +
      ggsignif::geom_signif(
        comparisons = comparisons, map_signif_level = TRUE,
        step_increase = 0.1
      )
    py <- py +
      ggsignif::geom_signif(
        comparisons = comparisons, map_signif_level = TRUE,
        step_increase = 0.1
      )
  } else if (is.double(static_mark) && 0 < static_mark && static_mark < 1) {
    px <- px +
      geom_text(
        data = group_signif_mark(pdata, xname, fillname,
          signif.alpha = static_mark
        ),
        aes_string(x = "locat.mark", y = fillname, label = "char")
      )
    py <- py +
      geom_text(
        data = group_signif_mark(pdata, yname, fillname,
          signif.alpha = static_mark
        ),
        aes_string(y = "locat.mark", x = fillname, label = "char")
      )
  } else {
    warning("unknown static marker: ", static_mark)
  }
  # draw plot at side <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

  # reset limits >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
  xlimits <- c(
    pbuild$layout$panel_scales_x[[1]]$get_limits(),
    ggplot_build(px)$layout$panel_scales_x[[1]]$get_limits()
  )
  xlimits <- c(min(xlimits), max(xlimits))
  ylimits <- c(
    pbuild$layout$panel_scales_y[[1]]$get_limits(),
    ggplot_build(py)$layout$panel_scales_y[[1]]$get_limits()
  )
  ylimits <- c(min(ylimits), max(ylimits))
  # reset limits <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

  subadjust <- function(x) {
    x + labs(x = "", y = "") +
      guides(fill = "none") +
      theme_void() +
      theme(axis.text = element_blank())
  }

  # use patchwork to organize picture >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
  p %>%
    {
      if (is.null(p_theme)) . else p_theme(.) # nolint: object_usage_linter.
    } +
    subadjust(px) + scale_x(limits = xlimits) +
    subadjust(py) + scale_y(limits = ylimits) +
    plot_layout(
      design = c(
        patchwork::area(t = 2, l = 1, b = 1 + p_height, r = p_width),
        patchwork::area(t = 1, l = 1, b = 1, r = p_width),
        patchwork::area(
          t = 2, l = 1 + p_width,
          b = 1 + p_height, r = 1 + p_width
        )
      ),
      guides = "collect"
    )
  # use patchwork to organize picture <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
}

group.signif.mark <- function(data, value, group, signif.alpha = 0.05, max.locat.scale = 0.1, glht.mcp.test = "Tukey") { # nolint.
  group_signif_mark(
    data, value, group, signif.alpha, max.locat.scale, glht.mcp.test
  )
}
ggpoint.signif <- function(p, igroup, x.geom = geom_boxplot, y.geom = geom_boxplot, static.mark = "*", p.width = 3, p.height = 3, p_theme = NULL, scale_x = scale_x_continuous, scale_y = scale_y_continuous) { # nolint
  ggpoint_signif(
    p, igroup, x.geom, y.geom, static.mark,
    p.width, p.height, p_theme, scale_x, scale_y
  )
}
