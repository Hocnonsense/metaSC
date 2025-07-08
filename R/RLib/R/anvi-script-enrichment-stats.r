###
#' @Date: 2025-07-08 10:30:04
#' @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
#' @LastEditTime: 2025-07-08 10:44:49
#' @FilePath: /metaSC/R/RLib/R/anvi-script-enrichment-stats.r
#' @Description:
###

# !/usr/bin/env Rscript
## Original Author: Amy Willis, March 2020
## Generalized by: Iva Veseli, Sept 2020
## Last updated: Aoran Hu, July 2025

suppressPackageStartupMessages(library(tidyverse))

#' @title enrichment analysis on multiple types of per-group inputs
#'
#' @param col_gm_row_gn tibble/dataframe of feature counts in
#'                      pangenome. The first column should be
#'                      the identifier of the feature (named whatever you like),
#'                      and the rest of the columns should be the
#'                      presence/absence of the feature in each genome
#'                      (named by the genome accession).
#'
#' @param gm_label named vector of group labels for each genome.
#'                 The names of the vector should be the genome accessions,
#'                 and the values should be the group labels.
#'
#' @return tibble with same number of rows as the number of features /
#'         input (col_gm_row_gn), and the following columns:
#'         - <first column name of col_gm_row_gn>: the identifier of the feature
#'         - `enrichment_score`: the Rao test statistic for the feature
#'         - `unadjusted_p_value`: the unadjusted p-value for the Rao test
#'         - `adjusted_q_value`: the adjusted q-value for the Rao test
#'         - <group_label>_<group_size>: name and size of each group,
#'           joined by an underscore. The value is the number of genomes
#'           in the group that contain the feature.
#'
#' @description outlinks:
#'  modified from: https://github.com/merenlab/anvio/blob/9f079b3952bf8add49bbf6312bf1b5fca5c255b5/sandbox/anvi-script-enrichment-stats  # nolint: line_length_linter.
#'  description: https://merenlab.org/2016/11/08/pangenomics-v2/
#'
#' @export
anvio_enrichment_stats <- (\(col_gm_row_gn, gm_label = c()) {
  # Set up functions to run
  run_test_no_spread <- function(df) {
    glm(cbind(x, N - x) ~ group, df,
      family = binomial(link = "logit")
    ) |>
      anova(test = "Rao")
  }
  get_enrichment_score <- function(anova_output) {
    anova_output$Rao[2]
  }
  get_unadjusted_p_value <- function(anova_output) {
    if (
      any(is.na(anova_output$`Pr(>Chi)`)[2] & (anova_output$Rao[2] >= 1e-3))
    ) {
      paste0(
        "anova output contains NA values. This shouldn't happen! ",
        "Please submit an issue on our GitHub and tag @adw96 :)"
      ) |>
        stop()
    }
    tidyr::replace_na(anova_output$`Pr(>Chi)`[2], 1)
  }

  col1_name <- colnames(col_gm_row_gn)[1]
  groups_size <- gm_label[colnames(col_gm_row_gn)[-1]] |>
    enframe("genome", "group") |>
    group_by(.data$group) |>
    summarise(N = n())
  group_count <- col_gm_row_gn |>
    pivot_longer(
      cols = !c(col1_name),
      names_to = "genome",
      values_to = "presence"
    ) |>
    group_by(accession = .data[[col1_name]], group = gm_label[.data$genome]) |>
    summarise(x = sum(.data$presence > 0)) |>
    ungroup() |>
    inner_join(groups_size)

  group_count_unique <- group_count |>
    pivot_wider(
      id_cols = "accession",
      names_from = c("group", "N"),
      values_from = "x",
      values_fill = 0
    ) |>
    group_by(across(!"accession")) |>
    summarise(
      accessions = list(.data$accession),
      accession = first(.data$accession),
      .groups = "drop"
    )

  pvalues_df_unique <- group_count |>
    filter(.data$accession %in% group_count_unique$accession) |>
    nest(data = c("group", "N", "x")) |>
    mutate(
      model = map(data, run_test_no_spread),
      unadjusted_p_value = map_dbl(.data$model, get_unadjusted_p_value),
      enrichment_score = map_dbl(.data$model, get_enrichment_score)
    )

  # Rao test statistic should always be positive. However,
  # it happens that Rao test statistics are numerically negative
  # but close to zero. If they are large and negative (say, smaller than
  # 0.00001), we should throw an error:
  pvalues_df_unique |>
    filter(.data$enrichment_score < -1e-5) |>
    (\(df) {
      if (nrow(df) == 0) {
        return(invisible())
      }
      message("A Rao test statistic is large and negative, oh my!")
      message("Here are the rows where the test statistic was negative:")
      df |>
        arrange(desc(enrichment_score)) |>
        print(n = Inf)
      paste0(
        "Something has gone terribly wrong, and it's Amy's fault. ",
        "Please let her know ASAP so she can fix it!"
      ) |>
        stop()
    })()

  pvalues_df <- pvalues_df_unique |>
    inner_join(group_count_unique) |>
    unnest("accessions")
  # find the maximum lambda for qvalue
  # reasons for this change documented at
  # https://github.com/merenlab/anvio/issues/1383
  lambdas <- pvalues_df$unadjusted_p_value |>
    (\(x) {
      lambdas <- seq(0.05, 0.95, 0.05)
      lambdas <- lambdas[min(x, na.rm = TRUE) < lambdas]
      lambdas <- lambdas[lambdas < max(x, na.rm = TRUE)]

      pi0_est <- \(lambdas) {
        est <- try(
          qvalue::pi0est(x, lambda = lambdas, pi0.method = "smoother")$pi0,
          silent = TRUE
        )
        !"try-error" %in% class(est)
      }
      if (pi0_est(lambdas)) {
        return(lambdas)
      }
      # Sometimes, especially when there are many of p-values near 0.5 and fewer
      # near 1 (documented at https://github.com/merenlab/anvio/issues/1828),
      # we can't reliably do the extrapolation we want. Let's instead just say
      # that features with p-values above 0.2 are likely null, and use that to
      # estimate the proportion of features that are truly null.
      lambdas <- 0.2
      if (pi0_est(lambdas)) {
        return(lambdas)
      }
      paste0(
        "Doh! We still can't estimate the proportion of features ",
        "that are truly null. It's Amy's fault, ",
        "so please let her know ASAP so she can fix it!"
      ) |>
        stop()
    })()

  # get the q-values, then format the data and save
  # If they are small and negative (checked earlier), just set them to zero.
  # Then, order the columns from largest to smallest enrichment score:
  qvalues_df_adjust <- pvalues_df |>
    mutate(
      adjusted_q_value = qvalue::qvalue(
        p = .data$unadjusted_p_value, lambda = lambdas, pi0.method = "smoother"
      )$qvalues,
      enrichment_score = pmax(0, enrichment_score)
    ) |>
    group_by(.data$accession) |>
    summarise(
      across(
        c(
          .data$enrichment_score,
          .data$unadjusted_p_value,
          .data$adjusted_q_value
        ),
        first
      ),
      .groups = "drop"
    ) |>
    arrange(desc(.data$enrichment_score))

  qvalues_adjust_report <- qvalues_df_adjust |>
    inner_join(group_count_unique) |>
    unnest(.data$accessions) |>
    mutate(accession = .data$accessions, accessions = NULL)

  names(qvalues_adjust_report)[1] <- col1_name
  qvalues_adjust_report
})


anvio_enrichment_stats.example <- (\() { # nolint: object_name_linter
  col_gm_row_gn <- tibble(
    gn = c("gn1", "gn2", "gn3"),
    set1 = c(0, 1, 1),
    set2 = c(0, 0, 1),
    set3 = c(1, 0, 0),
    set4 = c(1, 1, 0)
  )
  gm_label <- c(
    set1 = "g1", set2 = "g1", set3 = "g2", set4 = "g2"
  )
  anvio_enrichment_stats(col_gm_row_gn, gm_label)
  "# A tibble: 3 × 6"
  "  gn    enrichment_score unadjusted_p_value adjusted_q_value  g1_2  g2_2"
  "  <chr>            <dbl>              <dbl>            <dbl> <int> <int>"
  "1 gn1               4.00             0.0455           0.0683     0     2"
  "2 gn3               4.00             0.0455           0.0683     2     0"
  "3 gn2               0                1                1          1     1"
})
