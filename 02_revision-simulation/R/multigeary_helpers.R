snu_multigeary <- function(multigeary_res, alternative = "two.sided",
                           p.adj.type = "fdr", q.value = 0.05) {
  p.adj.type <- match.arg(p.adj.type, c("bonferroni", "sidak", "fdr", "none"))

  if (alternative == "two.sided") {
    multigeary_res[["p_vals"]] <- pmin(multigeary_res[["p_vals"]] * 2, 1)
    multigeary_res[["c_vals"]] <- dplyr::if_else(
      multigeary_res[["p_vals"]] > q.value,
      0L,
      as.integer(multigeary_res[["c_vals"]])
    )
  }

  res_df <- tibble::tibble(
    mci = multigeary_res[["lisa_vals"]],
    P.mci = multigeary_res[["p_vals"]],
    class = multigeary_res[["c_vals"]]
  )

  lab_col_df <- tibble::tibble(
    id = 0:(length(multigeary_res[["labels"]]) - 1),
    labels = multigeary_res[["labels"]],
    cols = multigeary_res[["colors"]]
  )

  m <- nrow(res_df)
  class_sim <- if (p.adj.type == "bonferroni") {
    rgeoda::lisa_clusters(multigeary_res, cutoff = q.value / m)
  } else if (p.adj.type == "sidak") {
    rgeoda::lisa_clusters(multigeary_res, cutoff = 1 - (1 - q.value)^(1 / m))
  } else if (p.adj.type == "fdr") {
    rgeoda::lisa_clusters(
      multigeary_res,
      cutoff = rgeoda::lisa_fdr(multigeary_res, current_p = q.value)
    )
  } else {
    res_df$class
  }

  res_df |>
    dplyr::left_join(lab_col_df, dplyr::join_by(class == id)) |>
    dplyr::mutate(
      class.sim = class_sim
    ) |>
    dplyr::left_join(lab_col_df, dplyr::join_by(class.sim == id)) |>
    dplyr::rename(
      labels = labels.x,
      cols = cols.x,
      labels.sim = labels.y,
      cols.sim = cols.y
    ) |>
    dplyr::mutate(
      labels.sim = forcats::fct(
        labels.sim,
        levels = c("Positive", "Negative", "Undefined", "Isolated", "Not significant")
      )
    )
}

classify_multigeary_clusters <- function(multigeary_df, alpha = 0.05) {
  if (!"P.mci.adj" %in% names(multigeary_df)) {
    multigeary_df <- multigeary_df |>
      dplyr::mutate(P.mci.adj = stats::p.adjust(P.mci, method = "fdr"))
  }

  multigeary_df |>
    dplyr::mutate(
      # ML-SAS와 비교 가능한 local similarity cluster만 남기되,
      # 유의성 통제는 별도로 계산한 FDR 보정 p-value를 사용한다.
      cluster_flag = dplyr::if_else(
        labels == "Positive" & P.mci.adj <= alpha,
        "Significant spatial cluster",
        "Not significant"
      ),
      cluster_flag = factor(
        cluster_flag,
        levels = c("Significant spatial cluster", "Not significant")
      )
    )
}
