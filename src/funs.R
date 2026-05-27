# Helper Functions --------------------------------------------------------

# Select samples from CodeSets and transpose to standard format:
# Rows: Samples, Columns: Genes, Rownames: Sample FileName
select_samples <- function(cs_norm, samples) {
  cs_norm %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(samples),
      names_to = "FileName",
      names_prefix = "X",
      values_to = "value"
    ) %>%
    tidyr::pivot_wider(id_cols = FileName,
                       names_from = Name,
                       values_from = value) %>%
    tibble::column_to_rownames("FileName")
}

# Join histotypes to codeset and keep or discard the common ids and average
# the duplicate samples
join_avg <- function(cs, hist, id, type = c("keep", "discard")) {
  type <- match.arg(type)
  join_fun <- switch(type,
                     keep = dplyr::semi_join,
                     discard = dplyr::anti_join)
  cs %>%
    join_fun(hist, by = id) %>%
    dplyr::group_by(!!rlang::sym(id)) %>%
    dplyr::summarize_if(is.double, mean) %>%
    dplyr::ungroup() %>%
    tibble::column_to_rownames(id)
}

# Calculate squared correlation and concordance correlation
cor_stats <- function(x, y) {
  R2 <- cor(x, y) ^ 2
  ccc <- epiR::epi.ccc(x, y)
  Ca <- purrr::pluck(ccc, "C.b")
  Rc <- purrr::pluck(ccc, "rho.c", "est")
  tibble::lst(R2, Ca, Rc)
}

# Split gene expression data by histotype
split_hist <- function(data, hist_df) {
  data %>%
    tibble::rownames_to_column("ottaID") %>%
    dplyr::inner_join(hist_df, by = "ottaID") %>%
    tibble::column_to_rownames("ottaID")%>%
    split(.$revHist) %>%
    purrr::map(dplyr::select, -"revHist")
}

# Geometric mean, implemented in yardstick format
gmean_vec <- function(truth,
                      estimate,
                      estimator = NULL,
                      na_rm = TRUE,
                      case_weights = NULL,
                      event_level = "first",
                      ...) {
  estimator <-
    yardstick::finalize_estimator(truth, estimator, metric_class = "gmean")
  yardstick::check_class_metric(truth, estimate, case_weights, estimator)

  if (na_rm) {
    result <-
      yardstick::yardstick_remove_missing(truth, estimate, case_weights)

    truth <- result$truth
    estimate <- result$estimate
    case_weights <- result$case_weights
  } else if (yardstick::yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }
  gmean_impl(truth, estimate, estimator, event_level)
}

finalize_estimator_internal.gmean <- function(metric_dispatcher, x, estimator, call) {

  validate_estimator(estimator, estimator_override = c("binary", "multiclass"))
  if (!is.null(estimator)) {
    return(estimator)
  }

  lvls <- levels(x)
  if (length(lvls) > 2) {
    "multiclass"
  } else {
    "binary"
  }
}

gmean_impl <- function(truth, estimate, estimator, event_level) {
  xtab <- table(estimate, truth)
  p <- colSums(xtab)
  tp <- diag(xtab)
  sens <- tp / p
  purrr::reduce(sens, `*`) ^ (1 / length(sens))
}

gmean <- function(data, ...) {
  UseMethod("gmean")
}

gmean <- yardstick::new_class_metric(gmean, direction = "maximize")

gmean.data.frame <- function(data,
                             truth,
                             estimate,
                             estimator = NULL,
                             na_rm = TRUE,
                             case_weights = NULL,
                             event_level = "first",
                             ...) {
  yardstick::class_metric_summarizer(
    name = "gmean",
    fn = gmean_vec,
    data = data,
    truth = !!rlang::enquo(truth),
    estimate = !!rlang::enquo(estimate),
    estimator = estimator,
    na_rm = na_rm,
    case_weights = !!rlang::enquo(case_weights),
    event_level = event_level
  )
}

# Compute one-vs-all metrics
ova_metrics <- function(x, truth, estimate, metric_set) {
  x %>%
    dplyr::mutate(
      truth_ova = purrr::map({{ truth }}, ~ {
        case_when(
          levels({{ truth }}) %in% .x ~ as.character(.x),
          is.na(.x) ~ NA_character_,
          .default = "class_0"
        ) %>%
          rlang::set_names(levels({{ truth }}))
      }),
      estimate_ova = purrr::map({{ estimate }}, ~ {
        case_when(
          levels({{ estimate }}) %in% .x ~ as.character(.x),
          is.na(.x) ~ NA_character_,
          .default = "class_0"
        ) %>%
          rlang::set_names(levels({{ estimate }}))
      })
    ) %>%
    tidyr::unnest_longer(col = c(truth_ova, estimate_ova)) %>%
    dplyr::mutate(class_group = purrr::map2_chr(truth_ova_id, estimate_ova_id, unique)) %>%
    tidyr::nest(.by = class_group) %>%
    dplyr::mutate(
      data =
        purrr::map2(data, class_group, \(x, y) dplyr::mutate(x, dplyr::across(
          dplyr::matches("ova"),
          ~ .x %>%
            forcats::fct_expand(y, "class_0") %>%
            forcats::fct_relevel("class_0", after = Inf)
        ))) %>%
        purrr::map(metric_set, truth = truth_ova, estimate = estimate_ova) %>%
        suppressWarnings()
    ) %>%
    tidyr::unnest(cols = data)
}

# Summarize mean metrics across CV folds
summarize_metrics <- function(x, metric, highlight = TRUE, digits = 3) {
  df <- x |>
    dplyr::filter(.metric == metric) |>
    dplyr::distinct(dplyr::pick(-fold_id, -.estimate)) |>
    dplyr::mutate(mean_estimate = round(mean_estimate, digits = digits))
  if (highlight) {
    df <- df |>
      dplyr::mutate(
        mean_estimate = dplyr::case_when(
          undefined == "all" ~ kableExtra::cell_spec(mean_estimate, background = "#FF0000"),
          undefined == "some" ~ kableExtra::cell_spec(mean_estimate, background = "#FFD700"),
          mean_estimate == max(mean_estimate[undefined == "none"], na.rm = TRUE) ~ kableExtra::cell_spec(mean_estimate, background = "#90ee90"),
          .default = as.character(mean_estimate)
        ),
        .by = .estimator
      )
  }
  df <- df |>
    dplyr::arrange(Subsampling, class_group, Algorithms) |>
    tidyr::pivot_wider(
      id_cols = c(Subsampling, Algorithms),
      names_from = class_group,
      values_from = mean_estimate
    )
  return(df)
}

#' Gene ranks from differential expression of a predicted class
#' using a particular method
#'
#' @param data expression data
#' @param preds prediction data frame
#' @param method prediction workflow
#' @param class class of interest to compare differential expression
#' @param candidates optional; smaller gene list of candidates to consider
#' for gene ranks
gene_rank_de <- function(data, preds, method, class, candidates = NULL) {
  condition <- preds |>
    dplyr::filter(.data$Method == method, .data$Truth == class) |>
    dplyr::mutate(dplyr::across(
      "Prediction",
      ~ forcats::fct_other(
        .data$Prediction,
        keep = class,
        other_level = paste0("non_", class)
      ),
      .names = "{.col}_{class}"
    ))

  exp <- data |>
    dplyr::rename_all(~ gsub("^X", "", .)) |>
    dplyr::select("Name", dplyr::all_of(condition[["FileName"]])) |>
    tibble::column_to_rownames("Name")

  desds <- DESeq2::DESeqDataSetFromMatrix(
    countData = exp,
    colData = condition,
    design = rlang::new_formula(NULL, rlang::sym(paste0(
      "Prediction_", class
    )))
  )
  desds <- desds[rowSums(BiocGenerics::counts(desds)) > 1, ]
  desds <- DESeq2::DESeq(desds)

  res <- desds |>
    DESeq2::results() |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "Gene") |>
    tibble::as_tibble() |>
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$log2FoldChange))) |>
    dplyr::mutate(Rank = seq_along(.data$Gene), .keep = "used") |>
    dplyr::mutate(wflow = gsub("-", "_", tolower(method)), .before = 1)

  if (!is.null(candidates)) {
    res |>
      dplyr::filter(.data$Gene %in% candidates) |>
      dplyr::mutate(Rank = seq_along(.data$Gene))
  } else {
    res
  }
}


# Plotting Functions ------------------------------------------------------

#' Plot confusion matrix from observed and predicted classes
#'
#' @param data data with columns `Truth` and `Prediction`
#' @param title plot title
#' @param cols vector of colours to use for classes
plot_conf_mat <- function(data, title, cols) {
  data |>
    yardstick::conf_mat(Truth, Prediction) |>
    ggplot2::autoplot(type = "heatmap") +
    ggplot2::scale_fill_distiller(palette = "YlOrBr", direction = 1) +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggtext::element_markdown(colour = rev(cols), margin = ggplot2::margin(r = -5)),
      axis.text.x = ggtext::element_markdown(colour = cols, margin = ggplot2::margin(t = -3))
    )
}

#' Plot ROC curve and report AUC
#'
#' @param model model object
#' @param data data used to evaluate `model`
#' @param class class labels for `data`
#' @param title plot title
#' @param cols optional vector of colours to use for classes
plot_roc_curve <- function(model, data, class, title, cols = NULL) {
  preds <- stats::predict(model, data, type = "prob") |>
    tibble::add_column(truth = factor(class))
  n_class <- length(model$fit$fit$lvl)
  prob_cols <- names(dplyr::select(preds, dplyr::matches(".pred")))
  if (n_class == 2) {
    prob_cols <- prob_cols[1]
  }
  roc_df <- preds |>
    yardstick::roc_curve(truth, dplyr::all_of(prob_cols))
  if (n_class == 4) {
    roc_df <- roc_df |>
      dplyr::mutate(.level = factor(
        .level,
        levels = c("CCOC", "ENOC", "MUOC", "LGSOC")
      ))
  } else if (n_class == 5) {
    roc_df <- roc_df |>
      dplyr::mutate(.level = factor(
        .level,
        levels = c("HGSOC", "CCOC", "ENOC", "MUOC", "LGSOC")
      ))
  }
  auc <- preds |>
    yardstick::roc_auc(truth, dplyr::all_of(prob_cols)) |>
    dplyr::pull(.estimate) |>
    scales::number(accuracy = 0.001, prefix = "AUC = ")

  p <- ggplot2::autoplot(roc_df) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      axis.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::labs(title = title)

  if (any(grepl("span style", title))) {
    p <- p +
      ggplot2::theme(plot.title = ggtext::element_markdown(face = "bold"))
  }

  if (n_class > 2) {
    if (is.null(cols)) {
      cols <- title |>
        stringr::str_split_1(pattern = " vs. ") |>
        gsub(".*color:(.*);.*", "\\1", x = _)
    }
    p +
      ggplot2::labs(subtitle = auc) +
      ggh4x::facet_wrap2(~ .level,
                         strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = cols)))
  } else {
    p +
      ggplot2::geom_label(
        x = 1,
        y = 0,
        hjust = 1,
        vjust = 0,
        label = auc,
        size = 3
      )
  }
}

#' Plot calibration curve
#'
#' @inheritParams plot_roc_curve
plot_calib_curve <- function(model, data, class, title, cols = NULL) {
  preds <- model |>
    stats::predict(data, type = "prob") |>
    tibble::add_column(truth = factor(class))

  preds |>
    probably::cal_plot_logistic(truth = truth) +
    ggh4x::facet_wrap2(~ factor(
      truth,
      levels = c("HGSOC", "CCOC", "ENOC", "MUOC", "LGSOC")
    ),
    strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = cols))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      axis.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::labs(title = title)
}

#' Get plot title to use as figure caption
#'
#' Extract the plot title from `ggplot2` or `patchwork` object and
#' use as figure caption by passing dynamic expression to `fig-cap`.
get_plot_title <- function(p = ggplot2::last_plot()) {
  if (inherits(p, "patchwork")) {
    p$patches$annotation$title
  } else {
    p$labels$title
  }
}
