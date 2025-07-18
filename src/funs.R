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

# Plot internal validitiy measures
plot_measure <- function(data) {
  p <-
    ggplot(data, aes(x = algorithm, y = percentile_50, color = sampling)) +
    geom_pointrange(aes(ymin = percentile_5, ymax = percentile_95),
                    position = position_dodge(width = 0.4)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold")) +
    xlab("Algorithm")
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
