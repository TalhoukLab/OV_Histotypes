# Helper Functions --------------------------------------------------------

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
