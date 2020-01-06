# Load datasets
library(ggplot2)
library(magrittr)
cs3 <- readr::read_csv(here::here("data-raw/cs3.csv"))
annot <- readr::read_csv(here::here("data-raw/annot.csv"))

# Format expression data
exp0 <- annot %>%
  purrr::set_names(gsub("RCC\\.", "", names(.))) %>%
  magrittr::inset(, "geneRLF", as.character(factor(
    x = .[["geneRLF"]],
    labels = c("HL1", "HL2", "HL3", "HuRef", "CS3", "mini", "CS1", "CS2")
  )))

# Run NanoString QC on CS3
exp_CS3_all <- exp0 %>%
  dplyr::filter(geneRLF == "CS3") %>%
  magrittr::inset(, "File.Name", paste0("X", .[["File.Name"]])) %>%
  nanostringr::NanoStringQC(raw = as.data.frame(cs3), exp = ., detect = 50, sn = 170)

# Normalize CS3 to housekeeping genes
cs3.norm <- nanostringr::HKnorm(cs3)

# Pools data
pools <- exp_CS3_all %>%
  dplyr::filter(sample.type %in% c("XSITE.ORIG", "XSITE.FOR", "POOL")) %>%
  dplyr::mutate(ref.type = factor(dplyr::case_when(
    POOL1 == "Yes" ~ "POOL 1",
    POOL2 == "Yes" ~ "POOL 2",
    POOL3 == "Yes" ~ "POOL 3",
    cross.site.control == "Yes" ~ "Cross Site",
    TRUE ~ NA_character_
  ))) %>%
  dplyr::filter(ref.type != "Cross Site") %>%
  dplyr::mutate(
    nanostring.date = as.Date(nanostring.date),
    averageHK = log2(averageHK),
    nanostring.site = dplyr::case_when(
      nanostring.site == "AOC" ~ "Melbourne",
      nanostring.site == "USC" ~ "San Francisco",
      TRUE ~ "Vancouver"
    )
  )

# Pools data for each site
vanpools <- pools %>%
  dplyr::filter(nanostring.site == "Vancouver") %>%
  dplyr::pull(File.Name) %>%
  dplyr::select(cs3.norm, .)

AOCpools <- pools %>%
  dplyr::filter(nanostring.site == "Melbourne") %>%
  dplyr::pull(File.Name) %>%
  dplyr::select(cs3.norm, .)

USCpools <- pools %>%
  dplyr::filter(nanostring.site == "San Francisco") %>%
  dplyr::pull(File.Name) %>%
  dplyr::select(cs3.norm, .)

# Cross-site data
cross.site <- exp_CS3_all %>%
  dplyr::filter(sample.type %in%
                  c("Cross site", "XSITE.ORIG", "XSITE.FOR", "POOL")) %>%
  dplyr::mutate(ref.type = factor(dplyr::case_when(
    POOL1 == "Yes" ~ "POOL 1",
    POOL2 == "Yes" ~ "POOL 2",
    POOL3 == "Yes" ~ "POOL 3",
    cross.site.control == "Yes" ~ "Cross Site",
    TRUE ~ NA_character_
  ))) %>%
  dplyr::filter(ref.type == "Cross Site",
                summaryID != "TVAN20681") %>%
  dplyr::mutate(nanostring.site = factor(dplyr::case_when(
    nanostring.site == "AOC" ~ "Melbourne",
    nanostring.site == "USC" ~ "San Francisco",
    TRUE ~ "Vancouver"
  )))

# Cross-site data for each site
VanRefs.gx <- cross.site %>%
  dplyr::filter(nanostring.site == "Vancouver") %>%
  dplyr::arrange(summaryID) %>%
  dplyr::pull(File.Name) %>%
  dplyr::select(cs3.norm, .)

AOCRefs.gx <- cross.site %>%
  dplyr::filter(nanostring.site == "Melbourne") %>%
  dplyr::arrange(summaryID) %>%
  dplyr::pull(File.Name) %>%
  dplyr::select(cs3.norm, .)

USCRefs.gx <- cross.site %>%
  dplyr::filter(nanostring.site == "San Francisco") %>%
  dplyr::arrange(summaryID) %>%
  dplyr::pull(File.Name) %>%
  dplyr::select(cs3.norm, .)

# Everything is calibrated to Vancouver
AOCRefs.gx2 <-
  t(nanostringr::refMethod(t(AOCRefs.gx), t(vanpools), t(AOCpools)))
USCRefs.gx2 <-
  t(nanostringr::refMethod(t(USCRefs.gx), t(vanpools), t(USCpools)))

# Combinations of cross-site gene expression
sites <- c("USC", "AOC", "VAN")
all_xsites <- combn(sites, 2) %>%
  tibble::as_tibble() %>%
  purrr::set_names(purrr::map_chr(., paste, collapse = "_vs_"))

# Combined gene expression
all_gx <- list(USCRefs.gx2, AOCRefs.gx2, VanRefs.gx) %>%
  purrr::set_names(sites) %>%
  purrr::map(~ as.data.frame(t(.)))

# Concordance measures for all genes averaged across samples
all_metrics <- all_xsites %>%
  purrr::imap_dfr(~ {
    purrr::pmap_dfr(all_gx[.x], ~ {
      R2 <- cor(.x, .y) ^ 2
      ccc <- epiR::epi.ccc(.x, .y)
      Ca <- purrr::pluck(ccc, "C.b")
      Rc <- purrr::pluck(ccc, "rho.c", "est")
      tibble::lst(R2, Ca, Rc)
    }) %>%
      dplyr::mutate(Sites = .y)
  }) %>%
  tidyr::gather(key = "Metric", value = "Expression", -Sites)

# Plot all combinations of cross-site concordance measure histograms
p <- ggplot(all_metrics, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Cross-Site Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p)
