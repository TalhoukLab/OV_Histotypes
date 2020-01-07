library(ggplot2)
library(magrittr)


# Data Processing ---------------------------------------------------------

# Read in raw data
annot <- readr::read_csv(here::here("data-raw/annot.csv"), col_types = readr::cols())
cs1 <- readr::read_csv(here::here("data-raw/cs1.csv"), col_types = readr::cols())
cs2 <- readr::read_csv(here::here("data-raw/cs2.csv"), col_types = readr::cols())
cs3 <- readr::read_csv(here::here("data-raw/cs3.csv"), col_types = readr::cols())
pools <- readxl::read_excel(here::here("data-raw/RNA-Pools-Source_CS1-2-3.xlsx"))
ref_pools <- readRDS(here::here("data/van_pools_cs3.rds"))

# Normalize to housekeeping genes
cs1_norm <- nanostringr::HKnorm(cs1)
cs2_norm <- nanostringr::HKnorm(cs2)
cs3_norm <- nanostringr::HKnorm(cs3)

# Filter annotations for CS 1, 2, 3
annot_cs_all <- annot %>%
  dplyr::filter(RCC.geneRLF %in% c("OvCa2103_C953", "PrOTYPE2_v2_C1645", "OTTA2014_C2822")) %>%
  dplyr::mutate(
    CodeSet = dplyr::recode_factor(
      RCC.geneRLF,
      `OvCa2103_C953` = "CS1",
      `PrOTYPE2_v2_C1645` = "CS2",
      `OTTA2014_C2822` = "CS3"
    )
  )

# Histotypes
hist <- annot_cs_all %>%
  dplyr::transmute(
    FileName = RCC.File.Name,
    CodeSet,
    revHist,
    hist_gr = ifelse(revHist == "HGSC", "HGSC", "non-HGSC")
  )

# Find summaryID common to all CS
cs_common <- annot_cs_all %>%
  dplyr::count(CodeSet, summaryID) %>%
  tidyr::spread(CodeSet, n, fill = 0) %>%
  dplyr::mutate(
    all_codesets = rowSums(.[, -1]),
    in_all_cs = dplyr::select(., 2:4) %>%
      purrr::pmap_lgl(~ purrr::every(list(..1, ..2, ..3), ~ . != 0))
  )

# Find common summaryID
common_id <- with(cs_common, summaryID[in_all_cs])

# Find common samples
common_samples <- annot_cs_all %>%
  dplyr::filter(summaryID %in% common_id) %>%
  dplyr::pull(RCC.File.Name)

# Find common genes
common_genes <- list(cs1_norm, cs2_norm, cs3_norm) %>%
  purrr::map(dplyr::filter, Code.Class == "Endogenous") %>%
  purrr::map("Name") %>%
  purrr::reduce(intersect)

# Clean data by keeping common samples and genes, add ottaID
cs1_clean <- cs1_norm %>%
  dplyr::rename_all(~ gsub("^X", "", .)) %>%
  dplyr::filter(Name %in% common_genes) %>%
  dplyr::select_if(names(.) %in% c("Name", common_samples)) %>%
  dplyr::mutate(Name = forcats::fct_inorder(Name)) %>%
  tidyr::gather(FileName, value, -Name) %>%
  dplyr::mutate(ottaID = annot_cs_all$ottaID[match(FileName, annot_cs_all$RCC.File.Name)]) %>%
  tidyr::spread(Name, value)

cs2_clean <- cs2_norm %>%
  dplyr::rename_all(~ gsub("^X", "", .)) %>%
  dplyr::filter(Name %in% common_genes) %>%
  dplyr::select_if(names(.) %in% c("Name", common_samples)) %>%
  dplyr::mutate(Name = forcats::fct_inorder(Name)) %>%
  tidyr::gather(FileName, value, -Name) %>%
  dplyr::mutate(ottaID = annot_cs_all$ottaID[match(FileName, annot_cs_all$RCC.File.Name)]) %>%
  tidyr::spread(Name, value)

cs3_clean <- cs3_norm %>%
  dplyr::rename_all(~ gsub("^X", "", .)) %>%
  dplyr::filter(Name %in% common_genes) %>%
  dplyr::select_if(names(.) %in% c("Name", common_samples)) %>%
  dplyr::mutate(Name = forcats::fct_inorder(Name)) %>%
  tidyr::gather(FileName, value, -Name) %>%
  dplyr::mutate(ottaID = annot_cs_all$ottaID[match(FileName, annot_cs_all$RCC.File.Name)]) %>%
  tidyr::spread(Name, value)
