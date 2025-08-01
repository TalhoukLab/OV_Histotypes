# Load Packages -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(naniar)
library(nanostringr)
library(ottaOvca)
library(knitr)
library(gtsummary)
library(here)


# Data Processing ---------------------------------------------------------

# Read in raw OTTA data
data("annot", "annotNEW", "rawOVCA2", "rawPROT", "rawOTTA",
     "rawOTTA2", "rawOTTA3", package = "otta")
cs1 <- rawOVCA2
cs2 <- rawPROT
cs3 <- rawOTTA
cs4 <- rawOTTA2
cs5 <- rawOTTA3

data("od.otta")
od.otta <- to_native_type(od.otta)

cohorts <- read_excel(here("data-raw/bc-842_avail_cosp_2020-09-15.xlsx"), sheet = "data")
pools <- read_excel(here("data-raw/RNA-Pools-Source_CS1-2-3.xlsx"))
ref_pools <- readRDS(here("data/van_pools_cs3.rds"))

# Pairwise CodeSet comparisons
codesets <- c("CS1", "CS2", "CS3")
all_codesets <- combn(codesets, 2) %>%
  as.data.frame() %>%
  set_names(map_chr(., paste, collapse = "_vs_"))

# Pairwise site comparisons
sites <- c("USC", "AOC", "VAN")
all_xsites <- combn(sites, 2) %>%
  as.data.frame() %>%
  set_names(map_chr(., paste, collapse = "_vs_"))

# Recode annotation geneRLF and rename file name and site columns
annot_all <- annot %>%
  replace_with_na(list(ottaID = c("", "N/A"))) %>%
  mutate(
    CodeSet = recode_factor(
      RCC.geneRLF,
      `OvCa2103_C953` = "CS1",
      `PrOTYPE2_v2_C1645` = "CS2",
      `OTTA2014_C2822` = "CS3"
    )
  ) %>%
  rename_with(~ gsub("^RCC\\.", "", .)) %>%
  rename(FileName = File.Name, site = nanostring.site)

# Gold standard histotypes for OVAR3 and ICON7 cohorts
hist_stan <- od.otta %>%
  filter(cohort %in% c("OVAR3", "ICON7")) %>%
  transmute(
    ottaID = as.character(otta_id),
    cohort,
    hist_rev = case_when(
      hist_rev == "high-grade serous" ~ "HGSC",
      hist_rev == "low-grade serous" ~ "LGSC",
      hist_rev == "mucinous" ~ "MUC",
      hist_rev == "endometrioid" ~ "ENOC",
      hist_rev == "clear cell" ~ "CCOC",
      hist_rev == "serous borderline tumour" ~ "SBOT",
      TRUE ~ NA_character_
    )
  )

# histology_mol_v3 histotypes (from Susan)
histology_mol_v3_df <- read_excel(
  here(
    "data-raw",
    "revHist_selected_cohorts compare AI mol his v3_editAB.xlsx"
  ),
  sheet = 1,
  na = c("", "NA")
) |>
  select(FileName, histology_mol_v3, histology_source_v3)

# IHC histotypes (includes COSP, from hist_rev_V2)
bc_1965_ihc_data <-
  read_excel(here("data-raw", "bc-1965_ihc_data_2025-06-04.xlsx"), sheet = 1) |>
  distinct(FileName, hist_rev_v2, hist_rev_v2_source)

# COSP histotypes
cosp_df <- cohorts |>
  mutate(FileName = gsub("^X", "", col_name)) |>
  select(FileName, has_cosp, hist_cosp, hist_cosp_details)

# Original reviewed histotypes
hist_all <- annot_all %>%
  left_join(hist_stan, by = "ottaID") %>%
  transmute(
    FileName,
    ottaID,
    CodeSet,
    revHist = case_when(
      !is.na(cohort) ~ hist_rev,
      is.na(cohort) & revHist == "CCC" ~ "CCOC",
      is.na(cohort) & revHist == "ENOCa" ~ "ENOC",
      revHist %in% c("", "UNK") ~ NA_character_,
      TRUE ~ revHist
    ),
    site
  )

# Add molecular-based histotypes on top of revHist
# Priority: histology_mol_v3 > hist_rev_v2 > hist_cosp
hist_all_mol <- hist_all |>
  left_join(histology_mol_v3_df, by = "FileName") |>
  left_join(bc_1965_ihc_data, by = "FileName") |>
  left_join(cosp_df, by = "FileName") |>
  mutate(
    histology_mol_v3_mapped = case_match(
      histology_mol_v3,
      "1H" ~ "HGSC",
      "1L" ~ "LGSC",
      "1" ~ "SC",
      "2" ~ "MUC",
      "3" ~ "ENOC",
      "4" ~ "CCOC",
      "6" ~ "Other specified epithelial ovarian cancer",
      "0" ~ "OTHER",
      "20" ~ "MUC LMP",
      "99" ~ "SC LMP",
      .default = histology_mol_v3
    ),
    hist_rev_v2_mapped = case_when(
      hist_rev_v2 == "clear cell" ~ "CCOC",
      hist_rev_v2 == "endometrioid" ~ "ENOC",
      hist_rev_v2 == "high-grade serous" ~ "HGSC",
      hist_rev_v2 == "low-grade serous" ~ "LGSC",
      hist_rev_v2 == "mucinous" ~ "MUC",
      hist_rev_v2 == "serous" & hist_rev_v2_source == "CB-81" ~ revHist,
      .default = hist_rev_v2
    ),
    hist_cosp_mapped = case_when(
      hist_cosp %in% c("HGSOC", "high-grade serous") ~ "HGSC",
      hist_cosp %in% c("LGSOC", "low-grade serous") ~ "LGSC",
      hist_cosp == "clear cell" ~ "CCOC",
      hist_cosp == "endometrioid" ~ "ENOC",
      hist_cosp == "mucinous" ~ "MUC",
      hist_cosp == "serous" & hist_cosp_details == "CB-81" ~ revHist,
      .default = hist_cosp
    ),
    across(
      c(
        histology_mol_v3_mapped,
        hist_rev_v2_mapped,
        hist_cosp_mapped
      ),
      ~ if_else(revHist != ., "discordant", "concordant"),
      .names = "revHist_vs_{gsub('_mapped', '', {col}, fixed = TRUE)}"
    ),
    hist_final_source = case_when(
      !is.na(histology_mol_v3_mapped) ~ "histology_mol_v3",
      !is.na(hist_rev_v2_mapped) ~ "hist_rev_v2",
      !is.na(hist_cosp_mapped) ~ "hist_cosp",
      !is.na(revHist) ~ "revHist"
    ),
    hist_final = case_match(
      hist_final_source,
      "histology_mol_v3" ~ histology_mol_v3_mapped,
      "hist_rev_v2" ~ hist_rev_v2_mapped,
      "hist_cosp" ~ hist_cosp_mapped,
      "revHist" ~ revHist
    ),
    hist_gr = if_else(hist_final == "HGSC", "HGSC", "non-HGSC"),
    hist_gr2 = hist_final |>
      factor(levels = c("HGSC", "CCOC", "ENOC", "MUC", "LGSC")) |>
      fct_na_value_to_level(level = "Other")
  ) |>
  relocate(hist_final_source, .after = hist_gr2) |>
  select(-matches("mapped"))

# Main Histotypes: "CCOC", "ENOC", "HGSC", "LGSC", "MUC"
hist_main <- hist_all_mol |>
  filter(hist_final %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC"))

# CS1/2/3 annotations
cs1_exp <- filter(annot_all, CodeSet == "CS1")
cs2_exp <- filter(annot_all, CodeSet == "CS2")
cs3_exp <- filter(annot_all, CodeSet == "CS3")

# Select specific cohorts and extract samples
cs1_samples_coh <- cohorts %>%
  filter(file_source == "cs1",
         cohort %in% c("OOU", "OOUE", "VOA", "MAYO", "MTL")) %>%
  pull(col_name)

cs2_samples_coh <- cohorts %>%
  filter(
    file_source == "cs2",
    cohort %in% c(
      "OOU",
      "OOUE",
      "VOA",
      "MAYO",
      "OVAR3",
      "ICON7",
      "JAPAN",
      "MTL",
      "POOL-CTRL"
    ),
  ) %>%
  pull(col_name)

cs3_samples_coh <- cohorts %>%
  filter(
    file_source == "cs3",
    cohort %in% c(
      "OOU",
      "OOUE",
      "VOA",
      "TNCO",
      "DOVE4",
      "POOL-1",
      "POOL-2",
      "POOL-3"
    )
  ) %>%
  pull(col_name)

# Select only samples from these cohorts in the data and expression
cs1_coh <- cs1 %>% select(Code.Class, Name, Accession, all_of(cs1_samples_coh))
cs2_coh <- cs2 %>% select(Code.Class, Name, Accession, all_of(cs2_samples_coh))
cs3_coh <- cs3 %>% select(Code.Class, Name, Accession, all_of(cs3_samples_coh))

cs1_exp_coh <- cs1_exp %>% filter(FileName %in% gsub("^X", "", cs1_samples_coh))
cs2_exp_coh <- cs2_exp %>% filter(FileName %in% gsub("^X", "", cs2_samples_coh))
cs3_exp_coh <- cs3_exp %>% filter(FileName %in% gsub("^X", "", cs3_samples_coh))

# Samples that failed QC from selected cohorts on CS1/2/3
cs1_qc <-
  NanoStringQC(raw = cs1_coh, exp = cs1_exp_coh, detect = 50, sn = 100)
cs1_qc_failed <-
  cs1_qc %>%
  filter(QCFlag == "Failed") %>%
  pull(FileName) %>%
  paste0("X", .)

cs2_qc <-
  NanoStringQC(raw = cs2_coh, exp = cs2_exp_coh, detect = 50, sn = 100)
cs2_qc_failed <-
  cs2_qc %>%
  filter(QCFlag == "Failed") %>%
  pull(FileName) %>%
  paste0("X", .)

cs3_qc <-
  NanoStringQC(raw = cs3_coh, exp = cs3_exp_coh, detect = 50, sn = 100)
cs3_qc_failed <-
  cs3_qc %>%
  filter(QCFlag == "Failed") %>%
  pull(FileName) %>%
  paste0("X", .)

cs1_samples_coh_qc <- setdiff(cs1_samples_coh, cs1_qc_failed)
cs2_samples_coh_qc <- setdiff(cs2_samples_coh, cs2_qc_failed)
cs3_samples_coh_qc <- setdiff(cs3_samples_coh, cs3_qc_failed)

cs1_samples <- hist_main %>%
  mutate(FileName = paste0("X", FileName)) %>%
  filter(FileName %in% cs1_samples_coh_qc) %>%
  pull(FileName)
cs2_samples <- hist_main %>%
  mutate(FileName = paste0("X", FileName)) %>%
  filter(FileName %in% cs2_samples_coh_qc) %>%
  pull(FileName)
cs3_samples <- hist_main %>%
  mutate(FileName = paste0("X", FileName)) %>%
  filter(FileName %in% cs3_samples_coh_qc) %>%
  pull(FileName)

cs123_samples <- gsub("^X", "", c(cs1_samples, cs2_samples, cs3_samples))

# Cohorts selected, QC fails removed
cs1_coh_qc <- cs1_coh %>%
  select(Code.Class, Name, Accession, all_of(cs1_samples))
cs2_coh_qc <- cs2_coh %>%
  select(Code.Class, Name, Accession, all_of(cs2_samples))
cs3_coh_qc <- cs3_coh %>%
  select(Code.Class, Name, Accession, all_of(cs3_samples))

# Normalize to housekeeping genes
cs1_norm <- HKnorm(cs1_coh_qc)
cs2_norm <- HKnorm(cs2_coh_qc)
cs3_norm <- HKnorm(cs3_coh_qc)

# Histotypes for CS1/CS2/CS3 samples
hist <- hist_main %>%
  filter(FileName %in% cs123_samples) %>%
  droplevels()

# Histotypes for each distinct ottaID
hist_df <- distinct(hist, ottaID, hist_final)

# Add ottaID and CodeSet to annotNEW
annotNEW <- annotNEW %>%
  mutate(
    ottaID = gsub(".*_(.*)_[0-9]{2}$", "\\1", File.Name) %>%
      gsub("-(N1|R)$", "", .),
    CodeSet = recode_factor(
      geneRLF,
      OTTA2017_C6082 = "CS4",
      OTTA2018_C6830 = "CS5",
      OTTA2018_C8440 = "CS6"
    )
  )

# CS4 data excluding samples that failed QC
cs4_exp <- filter(annotNEW, CodeSet == "CS4")
cs4_qc <- NanoStringQC(raw = cs4, exp = cs4_exp, detect = 50, sn = 100)
cs4_qc_failed <- filter(cs4_qc, QCFlag == "Failed")[["File.Name"]]
cs4_dat <- HKnorm(cs4) %>% select(-any_of(cs4_qc_failed))

# CS5 data excluding samples that failed QC
cs5_exp <- filter(annotNEW, CodeSet == "CS5")
cs5_qc <- NanoStringQC(raw = cs5, exp = cs5_exp, detect = 50, sn = 100)
cs5_qc_failed <- filter(cs5_qc, QCFlag == "Failed", !grepl("POOL", File.Name))[["File.Name"]]
cs5_dat <- HKnorm(cs5) %>% select(-any_of(cs5_qc_failed))

# Common samples amongst all 5 CodeSets
annot_cs123 <- annot_all %>%
  filter(FileName %in% cs123_samples) %>%
  droplevels()

annot_cs45 <- annotNEW %>%
  filter(CodeSet %in% c("CS4", "CS5"),
         !File.Name %in% c(cs4_qc_failed, cs5_qc_failed))

annot_full <- bind_rows(
  annot_cs123 %>% select(FileName, ottaID, CodeSet),
  annot_cs45 %>% select(FileName = File.Name, ottaID, CodeSet)
)

common_cs_full <- annot_full %>%
  count(CodeSet, ottaID) %>%
  spread(CodeSet, n, fill = 0) %>%
  mutate(
    all_codesets = rowSums(.[, -1]),
    in_cs123 = pmap_lgl(select(., matches("CS")),
                        ~ ..1 > 0 & ..2 > 0 & ..3 > 0 & ..4 == 0 & ..5 == 0),
    in_cs345 = pmap_lgl(select(., matches("CS")),
                        ~ ..1 == 0 & ..2 == 0 & ..3 > 0 & ..4 > 0 & ..5 > 0)
  ) %>%
  filter(!is.na(ottaID))

# Samples removed from multiple chains for downstrea validation
# common_cs_full %>% filter(CS1 > 0, CS2 > 0, CS3 > 0, CS4 > 0 | CS5 > 0)
# common_cs_full %>% filter(CS1 > 0 | CS2 > 0, CS3 > 0, CS4 > 0, CS5 > 0)

# CS3 site-specific samples
## Vancouver
van_samples <- hist %>%
  filter(CodeSet == "CS3", site == "Vancouver") %>%
  pull(FileName)

cs3_norm_van <- select(cs3_norm, 1:3, any_of(paste0("X", van_samples)))

## AOC
aoc_samples <- hist_all %>%
  filter(CodeSet == "CS3" &
           site == "AOC" & (FileName %in% cs123_samples |
                              grepl("POOL", FileName))) |>
  pull(FileName) |>
  paste0("X", ... = _)

cs3_norm_aoc <- cs3_coh %>%
  select(where(is.character), all_of(cs3_samples_coh_qc)) |>
  HKnorm() |>
  select(where(is.character), all_of(aoc_samples))

## USC
usc_samples <- hist_all %>%
  filter(CodeSet == "CS3" &
           site == "USC" & (FileName %in% cs123_samples |
                              grepl("POOL", FileName))) |>
  pull(FileName) |>
  paste0("X", ... = _)

cs3_norm_usc <- cs3_coh %>%
  select(where(is.character), all_of(cs3_samples_coh_qc)) |>
  HKnorm() |>
  select(where(is.character), all_of(usc_samples))

# Find summaryID common to CS1/CS2/CS3
common_cs <- annot_cs123 %>%
  count(CodeSet, summaryID) %>%
  spread(CodeSet, n, fill = 0) %>%
  mutate(
    all_codesets = rowSums(.[, -1]),
    in_all_cs = pmap_lgl(select(., matches("CS")),
                         ~ every(list(..1, ..2, ..3), ~ . != 0))
  )

# Find common summaryID
common_id123 <- with(common_cs_full, ottaID[in_cs123])
common_id345 <- with(common_cs_full, ottaID[in_cs345])

# Find common samples
common_samples123 <- hist %>%
  filter(ottaID %in% common_id123, site == "Vancouver") %>%
  pull(FileName)

# Find common genes
common_genes123 <- list(cs1_norm, cs2_norm, cs3_norm_van) %>%
  map(filter, Code.Class == "Endogenous") %>%
  map("Name") %>%
  reduce(intersect)
common_genes345 <- list(cs3_norm_van, cs4_dat, cs5_dat) %>%
  map(filter, Code.Class == "Endogenous") %>%
  map("Name") %>%
  reduce(intersect)

# Clean data by keeping common samples and genes, add ottaID
cs1_clean <- cs1_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  filter(Name %in% common_genes123) %>%
  select_if(names(.) %in% c("Name", common_samples123)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(FileName, ottaID, all_of(common_genes123))

cs2_clean <- cs2_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  filter(Name %in% common_genes123) %>%
  select_if(names(.) %in% c("Name", common_samples123)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(FileName, ottaID, all_of(common_genes123))

cs3_clean <- cs3_norm_van %>%
  rename_all(~ gsub("^X", "", .)) %>%
  filter(Name %in% common_genes123) %>%
  select_if(names(.) %in% c("Name", common_samples123)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(FileName, ottaID, all_of(common_genes123))

# Random selection of common samples with equal number of histotypes
# Ensure only technical (not biological) replicates are used for normalization
set.seed(1776)
hist_rand1 <- hist |>
  filter(FileName %in% c(cs1_clean$FileName, cs2_clean$FileName, cs3_clean$FileName)) |>
  inner_join(annot_all, by = c("FileName", "CodeSet", "ottaID")) |>
  select(CodeSet, ottaID, hist_final, tissue.source, sample.type) |>
  filter(sample.type != "Rep.BIO", .by = c(ottaID, CodeSet)) |>
  filter(n_distinct(CodeSet) == 3, .by = ottaID) |>
  distinct(ottaID, hist_final) |>
  slice_sample(n = 1, by = hist_final)

# CS3-VAN site mapping used for ensuring CS3 samples are from Vancouver
hist_cs3_van <- hist %>%
  filter(site == "Vancouver") %>%
  mutate(col_name = paste0("X", FileName), .keep = "none")

# Find summaryID common to all site
common_site <- annot_cs123 %>%
  count(site, summaryID) %>%
  spread(site, n, fill = 0) %>%
  mutate(
    all_sites = rowSums(.[, -1]),
    in_all_sites = select(., 2:4) %>%
      pmap_lgl(~ every(list(..1, ..2, ..3), ~ . != 0))
  )

# Find common summaryID
common_site_id <- with(common_site, summaryID[in_all_sites])

# Find common samples
common_site_samples <- hist %>%
  filter(ottaID %in% common_site_id) %>%
  pull(FileName)

# Clean data by keeping common samples and genes, add ottaID
cs3_clean_van <- cs3_norm_van %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select_if(names(.) %in% c("Name", common_site_samples)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-c(CodeSet, hist_final, hist_gr, site))

cs3_clean_aoc <- cs3_norm_aoc %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select_if(names(.) %in% c("Name", common_site_samples)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-c(CodeSet, hist_final, hist_gr, site))

cs3_clean_usc <- cs3_norm_usc %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select_if(names(.) %in% c("Name", common_site_samples)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-c(CodeSet, hist_final, hist_gr, site))

# CS3 site-specific expression and reference data
## Vancouver
cs3_van <- cs3_norm_van %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select(-c(Code.Class, Accession)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  anti_join(hist_rand1, by = c("ottaID", "hist_final")) %>%
  select(-c(CodeSet, hist_final, hist_gr, site))

cs3_van_X <- cs3_norm_van %>% select(Name, !matches("POOL"))
cs3_van_R <- cs3_norm_van %>% select(Name, matches("POOL"))

## AOC
cs3_aoc <- cs3_norm_aoc %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select(-c(Code.Class, Accession)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-c(CodeSet, hist_final, hist_gr, site))

cs3_aoc_X <- cs3_norm_aoc %>%
  select(Name, !matches(paste0("POOL|", paste(hist_rand1$ottaID, collapse = "|"))))
cs3_aoc_R <- cs3_norm_aoc %>%
  select(Name, matches("POOL")) %>%
  `rownames<-`(NULL)

## USC
cs3_usc <- cs3_norm_usc %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select(-c(Code.Class, Accession)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-c(CodeSet, hist_final, hist_gr, site))

cs3_usc_X <- cs3_norm_usc %>%
  select(Name, !matches(paste0("POOL|", paste(hist_rand1$ottaID, collapse = "|"))))
cs3_usc_R <- cs3_norm_usc %>%
  select(Name, matches("POOL")) %>%
  `rownames<-`(NULL)
