# Load Packages -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(naniar)
library(nanostringr)
library(ottaOvca)
library(here)


# Data Processing ---------------------------------------------------------

# Read in raw OTTA data
data("annot", "rawOVCA2", "rawPROT", "rawOTTA", package = "otta")
cs1 <- rawOVCA2
cs2 <- rawPROT
cs3 <- rawOTTA

data("od.otta")
od.otta <- to_native_type(od.otta)

cohorts <- read_excel(here("data-raw/bc-842_avail_cosp_2020-08-31.xlsx"),
                      sheet = "data")
pools <- read_excel(here("data-raw/RNA-Pools-Source_CS1-2-3.xlsx"))
ref_pools <- readRDS(here("data/van_pools_cs3.rds"))

# Pairwise CodeSet comparisons
codesets <- c("CS1", "CS2", "CS3")
all_codesets <- combn(codesets, 2) %>%
  as_tibble(.name_repair = "unique") %>%
  set_names(map_chr(., paste, collapse = "_vs_"))

# Pairwise site comparisons
sites <- c("USC", "AOC", "VAN")
all_xsites <- combn(sites, 2) %>%
  as_tibble(.name_repair = "unique") %>%
  set_names(map_chr(., paste, collapse = "_vs_"))

# Filter for specific cohorts and extract samples
cs1_samples <- cohorts %>%
  filter(file_source == "cs1",
         cohort %in% c("MAYO", "OOU", "OOUE", "VOA", "MTL")) %>%
  pull(col_name)

cs2_samples <- cohorts %>%
  filter(file_source == "cs2",
         cohort %in% c("MAYO", "OOU", "OOUE", "OVAR3", "VOA", "ICON7", "JAPAN", "MTL", "POOL-CTRL")) %>%
  pull(col_name)

cs3_samples <- cohorts %>%
  filter(file_source == "cs3",
         cohort %in% c("DOVE4", "OOU", "OOUE", "TNCO", "VOA", "POOL-1", "POOL-2", "POOL-3")) %>%
  pull(col_name)

cs123_samples <- gsub("^X", "", c(cs1_samples, cs2_samples, cs3_samples))

# Select only samples from these cohorts in the data
cs1_coh <- select(cs1, Code.Class, Name, Accession, all_of(cs1_samples))
cs2_coh <- select(cs2, Code.Class, Name, Accession, all_of(cs2_samples))
cs3_coh <- select(cs3, Code.Class, Name, Accession, all_of(cs3_samples))

# Normalize to housekeeping genes
cs1_norm <- HKnorm(cs1_coh)
cs2_norm <- HKnorm(cs2_coh)
cs3_norm <- HKnorm(cs3_coh)

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
  rename(FileName = RCC.File.Name, site = RCC.nanostring.site)

# Filter annotations for CS1/CS2/CS3 samples
annot_cs <- annot_all %>%
  filter(FileName %in% cs123_samples) %>%
  droplevels()

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

# Histotypes
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
    hist_gr = ifelse(revHist == "HGSC", "HGSC", "non-HGSC"),
    site
  )

# Histotypes for CS1/CS2/CS3 samples
hist <- hist_all %>%
  filter(FileName %in% cs123_samples) %>%
  droplevels()

# Histotypes for each distinct ottaID
hist_df <- distinct(hist, ottaID, revHist)

# CS3 site-specific samples
## Van
van_samples <- hist %>%
  filter(CodeSet == "CS3", site == "Vancouver") %>%
  pull(FileName)
full_samples <- hist_all %>%
  filter(CodeSet == "CS3") %>%
  pull(FileName)

cs3_norm_van <- select(cs3_norm, 1:3, any_of(paste0("X", van_samples)))
cs3_norm_full <- cs3 %>%
  HKnorm() %>%
  select(1:3, any_of(paste0("X", full_samples)))

cs3_van <- cs3_norm_van %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select(-c(Code.Class, Accession)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-c(CodeSet, revHist, hist_gr, site))

cs3_van_X <- cs3_norm_van %>% select(Name, !matches("POOL"))
cs3_van_R <- cs3_norm_van %>% select(Name, matches("POOL"))

## AOC
aoc_samples <- hist %>%
  filter(CodeSet == "CS3", site == "AOC") %>%
  pull(FileName)

cs3_norm_aoc <- select(cs3_norm, 1:3, any_of(paste0("X", aoc_samples)))

cs3_aoc <- cs3_norm_aoc %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select(-c(Code.Class, Accession)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-c(CodeSet, revHist, hist_gr, site))

cs3_aoc_X <- cs3_norm_aoc %>% select(Name, !matches("POOL"))
cs3_aoc_R <- cs3_norm_aoc %>% select(Name, matches("POOL"))

## USC
usc_samples <- hist %>%
  filter(CodeSet == "CS3", site == "USC") %>%
  pull(FileName)

cs3_norm_usc <- select(cs3_norm, 1:3, any_of(paste0("X", usc_samples)))

cs3_usc <- cs3_norm_usc %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select(-c(Code.Class, Accession)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-c(CodeSet, revHist, hist_gr, site))

cs3_usc_X <- cs3_norm_usc %>% select(Name, !matches("POOL"))
cs3_usc_R <- cs3_norm_usc %>% select(Name, matches("POOL"))

# Find summaryID common to all CS
common_cs <- annot_cs %>%
  count(CodeSet, summaryID) %>%
  spread(CodeSet, n, fill = 0) %>%
  mutate(
    all_codesets = rowSums(.[, -1]),
    in_all_cs = select(., 2:4) %>%
      pmap_lgl(~ every(list(..1, ..2, ..3), ~ . != 0))
  )

# Find common summaryID
common_id <- with(common_cs, summaryID[in_all_cs])

# Find common samples
common_samples <- hist %>%
  filter(ottaID %in% common_id, site == "Vancouver") %>%
  pull(FileName)

# Find common genes
common_genes <- list(cs1_norm, cs2_norm, cs3_norm_van) %>%
  map(filter, Code.Class == "Endogenous") %>%
  map("Name") %>%
  reduce(intersect)

# Clean data by keeping common samples and genes, add ottaID
cs1_clean <- cs1_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  filter(Name %in% common_genes) %>%
  select_if(names(.) %in% c("Name", common_samples)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(FileName, ottaID, all_of(common_genes))

cs2_clean <- cs2_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  filter(Name %in% common_genes) %>%
  select_if(names(.) %in% c("Name", common_samples)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(FileName, ottaID, all_of(common_genes))

cs3_clean <- cs3_norm_van %>%
  rename_all(~ gsub("^X", "", .)) %>%
  filter(Name %in% common_genes) %>%
  select_if(names(.) %in% c("Name", common_samples)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(FileName, ottaID, all_of(common_genes))

# Find summaryID common to all site
common_site <- annot_cs %>%
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
  select(-c(CodeSet, revHist, hist_gr, site))

cs3_clean_aoc <- cs3_norm_aoc %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select_if(names(.) %in% c("Name", common_site_samples)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-c(CodeSet, revHist, hist_gr, site))

cs3_clean_usc <- cs3_norm_usc %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select_if(names(.) %in% c("Name", common_site_samples)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-c(CodeSet, revHist, hist_gr, site))
