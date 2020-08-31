# Load Packages -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(nanostringr)
library(here)


# Data Processing ---------------------------------------------------------

# Read in raw data
annot <- read_csv(here("data-raw/annot.csv"), col_types = cols())
cs1 <- read_csv(here("data-raw/cs1.csv"), col_types = cols())
cs2 <- read_csv(here("data-raw/cs2.csv"), col_types = cols())
cs3 <- read_csv(here("data-raw/cs3.csv"), col_types = cols())
pools <- read_excel(here("data-raw/RNA-Pools-Source_CS1-2-3.xlsx"))
ref_pools <- readRDS(here("data/van_pools_cs3.rds"))

# Normalize to housekeeping genes
cs1_norm <- HKnorm(cs1)
cs2_norm <- HKnorm(cs2)
cs3_norm <- HKnorm(cs3)

# Filter annotations for CS 1, 2, 3
annot_cs_all <- annot %>%
  filter(RCC.geneRLF %in% c("OvCa2103_C953", "PrOTYPE2_v2_C1645", "OTTA2014_C2822")) %>%
  mutate(
    CodeSet = recode_factor(
      RCC.geneRLF,
      `OvCa2103_C953` = "CS1",
      `PrOTYPE2_v2_C1645` = "CS2",
      `OTTA2014_C2822` = "CS3"
    )
  )

# Histotypes
hist <- annot_cs_all %>%
  transmute(
    FileName = RCC.File.Name,
    CodeSet,
    revHist = case_when(
      revHist == "CCC" ~ "CCOC",
      revHist == "ENOCa" ~ "ENOC",
      TRUE ~ revHist
    ),
    hist_gr = ifelse(revHist == "HGSC", "HGSC", "non-HGSC")
  )

# Find summaryID common to all CS
cs_common <- annot_cs_all %>%
  count(CodeSet, summaryID) %>%
  spread(CodeSet, n, fill = 0) %>%
  mutate(
    all_codesets = rowSums(.[, -1]),
    in_all_cs = select(., 2:4) %>%
      pmap_lgl(~ every(list(..1, ..2, ..3), ~ . != 0))
  )

# Find common summaryID
common_id <- with(cs_common, summaryID[in_all_cs])

# Find common samples
common_samples <- annot_cs_all %>%
  filter(summaryID %in% common_id) %>%
  pull(RCC.File.Name)

# Find common genes
common_genes <- list(cs1_norm, cs2_norm, cs3_norm) %>%
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
  mutate(ottaID = annot_cs_all$ottaID[match(FileName, annot_cs_all$RCC.File.Name)]) %>%
  spread(Name, value)

cs2_clean <- cs2_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  filter(Name %in% common_genes) %>%
  select_if(names(.) %in% c("Name", common_samples)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  mutate(ottaID = annot_cs_all$ottaID[match(FileName, annot_cs_all$RCC.File.Name)]) %>%
  spread(Name, value)

cs3_clean <- cs3_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  filter(Name %in% common_genes) %>%
  select_if(names(.) %in% c("Name", common_samples)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  mutate(ottaID = annot_cs_all$ottaID[match(FileName, annot_cs_all$RCC.File.Name)]) %>%
  spread(Name, value)
