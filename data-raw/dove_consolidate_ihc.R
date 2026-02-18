## Consolidate Mike's IHC data

# Load packages
library(TMAtools)
library(here)
library(readxl)
library(writexl)
library(forcats)
tma_dir <- here("data-raw/CT14-05_Scores-Deconvoluted-Summarized_MA27SEPT2016/")
output_dir <- here(tma_dir, "outputs")

# Keep only relevant biomarkers from ENOC rules
biomarker_translation <-
  system.file("extdata", "biomarker_rules_enoc.xlsx", package = "TMAtools") |>
  read_excel(sheet = "translation") |>
  filter(biomarker %in% c("ARID1A", "Napsin A", "PR", "WT1", "p16", "p53"))

biomarker_consolidation <-
  system.file("extdata", "biomarker_rules_enoc.xlsx", package = "TMAtools") |>
  read_excel(sheet = "consolidation") |>
  filter(biomarker %in% c("ARID1A", "Napsin A", "PR", "WT1", "p16", "p53"),
         !(biomarker == "p16" & rule_type == "any;any"))

biomarker_rules_file <- file.path(output_dir, "biomarker_rules.xlsx")
write_xlsx(
  list(translation = biomarker_translation, consolidation = biomarker_consolidation),
  biomarker_rules_file
)

# Mastertable of deconvoluted IHC scores
mastertable_scores <- read_excel(file.path(tma_dir, "DOV_CT14-05 staining overview.xlsx"),
                                 sheet = "DOVE IHC score mastertable")

# Additional processing: rename columns, fix score values
deconvoluted_scores <- mastertable_scores |>
  mutate(across("WT1_GPEC_PR_score 1", ~ as.numeric(ifelse(.x ==  "`", "9", .x)))) |>
  mutate(
    `PR_GPEC_PR_score 2` = case_when(
      `DOVE ID` == "DOVE42315" ~ 1,
      `DOVE ID` == "DOVE42167" ~ 1,
      .default = `PR_GPEC_PR_score 2`
    ),
    `WT1_DAKO_MK_score 2` = case_when(
      `DOVE ID` == "DOVE42315" ~ 0,
      `DOVE ID` == "DOVE42167" ~ 2,
      .default = `WT1_DAKO_MK_score 2`
    ),
    `WT1_GPEC_PR_score 2` = case_when(
      `DOVE ID` == "DOVE42315" ~ 0,
      `DOVE ID` == "DOVE42167" ~ 2,
      .default = `WT1_GPEC_PR_score 2`
    ),
    `NAPSA_DAKO_MK_score 2` = case_when(
      `DOVE ID` == "DOVE42315" ~ 0,
      `DOVE ID` == "DOVE42167" ~ 0,
      .default = `NAPSA_DAKO_MK_score 2`
    ),
    `p16_GPEC_PR_score 2` = case_when(
      `DOVE ID` == "DOVE42315" ~ 1,
      `DOVE ID` == "DOVE42167" ~ 2,
      .default = `p16_GPEC_PR_score 2`
    )
  ) |>
  rename_with(~ gsub("TP53", "p53", .x), matches("score")) |>
  rename_with(~ gsub("NAPSA", "Napsin A", .x), matches("score")) |>
  rename_with(~ gsub("_score ", ".c", .x), matches("score")) |>
  rename(WT1.c3 = WT1_GPEC_PR.c1, WT1.c4 = WT1_GPEC_PR.c2) |>
  rename_with(~ gsub("_DAKO_MK|_GPEC_PR", "", .x)) |>
  select(`DOVE ID`, `CORE#`, matches("p53|ARID1A|PR|WT1|Napsin A|p16"))
deconvoluted_scores_file <- file.path(output_dir, "deconvoluted_scores.xlsx")
write_xlsx(deconvoluted_scores, deconvoluted_scores_file)

# Numerical scores translated to nominal scores
translated_scores <- translate_scores(
  biomarkers_file = deconvoluted_scores_file,
  biomarker_rules_file = biomarker_rules_file,
  output_file = file.path(output_dir, "translated_scores.xlsx")
)

# Consolidate biomarker scores into a single value per patient
consolidated_scores <- consolidate_scores(
  biomarkers_file = file.path(output_dir, "translated_scores.xlsx"),
  biomarker_rules_file = biomarker_rules_file
)

# Final biomarker scores with reordered levels and missing cases excluded
final_scores <- consolidated_scores |>
  mutate(
    ottaID = `DOVE ID`,
    core_id = as.numeric(`CORE#`),
    across(
      c("p53", "arid1a", "pr", "wt1", "napsin a", "p16"),
      ~ fct_na_level_to_value(.x, extra_levels = "Unk")
    ),
    p53 = fct_recode(p53, "wild type" = "__check__"),
    pr = fct_rev(pr),
    wt1 = fct_rev(wt1),
    p16 = fct_rev(p16),
    .keep = "none"
  ) |>
  relocate(ottaID, core_id) |>
  rename_with(toupper, where(is.factor)) |>
  filter(core_id == min(core_id), .by = ottaID) |>
  select(-core_id)
write_xlsx(final_scores, file.path(output_dir, "final_scores.xlsx"))
