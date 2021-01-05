# Load Packages -----------------------------------------------------------

library(tidyverse)
library(nanostringr)
library(magrittr)
library(here)
source(here("src/funs.R"))
cs3 <- read_csv(here("data-raw/cs3.csv"))
annot <- read_csv(here("data-raw/annot.csv"))

# Format expression data
exp0 <- annot %>%
  set_names(gsub("RCC\\.", "", names(.))) %>%
  inset(, "geneRLF", as.character(factor(
    x = .[["geneRLF"]],
    labels = c("HL1", "HL2", "HL3", "HuRef", "CS3", "mini", "CS1", "CS2")
  )))

# Run NanoString QC on CS3
exp_CS3_all <- exp0 %>%
  filter(geneRLF == "CS3") %>%
  inset(, "File.Name", paste0("X", .[["File.Name"]])) %>%
  NanoStringQC(raw = as.data.frame(cs3), exp = ., detect = 50, sn = 170)

# Normalize CS3 to housekeeping genes
cs3.norm <- HKnorm(cs3)

# Pools data
pools <- exp_CS3_all %>%
  filter(sample.type %in% c("XSITE.ORIG", "XSITE.FOR", "POOL")) %>%
  mutate(ref.type = factor(case_when(
    POOL1 == "Yes" ~ "POOL 1",
    POOL2 == "Yes" ~ "POOL 2",
    POOL3 == "Yes" ~ "POOL 3",
    cross.site.control == "Yes" ~ "Cross Site",
    TRUE ~ NA_character_
  ))) %>%
  filter(ref.type != "Cross Site") %>%
  mutate(
    nanostring.date = as.Date(nanostring.date),
    averageHK = log2(averageHK),
    nanostring.site = case_when(
      nanostring.site == "AOC" ~ "Melbourne",
      nanostring.site == "USC" ~ "San Francisco",
      TRUE ~ "Vancouver"
    )
  )

# Pools data for each site
vanpools <- pools %>%
  filter(nanostring.site == "Vancouver") %>%
  pull(File.Name) %>%
  select(cs3.norm, .)

AOCpools <- pools %>%
  filter(nanostring.site == "Melbourne") %>%
  pull(File.Name) %>%
  select(cs3.norm, .)

USCpools <- pools %>%
  filter(nanostring.site == "San Francisco") %>%
  pull(File.Name) %>%
  select(cs3.norm, .)

# Cross-site data
cross.site <- exp_CS3_all %>%
  filter(sample.type %in%
                  c("Cross site", "XSITE.ORIG", "XSITE.FOR", "POOL")) %>%
  mutate(ref.type = factor(case_when(
    POOL1 == "Yes" ~ "POOL 1",
    POOL2 == "Yes" ~ "POOL 2",
    POOL3 == "Yes" ~ "POOL 3",
    cross.site.control == "Yes" ~ "Cross Site",
    TRUE ~ NA_character_
  ))) %>%
  filter(ref.type == "Cross Site",
                summaryID != "TVAN20681") %>%
  mutate(nanostring.site = factor(case_when(
    nanostring.site == "AOC" ~ "Melbourne",
    nanostring.site == "USC" ~ "San Francisco",
    TRUE ~ "Vancouver"
  )))

# Cross-site data for each site
VanRefs.gx <- cross.site %>%
  filter(nanostring.site == "Vancouver") %>%
  arrange(summaryID) %>%
  pull(File.Name) %>%
  select(cs3.norm, .)

AOCRefs.gx <- cross.site %>%
  filter(nanostring.site == "Melbourne") %>%
  arrange(summaryID) %>%
  pull(File.Name) %>%
  select(cs3.norm, .)

USCRefs.gx <- cross.site %>%
  filter(nanostring.site == "San Francisco") %>%
  arrange(summaryID) %>%
  pull(File.Name) %>%
  select(cs3.norm, .)

# Everything is calibrated to Vancouver
AOCRefs.gx2 <-
  t(refMethod(t(AOCRefs.gx), t(vanpools), t(AOCpools)))
USCRefs.gx2 <-
  t(refMethod(t(USCRefs.gx), t(vanpools), t(USCpools)))

# Combinations of cross-site gene expression
sites <- c("USC", "AOC", "VAN")
all_xsites <- combn(sites, 2) %>%
  as_tibble(.name_repair = "unique") %>%
  set_names(map_chr(., paste, collapse = "_vs_"))

# Combined gene expression
all_gx <- list(USCRefs.gx2, AOCRefs.gx2, VanRefs.gx) %>%
  set_names(sites) %>%
  map(~ as.data.frame(t(.)))

# Concordance measures for all genes averaged across samples
all_metrics <- all_xsites %>%
  imap_dfr(~ {
    pmap_dfr(all_gx[.x], ~ cor_stats(.x, .y)) %>%
      mutate(Sites = .y)
  }) %>%
  gather(key = "Metric", value = "Expression", -Sites)

# Plot all combinations of cross-site concordance measure histograms
p <- ggplot(all_metrics, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Cross-Site Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p)
