# Prepare data for R package

#  accession numbers ---------------------------

library(blantyreESBL)
library(tidyverse)
library(devtools)

# load virulence determinents

vf <- read_csv(here("data-raw/all_dassim_esco_vf_ariba.csv"))


vf %>% 
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name)) %>% 
  pivot_longer(-name, 
               names_to = c( "cluster", ".value"), 
               names_sep = "\\.") %>%
  filter(match == "yes") %>% 
  select(name, ref_seq) ->
  dassimEcoli_BTEcoli.virulence
  
use_data(dassimEcoli_BTEcoli.virulence, overwrite = TRUE)


dassimEcoli_BTEcoli.virulence %>% 
  filter(grepl("stx|ltcA|sta|eae|aatA|aggR|aaiC|ipaH|ipaD", ref_seq)) %>% 
  mutate(ref_seq = gsub("_[0-9|A-Z|a-z]*$","", ref_seq)) %>% 
  mutate(Pathotype = case_when(
    grepl("stx", ref_seq) & grepl("eae", ref_seq) ~ "EHEC",
    grepl("stx", ref_seq) ~ "STEC",
    grepl("eae", ref_seq) ~ "aEPEC/EPEC",
    grepl("aatA|aggR|aaiC", ref_seq) ~ "EAEC",
    grepl("ltcA|sta", ref_seq) ~ "ETEC",
    grepl("ipaH|ipaD", ref_seq) ~ "EIEC",
    TRUE ~ NA_character_ )) %>% 
  select(name, Pathotype) %>% 
  unique() %>% 
  #pivot_wider(id_cols = name, values_from = pathotype, names_from = pathotype,
  #            values_fn = length, values_fill = NA_integer_) %>% 
  mutate(name = gsub("#", "_", name)) %>%  
  as.data.frame() ->
  vir

rownames(vir) <- vir$name           

vir



btESBL_sequence_sample_metadata %>%
  filter(lane %in% btESBL_coregene_tree_esco$tip.label) %>%
  left_join(
    rbind(read_tsv(#"~/Documents/PhD/Thesis/bookdown/chapter_7/
      # checkm_quast/D1/transposed_report.tsv"
      here(
        "data-raw/QUAST_report1.tsv"
      )) ,
      read_tsv(#"~/Documents/PhD/Thesis/bookdown/chapter_7/
        #checkm_quast/D220190318/transposed_report.tsv"
        here(
          "data-raw/QUAST_report2.tsv"
        )),
      read_tsv(# "~/Documents/PhD/Thesis/bookdown/chapter_7/
        #  checkm_quast/D220190503/transposed_report.tsv"
        here(
          "data-raw/QUAST_report3.tsv"
        ))) %>%
      transmute(
        lane = gsub("\\.contigs_spades", "",  Assembly),
        number_of_contigs = `# contigs`,
        N50 = N50
      ),
    by = "lane"
  ) %>%
  rename(date_of_collection = data_date) %>%
  select(accession,
         lane,
         supplier_name,
         pid,
         date_of_collection,
         number_of_contigs,
         N50) %>%
  # add in phylogroups
  left_join(# rownames(df_hclusts) <- df_hclusts$Taxon
    left_join(read_csv(#"~/Documents/PhD/Thesis/bookdown/chapter_7/
      #phylogroup_and_mlst/mlst.csv"
      here(
        "data-raw/mlst.csv"
      )),
      read_csv(#"~/Documents/PhD/Thesis/bookdown/chapter_7/
        #phylogroup_and_mlst/phylogroups.csv"
        here(
          "data-raw/phylogroups.csv"
        )),
      by = c("lane" = "Lane")) %>%
      mutate(ST = case_when(
        lane %in% c(
  # update those that were novel but now not
          "28099_1_144",
          "28099_1_18",
          "28099_1_159",
          "28099_1_249",
          "28099_1_23",
          "28099_1_249",
          "28099_1_330",
          "26141_1_141"
        ) ~ "9847",
        TRUE ~ ST
      )),
    by = "lane") %>% 
# add in pathotypes
  left_join(
    dassimEcoli_BTEcoli.virulence %>% 
      filter(grepl("stx|ltcA|sta|eae|aatA|aggR|aaiC|ipaH|ipaD", ref_seq)) %>% 
      mutate(ref_seq = gsub("_[0-9|A-Z|a-z]*$","", ref_seq)) %>% 
      mutate(Pathotype = case_when(
        grepl("stx", ref_seq) & grepl("eae", ref_seq) ~ "EHEC",
        grepl("stx", ref_seq) ~ "STEC",
        grepl("eae", ref_seq) ~ "aEPEC/EPEC",
        grepl("aatA|aggR|aaiC", ref_seq) ~ "EAEC",
        grepl("ltcA|sta", ref_seq) ~ "ETEC",
        grepl("ipaH|ipaD", ref_seq) ~ "EIEC",
        TRUE ~ NA_character_ )) %>% 
      select(name, Pathotype) %>% 
      unique() %>% 
      #pivot_wider(id_cols = name, values_from = pathotype, names_from = pathotype,
      #            values_fn = length, values_fill = NA_integer_) %>% 
      mutate(name = gsub("#", "_", name)),
    by = c("lane" = "name")
  ) -> 
  dassimEcoli_BTEcoli.accession

use_data(dassimEcoli_BTEcoli.accession, overwrite = TRUE)
