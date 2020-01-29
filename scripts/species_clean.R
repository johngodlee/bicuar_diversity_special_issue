# Cleaning the species list on the chosen plots
# John Godlee (johngodlee@gmail.com)
# 2019_11_28

# Preamble ----

# Remove old crap
rm(list=ls())
#dev.off()

# Set working directory to the location of the source file
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
library(dplyr)
library(BIOMASS)
library(taxize)

source("scripts/functions.R")

# Import data ----
stems_list <- readRDS("data/stems_list.rds")

# Inspect and correct species names ----

taxo_list <- lapply(stems_list, function(x){
  df <- x %>% mutate(genus = gsub("([A-z]+).*", "\\1", .$species_binomial),
    species = gsub("^([A-z]+) ", "", .$species_binomial))
  taxo <- correctTaxo(genus = df$genus, 
    species = df$species)
  taxo[taxo$nameModified %in% c("TRUE", "SpNotFound"),] %>% distinct()
})

# Bicuar big plots 
stems_list[[1]] <- stems_list[[1]] %>%
  mutate(species_binomial = case_when(
    species_binomial == "Brachystegia sp. plot_5" ~ "Brachystegia spiciformis",
    species_binomial == "Brachystegia sp. plot_6" ~ "Brachystegia spiciformis",
    species_binomial == "Brachystegia sp. plot_7" ~ "Brachystegia spiciformis",
    species_binomial == "Brachystegia sp. plot_11" ~ "Brachystegia spiciformis",
    species_binomial == "Combretum sp. plot_1" ~ "Combretum collinum",
    species_binomial == "Combretum sp. plot_10" ~ "Combretum psidioides",
    species_binomial == "Combretum sp. plot_11" ~ "Combretum psidioides",
    species_binomial == "Combretum sp. plot_13" ~ "Combretum hereroense",
    species_binomial == "Combretum sp. plot_14" ~ "Combretum celastroides",
    species_binomial == "Combretum sp. plot_2" ~ "Combretum psidioides",
    species_binomial == "Combretum sp. plot_9" ~ "Combretum collinum",
    species_binomial == "Combretum sellastroides" ~ "Combretum celastroides",
    species_binomial == "Monotes sp." ~ "Monotes angolensis",
    species_binomial == "Pterocarpus lucens subsp. antunesi" ~ "Pterocarpus lucens",
    species_binomial == "Baphia bequaerti" ~ "Baphia bequaertii",
    species_binomial == "Hippocratea parvifolia" ~ "Elachyptera parvifolia",
    TRUE ~ as.character(species_binomial)
  )) 

# Bicuar degrad plots
stems_list[[2]] <- stems_list[[2]] %>%
  mutate(species_binomial = case_when(
    species_binomial == "Combretum sp. degrad_plot_3" ~ "Combretum celastroides",
    species_binomial == "Combretum sp. degrad_plot_3_b" ~ "Combretum collinum",
    species_binomial == "Combretum sp. degrad_plot_2" ~ "Combretum celastroides",
    species_binomial == "Ekebergia sp." ~ "Ekebergia benguelensis",
    species_binomial == "Acacia sp." ~ "Acacia reficiens",
    species_binomial == "Vachellia reficiens" ~ "Acacia reficiens",
    species_binomial == "Baphia bequaerti" ~ "Baphia bequaertii",
    TRUE ~ as.character(species_binomial))) 

# DRC plots
stems_list[[3]] <- stems_list[[3]] %>%
  mutate(species_binomial = case_when(
    species_binomial == "Trillesanthus macroura" ~ "Trillesanthus macrourus", 
    TRUE ~ as.character(species_binomial)
  ))

# Kilwa plots
stems_list[[4]] <- stems_list[[4]] %>%
  mutate(species_binomial = case_when(
    species_binomial == "Securidaca longipedunculata" ~ "Securidaca longepedunculata",
    species_binomial == "Vachellia reficiens" ~ "Acacia reficiens",
    TRUE ~ as.character(species_binomial)))

# Nhambita plots
stems_list[[5]] <- stems_list[[5]] %>%
  mutate(species_binomial = case_when(
    species_binomial == "Searsia chiridensis" ~ "Rhus chirindensis",
    species_binomial == "Securidaca longipedunculata" ~ "Securidaca longepedunculata",
      TRUE ~ as.character(species_binomial)))

# families_list <- lapply(stems_list, function(x){
#   classification(unique(as.character(x$species_binomial)),
#     db = "gbif", accepted = TRUE)})
# saveRDS(families_list, "data/families_list.rds")
families_list <- readRDS("data/families_list.rds")

family <- unlist(sapply(families_list, function(x){
  sapply(x, function(y){
    c(as.character(taxon_extract(y, "family")))
  })
}))

species <- unlist(sapply(families_list, function(x){
  sapply(x, function(y){
    c(as.character(taxon_extract(y, "species")))
  })
}))

family_df <- data.frame(family, species, stringsAsFactors = FALSE)

# Manually fix some rows
family_df[row.names(family_df) == "kilwa.Margaritaria indet",1] <- "Phyllanthaceae"
family_df[row.names(family_df) == "kilwa.Margaritaria indet",2] <- "Margaritaria indet"
family_df[row.names(family_df) == "drc.Vangueriopsis indet",2] <- "Vangueriopsis indet"
family_df[row.names(family_df) == "kilwa.Diospyros indet",2] <- "Diospyros indet"
family_df[row.names(family_df) == "kilwa.Ehretia indet",2] <- "Ehretia indet"
family_df[row.names(family_df) == "kilwa.Markhamia indet",2] <- "Markhamia indet"
family_df[row.names(family_df) == "kilwa.Strychnos indet",2] <- "Strychnos indet"
family_df[row.names(family_df) == "kilwa.Sorindeia indet",2] <- "Sorindeia indet"
family_df[row.names(family_df) == "kilwa.Psidium indet",2] <- "Psidium indet"
family_df[row.names(family_df) == "nham.Combretum indet",2] <- "Combretum indet"
family_df[row.names(family_df) == "nham.Terminalia indet",2] <- "Terminalia indet"
family_df[row.names(family_df) == "nham.Diospyros indet",2] <- "Diospyros indet"
family_df[row.names(family_df) == "kilwa.Zanthoxylum indet",2] <- "Zanthoxylum indet"
family_df[row.names(family_df) == "kilwa.Huberantha stuhlmannii",2] <- "Huberantha stuhlmannii"
family_df <- rbind(family_df, 
  data.frame(family = "Combretaceae", species = "Pteleopsis anisoptera"),
  data.frame(family = "Dipterocarpaceae", species = "Monotes angolensis"),
  data.frame(family = "Rutaceae", species = "Zanthoxylum chalybeum"),
  data.frame(family = "Dipterocarpaceae", species = "Trillesanthus macrourus"),
  data.frame(family = "Fabaceae", species = "Aganope stuhlmannii"),
  data.frame(family = "Combretaceae", species = "Pteleopsis myrtifolia"),
  data.frame(family = "Fabaceae", species = "Piliostigma thonningii"),
  data.frame(family = "Fabaceae", species = "Philenoptera violacea"),
  data.frame(family = "Malvaceae", species = "Dombeya mupangae"),
  data.frame(family = "Malvaceae", species = "Rhodognaphalon schumannianum")
  )

rownames(family_df) <- NULL

names(family_df) <- c("family", "species_binomial")

family_df_clean <- family_df[!duplicated(family_df),]

saveRDS(family_df_clean, "data/tree_family_lookup.rds")

# Add family names to stems_list
stems_list_clean <- lapply(stems_list, function(x){
  left_join(x, family_df_clean, by = c("species_binomial" = "species_binomial"))
})

stems_list_clean <- lapply(stems_list_clean, function(x){
  x[x$species_binomial == "Indet indet" & !is.na(x$species_binomial),'species_binomial'] <- NA_character_
  return(x)
})

saveRDS(stems_list_clean, "data/stems_list.rds")
