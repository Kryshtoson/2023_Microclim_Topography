#-------------------------------------------------------------------------
# START OF THE CODE S2
#-------------------------------------------------------------------------
#' MANIPULATION OF SPECIES DATA
#' author Krystof Chytry, krystof.chytry@gmail.com
# R version 4.3.0
# IDE Pycharm
library(vegan) # version 2.6-4
library(tidyverse)
library(readxl)
library(DSSAT) # version 0.0.7

#------------------------------------------------------------------------
# species data into wide form
#------------------------------------------------------------------------
brbl2num <- function(v){
  factor(v, labels = c(0.01, 1, 3, 5, 10, 20, 38, 68, 88),
         levels = c('r', '+', '1', '2m', '2a', '2b', '3', '4', '5')) |>
    as.character() |>
    as.numeric()
}

df <- read_xlsx('raw_input\\T5_species_data_v2-0.xlsx') |>
  mutate(value = brbl2num(COVER)) |>
  select(logger_ID, name = TAXON, value) |>
  pivot_wider(values_fill = 0, values_fn = max)

#------------------------------------------------------------------------
# Traits availability
# We considered only these plots where total area covered by species consisted
# by 70% of species with traits available
#------------------------------------------------------------------------
traits_ava <- df |> pivot_longer(-1) |>
  left_join(
    read_excel('raw_input\\Species_data.xlsx') |>
      select(name = species,
             data_ava = trait_height)) |>
  mutate(data_ava = !is.na(data_ava)) |>
  group_by(logger_ID) |>
  mutate(cover_sum = sum(value)) |>
  group_by(logger_ID, data_ava, cover_sum) |>
  summarise(prop_ava = sum(value)) |>
  ungroup() |>
  mutate(prop_ava = prop_ava/cover_sum) |>
  filter(data_ava == T) |>
  mutate(traits_available = prop_ava > .7) |>
  select(logger_ID, traits_available)

#-------------------------------------------------------------------------
# Computation of mean EIVs, community means of traits, and diversity indices
#------------------------------------------------------------------------

mean_eivs <- read_xlsx('raw_input\\T5_species_data_v2-0.xlsx') |>
  filter(TAXON != 'no vegetation') |>
  left_join(read_excel('raw_input\\Species_data.xlsx') |>
              rename(TAXON = species)) |>
  group_by(logger_ID) |>
  mutate(across(c('Temperature', 'Moisture', 'Reaction', 'Nutrients'),
                as.numeric)) |>
  mutate('sp_richness' = n()) |>
  summarise_at(c('Temperature',
                 'Moisture',
                 'Reaction',
                 'Nutrients',
                 'sp_richness',
                 'trait_height',
                 'trait_sla',
                 'trait_N',
                 'trait_suc'),
               function(x) mean(x, na.rm = T)) |>
  left_join(tibble(logger_ID = df[[1]], Shannon = diversity(df[-1]))) |>
  left_join(traits_ava) |>
  mutate_cond(!traits_available,
              trait_height = NA,
              trait_sla = NA,
              trait_N = NA,
              trait_suc = NA)

mean_eivs |> pivot_longer(-1) |> group_by(name) |>
  summarise(mean = mean(value, na.rm = T),sd = sd(value, na.rm = T),
            max = max(value, na.rm = T), min = min(value, na.rm = T))

mean_eivs |> write_csv('data\\Comprop_means.csv')

#-------------------------------------------------------------------------
# END OF THE CODE S2
#-------------------------------------------------------------------------