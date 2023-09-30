#-------------------------------------------------------------------------
# START OF THE CODE S5
#-------------------------------------------------------------------------
library(tidyverse)
species_resolution_models <- read_csv('results\\meta_results\\models_resolution_species.csv')
comprops_resolution_models <- read_csv('results\\meta_results\\models_resolution_comprops.csv')

#' "Across all species, the variation in their distribution explained by one morphometric variable
#' alone at one scale varied between 0.0% and 31.2% (pseudo-R² values)"
round(range(species_resolution_models$dev_explained * 100), 2)
#> 0.00, 30.98

#' while the variation explained by elevation varied between
round(range(species_resolution_models$null_mod * 100), 2)

#' TWI had the highest pseudo-R2 values averaged over all species (Figure 1A)
#' at resolutions from 71 to 301 m (71: 5.8%, 101: 6.2%, 151: 6.6%, 201: 5.7% and 301: 6.0%
#' followed by TP5 (301: 5.7%), SLP (101: 5.7%)
#' and VR5 (151 and 201: 5.6%).
species_resolution_models |>
  group_by(index, resolution) |>
  summarise_at('dev_explained', mean) |>
  arrange(-dev_explained) |>
  mutate(dev_explained = round(dev_explained*100, 1)) |>
  print(n = 20)

#' Fine-scale morphometric parameters with the highest pseudo-R² values were
#' TR5 (19: 4.8% and 15: 4.7%),
#' SLP (19: 4.8% and 15: 4.6%)
#' and TR3 (19: 4.6% and 15: 4.3).
species_resolution_models |>
  group_by(index, resolution) |>
  summarise_at('dev_explained', mean) |>
  arrange(-dev_explained) |>
  filter(resolution < 25) |>
  mutate(dev_explained = round(dev_explained*100, 1)) |>
  print(n = 20)

#' For community attributes, the variance explained by models
#' including a single morphometric variable at one resolution ranged from 0.01% to 19.48%
round(range(comprops_resolution_models$dev_explained * 100), 2)

#' while the variation explained by elevation varied between 7.52% and 78.99%
round(range(comprops_resolution_models$null_mod * 100), 2)

#' In general, community attributes were on average (Figure 1B) best explained by
#' VR5 (6.6% at 151 m) and VR3 (6.4% at 201 m and 6.0% at 151 m).
comprops_resolution_models |>
  group_by(index, resolution) |>
  summarise_at('dev_explained', mean) |>
  arrange(-dev_explained) |>
  mutate(dev_explained = round(dev_explained*100, 1)) |>
  print(n = 20)

#' Best performing variables at fine resolutions were TR5 (15: 4.0%, 19: 4.0%), TR3 (19: 3.9%) and HLI (7: 3.8%)
comprops_resolution_models |>
  group_by(index, resolution) |>
  summarise_at('dev_explained', mean) |>
  arrange(-dev_explained) |>
  filter(resolution < 25) |>
  mutate(dev_explained = round(dev_explained*100, 1)) |>
  print(n = 20)

#' ========================================================================================
#' random models
#' ========================================================================================
random_mods_species <- read_csv('results\\SPECIES_random_models.csv') |>
  left_join(read_csv(paste0('results\\SPECIES_elevation_models.csv'))) |>
    mutate(r2_elevation_pure = rsq_full - rsq,
           r2_topography_pure = rsq_full - elevation_only) |>
  arrange(-r2_topography_pure) %>%
    group_by(species, group) %>%
    slice(1)

random_mods_comprops <- read_csv('results\\COMPROPS_random_models.csv') |>
  left_join(read_csv(paste0('results\\COMPROPS_elevation_models.csv'))) |>
    mutate(r2_elevation_pure = rsq_full - rsq,
           r2_topography_pure = rsq_full - elevation_only) |>
  arrange(-r2_topography_pure) %>%
    group_by(species, group) %>%
    slice(1)

#' The best of the 1,000 models with combinations of five predictors per response variable explained
#' the distribution of 56 (out of 79) species better with coarse than with fine-scale predictor variables
#' (Figure 2). Combining Mixing fine and coarse resolutions only marginally increased the explained
#' variation (of the better model) by 1.2% on average across all species. Models of species distribution
#' using only fine-resolution morphometric variables were most effective in (...)
#' Models with coarse resolution morphometric variables best predicted the distributions of (...)
#' Models with a combination mix of different all resolutions best explained the distributions of (...)
random_mods_species |>
  select(species, r2_topography_pure) |>
  mutate(r2_topography_pure = round(r2_topography_pure * 100, 2)) |>
  ungroup() |>
  arrange(-r2_topography_pure) |>
  mutate(group = factor(group, levels = c('fine', 'coarse', 'all'))) |>
  group_by(group) |>
  slice(1:10) |>
  mutate(label = paste0(species, ' (', r2_topography_pure, '%)')) |>
  select(label) |>
  group_split() |>
  map(~.x |> pull() |> paste(collapse = ', '))

#' On average, elevation alone explained more of the variation in species distribution
#' than the combination of five morphometric predictors in 43 species,
#' i.e., approximately half of all species (Fig. 4, Appendix S2).
#' On average, the best model using five morphometric predictors explained 15.2% of the variation,
#' vs. 14.1% explained by elevation alone.
random_mods_species |>
  group_by(species) |>
  arrange(-r2_topography_pure) |>
  slice(1) |>
  group_by(species) |>
  mutate(diff = round((r2_elevation_pure - r2_topography_pure) * 100, 2)) |>
  filter(0<diff) |> nrow()

#' On average across all species, the best model using five morphometric predictors
#' explained 15.7% of the variation, vs. 12.5% explained by elevation alone.
random_mods_species |>
  group_by(species) |>
  arrange(-r2_topography_pure) |>
  slice(1) |>
  ungroup() |>
  summarise_at(c('r2_topography_pure', 'r2_elevation_pure'), ~mean(round(.x*100, 2)))

#' Species where morphometric variables explained more than elevation included:
#' "Achillea erba-rotta subsp. moschata (25.93%), Oreochloa disticha (21.47%),
#' Trifolium pallescens (19.64%), Geum montanum (18.98%), Carex curvula (18.8%),
#' Jacobaea carniolica (18.42%), Primula glutinosa (17.8%), Salix helvetica (17.59%)"
random_mods_species |>
  group_by(species) |>
  arrange(-r2_topography_pure) |>
  slice(1) |>
  ungroup() |>
  mutate(diff = round((r2_topography_pure-r2_elevation_pure) * 100, 2)) |>
  select(species, diff, r2_elevation_pure, r2_topography_pure) |>
  arrange(-diff) |>
  slice(1:8) |>
  mutate(label = paste0(species, ' (', diff, '%)')) |> pull(label) |>
  paste(collapse = ', ')

#' In contrast species for which the distribution was better explained by elevation than topography
#' were mostly subalpine and nival species that dominate either
#' at very low or very high elevation in the study area, such as
#' "Rhododendron ferrugineum (23.2%), Avenella flexuosa (19.47%),
#' Poa laxa (19.33%), Festuca nigrescens (18.84%),
#' Leucanthemopsis alpina (11.18%), Saxifraga bryoides (10.69%),
#' Ranunculus glacialis (9.49%), Leontodon hispidus (8.54%)"
random_mods_species |>
  group_by(species) |>
  arrange(-r2_topography_pure) |>
  slice(1) |>
  ungroup() |>
  mutate(diff = round((r2_elevation_pure-r2_topography_pure) * 100, 2)) |>
  select(species, diff, r2_elevation_pure, r2_topography_pure) |>
  arrange(-diff) |>
  slice(1:8) |>
  mutate(label = paste0(species, ' (', diff, '%)')) |> pull(label) |>
  paste(collapse = ', ')

#' =======================================================================
random_mods_comprops |>
  select(species, r2_topography_pure) |>
  mutate(r2_topography_pure = round(r2_topography_pure * 100, 2)) |>
  ungroup() |>
  arrange(-r2_topography_pure) |>
  mutate(group = factor(group, levels = c('fine', 'coarse', 'all'))) |>
  group_by(group) |>
  slice(1:10) |>
  mutate(label = paste0(species, ' (', r2_topography_pure, '%)')) |>
  select(label) |>
  group_split() |>
  map(~.x |> pull() |> paste(collapse = ', '))

#'However, community attributes were distinct in their relation to morphometric variables
random_mods_comprops |>
  group_by(species) |>
  arrange(-r2_topography_pure) |>
  slice(1) |>
  ungroup() |>
  mutate(diff = round((r2_elevation_pure-r2_topography_pure) * 100, 2),
         r2_elevation_pure = round(r2_elevation_pure*100, 2),
  r2_topography_pure = round(r2_topography_pure*100, 2),
         rsq_full = round(rsq_full*100, 2),
  shared = rsq_full - (r2_elevation_pure + r2_topography_pure)) |>
  arrange(-diff) |>
  select(species, diff, r2_elevation_pure, r2_topography_pure, shared)

#' Values for table 1
read_csv('data\\Comprop_means.csv') |>
  rename_with(~gsub("_", ".", .x, fixed = TRUE)) |>
  summarise_all(list(mean = mean, sd = sd, min = min, max = max), na.rm = T) |>
  pivot_longer(everything()) |>
  separate(name, c('name', 'fn'), sep = '_') |>
  pivot_wider(names_from = fn)

#-------------------------------------------------------------------------
# END OF THE CODE S5
#-------------------------------------------------------------------------