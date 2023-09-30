#-------------------------------------------------------------------------
# START OF SCRIPT S4
#-------------------------------------------------------------------------
library(ggforce)
library(DSSAT)
library(tidyverse)
library(ggrepel)

# -------------------------------------------------------------------------
# DATA TO PLOT
# -------------------------------------------------------------------------
# SPECIES
# -------------------------------------------------------------------------
read_csv('results\\SPECIES_random_models.csv') |>
  left_join(read_csv(paste0('results\\SPECIES_elevation_models.csv'))) |>
  mutate(r2_elevation_pure = rsq_full - rsq,
         r2_topography_pure = rsq_full - elevation_only) |>
  arrange(-r2_topography_pure) |>
  group_by(species, group) |>
  slice(1) |>
  select(species, name = group, elevation_only, value = r2_topography_pure) |>
  pivot_wider() |>
  mutate(gr = fine > coarse) |>
  mutate(gr = factor(gr, levels = c('FALSE', 'TRUE'),
                     labels = c('(A)', '(B)'))) |>
  pivot_longer(c(all, coarse, fine),
               names_to = 'group', values_to = 'r2_topography_pure') |>
  ungroup() |>
  mutate(species = fct_reorder(.f = species, .x = r2_topography_pure, .fun = max)) -> SPECIES_step
list('all' = SPECIES_step |> filter(group == 'all'),
     'butall' = SPECIES_step |> filter(group != 'all')) -> SPECIES_out

# -------------------------------------------------------------------------
# DATA TO PLOT
# -------------------------------------------------------------------------
# COMPROPS
# -------------------------------------------------------------------------
read_csv('results\\COMPROPS_random_models.csv') |>
  left_join(read_csv(paste0('results\\COMPROPS_elevation_models.csv'))) |>
  mutate(species = fct_recode(species,
                              'EIV Moisture' = 'Moisture',
                              'EIV Reaction' = 'Reaction',
                              'EIV Nutrients' = 'Nutrients',
                              'EIV Temperature' = 'Temperature',
                              'CM Succulency' = 'trait_suc',
                              'CM Vegetative Height' = 'trait_height',
                              'CM Leaf Nitrogen' = 'trait_N',
                              'CM Specific Leaf Area' = 'trait_sla',
                              'Species Richness' = 'sp_richness',
                              "Shannon's Diversity" = 'Shannon')) |>
  mutate(r2_elevation_pure = rsq_full - rsq,
         r2_topography_pure = rsq_full - elevation_only) |>
  arrange(-r2_topography_pure) |>
  group_by(species, group) |>
  slice(1) |>
  select(species, name = group, value = r2_topography_pure) |>
  pivot_wider() |>
  mutate(gr = fine > coarse) |>
  mutate(gr = factor(gr, levels = c('FALSE', 'TRUE'),
                     labels = c('(A)', '(B)'))) |>
  pivot_longer(c(all, coarse, fine),
               names_to = 'group', values_to = 'r2_topography_pure') |>
  ungroup() |>
  mutate(species = fct_reorder(.f = species, .x = r2_topography_pure, .fun = max)) -> COMPROPS_step

list('all' = COMPROPS_step |> filter(group == 'all'),
     'butall' = COMPROPS_step |> filter(group != 'all')) -> COMPROPS_out

dem <- theme(
  axis.text.y = element_text(face = 'italic'),
  axis.title.y = element_blank(),
  axis.text = element_text(size = 12),
  axis.text.x = element_text(hjust = .3),
  plot.title = element_text(size = 18),
  strip.background = element_blank(),
  legend.position = c(1, 1),
  legend.background = element_blank(),
  legend.justification = c(1, 1),
  strip.placement = 'above',
  strip.text = element_text(hjust = 0, size = 12),
  legend.title = element_blank(),
  legend.text = element_text(size = 14))

# -------------------------------------------------------------------------
# PLOTTING
# -------------------------------------------------------------------------
# SPECIES
# -------------------------------------------------------------------------
SPECIES_out$butall |>
  mutate_cond(r2_topography_pure < 0, r2_topography_pure = 0) |>
  mutate(group = factor(group, levels = c('coarse', 'fine', 'combined'))) |>
  ggplot(aes(species, r2_topography_pure)) +
  facet_col(gr ~ ., scales = 'free_y', space = 'free', strip.position = "top") +
  coord_flip() +
  geom_bar(aes(fill = group),
           stat = 'identity',
           position = position_dodge(width = .8)) +
  geom_point(data = SPECIES_out$all, shape = 23, fill = 'red', size = 4) +
  geom_hline(yintercept = 0, colour = 'blue') +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(breaks = c(0, .05, .1, 0.15, 0.2, 0.25, 0.3, .35, .4, .45, .5, .55),
                     labels = c('0%', '5%', '10%', '15%', '20%', '25%', '30%', '35%', '40%', '45%', '50%', '55%'),
                     expand = c(0, 0, .03, 0)) +
  scale_fill_manual(values = c('#51718C', '#F2B33D', 'red'), drop = F) +
  labs(y = 'Variation explained by topography') +
  theme_bw() +
  dem

ggsave('results//figures//Figure2.png', height = 17, width = 9.5)

# -------------------------------------------------------------------------
# PLOTTING
# -------------------------------------------------------------------------
# COMPROPS
# -------------------------------------------------------------------------
COMPROPS_out$butall |>
  mutate_cond(r2_topography_pure < 0, r2_topography_pure = 0) |>
  mutate(group = factor(group, levels = c('coarse', 'fine', 'combined'))) |>
  filter(group != 'all') |>
  ggplot(aes(species, r2_topography_pure)) +
  coord_flip() +
  geom_bar(aes(fill = group), stat = 'identity', position = position_dodge(width = .8)) +
  geom_point(data = COMPROPS_out$all, shape = 23, fill = 'red', size = 4) +
  geom_hline(yintercept = 0, colour = 'blue') +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(breaks = c(0, .05, .1, 0.15, 0.2, 0.25, 0.3, .35, .4, .45, .5, .55),
                     labels = c('0%', '5%', '10%', '15%', '20%', '25%', '30%', '35%', '40%', '45%', '50%', '55%'),
                     expand = c(0, 0, .03, 0)) +
  scale_fill_manual(values = c('#51718C', '#F2B33D', 'red'), drop = F) +
  # expand_limits(y = c(-0.022, .34)) +
  labs(y = 'Variation explained by topography') +
  theme_bw() +
  dem

ggsave('results//figures//Figure3.png', height = 5, width = 9.5)

# -------------------------------------------------------------------------
# PLOTTING
# -------------------------------------------------------------------------
# USING THE BEST MODELS TO PLOT FIGURE 4
# -------------------------------------------------------------------------
shcut <- function(names) {
  names <- lapply(strsplit(names, "\\ "), function(x) substring(x, 1, 1))
  names <- lapply(names,
                  function(x) { paste(unlist(x)[1], unlist(x)[2], sep = "") })
  unlist(names)
}


read_csv('results\\SPECIES_random_models.csv') |>
  left_join(read_csv(paste0('results\\SPECIES_elevation_models.csv'))) |>
  mutate(r2_elevation_pure = rsq_full - rsq,
         r2_topography_pure = rsq_full - elevation_only) |>
  arrange(-r2_topography_pure) |>
  group_by(species, group) |>
  mutate(macro = 'Distribution of species') |>
  slice(1) -> best_models_species

read_csv('results\\COMPROPS_random_models.csv') |>
  left_join(read_csv(paste0('results\\COMPROPS_elevation_models.csv'))) |>
  mutate(r2_elevation_pure = rsq_full - rsq,
         r2_topography_pure = rsq_full - elevation_only) |>
  arrange(-r2_topography_pure) |>
  group_by(species, group) |>
  mutate(macro = 'Community attributes') |>
  slice(1) -> best_models_comprops

bind_rows(best_models_species, best_models_comprops) |>
  mutate(group = factor(group,
                        levels = c('fine', 'coarse'),
                        labels = c('Fine resolution', 'Coarse resolution'))) |>
  filter(group %in% c('Coarse resolution', 'Fine resolution')) |>
  mutate(species_abb = shcut(species)) -> step

step |>
  ggplot(aes(r2_topography_pure, r2_elevation_pure)) +
  geom_abline(slope = 1, linetype = "dashed", color = "Red") +
  geom_point() +
  geom_text_repel(data = step |>
    filter(macro != 'Distribution of species') |>
    mutate(species = fct_recode(species,
                                'EIV Moisture' = 'Moisture',
                                'EIV Reaction' = 'Reaction',
                                'EIV Nutrients' = 'Nutrients',
                                'EIV Temperature' = 'Temperature',
                                'CM Suculency' = 'trait_suc',
                                'CM Vegetative Height' = 'trait_height',
                                'CM Leaft Nitrogen' = 'trait_N',
                                'CM Specific Leaf Area' = 'trait_sla',
                                'Species Richness' = 'sp_richness',
                                "Shannon's Diversity" = 'Shannon')),
                  aes(label = species)) +
  geom_text_repel(data = step |>
    filter(macro == 'Distribution of species') |>
    filter(species %in% c('Gnaphalium supinum',
                          'Deschampsia cespitosa',
                          'Geum reptans',
                          'Vaccinium vitis-idaea',
                          'Cirsium spinossisimum',
                          'Festuca melanopsis',
                          'Helictochloa versicolor',
                          'Primula minima',
                          'Agrostis schraderiana',
                          'Kobresia myosuroides',
                          'Oreochloa disticha',
                          'Poa laxa',
                          'Ranunculus glacialis',
                          'Avenula flexuosa',
                          'Saxifraga bryoides',
                          'Achillea erba-rotta subsp. moschata',
                          'Festuca nigrescens',
                          'Avenulla flexuosa',
                          'Rhododendron ferrugineum',
                          'Carex curvula')),
                  colour = 'red', fontface = 'italic',
                  max.overlaps = Inf,
                  aes(label = species_abb)) +
  theme_bw() +
  expand_limits(x = c(0, .65), y = c(0, .65)) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = c(0:7 / 10),
                     labels = paste0((0:7) * 10, '%')) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(0:7 / 10),
                     labels = paste0((0:7) * 10, '%')) +
  #geom_text_repel(aes(label = species)) +
  facet_grid(macro ~ group) +
  labs(y = 'Variation explained by elevation',
       x = 'Variation explained by topography') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        strip.text = element_text(hjust = 0, size = 14),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14))

ggsave('results//figures//Figure4.png', width = 10.5, height = 9)

# -------------------------------------------------------------------------
# BOOTSTRAPING
# -------------------------------------------------------------------------
# species
# -------------------------------------------------------------------------
set.seed(1)
replicate(999, SPECIES_out$butall |>
  select(species, name = group, value = r2_topography_pure) |>
  pivot_wider(), simplify = F) |>
  map(~.x |>
    sample_n(10, replace = T) |>
    summarise_at(c('coarse', 'fine'), mean) |>
    mutate(result = fine / coarse) |>
    pull(result)
  ) |>
  unlist() -> SPECIES_ratios

set.seed(1)
replicate(999, COMPROPS_out$butall |>
  select(species, name = group, value = r2_topography_pure) |>
  pivot_wider(), simplify = F) |>
  map(~.x |>
    sample_n(10, replace = T) |>
    summarise_at(c('coarse', 'fine'), mean) |>
    mutate(result = fine / coarse) |>
    pull(result)
  ) |>
  unlist() -> COMPROPS_ratios

t.test(COMPROPS_ratios - SPECIES_ratios, mu = 0, alternative = 'greater')
quantile(COMPROPS_ratios - SPECIES_ratios, .05)
quantile(COMPROPS_ratios - SPECIES_ratios, .13085)
# t = 37.131, df = 998, p-value <0.001, 95% confidence interval > 0.104

#------------------------------------------------------------------------
# END OF SCRIPT S5
#------------------------------------------------------------------------
