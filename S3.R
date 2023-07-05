#-------------------------------------------------------------------------
# START OF THE CODE S3A
#-------------------------------------------------------------------------
#' RESOLUTION-WISE MODELS (M.V.)
#' author Krystof Chytry, krystof.chytry@gmail.com
# R version 4.3.0
# IDE Pycharm
library(tictoc) # version 1.2
tic() # ~84 sec
library(ggrepel) # version 0.9.3
library(broom) # version 1.0.4
library(sf)
library(furrr) # version 0.3.1
library(tidyverse)
library(readxl)
library(DSSAT)
library(gridExtra) # version 2.3

sites <- read_sf('raw_input\\MCdata_v2-0_GEO_vegetation.gpkg') %>%
  arrange(logger_ID)
veg_data <- read_csv('data\\species_data_wide.csv')
fairly_common <- names(veg_data)[colSums(veg_data) > 50][-1]
bank <- read_csv('data\\Morphometric_parameters_bank.csv')
comprop <- read_csv('data\\Comprop_means.csv')

indexbank <- tibble(index = names(bank)) %>%
  mutate(path = index) %>%
  filter(index != 'logger_ID' & index != 'elevation') %>%
  separate(index, c('index', 'resolution'), sep = '_') %>%
  mutate(resolution = as.numeric(resolution)) %>%
  arrange(index, resolution) %>%
  nest(data = path)

#-------------------------------------------------------------------------
# multithread implementation of species specific models
#-------------------------------------------------------------------------
d2 <- function(m) {
  devs <- (broom::glance(m))
  1 - (devs$deviance[1] / devs$null.deviance[1])
}

plan(multisession, workers = availableCores() - 1)
ls <- list()
for (i in fairly_common) {
  print(i)

  if (!is.null(bank$value)) {
    bank <- bank %>%
      dplyr::select(-value)
  }

  bank <- bank %>%
    left_join(veg_data %>%
                dplyr::select(logger_ID, value = all_of(i))) %>%
    mutate(value = replace_na(value, 0))

  elevation <- glm(value ~ elevation + I(elevation^2), data = bank, family = binomial()) %>% d2

  ls[[i]] <- indexbank %>%
    mutate(out = data %>%
      future_map(function(df) {
        m1 <- glm(formula(
          paste('value', '~', df[[1, 1]], '* elevation + I(',
                df[[1, 1]], '^2) + I(elevation^2)')),
                  data = bank,
                  family = binomial())
        tibble(species = i,
               deviance = d2(m1),
               null_mod = elevation) })) %>%
    unnest(out) %>%
    dplyr::select(species, deviance, null_mod, index, resolution)
}

plan(sequential)
species_out <- bind_rows(ls) %>%
  mutate(dev_explained = deviance - null_mod) %>%
  mutate(kind = 'species')
species_out %>% write_csv('results\\meta_results\\models_resolution_species.csv')
species_out$index |> unique()
# -------------------------------------------------------------------------
# multithread implementation of trait specific models
#-------------------------------------------------------------------------
d2 <- function(m) {
  devs <- (broom::glance(m))
  1 - (devs$deviance[1] / devs$null.deviance[1])
}

plan(multisession, workers = availableCores() - 1)
ls <- list()
for (i in c('sp_richness', 'Shannon', 'Temperature', 'Moisture', 'Reaction',
            'Nutrients', 'trait_height', 'trait_sla', 'trait_N', 'trait_suc')) {
  print(i)
  if (!is.null(bank$value)) {
    bank <- bank %>%
      dplyr::select(-value)
  }

  bank <- bank %>%
    left_join(comprop %>%
                dplyr::select(logger_ID, value = all_of(i))) %>%
    filter(!is.na(value))

  elevation <- glm(value ~ elevation + I(elevation^2),
                   data = bank, family = gaussian()) %>% d2

  ls[[i]] <- indexbank %>%
    mutate(out = data %>%
      future_map(function(df) {
        m1 <- glm(formula(
          paste('value', '~', df[[1, 1]], '* elevation + I(',
                df[[1, 1]], '^2) + I(elevation^2)')),
                  data = bank,
                  family = gaussian())
        tibble(species = i,
               deviance = d2(m1),
               null_mod = elevation) })) %>%
    unnest(out) %>%
    dplyr::select(species, deviance, null_mod, index, resolution)
}

plan(sequential)

species_out <- bind_rows(ls) %>%
  mutate(dev_explained = deviance - null_mod) %>%
  mutate(kind = 'species')
species_out %>% write_csv('results\\meta_results\\models_resolution_comprops.csv')

# -------------------------------------------------------------------------
# Species plotting
# -------------------------------------------------------------------------
df <- read_csv('results\\meta_results\\models_resolution_species.csv') %>%
  mutate_cond(dev_explained < 0, dev_explained = NA) %>%
  mutate(index2 = fct_recode(factor(index, levels = c('SLP', 'PLC', 'PRC', 'TLC', 'TP3',
                                                      'TP5', 'TR3', 'TR5', 'VR3', 'VR5',
                                                      'HLI', 'MSVF', 'MVS', 'MBI', 'MDI',
                                                      'HSD', 'TWI')),
                             Slope = 'SLP',
                             `Planform Curvature` = 'PLC',
                             `Profile Curvature` = 'PRC',
                             `Total Curvature` = 'TLC',
                             `Topographic position index 3×3` = 'TP3',
                             `Topographic position index 5×5` = 'TP5',
                             `Terrain Ruggedness Index 3×3` = 'TR3',
                             `Terrain Ruggedness Index 5×5` = 'TR5',
                             `Vector Ruggedenss Measure 3×3` = 'VR3',
                             `Vector Ruggedenss Measure 5×5` = 'VR5',
                             `Head Load Index` = 'HLI',
                             `Sky View Factor` = 'MSVF',
                             `Visible Sky` = 'MVS',
                             `Diffuse Solar radiation` = 'MBI',
                             `Direct Solar radiation` = 'MDI',
                             `Hill Shade 180°` = 'HSD',
                             `Topographic Wetness Index` = 'TWI')) %>%
  filter(index2 != index)

df %>% write_csv('results\\SPECIES_resolution_models.csv')

df %>%
  group_by(index2, resolution) %>%
  summarise(r2 = mean(dev_explained)) -> df_m

df %>%
  distinct(species, null_mod) %>%
  summarise(mean(null_mod)) %>%
  as.numeric() %>%
  round(., digits = 4) -> expl_elevation

df_m %>% ggplot(aes(resolution, r2)) +
  geom_path(aes(colour = index2), size = 1.2) +
  geom_text_repel(data = df_m %>% filter(resolution == 301),
                  aes(label = index2,
                      colour = index2),
                  force = 5,
                  nudge_x = 300,
                  direction = "y",
                  hjust = 0,
                  segment.size = 0.5,
                  segment.colour = 'grey20',
                  segment.alpha = .5) +
  expand_limits(x = c(0, 450)) +
  guides(colour = 'none') +
  theme_bw() +
  scale_x_continuous(breaks = c(1, 3, 7, 15, 19, 25, 31, 51, 71, 101, 151, 201, 301),
                     labels = c(1, ' ', '   ', 15, ' ', ' ', 31, 51, 71, 101, 151, 201, 301),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.01),
                     labels = paste0(seq(0, 0.1, by = 0.01) * 100, '%'),
                     expand = c(0, 0, 0, .005)) +
  labs(y = 'Mean explained variation', x = 'Spatial resolution') +
  geom_hline(yintercept = 0, colour = 'red') +
  theme(text = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid = element_line(colour = 'grey85'))

ggsave('results\\figures\\Figure1A.png', height = 6, width = 9)

#-------------------------------------------------------------------------
# Community attributes plotting
#-------------------------------------------------------------------------

cf <- read_csv('results\\meta_results\\models_resolution_comprops.csv') %>%
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
                              "Shannon's Diversity" = 'Shannon')) %>%
  mutate_cond(dev_explained < 0, dev_explained = NA) %>%
  mutate(index2 = fct_recode(factor(index, levels = c('SLP', 'PLC', 'PRC', 'TLC', 'TP3',
                                                      'TP5', 'TR3', 'TR5', 'VR3', 'VR5',
                                                      'HLI', 'MSVF', 'MVS', 'MBI', 'MDI',
                                                      'HSD', 'TWI', 'MFA')),
                             Slope = 'SLP',
                             `Planform Curvature` = 'PLC',
                             `Profile Curvature` = 'PRC',
                             `Total Curvature` = 'TLC',
                             `Topographic position index 3×3` = 'TP3',
                             `Topographic position index 5×5` = 'TP5',
                             `Terrain Ruggedness Index 3×3` = 'TR3',
                             `Terrain Ruggedness Index 5×5` = 'TR5',
                             `Vector Ruggedenss Measure 3×3` = 'VR3',
                             `Vector Ruggedenss Measure 5×5` = 'VR5',
                             `Head Load Index` = 'HLI',
                             `Sky View Factor` = 'MSVF',
                             `Visible Sky` = 'MVS',
                             `Diffuse Solar radiation` = 'MBI',
                             `Direct Solar radiation` = 'MDI',
                             `Hill Shade 180°` = 'HSD',
                             `Topographic Wetness Index` = 'TWI')) %>%
  filter(index2 != index)
cf %>% write_csv('results\\COMPROPS_resolution_models.csv')
cf %>%
  group_by(index2, resolution) %>%
  summarise(r2 = mean(dev_explained)) -> cf_m

cf %>%
  distinct(species, null_mod) %>%
  summarise(mean(null_mod)) %>%
  as.numeric() %>%
  round(., digits = 4) -> cf_expl_elevation

cf_m %>% ggplot(aes(resolution, r2)) +
  geom_path(aes(colour = index2), size = 1.2) +
  geom_text_repel(data = cf_m %>% filter(resolution == 301),
                  aes(label = index2,
                      colour = index2),
                  force = 5,
                  nudge_x = 300,
                  direction = "y",
                  hjust = 0,
                  segment.size = 0.5,
                  segment.colour = 'grey20',
                  segment.alpha = .5) +
  expand_limits(x = c(0, 450)) +
  guides(colour = 'none') +
  theme_bw() +
  scale_x_continuous(breaks = c(1, 3, 7, 15, 19, 25, 31, 51, 71, 101, 151, 201, 301),
                     labels = c(1, ' ', '   ', 15, ' ', ' ', 31, 51, 71, 101, 151, 201, 301),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.01),
                     labels = paste0(seq(0, 0.1, by = 0.01) * 100, '%'),
                     expand = c(0, 0, 0, .005)) +
  labs(y = 'Mean explained variation', x = 'Spatial resolution') +
  geom_hline(yintercept = 0, colour = 'red') +
  theme(text = element_text(size = 18),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid = element_line(colour = 'grey85'))

ggsave('results\\figures\\Figure1B.png', height = 6, width = 9)

# -------------------------------------------------------------------------
# SPECIES SPECIFIC CURVERS (for appendix)
# -------------------------------------------------------------------------

ls <- list()
for (i in sort(unique(df$species))) { print(i)
  df_step <- df %>% filter(species == i)

  ls[[i]] <- df_step %>% ggplot(aes(resolution, dev_explained)) +
    geom_path(aes(colour = index2), size = 1.2) +
    geom_text_repel(data = df_step %>% filter(resolution == 301),
                    aes(label = index2,
                        colour = index2),
                    force = 5,
                    nudge_x = 300,
                    direction = "y",
                    hjust = 0,
                    segment.size = 0.5,
                    segment.colour = 'grey20',
                    segment.alpha = .5
    ) +
    expand_limits(x = c(0, 450)) +
    guides(colour = 'none') +
    theme_bw() +
    scale_x_continuous(breaks = c(1, 3, 7, 15, 19, 25, 31, 51, 71, 101, 151, 201, 301),
                       labels = c(1, ' ', '   ', 15, ' ', ' ', 31, 51, 71, 101, 151, 201, 301),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.01),
                       labels = paste0(seq(0, 1, by = 0.01) * 100, '%'),
                       expand = c(0, 0, 0, .005)) +
    labs(y = 'Explained deviation', x = 'Spatial resolution',
         title = i,
         caption = paste0('Deviation explained by elevation: ',
                          paste0(round(mean(df_step$null_mod) * 100, 2), '%'))) +
    geom_hline(yintercept = 0, colour = 'red') +
    theme(text = element_text(size = 18),
          plot.title = element_text(face = 'italic'),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid = element_line(colour = 'grey85'))
}

m <- marrangeGrob(ls, ncol = 1, nrow = 1)
ggsave('results\\appendix_figures\\Appendix_S3A.pdf', m, height = 6, width = 9)


# -------------------------------------------------------------------------
# COMMUNITY ATTRIBUTES SPECIFIC CURVERS (for appendix)
# -------------------------------------------------------------------------

ls <- list()
for(i in unique(cf$species)){ print(i)
df_step <- cf %>% filter(species == i)
ls[[i]] <- df_step %>%  ggplot(aes(resolution, dev_explained))+
geom_path(aes(colour = index2), size = 1.2) +
geom_text_repel(data = df_step %>% filter(resolution == 301),
aes(label = index2,
colour = index2),
force        = 5,
nudge_x      = 300,
direction    = "y",
hjust        = 0,
segment.size = 0.5,
segment.colour = 'grey20',
segment.alpha = .5
) +
expand_limits(x = c(0, 450)) +
guides(colour = 'none') +
theme_bw()+
scale_x_continuous(breaks = c(1, 3, 7, 15, 19, 25, 31, 51, 71, 101, 151, 201, 301),
labels = c(1, ' ', '   ', 15, ' ', ' ', 31, 51, 71, 101, 151, 201, 301),
expand = c(0, 0)) +
scale_y_continuous(breaks = seq(0, 1, by = 0.01),
labels = paste0(seq(0, 1, by = 0.01)*100, '%'),
expand = c(0, 0, 0, .005)) +
labs(y = 'Explained deviation', x = 'Spatial resolution',
title = i,
caption =paste0('Deviation explained by elevation: ',
paste0(round(mean(df_step$null_mod)*100, 2), '%'))) +
geom_hline(yintercept = 0, colour = 'red') +
theme(text = element_text(size = 18),
plot.title = element_text(face = 'italic'),
panel.grid.minor = element_blank(),
panel.grid.major.y = element_blank(),
panel.grid = element_line(colour = 'grey85'))
}
m <- marrangeGrob(ls, ncol = 1, nrow = 1)
ggsave('results\\appendix_figures\\Appendix_S3B.pdf', m, height = 6, width = 9)

toc()
#------------------------------------------------------------------------
# END OF SCRIPT S3
#------------------------------------------------------------------------
