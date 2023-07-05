#-------------------------------------------------------------------------
# START OF SCRIPT S4
#-------------------------------------------------------------------------
library(tictoc)
tic()
library(DSSAT)
library(tidyverse)
library(furrr)
library(readxl)
library(ggrepel)
options(warn = -1)

spe <- read_csv('data\\species_data_wide.csv')
head <- read_csv('data\\Morphometric_parameters_bank.csv')
comprop <- read_csv('data\\Comprop_means.csv')
dir.create('results\\meta_results\\random_models')
# -------------------------------------------------------------------------
# CREATING RANDOM PREDICTORS FOR MODELS
# -------------------------------------------------------------------------
selected_combinations <- gsub(' ', '', apply(
  expand.grid(
    c('SLP', 'PCR', 'PLC', 'TLC', 'TP5', 'TR5', 'VR3', 'HLI', 'SRI', 'TWI', 'MDI'),
    c(1, 3, 7, 15, 19, 25, 31, 51, 71, 101, 151, 201, 301)),
  1, paste, sep = '', collapse = '_'))

head <- head %>% select('logger_ID', 'elevation', contains(selected_combinations))

indexbank <- tibble(index = names(head)[-(1:2)]) %>%
  mutate(path = index) %>%
  separate(index, c('index', 'resolution'), sep = '_') %>%
  drop_na() %>%
  mutate(resolution = as.numeric(resolution)) %>%
  arrange(index, resolution)

indexbank_fine <- indexbank %>%
  filter(resolution %in% 1:19) %>%
  nest(data = path)
unique(indexbank_fine$resolution)

indexbank_coarse <- indexbank %>%
  filter(resolution %in% c(51, 71, 101, 151, 201)) %>%
  nest(data = path)
unique(indexbank_coarse$resolution)

predictors <- function(source) {
  bank <- list()
  for (i in 1:1000) {
    set.seed(i)
    x <- source %>% unnest() %>% sample_n(5)
    bank[[i]] <- formula(
      paste('value ~ ',
            paste(
              c(paste0('poly(', x$path, ', 2)')), collapse = ' + ')
      )
    )
  }
  bank <- bank |> discard(is.null)
}

bank_fine <- predictors(source = indexbank_fine)
bank_coarse <- predictors(source = indexbank_coarse)
bank_all <- predictors(source = indexbank)

# -------------------------------------------------------------------------
# SUBSETTING SPECIES USED FOR THE ANALYSIS
# -------------------------------------------------------------------------
`%!in%` <- Negate(`%in%`)
fairly_common <- sort(names(spe[-1])[colSums(spe[-1]) > 50])
fairly_common <- fairly_common[fairly_common %!in% (read_excel('raw_input\\Species_data.xlsx') %>% filter(!keep))$species_epm]

# -------------------------------------------------------------------------
# ACTUAL MODELS, SPECIES (MOST-TIME-CONSUMING CHUNK)
# -------------------------------------------------------------------------
ls_outter <- list()
ls <- list()
v <- c(1, (1:10) * 100)

d2 <- function(m) {
  devs <- (broom::glance(m))
  1 - (devs$deviance[1] / devs$null.deviance[1])
}

banks <- tibble(bank_all,
                bank_fine,
                bank_coarse)

for (y in 1:3) {
  cat(paste0('\ngroup = ', y, '/3:'))
  file <- gsub('bank_', '', names(banks)[y])
  bank_meta <- banks[[y]]
  plan(multisession, workers = availableCores() - 1)
  for (j in 2:length(v)) {
    cat(paste0('\n', j-1, '/', length(v), ': ', Sys.time()))
    bank_step <- bank_meta[v[j - 1]:v[j]]

    for (i in fairly_common) {
      # print(i)
      if (!is.null(head$value)) {
        head <- head %>%
          dplyr::select(-value)
      }

      head <- head %>%
        left_join(spe %>%
                    dplyr::select(logger_ID, value = all_of(i)),
                  by = 'logger_ID') %>%
        mutate(value = replace_na(value, 0))

      bank_step_full <- map(bank_step, function(f) {
        as.formula(paste(paste(deparse(f, width.cutoff = 500), collapse = ""), ' + poly(elevation, 2)')) })

      ls[[i]] <- tibble(index = v[j - 1]:v[j],
                        formulas = bank_step,
                        formulas_full = bank_step_full) %>%
        mutate(species = i) %>%
        mutate(mods = formulas %>%
          future_map(function(fm) {
            glm(fm, data = head, family = binomial())
          }),
               mods_full = formulas_full %>%
                 future_map(function(fm) {
                   glm(fm, data = head, family = binomial())
                 })) %>%
        mutate(rsq = mods %>%
          map_dbl(d2),
               rsq_full = mods_full %>%
                 map_dbl(d2))
    }
    bind_rows(ls) %>%
      mutate(formulas = formulas %>%
        map_chr(function(c)as.character(c)[3])) %>%
      select(index, formulas, species, rsq, rsq_full) %>%
      write_csv(paste0('results\\meta_results\\random_models\\SPECIES_randmods_', file, '_', j, '.csv'))
  }
}

map2(
  list.files('results\\meta_results\\random_models', pattern = 'SPECIES_randmods', full.names = T),
  gsub('SPECIES_randmods_', '', list.files('results\\meta_results\\random_models', pattern = 'SPECIES_randmods')),
  function(path, name) {
    read_csv(path) |>
      mutate(
        group = gsub('\\.csv', '', name),
        group = gsub('\\d', '', group),
        group = gsub('_', '', group))
  }
) |>
  bind_rows() |>
  write_csv('results\\SPECIES_random_models.csv')

# -------------------------------------------------------------------------
# ACTUAL MODELS, SPECIES (FASTER THAN SPECIES)
# -------------------------------------------------------------------------
plan(sequential)
plan(multisession, workers = availableCores() - 1)
ls <- list()
for (y in 1:3) {
  cat(paste0('\ngroup = ', y, '/3:'))
  file <- gsub('bank_', '', names(banks)[y])
  bank_step <- banks[[y]]

  for (i in c('sp_richness',
              'Shannon',
              'Temperature',
              'Moisture',
              'Reaction',
              'Nutrients',
              'trait_height',
              'trait_sla',
              'trait_N',
              'trait_suc')) {

    print(i)
    if (!is.null(head$value)) {
      head <- head %>%
        dplyr::select(-value)
    }

    bank_step_full <- map(bank_step, function(f) {
      as.formula(paste(paste(deparse(f, width.cutoff = 500), collapse = ""), ' + poly(elevation, 2)')) })

    head_meta <- head %>%
      left_join(comprop %>%
                  dplyr::select(logger_ID, value = all_of(i)),
                by = 'logger_ID') %>%
      filter(!is.na(value))

    ls[[i]] <- tibble(formulas = bank_step,
                      formulas_full = bank_step_full) %>%
      mutate(species = i) %>%
      mutate(mods = formulas %>%
        future_map(function(fm) {
          glm(fm, data = head_meta, family = gaussian())
        }),
             mods_full = formulas_full %>%
               future_map(function(fm) {
                 glm(fm, data = head_meta, family = gaussian())
               })) %>%
      mutate(rsq = mods %>%
        map_dbl(d2),
             rsq_full = mods_full %>%
               map_dbl(d2))
  }
  bind_rows(ls) %>%
    mutate(formulas = formulas %>%
      map_chr(function(c)as.character(c)[3])) %>%
    select(formulas, species, rsq, rsq_full) %>%
    write_csv(paste0('results\\meta_results\\random_models\\COMPROPS_randmods_', file, '.csv'))
}

map2(
  list.files('results\\meta_results\\random_models', pattern = 'COMPROPS_randmods', full.names = T),
  gsub('COMPROPS_randmods_', '', list.files('results\\meta_results\\random_models', pattern = 'COMPROPS_randmods')),
  function(path, name) {
    read_csv(path) |>
      mutate(
        group = gsub('\\.csv', '', name),
        group = gsub('\\d', '', group),
        group = gsub('_', '', group))
  }
) |>
  bind_rows() |>
  write_csv('results\\COMPROPS_random_models.csv')

# -------------------------------------------------------------------------
# Models with elevation as single predictor (NULL models)
# -------------------------------------------------------------------------
comprop %>%
  pivot_longer(-1) %>%
  left_join(head %>%
              select(logger_ID, elevation)) %>%
  group_by(name) %>%
  nest() %>%
  filter(name %in% c('sp_richness',
                     'Shannon',
                     'Temperature',
                     'Moisture',
                     'Reaction',
                     'Nutrients',
                     'trait_height',
                     'trait_sla',
                     'trait_N',
                     'trait_suc')) %>%
  mutate(mod = data %>%
    map(function(df) {
      glm(value ~ poly(elevation, 2),
          data = df,
          family = gaussian()) }),
         elevation_only = mod %>%
           map_dbl(d2),
         species = name) %>%
  ungroup() %>%
  dplyr::select(species, elevation_only) %>%
  write_csv('results\\COMPROPS_elevation_models.csv')

spe %>%
  pivot_longer(-1) %>%
  left_join(head %>%
              select(logger_ID, elevation),
            by = 'logger_ID') %>%
  group_by(name) %>%
  nest() %>%
  filter(name %in% fairly_common) %>%
  mutate(mod = data %>%
    map(function(df) {
      glm(value ~ poly(elevation, 2),
          data = df,
          family = binomial()) }),
         elevation_only = mod %>%
           map_dbl(d2),
         species = name) %>%
  ungroup() %>%
  select(species, elevation_only) %>%
  write_csv('results\\SPECIES_elevation_models.csv')
toc()
#-------------------------------------------------------------------------
# END OF SCRIPT S4
#-------------------------------------------------------------------------
