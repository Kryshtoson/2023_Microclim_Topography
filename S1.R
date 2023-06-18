#-------------------------------------------------------------------------
# START OF THE CODE S1
#-------------------------------------------------------------------------
#' CALCULATION OF MORPHOMETRIC VARIABLES (M.V.)
#' author Krystof Chytry, krystof.chytry@gmail.com
# R version 4.3.0
# IDE Pycharm
library(exactextractr) # version 0.9.1
library(tictoc) # version 1.2
library(spatialEco) # version 1.3-6
library(sf) # version 1.0-12
library(readxl) # version 1.4.2
library(tidyverse) # version 2.0.0
#' dplyr v. 1.1.2
#' ggplot2 v. 3.4.2
#' purrr v. 1.0.1
#' readr v. 2.1.4
#' tibble v. 3.2.1
#' tidyr v 1.3.0
library(terra) # version 1.7-29

#' -------------------------------------------------------------------------
#' IMPORT
#' -------------------------------------------------------------------------
sites <- read_sf('raw_input\\MCdata_v2-0_GEO_vegetation.gpkg') %>%
  arrange(logger_ID)
#' original DTM of the study area, not supplied in the Appendix, on request: krystof.chytry@gmail.com
r <- rast('raw_input\\DTM_2017.tif')

#' -------------------------------------------------------------------------
#' ELEVATION OF PLOTS
#' -------------------------------------------------------------------------
tibble(
  logger_ID = sites$logger_ID,
  elevation = round(terra::extract(r, sites)[, 2])) %>%
  write_csv('data\\Elevation.csv')

#' -------------------------------------------------------------------------
#' CALCULATION OF LOCAL MORPHOMETRIC VARIABLES (M.V.)
#' - crop to cell's neighbourhood
#' - aggregation to coarser resolution
#' - calculation of M.V.s
#' - extraction of the value of the target cell
#' ||
#' - two-level loop:
#' - - first iteration goes resolution-wise and in the end stores values into separate csv files
#' - - second iteration goes plot-wise
#' -------------------------------------------------------------------------
for (res in c(1, 3, 7, 15, 19, 25, 31, 51, 71, 101, 151, 201, 301)) {
  ls <- list()
  cat(paste0('\n\n\n', res, ': '))

  for (id in sites$logger_ID) {
    cat(paste0(id, ','))
    one_site <- sites[sites$logger_ID == id,]
    if (res != 1) { one_r <- aggregate(crop(r, ext(one_site) + res * 3.5), res) }
    else { one_r <- crop(r, ext(one_site) + 3.5) }

    ls[[id]] <- data.frame(
      logger_ID = id,
      res = res,
      `PLC` = curvature(one_r, 'planform')[4, 4] |> as.numeric(),
      `PRC` = curvature(one_r, 'profile')[4, 4] |> as.numeric(),
      `TLC` = curvature(one_r, 'total')[4, 4] |> as.numeric(),
      `HLI` = suppressMessages(hli(one_r))[4, 4] |> as.numeric(),
      `SLP` = terrain(one_r, 'slope', unit = 'degrees')[4, 4] |> as.numeric(),
      `TP3` = tpi(one_r, s = 3)[4, 4] |> as.numeric(),
      `TP5` = tpi(one_r, s = 5)[4, 4] |> as.numeric(),
      `TR3` = tri(one_r, s = 3)[4, 4] |> as.numeric(),
      `TR5` = tri(one_r, s = 5)[4, 4] |> as.numeric(),
      `HSD` = shade(terrain(one_r, 'slope', unit = 'radians'), terrain(one_r, 'aspect', unit = 'radians'), 40, 180)[4, 4] |> as.numeric(),
      `VR3` = vrm(one_r, s = c(3, 3))[4, 4] |> as.numeric(),
      `VR5` = vrm(one_r, s = c(5, 5))[4, 4] |> as.numeric(),
      `SLR` = terrain(one_r, 'slope', unit = 'radians')[4, 4] |> as.numeric(),
      `ASP` = terrain(one_r, 'aspect', unit = 'degrees')[4, 4] |> as.numeric())
  }
  bind_rows(ls) %>% write_csv(paste0('data\\', res, '_topoindex.csv'))
}

#' -------------------------------------------------------------------------
#' LOCAL M.V. ALL IN ONE FILE
#' -------------------------------------------------------------------------

list.files('data', pattern = 'topoindex', full.names = T) |>
  as.list() |>
  map(read_csv) |>
  bind_rows() |>
  write_csv('data/Morphometric_parameters_meta.csv')

#' -------------------------------------------------------------------------
#' AGGREGATING DTM INTO COARSER RESOLUTIONS FOR REGIONAL M.V.
#' -------------------------------------------------------------------------
for (res in c(1, 3, 7, 15, 19, 25, 31, 51, 71, 101, 151, 201, 301)) {
  one_res <- aggregate(crop(r, ext(sites) + 2000), res)
  names(one_res) <- paste0('DEM_', res)
  writeRaster(one_res, paste0('meta_coarse_resolution\\DEM_', res, '.tif'), overwrite = T)
}

#' -------------------------------------------------------------------------
#' CALCULATION OF REGIONAL M.V. IN SAGA GIS!
#' -------------------------------------------------------------------------

#' -------------------------------------------------------------------------
#' READING REGIONAL M.V.
#' these variables were calculated in SAGA and supplied as .sgrd files
#' so here we read them one by one and extract values for individual plots
#' -------------------------------------------------------------------------
coarse_indices_all <- list.files('meta_coarse_resolution\\', pattern = '.sdat$', full.names = T)
coarse_indices <- coarse_indices_all[grepl('Flow|Sky|Insolation', coarse_indices_all)]
ls_coarse <- list()

for (index in coarse_indices) {
  r_meta <- raster::raster(index)

  ls_coarse[[index]] <- tibble(
    logger_ID = sites$logger_ID,
    value = extract(r_meta, sites), #exact_extract(r_meta, st_buffer(sites, 1), 'max'),
    name = paste0('M', gsub('\\b(\\pL)\\pL{1,}|.', '\\U\\1',
                            gsub('_', ' ', gsub('Diffuse', 'Blurred', names(r_meta))),
                            perl = TRUE)),
    res = res(r_meta)[[1]])
}

#' -------------------------------------------------------------------------
#' JOINING REGIONAL M.V. TO LOCAL ONES
#' -------------------------------------------------------------------------
read_csv('data\\Morphometric_parameters_meta.csv') |>
  left_join(bind_rows(ls_coarse) |>
              pivot_wider()) |>
  mutate(TWI = log((MFA / MFW) / tan(SLR))) |>
  left_join(read_csv('data\\Elevation.csv')) |>
  write_csv('data\\Morphometric_parameters.csv')

#' -------------------------------------------------------------------------
#' SPECIES DATA WIDE
#' In addition, species data are converted to wide form
#' -------------------------------------------------------------------------
read_xlsx('raw_input\\T5_species_data_v2-0.xlsx') %>%
  filter(TAXON != 'no vegetation') %>%
  left_join(read_excel('raw_input\\Species_data.xlsx') %>%
              dplyr::select(TAXON = species, species_epm)) %>%
  mutate(TAXON = species_epm,
         value = 1) %>%
  dplyr::select(logger_ID, name = TAXON, value) %>%
  pivot_wider(values_fn = max, values_fill = 0) %>%
  arrange(logger_ID) %>%
  write_csv('data\\species_data_wide.csv')

#' -------------------------------------------------------------------------
#' MORPHOMETRIC VARIABLES WIDE
#' -------------------------------------------------------------------------
read_csv('data\\Morphometric_parameters.csv') %>%
  pivot_longer(-c(logger_ID, res, elevation)) %>%
  mutate(name = paste0(name, '_', res)) %>%
  dplyr::select(-res) %>%
  pivot_wider() %>%
  write_csv('data\\Morphometric_parameters_bank.csv')

#-------------------------------------------------------------------------
# END OF THE CODE S1
#-------------------------------------------------------------------------