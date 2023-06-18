# Chytrý et al. Ecography: *Limited impact of microtopography on alpine plant distribution*
## Appendix S1: R Code and the data

### Structure of the project folder

* `S1.R`: Preparation of morphometric variables.
* `S2.R`: Manipulation of the species data.
* `S3.R`: Models of each morphometric variable over spatial resolutions.
* `S4.R`: Models using combinations of morphometric variables as predictors.
* `S5.R`: Extracting information used in the paper itself.


* `raw_input`
  * Original data. 
  * `DTM.tif`: Digital terrain model with 1 m resolution that was used for calculation of morphometric variables. It is too large to store online, hence on request: *krystof.chytry@gmail.com*
  * `Species_data.xlsx`: Ecological Indicator Values and trait values of individual species in the dataset used in the models.
  * `T5_species_data_v2-0.xlsx`: Original vegetation plot data in the long form. **logger_ID** is the primary key.
  * `MCdata_v2-0_GEO_vegetation.gpkg`: Header data for the vegetation plots including spatial information from the differential GPS. **logger_ID** is the primary key.


* `data`
   * Files in this folder were created in scripts `S1.R` and `S2.R`.
   * `Elevation.csv`: Elevation of vegetation plots.
   * `(\\d)_topoindex.csv`: Calculated morphometric variables, separate csv files per resolution.
   * `Morphometric_parameters.csv`: Calculated morphometric variables. Resolutions in rows, morphometric variables in columns.
   * `Morphometric_parameters_meta.csv`: Meta file.
   * `Morphometric_parameters_bank.csv`: Calculated morphometric variables the way they are later used in the models.
   * `species_data_wide.csv`: Vegetation plot data in a wide form (number of rows equal the number of analysed vegetation plots, i.e., 900).
   * `Comprop_means.csv`: Community means of traits and Ecological Indicator Values. Created in `S2.R`.
* `meta_coarse_resolution`
  * All files in this folder too large, hence on request: *krystof.chytry@gmail.com*
  * Files in this folder were created in the script `S1.R` and in SAGA GIS.
  * In the `S1.R` script, the DTM was aggregated to coarser resolutions and stored as individual `.tif` files. These raster files were later used in SAGA to derive regional morphometric variables (such as Flow Accumulation, Direct Solar Radiation etc.). They were stored as `.sgrd` raster files in this folder and later, they were, again in script `S1.R`, integrated in `data\\Morphometric parameters.csv`.


* `results`
  * Files in this folder were created in scripts `S3.R` and `S4.R`.
  * `meta_results`
    * Files in this folder are temporarily stored and later used.
  * `figures`
    * Figures used in the paper.
  * `appendix_figures`
    * Figures used in the appendix of the paper.
  * `(SPECIES|COMPROPS)_models_resolution.csv`: Outputs of models of each morphometric variable over spatial resolutions. 
  * `(SPECIES|COMPROPS)_random_models.csv`: Outputs of models using combinations of morphometric variables as predictors.
  * `(SPECIES|COMPROPS)_elevation_models.csv`: Outputs of models using with elevation as the only predictor.

### Brief description of the workflow

#### `S1.R`

* imports gpkg with location of sites `raw_input\\MCdata_v2-0_GEO_vegetation.gpkg` and dtm of the study area `DTM_2017.tif`
* `for` loop over resolutions.
  * `for` loop for individual sites.
    * morphometric variables are derived
  * and stored into a list
* the list is stored into a csv file (one per resolution)
* these files are integrated into one `data\\Morphometric_parameters_meta.csv`


* original raster is aggregated into coarser resolutions and stored into `meta_coarse_resolution` folder
* all `.tif` files from the folder are imported into SAGA GIS
* regional morphometric parameters are derived in the SAGA GIS and the SAGA project is stored in the same folder. This stores also the individual raster files for every morphometric variable as `.sgrd` rasters.
* back in R, the rasters with regional morphometric variables are imported, values for individual plots are extracted and integrated into `Morphometric_parameters_meta.csv` together with elevation (`data\\Elevation.csv`) creating a new `data\\Morphometric_parameters.csv` file.


* converting species data into wide form and storing as `data\\species_data_wide.csv`


* converting morphometric variables into a shape suitable for subsequent models: `data\\Morphometric_parameters_bank.csv`



#### `S2.R`

* imports original vegetation plot data
* checks, which plots have at least 70% of species with available trait information
* for these, computes the unweighted community means of traits, for all plots computes also Ecological Indicator Values
* stores these community attributes in `data\\Comprop_means.csv`



#### `S3.R`

* reads response variables: vegetation plot data wide (`data\\species_data_wide.csv`) and community attributes (`data\\Comprop_means.csv`) 
* as well as predictor variables, i.e., morphometric variables (`data\\Morphometric_parameters_bank.csv`)
* out of column names of table with morphometric variables, creates `indexbank`
* `for` loop per (species|community attributes)
  * calls `indexbank` in every iteration to extract predictors for the models
  * all models per one species are run with `future_map` loop
  * outside this `future_map` loop a model with elevation as the only predictors are fitted (= null model)
  * the chunk with `future_map` returns a tibble per every species
* these tibbles are stored in a list
* variation partitioning is performed (using the null models calculated above)
* outputs of all models is stored in a `results\\(SPECIES|COMPROPS)_resolution_models.csv` file and plotted `results\\appendix_figures\\Appendix_S3(A|B).pdf`
* the models are also averaged across all species and all community attributes and plotted for the paper itself `results\\figures\\Figure1(A|B).png` 



#### `S4.R`

* reads response variables: vegetation plot data wide (`data\\species_data_wide.csv`) and community attributes (`data\\Comprop_means.csv`) 
* as well as predictor variables, i.e., morphometric variables (`data\\Morphometric_parameters_bank.csv`)
* creates a source for allowed combinations of morphometric variables and resolutions (some morphometric variables are omitted in this case)
* creates a 1000 random predictor sets of 5 morphometric variables with fine, coarse and combined resolutions
* `for` loops over fine, coarse and combined predictor sets
  * `for` loops over 10 chunks (1/10 of 1000 models) which it stores into separate csv files
    * `for` loops over all species in the dataset
      * uses the future map to fit the first 1000 of models
    * stores the models into a csv file
* reads all csv files and combine them into one file: `results\\SPECIES_random_models.csv`
* later uses this file for plotting the output: `results\\figures\\Figure(2|3).png`



#### `S5.R`

* uses the outputs of scripts `S3.R` and `S4.R`, i.e., results, to derive information referred to in the main paper.
