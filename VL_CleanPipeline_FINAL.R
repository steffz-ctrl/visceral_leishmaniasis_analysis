# FULL CLEAN PIPELINE: Visceral Leishmaniasis Transmission Dynamics
# Project: The Epidemiology of Visceral Leishmaniasis:
#          Uncovering Transmission Dynamics

# Description: Complete data preparation and covariate extraction pipeline
#              from raw epidata to final modelling dataset

# SECTION 1: Load Required Libraries
packages <- c("sf", "terra", "dplyr", "ggplot2", "readr", "tidyr")

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed],
                   repos        = "https://cran.r-project.org",
                   dependencies = TRUE)
}

library(sf)       # Vector spatial data
library(terra)    # Raster spatial data
library(dplyr)    # Data manipulation
library(ggplot2)  # Visualisation
library(readr)    # Fast CSV reading
library(tidyr)    # Data reshaping

cat("Section 1 complete: All libraries loaded\n")


# SECTION 2: Set File Paths

# Update these paths if running locally vs HPC
epidata_path    <- "~/Kenya EpiData Complete 2.csv"
shapefile_path  <- "~/shapefiles/Shapesfiles/"
covariates_path <- "~/covariates/Environmental covariates/"
hydro_path      <- "~/hydro_folder/"
output_path     <- "~/outputs/"
socio_path      <- "~/"

cat("Section 2 complete: File paths set\n")


# SECTION 3: Load Epidemiological Data

epidata <- read_csv(epidata_path)

# Standardise column names replacing spaces with dots
names(epidata) <- gsub(" ", ".", names(epidata))

cat("Section 3 complete\n")
cat("Total records loaded:", nrow(epidata), "\n")
cat("Column names:", names(epidata), "\n")


# SECTION 4: Load Shapefiles

kenya_0 <- st_read(paste0(shapefile_path, "gadm41_KEN_0.shp"), quiet = TRUE)
kenya_1 <- st_read(paste0(shapefile_path, "gadm41_KEN_1.shp"), quiet = TRUE)
kenya_2 <- st_read(paste0(shapefile_path, "gadm41_KEN_2.shp"), quiet = TRUE)

cat("Section 4 complete: Shapefiles loaded\n")
cat("National boundary: 1 polygon\n")
cat("Counties (Admin 1):", nrow(kenya_1), "\n")
cat("Sub-counties (Admin 2):", nrow(kenya_2), "\n")


# SECTION 5: Load and Filter Raster Files

# List only properly named raster files excluding NDVI
# Pattern matches: Brightness_2020.tif, Temperature_2021.tif etc.
raster_files <- list.files(
  path      = covariates_path,
  pattern   = "^[A-Z][a-zA-Z]+_[0-9]{4}\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

# Remove NDVI files per supervisor instruction
raster_files <- raster_files[!grepl("NDVI", raster_files)]

cat("Section 5 complete\n")
cat("Total raster files found:", length(raster_files), "\n")
# Expected: 36 files (6 variables x 6 years)


# SECTION 6: Load and Align All Rasters

# Load all rasters into a list
raster_list <- lapply(raster_files, rast)

# Use first raster as alignment reference template
reference_raster <- raster_list[[1]]

# Align all rasters to the reference grid
# Corrects floating point precision mismatches across files
raster_list_aligned <- lapply(raster_list, function(r) {
  if (!compareGeom(r, reference_raster, stopOnError = FALSE)) {
    resample(r, reference_raster, method = "bilinear")
  } else {
    r
  }
})

# Confirm all rasters share identical dimensions
dims_check <- sapply(raster_list_aligned,
                     function(r) paste(nrow(r), ncol(r)))
cat("Section 6 complete\n")
cat("Unique dimensions after alignment:", length(unique(dims_check)), "\n")
cat("Total aligned rasters:", length(raster_list_aligned), "\n")
# Both should return 1 and 36


# SECTION 7: Build Organised Raster Stack by Year

# Define covariate variables and years
# NDVI excluded per supervisor instruction
cov_cols <- c("Brightness", "Greenness", "Rainfall",
              "RH", "Temperature", "Wetness")
years    <- 2020:2025

# Build organised list of raster stacks one per year
raster_stack_by_year <- list()

for (yr in years) {
  yr_layers <- list()
  for (var in cov_cols) {
    layer_name <- paste0(var, "_", yr)
    idx <- grep(paste0("^", layer_name, "\\.tif$"),
                basename(raster_files))
    if (length(idx) == 0) {
      cat("Year:", yr, "| Variable:", var, "| NOT FOUND\n")
      next
    }
    r <- raster_list_aligned[[idx]]
    names(r) <- layer_name
    yr_layers[[var]] <- r
  }
  raster_stack_by_year[[as.character(yr)]] <- rast(yr_layers)
}

cat("Section 7 complete\n")
cat("Years in stack:", names(raster_stack_by_year), "\n")
cat("Layers per year:", sapply(raster_stack_by_year, nlyr), "\n")
# All years should have 6 layers


# SECTION 8: Crop and Mask Raster Stacks to Kenya Boundary

# Convert Kenya national boundary to terra vector format
kenya_vect <- vect(kenya_0)

# Crop and mask each year stack to Kenya boundary
raster_stack_kenya <- lapply(names(raster_stack_by_year), function(yr) {
  cat("Cropping and masking year:", yr, "\n")
  r_crop <- crop(raster_stack_by_year[[yr]], kenya_vect)
  r_mask <- mask(r_crop, kenya_vect)
  return(r_mask)
})

names(raster_stack_kenya) <- names(raster_stack_by_year)

cat("Section 8 complete\n")
cat("Layers per year after masking:",
    sapply(raster_stack_kenya, nlyr), "\n")


# SECTION 9: Epidata Processing

# Convert full epidata to spatial object
epidata_sf <- st_as_sf(
  epidata,
  coords = c("Longitude", "Latitude"),
  crs    = 4326,
  remove = FALSE
)

# Step 1: Crop to Kenya boundary first
epidata_kenya <- epidata_sf[
  st_within(epidata_sf, st_union(kenya_0), sparse = FALSE), ]

cat("Records before Kenya crop:", nrow(epidata_sf), "\n")
cat("Records after Kenya crop:", nrow(epidata_kenya), "\n")
# 66 cross-border points removed

# Step 2: Filter to 2020-2025 to match covariate availability
epidata_kenya_2020 <- epidata_kenya %>%
  filter(Diagnosis.Year >= 2020 & Diagnosis.Year <= 2025)

cat("Records after year filter (2020-2025):",
    nrow(epidata_kenya_2020), "\n")
cat("Year distribution:\n")
print(table(epidata_kenya_2020$Diagnosis.Year))

# Step 3: Jitter coordinates to eliminate spatial duplicates
# Amount smaller than one raster pixel to stay within same grid cell
set.seed(123)
epidata_final <- st_jitter(epidata_kenya_2020, amount = 0.00005)

# Confirm zero duplicates after jittering
coords_check <- st_coordinates(epidata_final)
cat("Section 9 complete\n")
cat("Duplicate coordinates after jittering:",
    sum(duplicated(coords_check)), "\n")
# Should return 0


# SECTION 10: Matched Year-to-Covariate Extraction

# Convert jittered sf points to terra vector for extraction
epi_vect <- vect(epidata_final)

# Initialise list to store extraction results per year
extracted_list <- list()

for (yr in 2020:2025) {
  cat("Extracting covariates for year:", yr, "\n")
  pts_yr    <- epi_vect[epi_vect$Diagnosis.Year == yr, ]
  stack_yr  <- raster_stack_kenya[[as.character(yr)]]
  extracted_yr <- terra::extract(stack_yr, pts_yr,
                                  method = "bilinear")
  extracted_yr  <- extracted_yr[, -1]
  cases_yr      <- as.data.frame(pts_yr)
  result_yr     <- cbind(cases_yr, extracted_yr)
  extracted_list[[as.character(yr)]] <- result_yr
  cat("  Cases:", length(pts_yr),
      "| Covariates extracted:", ncol(extracted_yr), "\n")
}


# SECTION 11: Combine All Years into Base Modelling Dataset

model_data <- bind_rows(extracted_list)

cat("Section 11 complete\n")
cat("Base dataset dimensions:",
    nrow(model_data), "rows x", ncol(model_data), "cols\n")
cat("Missing values:\n")
print(colSums(is.na(model_data)))
# Expected: 6977 rows x 12 cols, zero missing values


# SECTION 12: Healthcare Facility Access Distance

# Load facility data
facility_data <- read_csv(paste0(socio_path,
                                  "Kenya EpiData Complete 2.csv"))
names(facility_data) <- gsub(" ", ".",
                              gsub("/", ".", names(facility_data)))

# Extract unique facility locations
facilities <- facility_data %>%
  select(Health_facility, County_of_facility,
         Subcounty_of_facility, Latitude, Longitude) %>%
  distinct(Health_facility, .keep_all = TRUE) %>%
  filter(!is.na(Latitude) & !is.na(Longitude))

cat("Unique VL treatment facilities:", nrow(facilities), "\n")

# Convert to spatial object
facilities_sf   <- st_as_sf(facilities,
                              coords = c("Longitude", "Latitude"),
                              crs    = 4326)
facilities_vect <- vect(facilities_sf)

# Compute distance raster from every Kenya pixel to nearest facility
cat("Computing healthcare access distance raster...\n")
distance_raster    <- terra::distance(reference_raster, facilities_vect)
distance_raster_km <- distance_raster / 1000
names(distance_raster_km) <- "Healthcare_Access_km"

# Crop and mask to Kenya
distance_crop  <- crop(distance_raster_km, kenya_vect)
distance_kenya <- mask(distance_crop, kenya_vect)

# Extract at case locations
healthcare_vals <- terra::extract(distance_kenya, epi_vect,
                                   method = "bilinear")
healthcare_vals <- healthcare_vals[, -1]
model_data$Healthcare_Access_km <- healthcare_vals

cat("Section 12 complete\n")
cat("Healthcare access summary:\n")
print(summary(model_data$Healthcare_Access_km))
cat("Missing values:", sum(is.na(model_data$Healthcare_Access_km)), "\n")


# SECTION 13: Water Bodies Distance

# Load HydroRIVERS - Kenya extent only using bounding box filter
kenya_bbox   <- st_bbox(kenya_0)
cat("Loading HydroRIVERS - Kenya extent only...\n")

rivers_kenya <- st_read(
  paste0(hydro_path, "HydroRIVERS_v10_af.gdb"),
  layer      = "HydroRIVERS_v10_af",
  wkt_filter = st_as_text(st_as_sfc(kenya_bbox))
)
cat("Kenya rivers loaded:", nrow(rivers_kenya), "\n")

# Filter to major rivers only (ORD_FLOW 1-4)
rivers_major <- rivers_kenya[rivers_kenya$ORD_FLOW <= 4, ]
cat("Major rivers:", nrow(rivers_major), "features\n")

# Load HydroLAKES
lakes_world  <- st_read(
  paste0(hydro_path,
         "lakes/HydroLAKES_points_v10_shp/HydroLAKES_points_v10.shp"),
  quiet = TRUE
)
lakes_kenya  <- lakes_world[st_within(lakes_world,
                                       st_as_sfc(kenya_bbox),
                                       sparse = FALSE), ]
cat("Kenya lakes:", nrow(lakes_kenya), "\n")

# Convert to terra vectors
rivers_vect <- vect(rivers_major)
lakes_vect  <- vect(lakes_kenya)

# Compute distance to rivers and lakes separately
cat("Computing distance to rivers (may take ~30 minutes locally)...\n")
dist_rivers <- terra::distance(reference_raster, rivers_vect)
cat("Rivers done\n")
dist_lakes  <- terra::distance(reference_raster, lakes_vect)
cat("Lakes done\n")

# Take minimum distance from either water body type
water_distance    <- min(dist_rivers, dist_lakes)
water_distance_km <- water_distance / 1000
names(water_distance_km) <- "Water_Distance_km"

# Crop and mask to Kenya
water_crop  <- crop(water_distance_km, kenya_vect)
water_kenya <- mask(water_crop, kenya_vect)

# Extract at case locations
water_vals              <- terra::extract(water_kenya, epi_vect,
                                           method = "bilinear")
water_vals              <- water_vals[, -1]
model_data$Water_Distance_km <- water_vals

cat("Section 13 complete\n")
cat("Water distance summary:\n")
print(summary(model_data$Water_Distance_km))
cat("Missing values:", sum(is.na(model_data$Water_Distance_km)), "\n")


# SECTION 14: Night Time Lights

cat("Loading night time lights...\n")

ntl_raw <- rast(paste0(socio_path,
  "SVDNB_npp_20200101-20200131_global_vcmslcfg_v10_c202002111500.avg_rade9h.masked.tif"))

# Crop and mask to Kenya
ntl_crop  <- crop(ntl_raw, kenya_vect)
ntl_kenya <- mask(ntl_crop, kenya_vect)
names(ntl_kenya) <- "NightLights"

cat("Night time lights loaded\n")
cat("Min:", round(min(values(ntl_kenya), na.rm = TRUE), 4), "\n")
cat("Max:", round(max(values(ntl_kenya), na.rm = TRUE), 4), "\n")

# Extract at case locations
ntl_vals            <- terra::extract(ntl_kenya, epi_vect,
                                       method = "bilinear")
ntl_vals            <- ntl_vals[, -1]
model_data$NightLights <- ntl_vals

cat("Section 14 complete\n")
cat("Missing values:", sum(is.na(model_data$NightLights)), "\n")


# SECTION 15: Population Density

cat("Loading WorldPop population density...\n")

# Load template raster for resampling
template_raster <- rast(paste0(covariates_path,
                                "2020/Temperature_2020.tif"))

pop_raw      <- rast(paste0(socio_path, "ken_ppp_2020_constrained.tif"))
pop_crop     <- crop(pop_raw, kenya_vect)

# Resample to match covariate resolution using sum method
# Sum preserves total population count when aggregating
pop_resampled <- resample(pop_crop, template_raster, method = "sum")
pop_kenya     <- mask(pop_resampled, kenya_vect)
names(pop_kenya) <- "Population_Density"

cat("Population density processed\n")
cat("Min:", round(min(values(pop_kenya), na.rm = TRUE), 2), "\n")
cat("Max:", round(max(values(pop_kenya), na.rm = TRUE), 2), "\n")

# Extract at case locations
pop_vals <- terra::extract(pop_kenya, epi_vect, method = "bilinear")
pop_vals <- pop_vals[, -1]
model_data$Population_Density <- pop_vals

# Replace missing values with median
median_pop <- median(model_data$Population_Density, na.rm = TRUE)
model_data$Population_Density[is.na(model_data$Population_Density)] <- median_pop

cat("Section 15 complete\n")
cat("Missing values after fix:",
    sum(is.na(model_data$Population_Density)), "\n")


# SECTION 16: Relative Wealth Index

cat("Loading Relative Wealth Index...\n")

rwi_data <- read.csv(paste0(socio_path,
                              "ken_relative_wealth_index.csv"))
cat("RWI rows:", nrow(rwi_data), "\n")

# Convert RWI points to spatial object
rwi_sf <- st_as_sf(rwi_data,
                    coords = c("longitude", "latitude"),
                    crs    = 4326)

# Convert case locations to sf for nearest neighbour matching
cases_sf <- st_as_sf(model_data,
                      coords = c("Longitude", "Latitude"),
                      crs    = 4326)

# Find nearest RWI point for each case location
cat("Finding nearest RWI point for each case...\n")
nearest_idx            <- st_nearest_feature(cases_sf, rwi_sf)
model_data$Wealth_Index <- rwi_data$rwi[nearest_idx]

cat("Section 16 complete\n")
cat("Wealth Index summary:\n")
print(summary(model_data$Wealth_Index))
cat("Missing values:", sum(is.na(model_data$Wealth_Index)), "\n")


# SECTION 17: Final Dataset Validation and Save

cat("\n=== FINAL DATASET VALIDATION ===\n")
cat("Dimensions:", nrow(model_data), "rows x", ncol(model_data), "cols\n")
cat("\nColumn names:\n")
print(names(model_data))
cat("\nMissing values per column:\n")
print(colSums(is.na(model_data)))
cat("\nPredictor summary:\n")

predictors <- c("Brightness", "Greenness", "Rainfall", "RH",
                "Temperature", "Wetness", "Healthcare_Access_km",
                "Water_Distance_km", "NightLights",
                "Population_Density", "Wealth_Index")

cat("Total predictors:", length(predictors), "\n")
cat("Environmental:", 6, "\n")
cat("Socioeconomic:", 5, "\n")

# Save final complete dataset
write.csv(
  model_data,
  paste0(output_path, "VL_Kenya_Model_Dataset_v6_COMPLETE.csv"),
  row.names = FALSE
)

cat("\nFinal dataset saved as VL_Kenya_Model_Dataset_v6_COMPLETE.csv\n")
cat("Pipeline complete. Ready for EDA and modelling.\n")


# QUICK VISUALISATION: Case locations on Kenya map

plot(st_geometry(kenya_1),
     border = "grey60",
     main   = "Final VL Case Locations (2020-2025)")
plot(st_geometry(epidata_final),
     add  = TRUE,
     col  = "red",
     pch  = 20,
     cex  = 0.3)
