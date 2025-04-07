# 01_data_prep.R
# Data preparation script for LA Wildfire Spatial Analysis Project
# 
# This script:
# 1. Loads and standardizes datasets from different sources
# 2. Performs basic cleaning and transformation
# 3. Converts to consistent coordinate reference systems (CRS)
# 4. Creates spatial joins between different geographic datasets
# 5. Saves processed data for further analysis

# Load required libraries
library(sf)          # For spatial data handling
library(dplyr)       # For data manipulation
library(readr)       # For reading CSV files
library(stringr)     # For string manipulation
library(lubridate)   # For date handling
library(tidyr)       # For data reshaping
library(tigris)      # For census tract boundaries
options(tigris_use_cache = TRUE)  # Cache tigris downloads

# Set up directories
root_dir <- getwd()  # Assumes script is run from the repository root
data_dir <- file.path(root_dir, "data")
out_dir <- file.path(data_dir, "processed")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Set CRS - using California State Plane Zone 5 (NAD83) - EPSG:2229
target_crs <- 2229

# Define file paths for input data
# Fire Perimeter Data
eaton_fire_path <- file.path(data_dir, "fires", "Eaton_fire_perimeter.shp")
palisades_fire_path <- file.path(data_dir, "fires", "Palisades_fire_perimeter.shp")
historical_fires_path <- file.path(data_dir, "fires", "historical_ca_fires.shp")

# Social Vulnerability Index
svi_path <- file.path(data_dir, "svi", "SVI2020_CALIFORNIA_tract.csv")

# CalEnviroScreen
ces_path <- file.path(data_dir, "calenviroscreen", "calenviroscreen_4.0_tract_data.csv")

# Air Quality Data
aqi_path <- file.path(data_dir, "aqi", "purpleair_jan2025.csv")

# Building Damage Assessment
bldg_damage_path <- file.path(data_dir, "building_damage", "LA_fires_bldg_damage.shp")

# Load Administrative Boundaries
cat("Loading administrative boundaries...\n")
# Get LA County boundary
la_county <- tigris::counties(state = "CA") %>%
  filter(NAME == "Los Angeles") %>%
  st_transform(target_crs)

# Get census tracts for LA County
la_tracts <- tigris::tracts(state = "CA", county = "Los Angeles") %>%
  st_transform(target_crs)

# Function to standardize GEOID in different datasets
standardize_geoid <- function(geoid) {
  # Ensure GEOID is 11 characters (state(2) + county(3) + tract(6))
  geoid <- as.character(geoid)
  # If GEOID is missing leading zeros, add them
  if (nchar(geoid) < 11) {
    geoid <- str_pad(geoid, 11, pad = "0")
  }
  return(geoid)
}

# Load and process Social Vulnerability Index data
cat("Processing Social Vulnerability Index data...\n")
svi_data <- read_csv(svi_path) %>%
  # Select relevant columns - adjust based on actual data structure
  select(FIPS, RPL_THEMES, RPL_THEME1, RPL_THEME2, RPL_THEME3, RPL_THEME4,
         EP_POV, EP_UNEMP, EP_PCI, EP_NOHSDP, EP_AGE65, EP_AGE17,
         EP_DISABL, EP_SNGPNT, EP_MINRTY, EP_LIMENG, EP_MUNIT, EP_MOBILE,
         EP_CROWD, EP_NOVEH, EP_GROUPQ) %>%
  # Rename columns for clarity
  rename(
    geoid = FIPS,
    svi_score = RPL_THEMES,
    svi_socioeconomic = RPL_THEME1,
    svi_household = RPL_THEME2,
    svi_minority = RPL_THEME3,
    svi_housing_transport = RPL_THEME4
  ) %>%
  # Standardize GEOID
  mutate(geoid = standardize_geoid(geoid))

# Load and process CalEnviroScreen data
cat("Processing CalEnviroScreen data...\n")
ces_data <- read_csv(ces_path) %>%
  # Select relevant columns - adjust based on actual data structure
  select(Census.Tract, CES.4.0.Score, CES.4.0.Percentile, Asthma, Asthma.Pctl, 
         Poverty, Unemployment, Population, 
         PM2.5, Ozone, Diesel.PM, contains("Pctl")) %>%
  # Rename columns for clarity
  rename(
    geoid = Census.Tract,
    ces_score = CES.4.0.Score,
    ces_percentile = CES.4.0.Percentile,
    asthma_rate = Asthma,
    asthma_percentile = Asthma.Pctl
  ) %>%
  # Standardize GEOID
  mutate(geoid = standardize_geoid(geoid))

# Load fire perimeter data
cat("Processing fire perimeter data...\n")
# Check if files exist since they're hypothetical in this script
if (file.exists(eaton_fire_path)) {
  eaton_fire <- st_read(eaton_fire_path) %>%
    st_transform(target_crs)
} else {
  warning("Eaton fire perimeter file not found. Skipping.")
  eaton_fire <- NULL
}

if (file.exists(palisades_fire_path)) {
  palisades_fire <- st_read(palisades_fire_path) %>%
    st_transform(target_crs)
} else {
  warning("Palisades fire perimeter file not found. Skipping.")
  palisades_fire <- NULL
}

# Combine 2025 fire perimeters
if (!is.null(eaton_fire) && !is.null(palisades_fire)) {
  fires_2025 <- bind_rows(
    eaton_fire %>% mutate(FIRE_NAME = "Eaton"),
    palisades_fire %>% mutate(FIRE_NAME = "Palisades")
  )
  
  # Save combined fire perimeters
  st_write(fires_2025, file.path(out_dir, "fires_2025.shp"), append = FALSE)
} else {
  warning("Could not combine 2025 fire perimeters due to missing data.")
}

# Load historical fires if available
if (file.exists(historical_fires_path)) {
  historical_fires <- st_read(historical_fires_path) %>%
    st_transform(target_crs)
  
  # Filter for relevant historical fires (example - adjust criteria based on your research)
  relevant_hist_fires <- historical_fires %>%
    filter(
      YEAR_ >= 2010,
      STATE == "CA",
      GIS_ACRES > 1000  # Only large fires
    )
  
  # Save filtered historical fires
  st_write(relevant_hist_fires, file.path(out_dir, "relevant_historical_fires.shp"), append = FALSE)
} else {
  warning("Historical fires file not found. Skipping.")
}

# Process air quality data
cat("Processing air quality data...\n")
if (file.exists(aqi_path)) {
  aqi_data <- read_csv(aqi_path) %>%
    # Example processing - adjust based on actual data structure
    filter(!is.na(latitude), !is.na(longitude)) %>%
    mutate(
      date_time = ymd_hms(datetime),
      date = as.Date(date_time)
    )
  
  # Convert to SF object
  aqi_sf <- aqi_data %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(target_crs)
  
  # Save processed AQI data
  write_csv(aqi_data, file.path(out_dir, "processed_aqi.csv"))
  st_write(aqi_sf, file.path(out_dir, "aqi_points.shp"), append = FALSE)
} else {
  warning("AQI data file not found. Skipping.")
}

# Process building damage data
cat("Processing building damage data...\n")
if (file.exists(bldg_damage_path)) {
  bldg_damage <- st_read(bldg_damage_path) %>%
    st_transform(target_crs)
  
  # Example: Calculate damage statistics by census tract
  tracts_with_damage <- st_join(la_tracts, bldg_damage) %>%
    group_by(GEOID) %>%
    summarize(
      buildings_affected = n(),
      buildings_destroyed = sum(DAMAGE_CLS == "Destroyed", na.rm = TRUE),
      buildings_major_damage = sum(DAMAGE_CLS == "Major", na.rm = TRUE),
      buildings_minor_damage = sum(DAMAGE_CLS == "Minor", na.rm = TRUE)
    )
  
  # Save processed building damage data
  st_write(tracts_with_damage, file.path(out_dir, "tracts_with_damage.shp"), append = FALSE)
} else {
  warning("Building damage file not found. Skipping.")
}

# Join all datasets to census tracts
cat("Creating integrated spatial dataset...\n")

# Start with LA County tracts
integrated_data <- la_tracts %>%
  # Add a standardized GEOID
  mutate(geoid = standardize_geoid(GEOID))

# Join SVI data
integrated_data <- integrated_data %>%
  left_join(svi_data, by = "geoid")

# Join CalEnviroScreen data
integrated_data <- integrated_data %>%
  left_join(ces_data, by = "geoid")

# Identify tracts within fire buffer zones
if (exists("fires_2025") && !is.null(fires_2025)) {
  # Create buffers at multiple distances
  buffer_5km <- st_buffer(fires_2025, dist = 5000)
  buffer_10km <- st_buffer(fires_2025, dist = 10000)
  buffer_20km <- st_buffer(fires_2025, dist = 20000)
  
  # Identify tracts within each buffer zone
  integrated_data <- integrated_data %>%
    mutate(
      in_fire_5km = as.integer(st_intersects(geometry, buffer_5km, sparse = FALSE)[,1]),
      in_fire_10km = as.integer(st_intersects(geometry, buffer_10km, sparse = FALSE)[,1]),
      in_fire_20km = as.integer(st_intersects(geometry, buffer_20km, sparse = FALSE)[,1])
    )
  
  # Calculate distance to nearest fire perimeter
  tract_centroids <- st_centroid(integrated_data)
  distances <- st_distance(tract_centroids, fires_2025)
  integrated_data$dist_to_fire_m <- apply(distances, 1, min)
  
  # Convert to kilometers for easier interpretation
  integrated_data$dist_to_fire_km <- integrated_data$dist_to_fire_m / 1000
}

# Save final integrated dataset
cat("Saving final integrated dataset...\n")
st_write(integrated_data, file.path(out_dir, "la_integrated_data.shp"), append = FALSE)
# Also save as GeoJSON for web mapping and other uses
st_write(integrated_data, file.path(out_dir, "la_integrated_data.geojson"), append = FALSE)

# Create simplified version with key variables for easier analysis
integrated_data_key_vars <- integrated_data %>%
  select(geoid, NAME, svi_score, asthma_rate, ces_score, 
         dist_to_fire_km, in_fire_5km, in_fire_10km, in_fire_20km)

# Save simplified version
write_csv(st_drop_geometry(integrated_data_key_vars), 
          file.path(out_dir, "la_integrated_data_key_vars.csv"))

cat("Data preparation complete!\n")
cat("Processed data saved to:", out_dir, "\n")

# Next steps - additional processing that could be added:
# 1. Interpolate air quality data to create continuous PM2.5 surface
# 2. Calculate temporal aggregates of AQI data (daily averages, maximums)
# 3. Conduct wind direction analysis to model smoke dispersion
# 4. Extract additional variables from CalEnviroScreen/SVI for specific analysis needs
# 5. Prepare historical fire comparison datasets