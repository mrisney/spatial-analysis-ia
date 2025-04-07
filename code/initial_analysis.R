# Initial R Code for LA Wildfire Impact Analysis
# This script demonstrates loading data, creating visualizations,
# and conducting basic spatial analysis for the project

# Load required libraries
library(sf)         # For spatial data handling
library(dplyr)      # For data manipulation
library(ggplot2)    # For data visualization - or gspatial
library(tmap)       # For thematic maps
library(spdep)      # For spatial dependence analysis
library(spatstat)   # For point pattern analysis
library(geoR)       # For geostatistical analysis
library(raster)     # For raster data
library(RColorBrewer) # For color palettes
library(viridis)    # For colorblind-friendly palettes

# Set working directory (adjust as needed)
# setwd("/Users/marcrisney/Projects/jhu/mas/Spring2025/spatial-analysis-ia/")

# Data Loading ----------------------------------------------------------------

# 1. Load LA County boundary
la_county <- st_read("data/boundaries/la_county_boundary.shp")

# 2. Load census tracts with Social Vulnerability Index
tracts <- st_read("data/boundaries/la_census_tracts.shp")

# 3. Load fire perimeters
fire_perimeters <- st_read("data/fires/la_fires_jan2025.shp")

# 4. Load CalEnviroScreen asthma ED data
asthma_data <- read.csv("data/health/calenviroscreen_asthma.csv")

# 5. Load air quality data from PurpleAir
purpleair_data <- read.csv("data/air_quality/purpleair_jan2025.csv")

# Data Preparation ------------------------------------------------------------

# Ensure all data is in the same coordinate reference system (CRS)
la_county <- st_transform(la_county, 2229)  # California State Plane Zone 5 (NAD83)
tracts <- st_transform(tracts, 2229)
fire_perimeters <- st_transform(fire_perimeters, 2229)

# Join CalEnviroScreen asthma data to census tracts
tracts <- tracts %>%
  left_join(asthma_data, by = c("GEOID" = "census_tract"))

# Create spatial points object for PurpleAir sensors
purpleair_sf <- st_as_sf(purpleair_data, 
                        coords = c("longitude", "latitude"), 
                        crs = 4326) %>%  # WGS84
  st_transform(2229)  # Transform to match other layers

# Exploratory Data Analysis ---------------------------------------------------

# Basic map of fire perimeters
tm_shape(la_county) +
  tm_borders() +
  tm_shape(fire_perimeters) +
  tm_fill(col = "red", alpha = 0.5) +
  tm_borders(col = "red", lwd = 1) +
  tm_layout(title = "January 2025 LA County Fires")

# Map of asthma ED visits by census tract
tm_shape(tracts) +
  tm_fill("asthma_rate", 
          style = "quantile", 
          palette = "YlOrRd",
          title = "Asthma ED Visits\n(per 10,000)") +
  tm_borders(alpha = 0.2) +
  tm_shape(fire_perimeters) +
  tm_borders(col = "red", lwd = 2) +
  tm_layout(title = "Asthma ED Visits and Fire Perimeters")

# Social Vulnerability Index map
tm_shape(tracts) +
  tm_fill("svi_score", 
          style = "quantile", 
          palette = "PuBu",
          title = "Social Vulnerability\nIndex") +
  tm_borders(alpha = 0.2) +
  tm_shape(fire_perimeters) +
  tm_borders(col = "red", lwd = 2) +
  tm_layout(title = "Social Vulnerability and Fire Perimeters")

# Fire Exposure Analysis ------------------------------------------------------

# Create buffers around fire perimeters
buffer_5km <- st_buffer(fire_perimeters, dist = 5000)
buffer_10km <- st_buffer(fire_perimeters, dist = 10000)
buffer_20km <- st_buffer(fire_perimeters, dist = 20000)

# Identify tracts within each buffer zone
tracts$in_5km <- as.integer(st_intersects(tracts, buffer_5km, sparse = FALSE))
tracts$in_10km <- as.integer(st_intersects(tracts, buffer_10km, sparse = FALSE))
tracts$in_20km <- as.integer(st_intersects(tracts, buffer_20km, sparse = FALSE))

# Map showing buffer zones
tm_shape(tracts) +
  tm_borders() +
  tm_shape(buffer_20km) +
  tm_fill(col = "yellow", alpha = 0.2) +
  tm_shape(buffer_10km) +
  tm_fill(col = "orange", alpha = 0.3) +
  tm_shape(buffer_5km) +
  tm_fill(col = "red", alpha = 0.4) +
  tm_shape(fire_perimeters) +
  tm_fill(col = "darkred") +
  tm_layout(title = "Fire Impact Zones")

# Compute distance from each tract centroid to nearest fire perimeter
tracts_centroid <- st_centroid(tracts)
tracts$dist_to_fire <- st_distance(tracts_centroid, fire_perimeters)
tracts$dist_to_fire <- as.numeric(tracts$dist_to_fire) / 1000  # Convert to km

# Air Quality Analysis --------------------------------------------------------

# Interpolate air quality data (example using IDW)
# Create a grid over LA County
la_grid <- st_make_grid(la_county, cellsize = 1000, what = "centers") %>%
  st_as_sf() %>%
  st_intersection(la_county)

# Extract PM2.5 values from PurpleAir data
pm25_values <- purpleair_sf$pm25
coords <- st_coordinates(purpleair_sf)

# Use gstat package for IDW interpolation
library(gstat)
idw_model <- gstat(formula = z ~ 1, locations = coords, data = data.frame(z = pm25_values))
idw_pred <- predict(idw_model, newdata = la_grid)

# Convert to raster for visualization
pm25_raster <- rasterFromXYZ(data.frame(
  x = st_coordinates(la_grid)[,1],
  y = st_coordinates(la_grid)[,2],
  z = idw_pred$var1.pred
))

# Plot air quality raster
tm_shape(pm25_raster) +
  tm_raster(palette = "plasma", title = "PM2.5 (μg/m³)") +
  tm_shape(fire_perimeters) +
  tm_borders(col = "red", lwd = 2) +
  tm_layout(title = "Estimated PM2.5 Levels During Fires")

# Spatial Analysis of Health and Vulnerability --------------------------------

# Create spatial weights matrix for census tracts
tracts_nb <- poly2nb(tracts, queen = TRUE)
tracts_weights <- nb2listw(tracts_nb, style = "W")

# Compute global Moran's I for asthma rates
moran_asthma <- moran.test(tracts$asthma_rate, tracts_weights)
print(moran_asthma)

# Compute local Moran's I to identify clusters
local_moran <- localmoran(tracts$asthma_rate, tracts_weights)
tracts$local_moran_p <- local_moran[, 5]
tracts$local_moran <- local_moran[, 1]

# Identify significant clusters
tracts$cluster_type <- "Not Significant"
tracts$cluster_type[tracts$asthma_rate > mean(tracts$asthma_rate, na.rm = TRUE) & 
                   tracts$local_moran > 0 & 
                   tracts$local_moran_p < 0.05] <- "High-High"
tracts$cluster_type[tracts$asthma_rate < mean(tracts$asthma_rate, na.rm = TRUE) & 
                   tracts$local_moran > 0 & 
                   tracts$local_moran_p < 0.05] <- "Low-Low"
tracts$cluster_type[tracts$asthma_rate > mean(tracts$asthma_rate, na.rm = TRUE) & 
                   tracts$local_moran < 0 & 
                   tracts$local_moran_p < 0.05] <- "High-Low"
tracts$cluster_type[tracts$asthma_rate < mean(tracts$asthma_rate, na.rm = TRUE) & 
                   tracts$local_moran < 0 & 
                   tracts$local_moran_p < 0.05] <- "Low-High"

# Map LISA clusters
tm_shape(tracts) +
  tm_fill("cluster_type", 
          palette = c("High-High" = "red", 
                     "Low-Low" = "blue", 
                     "High-Low" = "pink", 
                     "Low-High" = "lightblue", 
                     "Not Significant" = "white"),
          title = "Asthma Rate\nClusters") +
  tm_borders(alpha = 0.2) +
  tm_shape(fire_perimeters) +
  tm_borders(col = "black", lwd = 2) +
  tm_layout(title = "Spatial Clusters of Asthma ED Visits")

# Regression Analysis ---------------------------------------------------------

# Basic linear model (non-spatial)
model1 <- lm(asthma_rate ~ dist_to_fire + svi_score + poverty_pct + median_income, data = tracts)
summary(model1)

# Create interaction term between fire distance and social vulnerability
tracts$dist_x_svi <- tracts$dist_to_fire * tracts$svi_score

# Model with interaction term
model2 <- lm(asthma_rate ~ dist_to_fire + svi_score + dist_x_svi + poverty_pct + median_income, 
            data = tracts)
summary(model2)

# Spatial lag model to account for spatial autocorrelation
library(spatialreg)
model_spatial <- lagsarlm(asthma_rate ~ dist_to_fire + svi_score + dist_x_svi + poverty_pct + median_income, 
                         data = tracts, 
                         listw = tracts_weights)
summary(model_spatial)

# Compare model results
anova(model1, model2)
AIC(model1, model2)

# Map residuals to check for spatial patterns
tracts$residuals <- residuals(model_spatial)
tm_shape(tracts) +
  tm_fill("residuals", 
          style = "quantile", 
          palette = "RdBu",
          title = "Model Residuals") +
  tm_borders(alpha = 0.2) +
  tm_layout(title = "Spatial Distribution of Model Residuals")

# Differential Impact Analysis ------------------------------------------------

# Compare asthma rates in fire-affected vs. non-affected areas across SVI categories
tracts$svi_category <- cut(tracts$svi_score, breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         labels = c("Low", "Medium-Low", "Medium-High", "High"))

# Calculate mean asthma rates for each combination
impact_summary <- tracts %>%
  st_drop_geometry() %>%
  group_by(in_10km, svi_category) %>%
  summarize(
    mean_asthma = mean(asthma_rate, na.rm = TRUE),
    sd_asthma = sd(asthma_rate, na.rm = TRUE),
    count = n()
  )

# Create comparative bar plot
ggplot(impact_summary, aes(x = svi_category, y = mean_asthma, fill = factor(in_10km))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_asthma - sd_asthma/sqrt(count), 
                   ymax = mean_asthma + sd_asthma/sqrt(count)),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "firebrick"),
                   labels = c("0" = "Not Fire-Affected", "1" = "Fire-Affected"),
                   name = "Fire Exposure") +
  labs(title = "Asthma ED Rates by Social Vulnerability and Fire Exposure",
       x = "Social Vulnerability Index Category",
       y = "Asthma ED Rate (per 10,000)") +
  theme_minimal()

# Save processed data for further analysis
st_write(tracts, "data/processed/la_tracts_analysis.shp")
saveRDS(model_spatial, "data/processed/spatial_model.rds")

# Next steps - Additional analyses to consider:
# 1. Incorporate temporal data to compare before/after fire impacts
# 2. Add building damage assessment data
# 3. Implement proper smoke dispersion modeling using wind data
# 4. Perform time series analysis on emergency department visits
# 5. Integrate Caltech soil testing data when available