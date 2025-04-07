# 02_eda.R
# Exploratory Data Analysis for LA Wildfire Spatial Analysis Project
#
# This script:
# 1. Loads processed data from 01_data_prep.R
# 2. Produces descriptive statistics for key variables
# 3. Creates exploratory maps and visualizations
# 4. Examines spatial patterns and relationships
# 5. Tests for spatial autocorrelation in key variables

# Load required libraries
library(sf)          # For spatial data handling
library(dplyr)       # For data manipulation
library(ggplot2)     # For data visualization
library(tmap)        # For thematic maps
library(spdep)       # For spatial dependence analysis
library(corrplot)    # For correlation matrices
library(RColorBrewer) # For color palettes
library(gridExtra)   # For arranging multiple plots
library(knitr)       # For formatted tables
library(scales)      # For formatting numbers and scales

# Set up directories
root_dir <- getwd()  # Assumes script is run from the repository root
data_dir <- file.path(root_dir, "data")
processed_dir <- file.path(data_dir, "processed")
results_dir <- file.path(root_dir, "results")
fig_dir <- file.path(results_dir, "figures")
map_dir <- file.path(results_dir, "maps")

# Create directories if they don't exist
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)

# Set the tmap mode to plotting
tmap_mode("plot")

# Load processed data
cat("Loading processed data...\n")
integrated_data <- st_read(file.path(processed_dir, "la_integrated_data.shp"))
fires_2025 <- st_read(file.path(processed_dir, "fires_2025.shp"))

# Basic dataset summary
cat("Dataset summary:\n")
print(paste("Number of census tracts:", nrow(integrated_data)))
print(paste("Number of variables:", ncol(integrated_data) - 1))  # Subtract geometry column

# Summary of key variables
cat("\nDescriptive statistics for key variables:\n")
key_vars <- integrated_data %>%
  st_drop_geometry() %>%
  select(svi_score, asthma_rate, ces_score, dist_to_fire_km) %>%
  summary()
print(key_vars)

# Fire exposure summary
exposure_summary <- integrated_data %>%
  st_drop_geometry() %>%
  summarize(
    tracts_within_5km = sum(in_fire_5km, na.rm = TRUE),
    pct_within_5km = mean(in_fire_5km, na.rm = TRUE) * 100,
    tracts_within_10km = sum(in_fire_10km, na.rm = TRUE),
    pct_within_10km = mean(in_fire_10km, na.rm = TRUE) * 100,
    tracts_within_20km = sum(in_fire_20km, na.rm = TRUE),
    pct_within_20km = mean(in_fire_20km, na.rm = TRUE) * 100
  )
print("Fire exposure summary:")
print(exposure_summary)

# Create a histogram of fire distances
p1 <- ggplot(integrated_data %>% st_drop_geometry(), aes(x = dist_to_fire_km)) +
  geom_histogram(binwidth = 2, fill = "firebrick", alpha = 0.7) +
  labs(
    title = "Distribution of Distances to Nearest Fire",
    x = "Distance (km)",
    y = "Number of Census Tracts"
  ) +
  theme_minimal()
ggsave(file.path(fig_dir, "fire_distance_histogram.png"), p1, width = 8, height = 6)

# Create a scatterplot of SVI vs. asthma rate
p2 <- ggplot(integrated_data %>% st_drop_geometry(), 
             aes(x = svi_score, y = asthma_rate)) +
  geom_point(aes(color = factor(in_fire_10km)), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "darkblue") +
  scale_color_manual(values = c("0" = "steelblue", "1" = "firebrick"),
                    labels = c("0" = "Not Fire-Affected", "1" = "Fire-Affected"),
                    name = "Within 10km of Fire") +
  labs(
    title = "Relationship Between Social Vulnerability and Asthma ED Rates",
    x = "Social Vulnerability Index Score",
    y = "Asthma ED Rate (per 10,000)"
  ) +
  theme_minimal()
ggsave(file.path(fig_dir, "svi_asthma_scatter.png"), p2, width = 8, height = 6)

# Create correlation matrix of key variables
corr_vars <- integrated_data %>%
  st_drop_geometry() %>%
  select(svi_score, svi_socioeconomic, svi_household, svi_minority, 
         svi_housing_transport, asthma_rate, ces_score, dist_to_fire_km)

corr_matrix <- cor(corr_vars, use = "pairwise.complete.obs")
png(file.path(fig_dir, "correlation_matrix.png"), width = 8, height = 8, units = "in", res = 300)
corrplot(corr_matrix, method = "circle", type = "upper", 
        tl.col = "black", tl.srt = 45, addCoef.col = "black", 
        number.cex = 0.7, mar = c(0,0,2,0),
        title = "Correlation Matrix of Key Variables")
dev.off()

# Create comparative boxplots for asthma rates by fire exposure
p3 <- ggplot(integrated_data %>% st_drop_geometry(), 
             aes(x = factor(in_fire_10km), y = asthma_rate)) +
  geom_boxplot(aes(fill = factor(in_fire_10km)), alpha = 0.7) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "firebrick"),
                   labels = c("0" = "Not Fire-Affected", "1" = "Fire-Affected"),
                   name = "Within 10km of Fire") +
  labs(
    title = "Asthma ED Rates by Fire Exposure",
    x = "Within 10km of Fire",
    y = "Asthma ED Rate (per 10,000)"
  ) +
  theme_minimal()
ggsave(file.path(fig_dir, "asthma_by_fire_exposure.png"), p3, width = 8, height = 6)

# Create maps of key variables ------------------------------------------------

# Map 1: Social Vulnerability Index
map1 <- tm_shape(integrated_data) +
  tm_fill("svi_score", 
          style = "quantile", 
          palette = "PuBu", 
          title = "SVI Score",
          n = 5) +
  tm_borders(alpha = 0.2) +
  tm_shape(fires_2025) +
  tm_borders(col = "red", lwd = 2) +
  tm_layout(
    title = "Social Vulnerability Index",
    legend.position = c("right", "bottom"),
    legend.text.size = 0.7,
    legend.title.size = 0.9
  )
tmap_save(map1, file.path(map_dir, "svi_map.png"), width = 8, height = 8)

# Map 2: Asthma ED Rates
map2 <- tm_shape(integrated_data) +
  tm_fill("asthma_rate", 
          style = "quantile", 
          palette = "YlOrRd", 
          title = "Asthma ED Rate\n(per 10,000)",
          n = 5) +
  tm_borders(alpha = 0.2) +
  tm_shape(fires_2025) +
  tm_borders(col = "red", lwd = 2) +
  tm_layout(
    title = "Asthma ED Rates",
    legend.position = c("right", "bottom"),
    legend.text.size = 0.7,
    legend.title.size = 0.9
  )
tmap_save(map2, file.path(map_dir, "asthma_map.png"), width = 8, height = 8)

# Map 3: Fire Proximity
map3 <- tm_shape(integrated_data) +
  tm_fill("dist_to_fire_km", 
          style = "fisher", 
          palette = "RdYlGn", 
          title = "Distance to Fire (km)",
          n = 6) +
  tm_borders(alpha = 0.2) +
  tm_shape(fires_2025) +
  tm_borders(col = "black", lwd = 2) +
  tm_layout(
    title = "Distance to Nearest Fire",
    legend.position = c("right", "bottom"),
    legend.text.size = 0.7,
    legend.title.size = 0.9
  )
tmap_save(map3, file.path(map_dir, "fire_distance_map.png"), width = 8, height = 8)

# Map 4: Combined Fire Buffers
buffer_colors <- c("#fee8c8", "#fdbb84", "#e34a33")
map4 <- tm_shape(integrated_data) +
  tm_borders(alpha = 0.2) +
  tm_shape(st_buffer(fires_2025, dist = 20000)) +
  tm_fill(col = buffer_colors[1], alpha = 0.3) +
  tm_shape(st_buffer(fires_2025, dist = 10000)) +
  tm_fill(col = buffer_colors[2], alpha = 0.4) +
  tm_shape(st_buffer(fires_2025, dist = 5000)) +
  tm_fill(col = buffer_colors[3], alpha = 0.5) +
  tm_shape(fires_2025) +
  tm_borders(col = "red", lwd = 2) +
  tm_fill(col = "red", alpha = 0.7) +
  tm_layout(
    title = "Fire Impact Zones",
    legend.position = c("right", "bottom"),
    legend.text.size = 0.7,
    legend.title.size = 0.9
  ) +
  tm_add_legend(type = "fill", 
                col = c("red", buffer_colors[3], buffer_colors[2], buffer_colors[1]),
                labels = c("Fire Perimeter", "5km Buffer", "10km Buffer", "20km Buffer"),
                alpha = c(0.7, 0.5, 0.4, 0.3))
tmap_save(map4, file.path(map_dir, "fire_buffer_map.png"), width = 8, height = 8)

# Spatial autocorrelation analysis --------------------------------------------

# Create spatial weights matrix based on census tract adjacency
cat("Creating spatial weights matrix...\n")
tracts_nb <- poly2nb(integrated_data, queen = TRUE)
tracts_weights <- nb2listw(tracts_nb, style = "W")

# Calculate global Moran's I for asthma rates
cat("Testing for spatial autocorrelation in asthma rates...\n")
asthma_moran <- moran.test(integrated_data$asthma_rate, tracts_weights)
print(asthma_moran)

# Calculate global Moran's I for SVI
cat("Testing for spatial autocorrelation in social vulnerability...\n")
svi_moran <- moran.test(integrated_data$svi_score, tracts_weights)
print(svi_moran)

# Calculate local Moran's I for asthma rates
cat("Calculating local indicators of spatial association...\n")
asthma_lisa <- localmoran(integrated_data$asthma_rate, tracts_weights)

# Add LISA results to integrated_data
integrated_data$asthma_local_i <- asthma_lisa[, 1]
integrated_data$asthma_lisa_p <- asthma_lisa[, 5]

# Identify significant LISA clusters
integrated_data$asthma_cluster <- "Not Significant"
mean_asthma <- mean(integrated_data$asthma_rate, na.rm = TRUE)

integrated_data$asthma_cluster[integrated_data$asthma_rate > mean_asthma & 
                              integrated_data$asthma_local_i > 0 & 
                              integrated_data$asthma_lisa_p < 0.05] <- "High-High"

integrated_data$asthma_cluster[integrated_data$asthma_rate < mean_asthma & 
                              integrated_data$asthma_local_i > 0 & 
                              integrated_data$asthma_lisa_p < 0.05] <- "Low-Low"

integrated_data$asthma_cluster[integrated_data$asthma_rate > mean_asthma & 
                              integrated_data$asthma_local_i < 0 & 
                              integrated_data$asthma_lisa_p < 0.05] <- "High-Low"

integrated_data$asthma_cluster[integrated_data$asthma_rate < mean_asthma & 
                              integrated_data$asthma_local_i < 0 & 
                              integrated_data$asthma_lisa_p < 0.05] <- "Low-High"

# Map LISA clusters
map5 <- tm_shape(integrated_data) +
  tm_fill("asthma_cluster", 
          palette = c("High-High" = "red", 
                     "Low-Low" = "blue", 
                     "High-Low" = "pink", 
                     "Low-High" = "lightblue", 
                     "Not Significant" = "white"),
          title = "Asthma Rate\nClusters") +
  tm_borders(alpha = 0.2) +
  tm_shape(fires_2025) +
  tm_borders(col = "black", lwd = 2) +
  tm_layout(
    title = "Spatial Clusters of Asthma ED Rates",
    legend.position = c("right", "bottom"),
    legend.text.size = 0.7,
    legend.title.size = 0.9
  )
tmap_save(map5, file.path(map_dir, "asthma_lisa_map.png"), width = 8, height = 8)

# Examine relationship between fire exposure and social vulnerability ----------

# Create categorical SVI variable for analysis
integrated_data$svi_cat <- cut(integrated_data$svi_score, 
                              breaks = c(0, 0.25, 0.5, 0.75, 1),
                              labels = c("Low", "Medium-Low", "Medium-High", "High"))

# Create summary table of SVI categories by fire exposure
svi_fire_table <- integrated_data %>%
  st_drop_geometry() %>%
  group_by(svi_cat, in_fire_10km) %>%
  summarize(
    n_tracts = n(),
    mean_asthma = mean(asthma_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_tracts = n_tracts / sum(n_tracts) * 100,
    fire_status = ifelse(in_fire_10km == 1, "Fire-Affected", "Not Fire-Affected")
  )

# Plot relationship between SVI category and fire exposure
p4 <- ggplot(svi_fire_table, aes(x = svi_cat, y = mean_asthma, fill = fire_status)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Not Fire-Affected" = "steelblue", "Fire-Affected" = "firebrick