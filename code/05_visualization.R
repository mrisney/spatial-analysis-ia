# 05_visualization.R
# Advanced Visualization for LA Wildfire Health Impact Study
#
# This script:
# 1. Creates publication-quality maps for key findings
# 2. Generates interactive visualizations for exploration
# 3. Produces composite visualizations combining multiple variables
# 4. Creates animated visualizations showing temporal patterns
# 5. Generates dashboards summarizing key results

# Load required libraries
library(sf)             # For spatial data handling
library(dplyr)          # For data manipulation
library(ggplot2)        # For data visualization
library(tmap)           # For thematic maps
library(leaflet)        # For interactive maps
library(plotly)         # For interactive plots
library(RColorBrewer)   # For color palettes
library(viridis)        # For colorblind-friendly palettes
library(gridExtra)      # For arranging multiple plots
library(cowplot)        # For plot composition
library(patchwork)      # For combining ggplots
library(gganimate)      # For animated visualizations
library(kableExtra)     # For nice tables
library(scales)         # For formatting numbers and scales
library(raster)         # For raster operations

# Set up directories
root_dir <- getwd()  # Assumes script is run from the repository root
data_dir <- file.path(root_dir, "data")
processed_dir <- file.path(data_dir, "processed")
results_dir <- file.path(root_dir, "results")
fig_dir <- file.path(results_dir, "figures")
map_dir <- file.path(results_dir, "maps")
viz_dir <- file.path(results_dir, "visualizations")
html_dir <- file.path(results_dir, "html")

# Create directories if they don't exist
dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(html_dir, showWarnings = FALSE, recursive = TRUE)

# Load data from previous analyses
cat("Loading processed data...\n")
model_data <- st_read(file.path(processed_dir, "la_model_results.shp"))
fires_2025 <- st_read(file.path(processed_dir, "fires_2025.shp"))

# Also load historical comparison data if available
if(file.exists(file.path(processed_dir, "historical_tracts_data.shp"))) {
  historical_data <- st_read(file.path(processed_dir, "historical_tracts_data.shp"))
  has_historical <- TRUE
} else {
  has_historical <- FALSE
}

# Set theme for all ggplot visualizations
theme_set(theme_minimal(base_size = 12))
custom_theme <- theme(
  plot.title = element_text(face = "bold", size = 14),
  plot.subtitle = element_text(size = 12, color = "gray30"),
  legend.title = element_text(face = "bold", size = 10),
  legend.text = element_text(size = 9),
  axis.title = element_text(face = "bold"),
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "gray90", color = NA),
  strip.text = element_text(face = "bold", size = 10)
)

# 1. Publication-quality maps for key findings --------------------------------

# Setup for publication-quality maps
tmap_mode("plot")  # Set tmap to static plotting mode

# Custom color palettes
asthma_palette <- colorRampPalette(c("#ffffcc", "#fed976", "#fd8d3c", "#bd0026"))(7)
svi_palette <- colorRampPalette(c("#f1eef6", "#bdc9e1", "#74a9cf", "#0570b0", "#034e7b"))(7)
interaction_palette <- colorRampPalette(c("#2166ac", "#d1e5f0", "#f7f7f7", "#fddbc7", "#b2182b"))(9)

# Function to create a map with consistent styling
create_styled_map <- function(data, variable, palette, title, legend_title, 
                             style = "quantile", n = 7, midpoint = NULL,
                             show_fires = TRUE) {
  # Base map with census tracts
  map <- tm_shape(data) +
    tm_fill(variable, 
            style = style, 
            palette = palette, 
            title = legend_title,
            n = n,
            midpoint = midpoint) +
    tm_borders(alpha = 0.2)
  
  # Add fire perimeters if requested
  if(show_fires) {
    map <- map +
      tm_shape(fires_2025) +
      tm_borders(col = "red", lwd = 2)
  }
  
  # Add layout elements
  map <- map +
    tm_layout(
      title = title,
      frame = FALSE,
      legend.outside = FALSE,
      legend.position = c("right", "bottom"),
      legend.title.size = 0.8,
      legend.text.size = 0.6,
      legend.frame = FALSE
    ) +
    tm_compass(position = c("left", "top")) +
    tm_scale_bar(position = c("right", "bottom"))
  
  return(map)
}

# Map 1: Social Vulnerability Index with 2025 Fire Overlay
map1 <- create_styled_map(
  data = model_data,
  variable = "svi_score",
  palette = svi_palette,
  title = "Social Vulnerability in LA County and 2025 Fires",
  legend_title = "SVI Score"
)
tmap_save(map1, file.path(map_dir, "svi_with_fires_map.png"), width = 8, height = 8, dpi = 300)

# Map 2: Asthma ED Rates with 2025 Fire Overlay
map2 <- create_styled_map(
  data = model_data,
  variable = "asthma_rate",
  palette = asthma_palette,
  title = "Asthma ED Rates in LA County and 2025 Fires",
  legend_title = "Asthma ED Rate\n(per 10,000)"
)
tmap_save(map2, file.path(map_dir, "asthma_with_fires_map.png"), width = 8, height = 8, dpi = 300)

# Map 3: Fire Impact Zones (multiple buffers)
buffer_colors <- c("#fee8c8", "#fdbb84", "#e34a33")
map3 <- tm_shape(model_data) +
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
    title = "2025 LA Fires Impact Zones",
    frame = FALSE,
    legend.outside = FALSE,
    legend.position = c("right", "bottom"),
    legend.title.size = 0.8,
    legend.text.size = 0.6,
    legend.frame = FALSE
  ) +
  tm_compass(position = c("left", "top")) +
  tm_scale_bar(position = c("right", "bottom")) +
  tm_add_legend(type = "fill", 
                col = c("red", buffer_colors[3], buffer_colors[2], buffer_colors[1]),
                labels = c("Fire Perimeter", "5km Buffer", "10km Buffer", "20km Buffer"),
                alpha = c(0.7, 0.5, 0.4, 0.3))
tmap_save(map3, file.path(map_dir, "fire_impact_zones_map.png"), width = 8, height = 8, dpi = 300)

# Map 4: Model Predicted Asthma Rates
map4 <- create_styled_map(
  data = model_data,
  variable = "predicted_",  # Adjust column name based on actual data
  palette = asthma_palette,
  title = "Predicted Asthma ED Rates from Spatial Model",
  legend_title = "Predicted Rate\n(per 10,000)"
)
tmap_save(map4, file.path(map_dir, "predicted_asthma_map.png"), width = 8, height = 8, dpi = 300)

# Map 5: Local Coefficients from Geographically Weighted Regression (if available)
if("gwr_svi" %in% names(model_data)) {
  map5 <- create_styled_map(
    data = model_data,
    variable = "gwr_intera",  # Adjust column name based on actual data
    palette = interaction_palette,
    title = "Local Coefficient: SVI × Fire Interaction",
    legend_title = "Coefficient",
    style = "quantile",
    midpoint = 0
  )
  tmap_save(map5, file.path(map_dir, "gwr_interaction_map.png"), width = 8, height = 8, dpi = 300)
}

# Map 6: Spatial Clustering (LISA) of Asthma Rates
if("asthma_clu" %in% names(model_data)) {
  # Create categorical legend colors
  lisa_colors <- c("High-High" = "red", "Low-Low" = "blue", 
                  "High-Low" = "pink", "Low-High" = "lightblue", 
                  "Not Significant" = "white")
  
  map6 <- tm_shape(model_data) +
    tm_fill("asthma_clu",
            palette = lisa_colors,
            title = "Spatial Clusters") +
    tm_borders(alpha = 0.2) +
    tm_shape(fires_2025) +
    tm_borders(col = "black", lwd = 1.5) +
    tm_layout(
      title = "Spatial Clusters of Asthma ED Rates",
      frame = FALSE,
      legend.outside = FALSE,
      legend.position = c("right", "bottom"),
      legend.title.size = 0.8,
      legend.text.size = 0.6,
      legend.frame = FALSE
    ) +
    tm_compass(position = c("left", "top")) +
    tm_scale_bar(position = c("right", "bottom"))
  
  tmap_save(map6, file.path(map_dir, "asthma_clusters_map.png"), width = 8, height = 8, dpi = 300)
}

# 2. Interactive visualizations for exploration -------------------------------

# Switch tmap to interactive mode
tmap_mode("view")

# Interactive map of fire impacts and health outcomes
interactive_map <- tm_shape(model_data) +
  tm_fill("asthma_rate", 
          style = "quantile",
          palette = "YlOrRd",
          alpha = 0.7,
          id = "NAME",
          popup.vars = c("Asthma ED Rate" = "asthma_rate",
                         "SVI Score" = "svi_score",
                         "In 10km of Fire" = "in_fire_10km",
                         "Distance to Fire (km)" = "dist_to_fi")) +
  tm_borders(alpha = 0.2) +
  tm_shape(fires_2025) +
  tm_borders(col = "red", lwd = 2) +
  tm_shape(st_buffer(fires_2025, dist = 10000)) +
  tm_borders(col = "orange", lwd = 1.5, alpha = 0.5)

# Save interactive map
tmap_save(interactive_map, file.path(html_dir, "interactive_asthma_map.html"))

# Create an interactive scatterplot of SVI vs asthma rates colored by fire proximity
# Convert to regular data frame for plotly
plot_data <- model_data %>%
  st_drop_geometry() %>%
  mutate(fire_status = ifelse(in_fire_10km == 1, "Within 10km of Fire", "Beyond 10km of Fire"))

p <- ggplot(plot_data, aes(x = svi_score, y = asthma_rate, color = fire_status)) +
  geom_point(alpha = 0.7, aes(text = paste0("Tract: ", NAME, 
                                           "<br>SVI: ", round(svi_score, 2),
                                           "<br>Asthma Rate: ", round(asthma_rate, 1),
                                           "<br>Distance to Fire: ", round(dist_to_fi, 1), " km"))) +
  geom_smooth(method = "lm", se = TRUE, aes(group = fire_status)) +
  scale_color_manual(values = c("Within 10km of Fire" = "firebrick", 
                              "Beyond 10km of Fire" = "steelblue")) +
  labs(
    title = "Relationship Between Social Vulnerability and Asthma ED Rates",
    x = "Social Vulnerability Index Score",
    y = "Asthma ED Rate (per 10,000)",
    color = "Fire Proximity"
  ) +
  custom_theme

# Convert to interactive plotly
p_interactive <- ggplotly(p, tooltip = "text")

# Save interactive plot
htmlwidgets::saveWidget(p_interactive, file.path(html_dir, "interactive_svi_asthma_scatter.html"))

# 3. Composite visualizations combining multiple variables --------------------

# Create a bivariate choropleth map of SVI and asthma rates
# First, classify both variables into categories
model_data$svi_cat <- cut(model_data$svi_score, 
                        breaks = quantile(model_data$svi_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE),
                        labels = c("Low", "Medium", "High"),
                        include.lowest = TRUE)

model_data$asthma_cat <- cut(model_data$asthma_rate, 
                           breaks = quantile(model_data$asthma_rate, probs = seq(0, 1, by = 1/3), na.rm = TRUE),
                           labels = c("Low", "Medium", "High"),
                           include.lowest = TRUE)

# Create bivariate categories
model_data$bivar_cat <- paste(as.character(model_data$svi_cat), 
                            as.character(model_data$asthma_cat))

# Define a 3x3 bivariate color scheme
bivar_colors <- c(
  "Low Low" = "#e8e8e8",       # Light gray
  "Low Medium" = "#abd9e9",    # Light blue
  "Low High" = "#2c7fb8",      # Dark blue
  "Medium Low" = "#fdae61",    # Light orange
  "Medium Medium" = "#9eabd3", # Purple-blue
  "Medium High" = "#225ea8",   # Darker blue
  "High Low" = "#d7191c",      # Red
  "High Medium" = "#976b82",   # Purple
  "High High" = "#542788"      # Dark purple
)

# Create bivariate map
bivar_map <- tm_shape(model_data) +
  tm_fill("bivar_cat", 
          palette = bivar_colors,
          title = "SVI & Asthma ED Rates") +
  tm_borders(alpha = 0.2) +
  tm_shape(fires_2025) +
  tm_borders(col = "black", lwd = 1.5) +
  tm_layout(
    title = "Social Vulnerability and Asthma ED Rates",
    frame = FALSE,
    legend.outside = TRUE,
    legend.position = c("right", "center"),
    legend.title.size = 0.8,
    legend.text.size = 0.6
  )

# Add custom bivariate legend
bivar_legend <- tm_add_legend(
  type = "fill",
  labels = c("Low SVI, Low Asthma", "Low SVI, Med Asthma", "Low SVI, High Asthma",
            "Med SVI, Low Asthma", "Med SVI, Med Asthma", "Med SVI, High Asthma",
            "High SVI, Low Asthma", "High SVI, Med Asthma", "High SVI, High Asthma"),
  col = bivar_colors,
  title = "SVI & Asthma ED Rates"
)

# Save bivariate map
tmap_save(bivar_map, file.path(map_dir, "bivariate_svi_asthma_map.png"), width = 8, height = 8, dpi = 300)

# Create a faceted map series showing fire exposure by SVI category
facet_data <- model_data %>%
  mutate(svi_quartile = cut(svi_score, 
                          breaks = quantile(svi_score, probs = seq(0, 1, 0.25), na.rm = TRUE),
                          labels = c("Low SVI (Q1)", "Medium-Low SVI (Q2)", 
                                   "Medium-High SVI (Q3)", "High SVI (Q4)"),
                          include.lowest = TRUE))

facet_map <- tm_shape(facet_data) +
  tm_fill("asthma_rate", 
          style = "quantile",
          palette = asthma_palette,
          title = "Asthma ED Rate") +
  tm_borders(alpha = 0.2) +
  tm_facets(by = "svi_quartile", ncol = 2) +
  tm_shape(fires_2025) +
  tm_borders(col = "red", lwd = 1.5) +
  tm_layout(
    title = "Asthma ED Rates by Social Vulnerability Quartile",
    frame = FALSE,
    legend.outside = TRUE,
    legend.position = c("right", "center"),
    panel.labels = c("Low SVI (Q1)", "Medium-Low SVI (Q2)", 
                    "Medium-High SVI (Q3)", "High SVI (Q4)"),
    panel.label.size = 1
  )

# Save faceted map
tmap_save(facet_map, file.path(map_dir, "faceted_svi_asthma_map.png"), width = 10, height = 10, dpi = 300)

# 4. Advanced visualizations showing model outputs ---------------------------

# Plot showing the interaction effect between fire proximity and SVI
# Use stratified coefficients if available

# Check if coefficients from stratified analysis exist (from 03_spatial_analysis.R)
if(file.exists(file.path(results_dir, "models", "model_comparison.csv"))) {
  # Read model comparison results
  model_comparison <- read.csv(file.path(results_dir, "models", "model_comparison.csv"))
  
  # Create visualization of model comparison
  p_models <- ggplot(model_comparison, aes(x = Description, y = AIC)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(
      title = "Model Comparison",
      subtitle = "Lower AIC indicates better model fit",
      x = "",
      y = "Akaike Information Criterion (AIC)"
    ) +
    custom_theme
  
  ggsave(file.path(fig_dir, "model_comparison.png"), p_models, width = 10, height = 6, dpi = 300)
}

# Check if stratified model results exist
if(file.exists(file.path(results_dir, "models", "stratified_fire_effects.csv"))) {
  # Read stratified model results
  stratified_effects <- read.csv(file.path(results_dir, "models", "stratified_fire_effects.csv"))
  
  # Create visualization of stratified effects
  p_strat <- ggplot(stratified_effects, aes(x = SVI_Category, y = Fire_Proximity_Coef)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    geom_errorbar(aes(ymin = Fire_Proximity_Coef - 1.96*SE, 
                     ymax = Fire_Proximity_Coef + 1.96*SE), 
                 width = 0.2) +
    labs(
      title = "Effect of Fire Proximity on Asthma Rates by SVI Category",
      subtitle = "Coefficient estimates with 95% confidence intervals",
      x = "Social Vulnerability Index Category",
      y = "Coefficient for 10km Fire Proximity"
    ) +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(fig_dir, "stratified_fire_effects.png"), p_strat, width = 10, height = 6, dpi = 300)
}

# 5. Animated visualization of changing impacts over time -------------------

# Check if historical comparison data exists
if(has_historical) {
  # Create a dataset for animation that combines current and historical data
  animation_data <- bind_rows(
    # Current data
    model_data %>%
      st_drop_geometry() %>%
      select(NAME, geoid, asthma_rate, svi_score, in_fire_10km) %>%
      mutate(time_period = "2025"),
    
    # Historical data
    historical_data %>%
      st_drop_geometry() %>%
      select(NAME, geoid, asthma_rate, svi_score, in_hist_fire_10km) %>%
      rename(in_fire_10km = in_hist_fire_10km) %>%
      mutate(time_period = "Historical")
  )
  
  # Create animated plot
  p_animated <- ggplot(animation_data, aes(x = svi_score, y = asthma_rate, color = factor(in_fire_10km))) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = c("0" = "steelblue", "1" = "firebrick"),
                      labels = c("0" = "Not Fire-Affected", "1" = "Fire-Affected"),
                      name = "Fire Exposure") +
    labs(
      title = "Relationship Between Social Vulnerability and Asthma Rates",
      subtitle = "Time Period: {closest_state}",
      x = "Social Vulnerability Index Score",
      y = "Asthma ED Rate (per 10,000)"
    ) +
    custom_theme +
    transition_states(
      time_period,
      transition_length = 2,
      state_length = 3
    ) +
    enter_fade() +
    exit_fade()
  
  # Save animation
  anim_save(file.path(viz_dir, "fire_impact_animation.gif"), p_animated, width = 8, height = 6, units = "in", res = 150)
}

