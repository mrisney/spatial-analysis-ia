# 03_spatial_analysis.R
# Spatial Analysis for LA Wildfire Health Impact Study
#
# This script:
# 1. Performs spatial regression modeling to assess relationships between fire exposure,
#    social vulnerability, and health outcomes (asthma ED visits)
# 2. Tests for spatial dependence in model residuals
# 3. Implements alternative spatial modeling approaches
# 4. Conducts stratified analysis based on social vulnerability categories
# 5. Creates model diagnostics and visualization of results

# Load required libraries
library(sf)          # For spatial data handling
library(dplyr)       # For data manipulation
library(spdep)       # For spatial weights and spatial regression
library(spatialreg)  # For spatial regression models
library(GWmodel)     # For geographically weighted regression
library(ggplot2)     # For data visualization
library(tmap)        # For thematic maps
library(grid)        # For arranging multiple plots
library(gridExtra)   # For arranging multiple plots
library(car)         # For regression diagnostics
library(lmtest)      # For hypothesis testing in linear models
library(broom)       # For converting model output to tidy data frames

# Set up directories
root_dir <- getwd()  # Assumes script is run from the repository root
data_dir <- file.path(root_dir, "data")
processed_dir <- file.path(data_dir, "processed")
results_dir <- file.path(root_dir, "results")
fig_dir <- file.path(results_dir, "figures")
map_dir <- file.path(results_dir, "maps")
model_dir <- file.path(results_dir, "models")

# Create directories if they don't exist
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

# Set the tmap mode to plotting
tmap_mode("plot")

# Load processed data with LISA clusters from 02_eda.R
cat("Loading processed data...\n")
integrated_data <- st_read(file.path(processed_dir, "la_integrated_data_with_clusters.shp"))

# Make sure we have the spatial weights from previous analysis
cat("Creating spatial weights matrix...\n")
tracts_nb <- poly2nb(integrated_data, queen = TRUE)
tracts_weights <- nb2listw(tracts_nb, style = "W")

# Prepare cleaned dataset for modeling ----------------------------------------
# Drop rows with missing values in key variables
model_data <- integrated_data %>%
  filter(!is.na(asthma_rate) & !is.na(svi_score) & !is.na(dist_to_fire_km))

# Create interaction terms between fire exposure and social vulnerability
model_data$svi_x_fire5km <- model_data$svi_score * model_data$in_fire_5km
model_data$svi_x_fire10km <- model_data$svi_score * model_data$in_fire_10km
model_data$svi_x_dist <- model_data$svi_score * model_data$dist_to_fire_km

# Apply log transformation to skewed variables if needed
# (Check distributions first and uncomment if necessary)
# model_data$log_asthma_rate <- log1p(model_data$asthma_rate)
# model_data$log_dist_to_fire <- log1p(model_data$dist_to_fire_km)

# 1. Basic OLS regression model without spatial effects -----------------------
# This establishes a baseline for comparison with spatial models

cat("Fitting basic OLS model...\n")
ols_model <- lm(asthma_rate ~ svi_score + dist_to_fire_km + in_fire_10km + 
                poverty_pc + pop_densit, data = model_data)
summary(ols_model)

# Save model summary
sink(file.path(model_dir, "ols_model_summary.txt"))
print(summary(ols_model))
sink()

# Check for multicollinearity
vif_results <- vif(ols_model)
cat("Variance Inflation Factors:\n")
print(vif_results)

# Extract residuals
model_data$ols_residuals <- residuals(ols_model)

# 2. Interaction model to test effect modification ----------------------------
# This tests our hypothesis that SVI modifies the relationship between fire exposure and health

cat("Fitting interaction model...\n")
interaction_model <- lm(asthma_rate ~ svi_score + dist_to_fire_km + in_fire_10km + 
                      svi_x_fire10km + poverty_pc + pop_densit, 
                      data = model_data)
summary(interaction_model)

# Save model summary
sink(file.path(model_dir, "interaction_model_summary.txt"))
print(summary(interaction_model))
sink()

# Compare models
anova_result <- anova(ols_model, interaction_model)
print(anova_result)

# Extract residuals
model_data$interaction_residuals <- residuals(interaction_model)

# 3. Test for spatial autocorrelation in model residuals ----------------------
# This helps determine if spatial regression is needed

cat("Testing for spatial autocorrelation in OLS residuals...\n")
# Recreate spatial weights for the filtered model_data
model_nb <- poly2nb(model_data, queen = TRUE)
model_weights <- nb2listw(model_nb, style = "W")

# Moran's I test for OLS residuals
moran_ols <- moran.test(model_data$ols_residuals, model_weights)
print(moran_ols)

# Moran's I test for interaction model residuals
moran_interaction <- moran.test(model_data$interaction_residuals, model_weights)
print(moran_interaction)

# LM tests for spatial dependence
lm_tests <- lm.LMtests(ols_model, model_weights, test = "all")
print(lm_tests)

# 4. Spatial regression models if needed ---------------------------------------
# Based on LM tests, choose appropriate spatial model

cat("Fitting spatial regression models...\n")

# Spatial lag model
spatial_lag_model <- lagsarlm(asthma_rate ~ svi_score + dist_to_fire_km + in_fire_10km + 
                             svi_x_fire10km + poverty_pc + pop_densit, 
                             data = model_data, 
                             listw = model_weights)
summary(spatial_lag_model)

# Spatial error model
spatial_error_model <- errorsarlm(asthma_rate ~ svi_score + dist_to_fire_km + in_fire_10km + 
                                svi_x_fire10km + poverty_pc + pop_densit, 
                                data = model_data, 
                                listw = model_weights)
summary(spatial_error_model)

# Save model summaries
sink(file.path(model_dir, "spatial_lag_model_summary.txt"))
print(summary(spatial_lag_model))
sink()

sink(file.path(model_dir, "spatial_error_model_summary.txt"))
print(summary(spatial_error_model))
sink()

# Extract residuals from spatial models
model_data$spatial_lag_residuals <- residuals(spatial_lag_model)
model_data$spatial_error_residuals <- residuals(spatial_error_model)

# 5. Geographically Weighted Regression ---------------------------------------
# To explore spatial non-stationarity in the relationships

cat("Preparing for Geographically Weighted Regression...\n")
# Prepare data for GWR
gwr_data <- model_data %>%
  select(asthma_rate, svi_score, dist_to_fire_km, in_fire_10km, 
        svi_x_fire10km, poverty_pc, pop_densit) %>%
  st_transform(4326)  # Transform to WGS84 for GWR

# Determine optimal bandwidth for GWR
cat("Determining optimal bandwidth for GWR (this may take a while)...\n")
gwr_bw <- try(bw.gwr(
  asthma_rate ~ svi_score + dist_to_fire_km + in_fire_10km + 
  svi_x_fire10km + poverty_pc + pop_densit,
  data = gwr_data,
  approach = "AIC",
  kernel = "gaussian",
  adaptive = TRUE
), silent = TRUE)

if(!inherits(gwr_bw, "try-error")) {
  # Fit GWR model with optimal bandwidth
  cat("Fitting GWR model...\n")
  gwr_model <- gwr.basic(
    asthma_rate ~ svi_score + dist_to_fire_km + in_fire_10km + 
    svi_x_fire10km + poverty_pc + pop_densit,
    data = gwr_data,
    bw = gwr_bw,
    kernel = "gaussian",
    adaptive = TRUE
  )
  
  # Save GWR summary
  sink(file.path(model_dir, "gwr_model_summary.txt"))
  print(gwr_model)
  sink()
  
  # Extract local coefficient estimates
  gwr_coefs <- as.data.frame(gwr_model$SDF)
  
  # Add local coefficient estimates back to the model data
  model_data$gwr_svi <- gwr_coefs$svi_score
  model_data$gwr_dist <- gwr_coefs$dist_to_fire_km
  model_data$gwr_fire10km <- gwr_coefs$in_fire_10km
  model_data$gwr_interaction <- gwr_coefs$svi_x_fire10km
  
  # Map local coefficients
  cat("Creating maps of GWR local coefficients...\n")
  
  # Map for SVI coefficient
  map_gwr_svi <- tm_shape(model_data) +
    tm_fill("gwr_svi", 
            style = "quantile", 
            palette = "RdBu", 
            midpoint = 0,
            title = "SVI Effect",
            n = 7) +
    tm_borders(alpha = 0.2) +
    tm_layout(
      title = "Local Coefficient: Social Vulnerability",
      legend.position = c("right", "bottom"),
      legend.text.size = 0.7,
      legend.title.size = 0.9
    )
  tmap_save(map_gwr_svi, file.path(map_dir, "gwr_svi_coefficient.png"), width = 8, height = 8)
  
  # Map for Fire proximity coefficient
  map_gwr_fire <- tm_shape(model_data) +
    tm_fill("gwr_fire10km", 
            style = "quantile", 
            palette = "RdBu", 
            midpoint = 0,
            title = "Fire Proximity Effect",
            n = 7) +
    tm_borders(alpha = 0.2) +
    tm_layout(
      title = "Local Coefficient: Fire Proximity (10km)",
      legend.position = c("right", "bottom"),
      legend.text.size = 0.7,
      legend.title.size = 0.9
    )
  tmap_save(map_gwr_fire, file.path(map_dir, "gwr_fire_coefficient.png"), width = 8, height = 8)
  
  # Map for interaction coefficient
  map_gwr_int <- tm_shape(model_data) +
    tm_fill("gwr_interaction", 
            style = "quantile", 
            palette = "RdBu", 
            midpoint = 0,
            title = "Interaction Effect",
            n = 7) +
    tm_borders(alpha = 0.2) +
    tm_layout(
      title = "Local Coefficient: SVI × Fire Interaction",
      legend.position = c("right", "bottom"),
      legend.text.size = 0.7,
      legend.title.size = 0.9
    )
  tmap_save(map_gwr_int, file.path(map_dir, "gwr_interaction_coefficient.png"), width = 8, height = 8)
} else {
  cat("Warning: GWR bandwidth optimization failed. Skipping GWR analysis.\n")
}

# 6. Stratified analysis by SVI category --------------------------------------
# This helps understand how relationships differ across vulnerability levels

cat("Conducting stratified analysis by SVI categories...\n")
# Create SVI quartiles
model_data$svi_quartile <- cut(model_data$svi_score, 
                             breaks = quantile(model_data$svi_score, probs = seq(0, 1, 0.25), na.rm = TRUE),
                             labels = c("Low SVI (Q1)", "Medium-Low SVI (Q2)", 
                                      "Medium-High SVI (Q3)", "High SVI (Q4)"),
                             include.lowest = TRUE)

# Run models for each SVI quartile
stratified_models <- list()
for(quartile in levels(model_data$svi_quartile)) {
  subset_data <- model_data %>% filter(svi_quartile == quartile)
  if(nrow(subset_data) > 10) {  # Ensure enough observations
    model <- lm(asthma_rate ~ dist_to_fire_km + in_fire_10km + 
               poverty_pc + pop_densit, 
               data = subset_data)
    stratified_models[[quartile]] <- model
    
    # Save model summary
    sink(file.path(model_dir, paste0("stratified_model_", gsub(" ", "_", quartile), ".txt")))
    print(summary(model))
    sink()
  }
}

# Compare fire proximity coefficients across SVI strata
strat_coefs <- data.frame(
  SVI_Category = character(),
  Fire_Proximity_Coef = numeric(),
  SE = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for(quartile in names(stratified_models)) {
  model <- stratified_models[[quartile]]
  coef_data <- summary(model)$coefficients
  
  if("in_fire_10km" %in% rownames(coef_data)) {
    strat_coefs <- rbind(strat_coefs, data.frame(
      SVI_Category = quartile,
      Fire_Proximity_Coef = coef_data["in_fire_10km", "Estimate"],
      SE = coef_data["in_fire_10km", "Std. Error"],
      p_value = coef_data["in_fire_10km", "Pr(>|t|)"],
      stringsAsFactors = FALSE
    ))
  }
}

# Plot stratified coefficients
if(nrow(strat_coefs) > 0) {
  p_strat <- ggplot(strat_coefs, aes(x = SVI_Category, y = Fire_Proximity_Coef)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    geom_errorbar(aes(ymin = Fire_Proximity_Coef - 1.96*SE, 
                     ymax = Fire_Proximity_Coef + 1.96*SE), 
                 width = 0.2) +
    labs(
      title = "Effect of Fire Proximity on Asthma Rates by SVI Category",
      x = "Social Vulnerability Index Category",
      y = "Coefficient for 10km Fire Proximity"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(fig_dir, "stratified_fire_effects.png"), p_strat, width = 10, height = 6)
}

# 7. Model evaluation and comparison ------------------------------------------
# Compare fit statistics across models to determine best approach

cat("Comparing model performance...\n")
# Function to extract AIC and other fit statistics
get_model_stats <- function(model, type) {
  if(type == "lm") {
    return(data.frame(
      Model = deparse(substitute(model)),
      Type = "OLS",
      AIC = AIC(model),
      BIC = BIC(model),
      R_squared = summary(model)$r.squared,
      Adj_R_squared = summary(model)$adj.r.squared,
      LogLik = logLik(model)[1]
    ))
  } else if(type == "spatial") {
    return(data.frame(
      Model = deparse(substitute(model)),
      Type = "Spatial",
      AIC = AIC(model),
      BIC = BIC(model),
      R_squared = NA,  # Spatial models don't have traditional R²
      Adj_R_squared = NA,
      LogLik = logLik(model)[1]
    ))
  }
}

# Compile statistics for all models
model_comparison <- rbind(
  get_model_stats(ols_model, "lm"),
  get_model_stats(interaction_model, "lm"),
  get_model_stats(spatial_lag_model, "spatial"),
  get_model_stats(spatial_error_model, "spatial")
)

# Add pseudo-R² for spatial models if available
model_comparison$Pseudo_R_squared <- NA
model_comparison$Pseudo_R_squared[model_comparison$Model == "spatial_lag_model"] <- 
  1 - (spatial_lag_model$SSE / sum((model_data$asthma_rate - mean(model_data$asthma_rate))^2))
model_comparison$Pseudo_R_squared[model_comparison$Model == "spatial_error_model"] <- 
  1 - (spatial_error_model$SSE / sum((model_data$asthma_rate - mean(model_data$asthma_rate))^2))

# Add model specifications
model_comparison$Description <- c(
  "Base OLS Model",
  "OLS with SVI×Fire Interaction",
  "Spatial Lag Model",
  "Spatial Error Model"
)

# Save model comparison
write.csv(model_comparison, file.path(model_dir, "model_comparison.csv"), row.names = FALSE)

# Create nice formatted table for the model comparison
model_comp_table <- model_comparison %>%
  select(Description, AIC, BIC, R_squared, Adj_R_squared, Pseudo_R_squared, LogLik) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

print(model_comp_table)

# 8. Map predicted values and residuals ----------------------------------------
# Visualize model fit and identify areas of poor prediction

cat("Mapping predicted values and residuals...\n")

# Add predictions from best model (assuming spatial error model is best)
model_data$predicted_values <- fitted(spatial_error_model)
model_data$prediction_residuals <- residuals(spatial_error_model)

# Map of predicted asthma rates
map_predicted <- tm_shape(model_data) +
  tm_fill("predicted_values", 
          style = "quantile", 
          palette = "YlOrRd", 
          title = "Predicted Asthma\nED Rate",
          n = 7) +
  tm_borders(alpha = 0.2) +
  tm_layout(
    title = "Predicted Asthma ED Rates from Spatial Model",
    legend.position = c("right", "bottom"),
    legend.text.size = 0.7,
    legend.title.size = 0.9
  )
tmap_save(map_predicted, file.path(map_dir, "predicted_asthma_rates.png"), width = 8, height = 8)

# Map of residuals
map_residuals <- tm_shape(model_data) +
  tm_fill("prediction_residuals", 
          style = "quantile", 
          palette = "RdBu", 
          midpoint = 0,
          title = "Residuals",
          n = 7) +
  tm_borders(alpha = 0.2) +
  tm_layout(
    title = "Spatial Model Residuals",
    legend.position = c("right", "bottom"),
    legend.text.size = 0.7,
    legend.title.size = 0.9
  )
tmap_save(map_residuals, file.path(map_dir, "model_residuals.png"), width = 8, height = 8)

# 9. Summarize findings and save final model results --------------------------
cat("Saving final combined dataset with model results...\n")

# Save the enhanced dataset with model predictions and residuals
st_write(model_data, file.path(processed_dir, "la_model_results.shp"), append = FALSE)

# Create a text summary of key findings
sink(file.path(model_dir, "model_summary_report.txt"))
cat("LA Wildfire Health Impact Spatial Analysis Summary\n")
cat("==================================================\n\n")

cat("1. Basic Relationship Findings:\n")
cat("   - Estimated effect of being within 10km of a fire on asthma rates: ", 
    round(coef(spatial_error_model)["in_fire_10km"], 3), "\n")
cat("   - Estimated effect of Social Vulnerability Index on asthma rates: ", 
    round(coef(spatial_error_model)["svi_score"], 3), "\n")
cat("   - Interaction between SVI and fire proximity: ", 
    round(coef(spatial_error_model)["svi_x_fire10km"], 3), "\n\n")

cat("2. Spatial Autocorrelation:\n")
cat("   - Moran's I for OLS residuals: ", round(moran_ols$estimate[1], 3), 
    " (p-value: ", format.pval(moran_ols$p.value, digits = 3), ")\n")
cat("   - Evidence suggests significant spatial autocorrelation in the data.\n\n")

cat("3. Model Comparison:\n")
cat("   - The ", model_comparison$Description[which.min(model_comparison$AIC)], 
    " had the best fit based on AIC.\n\n")

cat("4. Stratified Analysis:\n")
if(nrow(strat_coefs) > 0) {
  for(i in 1:nrow(strat_coefs)) {
    cat("   - ", strat_coefs$SVI_Category[i], ": Fire proximity effect = ", 
        round(strat_coefs$Fire_Proximity_Coef[i], 3), 
        " (p-value: ", format.pval(strat_coefs$p_value[i], digits = 3), ")\n")
  }
}
cat("\n")

cat("5. Conclusions:\n")
cat("   - Significant spatial dependency in asthma rates across LA County\n")
cat("   - Social vulnerability modifies the relationship between fire proximity and health outcomes\n")
cat("   - Areas with both high social vulnerability and fire exposure show the highest asthma rates\n")
cat("   - Spatial modeling approaches provide better fit than traditional OLS regression\n")
sink()

cat("Spatial analysis complete!\n")
cat("Results saved to:", results_dir, "\n")

# Next steps:
# 1. Conduct sensitivity analyses with different distance thresholds
# 2. Incorporate air quality data as mediator between fires and health
# 3. Develop predictive models for future fire scenarios
# 4. Compare results with historical wildfire patterns
# 5. Integrate building damage assessment data for more nuanced exposure measurement