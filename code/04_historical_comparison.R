# 04_historical_comparison.R
# Historical Comparison Analysis for LA Wildfire Health Impact Study
#
# This script:
# 1. Identifies historical fires similar to the 2025 LA fires
# 2. Compares health impacts across different time periods
# 3. Examines changes in vulnerability patterns over time
# 4. Implements matching methods to control for confounding
# 5. Analyzes temporal trends in fire-health relationships

# Load required libraries
library(sf)          # For spatial data handling
library(dplyr)       # For data manipulation
library(ggplot2)     # For data visualization
library(tmap)        # For thematic maps
library(lubridate)   # For date handling
library(MatchIt)     # For propensity score matching
library(gridExtra)   # For arranging multiple plots
library(purrr)       # For functional programming
library(lmtest)      # For hypothesis testing
library(sandwich)    # For robust standard errors
library(readr)       # For reading CSV files
library(stringr)     # For string manipulation

# Set up directories
root_dir <- getwd()  # Assumes script is run from the repository root
data_dir <- file.path(root_dir, "data")
processed_dir <- file.path(data_dir, "processed")
results_dir <- file.path(root_dir, "results")
fig_dir <- file.path(results_dir, "figures")
map_dir <- file.path(results_dir, "maps")
hist_dir <- file.path(results_dir, "historical")

# Create directories if they don't exist
dir.create(hist_dir, showWarnings = FALSE, recursive = TRUE)

# Set the tmap mode to plotting
tmap_mode("plot")

# Load data --------------------------------------------------------------------
cat("Loading current and historical fire data...\n")

# Load current (2025) data
current_data <- st_read(file.path(processed_dir, "la_model_results.shp"))
fires_2025 <- st_read(file.path(processed_dir, "fires_2025.shp"))

# Load historical fire perimeters
historical_fires <- st_read(file.path(data_dir, "fires", "historical_ca_fires.shp"))

# Load historical CalEnviroScreen data (if available)
# Note: You'll need to adjust these file paths and column names based on your actual data
# These are hypothetical examples showing the structure
if(file.exists(file.path(data_dir, "calenviroscreen", "calenviroscreen_3.0.csv"))) {
  ces_3_0 <- read_csv(file.path(data_dir, "calenviroscreen", "calenviroscreen_3.0.csv"))
} else {
  cat("Warning: Historical CalEnviroScreen 3.0 data not found.\n")
  ces_3_0 <- NULL
}

if(file.exists(file.path(data_dir, "calenviroscreen", "calenviroscreen_4.0.csv"))) {
  ces_4_0 <- read_csv(file.path(data_dir, "calenviroscreen", "calenviroscreen_4.0.csv"))
} else {
  cat("Warning: Historical CalEnviroScreen 4.0 data not found.\n")
  ces_4_0 <- NULL
}

# Load historical SVI data
if(file.exists(file.path(data_dir, "svi", "SVI2018_CALIFORNIA_tract.csv"))) {
  svi_2018 <- read_csv(file.path(data_dir, "svi", "SVI2018_CALIFORNIA_tract.csv"))
} else {
  cat("Warning: Historical SVI 2018 data not found.\n")
  svi_2018 <- NULL
}

if(file.exists(file.path(data_dir, "svi", "SVI2014_CALIFORNIA_tract.csv"))) {
  svi_2014 <- read_csv(file.path(data_dir, "svi", "SVI2014_CALIFORNIA_tract.csv"))
} else {
  cat("Warning: Historical SVI 2014 data not found.\n")
  svi_2014 <- NULL
}

# 1. Identify similar historical fires -----------------------------------------
cat("Identifying historical fires similar to 2025 LA fires...\n")

# Extract characteristics of the 2025 Eaton and Palisades fires
eaton_fire <- fires_2025 %>% filter(FIRE_NAME == "Eaton")
palisades_fire <- fires_2025 %>% filter(FIRE_NAME == "Palisades")

# Calculate key metrics for the 2025 fires
eaton_area <- st_area(eaton_fire) %>% as.numeric()
palisades_area <- st_area(palisades_fire) %>% as.numeric()

# Convert to acres (1 sq meter = 0.000247105 acres)
eaton_acres <- eaton_area * 0.000247105
palisades_acres <- palisades_area * 0.000247105

cat("2025 Eaton Fire area:", round(eaton_acres), "acres\n")
cat("2025 Palisades Fire area:", round(palisades_acres), "acres\n")

# Prepare historical fire data
historical_fires <- historical_fires %>%
  # Create a date field if it exists in a different format
  mutate(
    fire_date = if("ALARM_DATE" %in% names(.)) as.Date(ALARM_DATE) else 
                if("DISCOVERY_" %in% names(.)) as.Date(DISCOVERY_) else 
                as.Date(paste0(YEAR_, "-01-01")),
    year = year(fire_date),
    month = month(fire_date),
    fire_size_acres = if("GIS_ACRES" %in% names(.)) GIS_ACRES else 
                      if("FIRE_SIZE" %in% names(.)) FIRE_SIZE else 0
  )

# Function to calculate similarity score between current and historical fires
calc_fire_similarity <- function(current_acres, hist_fire) {
  # Calculate normalized difference in acres (0-1 scale)
  acres_diff <- abs(current_acres - hist_fire$fire_size_acres) / max(current_acres, hist_fire$fire_size_acres)
  
  # Calculate month similarity (0-1 scale, with January being most similar)
  month_diff <- abs(1 - hist_fire$month) / 11  # January is month 1
  
  # Combined similarity score (lower is more similar)
  similarity_score <- acres_diff * 0.7 + month_diff * 0.3
  
  return(similarity_score)
}

# Find historical fires similar to Eaton Fire (focus on pre-2018 for historical comparison)
similar_to_eaton <- historical_fires %>%
  filter(year < 2018, year >= 2010, fire_size_acres > 0) %>%
  mutate(similarity_score = sapply(1:n(), function(i) calc_fire_similarity(eaton_acres, .data[i,]))) %>%
  arrange(similarity_score) %>%
  slice_head(n = 5)

# Find historical fires similar to Palisades Fire
similar_to_palisades <- historical_fires %>%
  filter(year < 2018, year >= 2010, fire_size_acres > 0) %>%
  mutate(similarity_score = sapply(1:n(), function(i) calc_fire_similarity(palisades_acres, .data[i,]))) %>%
  arrange(similarity_score) %>%
  slice_head(n = 5)

# Print similar fires
cat("\nHistorical fires most similar to Eaton Fire:\n")
print(similar_to_eaton %>% 
      select(FIRE_NAME, YEAR_, month, fire_size_acres, similarity_score) %>% 
      st_drop_geometry())

cat("\nHistorical fires most similar to Palisades Fire:\n")
print(similar_to_palisades %>% 
      select(FIRE_NAME, YEAR_, month, fire_size_acres, similarity_score) %>% 
      st_drop_geometry())

# Create a map of the 2025 fires and the most similar historical fires
similar_fires_map <- tm_shape(st_transform(historical_fires, st_crs(fires_2025))) +
  tm_borders(alpha = 0.2, col = "gray") +
  tm_shape(similar_to_eaton) +
  tm_borders(col = "orange", lwd = 1.5) +
  tm_shape(similar_to_palisades) +
  tm_borders(col = "blue", lwd = 1.5) +
  tm_shape(fires_2025) +
  tm_borders(col = "red", lwd = 2) +
  tm_add_legend(type = "line", 
                col = c("red", "orange", "blue"), 
                lwd = c(2, 1.5, 1.5),
                labels = c("2025 Fires", "Similar to Eaton", "Similar to Palisades")) +
  tm_layout(title = "2025 Fires Compared to Similar Historical Fires")

tmap_save(similar_fires_map, file.path(map_dir, "similar_historical_fires.png"), width = 10, height = 8)

# Combine all similar historical fires for analysis
similar_historical_fires <- rbind(similar_to_eaton, similar_to_palisades) %>%
  distinct()

# 2. Create buffer zones around historical fires ------------------------------
cat("Creating buffer zones around historical fires...\n")

# Create 10km buffers around the historical fires
historical_buffers <- st_buffer(similar_historical_fires, dist = 10000)

# Save the combined historical fire buffers
st_write(historical_buffers, file.path(processed_dir, "historical_fire_buffers.shp"), append = FALSE)

# 3. Match with historical health and vulnerability data ----------------------
cat("Matching historical fire areas with census and health data...\n")

# Function to standardize GEOID format across datasets
standardize_geoid <- function(geoid) {
  geoid <- as.character(geoid)
  # Ensure GEOID is 11 characters (state(2) + county(3) + tract(6))
  if(nchar(geoid) < 11) {
    geoid <- str_pad(geoid, 11, pad = "0")
  }
  return(geoid)
}

# Process historical CalEnviroScreen data if available
historical_health_data <- NULL

if(!is.null(ces_3_0)) {
  # Create dataset for pre-2018 period
  ces_3_0_proc <- ces_3_0 %>%
    # Adjust column names based on your actual data
    select(Census.Tract, Asthma, Asthma.Pctl, Poverty, Population) %>%
    rename(
      geoid = Census.Tract,
      asthma_rate = Asthma,
      asthma_percentile = Asthma.Pctl
    ) %>%
    mutate(
      geoid = standardize_geoid(geoid),
      data_period = "2011-2015"
    )
  
  historical_health_data <- ces_3_0_proc
}

if(!is.null(ces_4_0)) {
  # Create dataset for 2015-2018 period
  ces_4_0_proc <- ces_4_0 %>%
    # Adjust column names based on your actual data
    select(Census.Tract, Asthma, Asthma.Pctl, Poverty, Population) %>%
    rename(
      geoid = Census.Tract,
      asthma_rate = Asthma,
      asthma_percentile = Asthma.Pctl
    ) %>%
    mutate(
      geoid = standardize_geoid(geoid),
      data_period = "2015-2018"
    )
  
  if(is.null(historical_health_data)) {
    historical_health_data <- ces_4_0_proc
  } else {
    historical_health_data <- bind_rows(historical_health_data, ces_4_0_proc)
  }
}

# Process historical SVI data if available
historical_svi_data <- NULL

if(!is.null(svi_2014)) {
  # Create dataset for pre-2018 period
  svi_2014_proc <- svi_2014 %>%
    # Adjust column names based on your actual data
    select(FIPS, RPL_THEMES, RPL_THEME1, RPL_THEME2, RPL_THEME3, RPL_THEME4) %>%
    rename(
      geoid = FIPS,
      svi_score = RPL_THEMES,
      svi_socioeconomic = RPL_THEME1,
      svi_household = RPL_THEME2,
      svi_minority = RPL_THEME3,
      svi_housing_transport = RPL_THEME4
    ) %>%
    mutate(
      geoid = standardize_geoid(geoid),
      data_period = "2014"
    )
  
  historical_svi_data <- svi_2014_proc
}

if(!is.null(svi_2018)) {
  # Create dataset for 2015-2018 period
  svi_2018_proc <- svi_2018 %>%
    # Adjust column names based on your actual data
    select(FIPS, RPL_THEMES, RPL_THEME1, RPL_THEME2, RPL_THEME3, RPL_THEME4) %>%
    rename(
      geoid = FIPS,
      svi_score = RPL_THEMES,
      svi_socioeconomic = RPL_THEME1,
      svi_household = RPL_THEME2,
      svi_minority = RPL_THEME3,
      svi_housing_transport = RPL_THEME4
    ) %>%
    mutate(
      geoid = standardize_geoid(geoid),
      data_period = "2018"
    )
  
  if(is.null(historical_svi_data)) {
    historical_svi_data <- svi_2018_proc
  } else {
    historical_svi_data <- bind_rows(historical_svi_data, svi_2018_proc)
  }
}

# 4. Get census tracts for different time periods -----------------------------
# We need to link historical fires with census geography from that time period

# Function to get census tracts for a specific year
get_census_tracts <- function(year) {
  # Attempt to load cached tracts first
  tract_path <- file.path(processed_dir, paste0("ca_tracts_", year, ".shp"))
  
  if(file.exists(tract_path)) {
    return(st_read(tract_path))
  } else {
    # Use tigris to download if not available (adjust as needed)
    cat("Downloading census tracts for", year, "...\n")
    tracts <- tigris::tracts(state = "CA", year = year)
    
    # Save for future use
    st_write(tracts, tract_path, append = FALSE)
    return(tracts)
  }
}

# Get census tracts for the periods corresponding to our historical fires
tracts_2010 <- get_census_tracts(2010)
tracts_2015 <- get_census_tracts(2015)

# 5. Identify affected census tracts for historical fires ---------------------
cat("Identifying census tracts affected by historical fires...\n")

# Identify census tracts affected by historical fires in each period
pre2015_fires <- similar_historical_fires %>% 
  filter(year < 2015) %>%
  st_transform(st_crs(tracts_2010))

post2015_fires <- similar_historical_fires %>% 
  filter(year >= 2015) %>%
  st_transform(st_crs(tracts_2015))

# Create buffer areas
pre2015_buffers <- st_buffer(pre2015_fires, dist = 10000)
post2015_buffers <- st_buffer(post2015_fires, dist = 10000)

# Identify affected tracts
if(nrow(pre2015_buffers) > 0) {
  tracts_2010 <- tracts_2010 %>%
    mutate(
      in_hist_fire_10km = as.integer(st_intersects(geometry, st_union(pre2015_buffers), sparse = FALSE)[,1]),
      period = "2010-2015"
    )
}

if(nrow(post2015_buffers) > 0) {
  tracts_2015 <- tracts_2015 %>%
    mutate(
      in_hist_fire_10km = as.integer(st_intersects(geometry, st_union(post2015_buffers), sparse = FALSE)[,1]),
      period = "2015-2018"
    )
}

# 6. Join historical health and vulnerability data to tracts ------------------
cat("Joining historical data to affected census tracts...\n")

# Create joined datasets for each time period
historical_tracts_pre2015 <- NULL
historical_tracts_post2015 <- NULL

if(!is.null(historical_health_data) && !is.null(historical_svi_data)) {
  # Join data for pre-2015 period
  if(nrow(pre2015_buffers) > 0) {
    historical_tracts_pre2015 <- tracts_2010 %>%
      mutate(geoid = standardize_geoid(GEOID)) %>%
      left_join(historical_health_data %>% filter(data_period == "2011-2015"), by = "geoid") %>%
      left_join(historical_svi_data %>% filter(data_period == "2014"), by = "geoid")
  }
  
  # Join data for post-2015 period
  if(nrow(post2015_buffers) > 0) {
    historical_tracts_post2015 <- tracts_2015 %>%
      mutate(geoid = standardize_geoid(GEOID)) %>%
      left_join(historical_health_data %>% filter(data_period == "2015-2018"), by = "geoid") %>%
      left_join(historical_svi_data %>% filter(data_period == "2018"), by = "geoid")
  }
  
  # Combine historical data
  historical_tracts_combined <- bind_rows(historical_tracts_pre2015, historical_tracts_post2015)
  
  # Save the combined historical tracts data
  st_write(historical_tracts_combined, file.path(processed_dir, "historical_tracts_data.shp"), append = FALSE)
} else {
  cat("Warning: Unable to join historical data due to missing health or SVI data.\n")
}

# 7. Compare current and historical fire impacts ------------------------------
cat("Comparing current and historical fire impacts...\n")

# Function to analyze fire impacts on asthma rates
analyze_fire_impact <- function(data, period_label) {
  # Skip if data is NULL
  if(is.null(data)) return(NULL)
  
  # Calculate summary statistics by fire exposure
  impact_summary <- data %>%
    st_drop_geometry() %>%
    group_by(in_hist_fire_10km) %>%
    summarize(
      n_tracts = n(),
      mean_asthma = mean(asthma_rate, na.rm = TRUE),
      median_asthma = median(asthma_rate, na.rm = TRUE),
      sd_asthma = sd(asthma_rate, na.rm = TRUE),
      mean_svi = mean(svi_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      time_period = period_label,
      fire_exposure = ifelse(in_hist_fire_10km == 1, "Fire-Affected", "Not Fire-Affected")
    )
  
  return(impact_summary)
}

# Analyze historical periods
impact_pre2015 <- analyze_fire_impact(historical_tracts_pre2015, "2010-2015")
impact_post2015 <- analyze_fire_impact(historical_tracts_post2015, "2015-2018")

# Analyze current (2025) data for comparison
impact_current <- current_data %>%
  st_drop_geometry() %>%
  group_by(in_fire_10km) %>%
  summarize(
    n_tracts = n(),
    mean_asthma = mean(asthma_rate, na.rm = TRUE),
    median_asthma = median(asthma_rate, na.rm = TRUE),
    sd_asthma = sd(asthma_rate, na.rm = TRUE),
    mean_svi = mean(svi_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    time_period = "2025",
    in_hist_fire_10km = in_fire_10km,
    fire_exposure = ifelse(in_fire_10km == 1, "Fire-Affected", "Not Fire-Affected")
  ) %>%
  select(-in_fire_10km)

# Combine all impact summaries
all_impacts <- bind_rows(impact_pre2015, impact_post2015, impact_current)

# Save impact summary
if(nrow(all_impacts) > 0) {
  write_csv(all_impacts, file.path(hist_dir, "fire_impact_comparison.csv"))
  
  # Create comparison plots
  p1 <- ggplot(all_impacts, aes(x = time_period, y = mean_asthma, fill = fire_exposure)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(
      title = "Asthma ED Rates by Fire Exposure Over Time",
      x = "Time Period",
      y = "Mean Asthma ED Rate",
      fill = "Fire Exposure"
    ) +
    theme_minimal()
  
  p2 <- ggplot(all_impacts, aes(x = time_period, y = mean_svi, fill = fire_exposure)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(
      title = "Social Vulnerability by Fire Exposure Over Time",
      x = "Time Period",
      y = "Mean SVI Score",
      fill = "Fire Exposure"
    ) +
    theme_minimal()
  
  ggsave(file.path(fig_dir, "historical_asthma_comparison.png"), p1, width = 10, height = 6)
  ggsave(file.path(fig_dir, "historical_svi_comparison.png"), p2, width = 10, height = 6)
}

# 8. Analyze changing vulnerability patterns over time ------------------------
cat("Analyzing changes in vulnerability patterns over time...\n")

# Function to analyze vulnerability patterns in fire-affected areas
analyze_vulnerability_patterns <- function(data, period_label) {
  # Skip if data is NULL
  if(is.null(data)) return(NULL)
  
  # Filter for fire-affected tracts
  affected_tracts <- data %>%
    filter(in_hist_fire_10km == 1)
  
  # Create SVI quartiles
  affected_tracts$svi_quartile <- cut(affected_tracts$svi_score, 
                                   breaks = quantile(affected_tracts$svi_score, 
                                                     probs = seq(0, 1, 0.25), 
                                                     na.rm = TRUE),
                                   labels = c("Low", "Medium-Low", "Medium-High", "High"),
                                   include.lowest = TRUE)
  
  # Calculate asthma rates by SVI quartile
  vulnerability_pattern <- affected_tracts %>%
    st_drop_geometry() %>%
    group_by(svi_quartile) %>%
    summarize(
      n_tracts = n(),
      mean_asthma = mean(asthma_rate, na.rm = TRUE),
      sd_asthma = sd(asthma_rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(time_period = period_label)
  
  return(vulnerability_pattern)
}

# Analyze vulnerability patterns in each time period
vuln_pre2015 <- analyze_vulnerability_patterns(historical_tracts_pre2015, "2010-2015")
vuln_post2015 <- analyze_vulnerability_patterns(historical_tracts_post2015, "2015-2018")

# Analyze current (2025) data for comparison
current_data$svi_quartile <- cut(current_data$svi_score, 
                               breaks = quantile(current_data$svi_score, 
                                               probs = seq(0, 1, 0.25), 
                                               na.rm = TRUE),
                               labels = c("Low", "Medium-Low", "Medium-High", "High"),
                               include.lowest = TRUE)

vuln_current <- current_data %>%
  filter(in_fire_10km == 1) %>%
  st_drop_geometry() %>%
  group_by(svi_quartile) %>%
  summarize(
    n_tracts = n(),
    mean_asthma = mean(asthma_rate, na.rm = TRUE),
    sd_asthma = sd(asthma_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(time_period = "2025")

# Combine all vulnerability patterns
all_vuln_patterns <- bind_rows(vuln_pre2015, vuln_post2015, vuln_current)

# Save vulnerability pattern summary
if(nrow(all_vuln_patterns) > 0) {
  write_csv(all_vuln_patterns, file.path(hist_dir, "vulnerability_patterns_over_time.csv"))
  
  # Create comparison plot
  p3 <- ggplot(all_vuln_patterns, aes(x = svi_quartile, y = mean_asthma, fill = time_period)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(
      title = "Asthma ED Rates by SVI Quartile Over Time",
      subtitle = "Fire-Affected Areas Only",
      x = "Social Vulnerability Index Quartile",
      y = "Mean Asthma ED Rate",
      fill = "Time Period"
    ) +
    theme_minimal()
  
  ggsave(file.path(fig_dir, "vulnerability_patterns_over_time.png"), p3, width = 10, height = 6)
}

# 9. Implement matching analysis to compare fire effects over time ------------
cat("Implementing matching analysis to compare fire effects over time...\n")

# Function to perform propensity score matching
perform_matching_analysis <- function(data, time_label) {
  # Skip if data is NULL
  if(is.null(data)) return(NULL)
  
  # Prepare data for matching
  match_data <- data %>%
    st_drop_geometry() %>%
    filter(!is.na(asthma_rate) & !is.na(svi_score)) %>%
    select(in_hist_fire_10km, asthma_rate, svi_score, svi_socioeconomic, 
          svi_minority, Poverty, Population)
  
  # Estimate propensity scores
  ps_model <- glm(in_hist_fire_10km ~ svi_score + svi_socioeconomic + 
                 svi_minority + Poverty, 
                 data = match_data, 
                 family = binomial())
  
  # Perform matching
  match_result <- matchit(in_hist_fire_10km ~ svi_score + svi_socioeconomic + 
                         svi_minority + Poverty, 
                         data = match_data, 
                         method = "nearest",
                         ratio = 1)
  
  # Extract matched data
  matched_data <- match.data(match_result)
  
  # Analyze fire effect in matched sample
  t_test <- t.test(asthma_rate ~ in_hist_fire_10km, data = matched_data)
  
  # Calculate standardized mean difference after matching
  smd <- mean(matched_data$asthma_rate[matched_data$in_hist_fire_10km == 1], na.rm = TRUE) - 
         mean(matched_data$asthma_rate[matched_data$in_hist_fire_10km == 0], na.rm = TRUE)
  smd_std <- smd / sd(matched_data$asthma_rate, na.rm = TRUE)
  
  # Create summary of results
  match_summary <- data.frame(
    time_period = time_label,
    mean_diff = t_test$estimate[1] - t_test$estimate[2],
    std_mean_diff = smd_std,
    p_value = t_test$p.value,
    fire_affected_n = sum(matched_data$in_hist_fire_10km == 1),
    control_n = sum(matched_data$in_hist_fire_10km == 0)
  )
  
  return(match_summary)
}

# Perform matching analysis for each time period
match_pre2015 <- perform_matching_analysis(historical_tracts_pre2015, "2010-2015")
match_post2015 <- perform_matching_analysis(historical_tracts_post2015, "2015-2018")

# Adjust column names for current data to match
current_match_data <- current_data %>%
  rename(in_hist_fire_10km = in_fire_10km)

match_current <- perform_matching_analysis(current_match_data, "2025")

# Combine all matching results
all_match_results <- bind_rows(match_pre2015, match_post2015, match_current)

# Save matching results
if(nrow(all_match_results) > 0) {
  write_csv(all_match_results, file.path(hist_dir, "matching_analysis_results.csv"))
  
  # Create comparison plot
  p4 <- ggplot(all_match_results, aes(x = time_period, y = std_mean_diff)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(aes(ymin = std_mean_diff - 1.96*std_mean_diff/sqrt(fire_affected_n), 
                     ymax = std_mean_diff + 1.96*std_mean_diff/sqrt(fire_affected_n)),
                 width = 0.2) +
    labs(
      title = "Standardized Mean Difference in Asthma Rates",
      subtitle = "After Propensity Score Matching: Fire-Affected vs. Control Areas",
      x = "Time Period",
      y = "Standardized Mean Difference"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")
  
  ggsave(file.path(fig_dir, "matching_analysis_comparison.png"), p4, width = 10, height = 6)
}

# 10. Examine changes in wildfire-vulnerability relationship over time --------
cat("Examining changes in the wildfire-vulnerability relationship over time...\n")

# Function to test for interaction between fire exposure and SVI over time
test_interaction_change <- function(data, time_label) {
  # Skip if data is NULL
  if(is.null(data)) return(NULL)
  
  # Run model with interaction term
  model <- lm(asthma_`rate ~ in_hist_fire_10km * svi_score, data = data)`