# Historical Fire Comparison Analysis
# This script demonstrates how to compare the January 2025 LA fires
# with historical fires to examine patterns across different social vulnerability profiles

# Load required libraries
library(sf)
library(dplyr)
library(ggplot2)
library(tmap)
library(spdep)
library(lubridate)
library(MatchIt)  # For propensity score matching
library(purrr)    # For functional programming tools

# Set working directory (adjust as needed)
# setwd("/Users/marcrisney/Projects/jhu/mas/Spring2025/spatial-analysis-ia")

# Data Loading ----------------------------------------------------------------

# 1. Load LA County boundary
la_county <- st_read("data/boundaries/la_county_boundary.shp")

# 2. Load census tracts with Social Vulnerability Index (current)
tracts_current <- st_read("data/boundaries/la_census_tracts.shp")

# 3. Load January 2025 fire perimeters
current_fires <- st_read("data/fires/la_fires_jan2025.shp")

# 4. Load historical fire perimeters from CAL FIRE FRAP
historical_fires <- st_read("data/fires/historical_fires.shp")

# 5. Load historical CalEnviroScreen data (versions 3.0 and 4.0)
calenviroscreen_3 <- read.csv("data/health/calenviroscreen_3.0.csv")
calenviroscreen_4 <- read.csv("data/health/calenviroscreen_4.0.csv")

# 6. Load historical census and SVI data
census_tracts_2014 <- st_read("data/boundaries/ca_census_tracts_2014.shp")
census_tracts_2018 <- st_read("data/boundaries/ca_census_tracts_2018.shp")
svi_2014 <- read.csv("data/social/svi_2014.csv")
svi_2018 <- read.csv("data/social/svi_2018.csv")

# Data Preparation ------------------------------------------------------------

# Ensure all spatial data is in the same CRS
la_county <- st_transform(la_county, 2229)
tracts_current <- st_transform(tracts_current, 2229)
current_fires <- st_transform(current_fires, 2229)
historical_fires <- st_transform(historical_fires, 2229)
census_tracts_2014 <- st_transform(census_tracts_2014, 2229)
census_tracts_2018 <- st_transform(census_tracts_2018, 2229)

# Add year and extract month information for historical fires
historical_fires$year <- year(as.Date(historical_fires$ALARM_DATE, format = "%Y-%m-%d"))
historical_fires$month <- month(as.Date(historical_fires$ALARM_DATE, format = "%Y-%m-%d"))

# Filter historical fires for the relevant time periods
fires_2011_2014 <- historical_fires %>% 
  filter(year >= 2011 & year <= 2014)

fires_2015_2018 <- historical_fires %>% 
  filter(year >= 2015 & year <= 2018)

# Join SVI data to historical census tracts
census_tracts_2014 <- census_tracts_2014 %>%
  left_join(svi_2014, by = c("GEOID" = "FIPS"))

census_tracts_2018 <- census_tracts_2018 %>%
  left_join(svi_2018, by = c("GEOID" = "FIPS"))

# Join CalEnviroScreen data to historical census tracts
census_tracts_2014 <- census_tracts_2014 %>%
  left_join(calenviroscreen_3, by = c("GEOID" = "Census.Tract"))

census_tracts_2018 <- census_tracts_2018 %>%
  left_join(calenviroscreen_4, by = c("GEOID" = "Census.Tract"))

# Identify the Current Fire Characteristics -------------------------------------

# Calculate attributes of the January 2025 fires
eaton_fire <- current_fires %>% filter(FIRE_NAME == "Eaton")
palisades_fire <- current_fires %>% filter(FIRE_NAME == "Palisades")

eaton_fire_attributes <- tibble(
  name = "Eaton",
  acres = as.numeric(st_area(eaton_fire) * 0.000247105), # convert to acres
  perimeter_length = as.numeric(st_length(st_cast(st_boundary(eaton_fire), "MULTILINESTRING")) * 0.000621371), # convert to miles
  month = 1, # January
  land_use = "Urban/Wildland Interface"
)

palisades_fire_attributes <- tibble(
  name = "Palisades",
  acres = as.numeric(st_area(palisades_fire) * 0.000247105),
  perimeter_length = as.numeric(st_length(st_cast(st_boundary(palisades_fire), "MULTILINESTRING")) * 0.000621371),
  month = 1, # January
  land_use = "Urban/Wildland Interface"
)

current_fire_attributes <- bind_rows(eaton_fire_attributes, palisades_fire_attributes)
print(current_fire_attributes)

# Find Similar Historical Fires ------------------------------------------------

# Function to calculate similarity score between fires
calculate_fire_similarity <- function(current_fire, historical_fire) {
  # Calculate normalized difference in acreage (0-1 scale)
  acres_diff <- abs(current_fire$acres - historical_fire$GIS_ACRES) / max(current_fire$acres, historical_fire$GIS_ACRES)
  
  # Calculate similarity based on month (0-1 scale)
  month_diff <- abs(current_fire$month - historical_fire$month) / 11
  
  # Combine scores (lower is more similar)
  similarity_score <- acres_diff * 0.7 + month_diff * 0.3
  
  return(similarity_score)
}

# Find similar fires to Eaton Fire
eaton_similar_fires <- fires_2011_2014 %>%
  rowwise() %>%
  mutate(similarity_score = calculate_fire_similarity(eaton_fire_attributes, .)) %>%
  ungroup() %>%
  arrange(similarity_score) %>%
  slice_head(n = 5)

# Find similar fires to Palisades Fire
palisades_similar_fires <- fires_2015_2018 %>%
  rowwise() %>%
  mutate(similarity_score = calculate_fire_similarity(palisades_fire_attributes, .)) %>%
  ungroup() %>%
  arrange(similarity_score) %>%
  slice_head(n = 5)

# Print the most similar historical fires
print("Most similar historical fires to Eaton Fire:")
print(eaton_similar_fires %>% select(FIRE_NAME, YEAR_, GIS_ACRES, similarity_score))

print("Most similar historical fires to Palisades Fire:")
print(palisades_similar_fires %>% select(FIRE_NAME, YEAR_, GIS_ACRES, similarity_score))

# Identify Census Tracts Affected by Historical Fires --------------------------

# Create buffers around historical fires (similar to current analysis)
eaton_similar_buffers <- eaton_similar_fires %>%
  st_buffer(dist = 10000) # 10km buffer

palisades_similar_buffers <- palisades_similar_fires %>%
  st_buffer(dist = 10000) # 10km buffer

# Identify census tracts within the buffer zones of similar historical fires
tracts_2014_affected <- census_tracts_2014 %>%
  mutate(affected_by_similar_fire = as.integer(
    st_intersects(., st_union(eaton_similar_buffers), sparse = FALSE)
  ))

tracts_2018_affected <- census_tracts_2018 %>%
  mutate(affected_by_similar_fire = as.integer(
    st_intersects(., st_union(palisades_similar_buffers), sparse = FALSE)
  ))

# Compare Social Vulnerability Profiles ----------------------------------------

# Compare SVI distributions in current fire-affected areas vs historical fire-affected areas
current_affected_tracts <- tracts_current %>%
  mutate(affected_by_fire = as.integer(
    st_intersects(., st_buffer(current_fires, 10000), sparse = FALSE)
  ))

# Create summary statistics for SVI in fire-affected areas
svi_current_summary <- current_affected_tracts %>%
  st_drop_geometry() %>%
  group_by(affected_by_fire) %>%
  summarize(
    mean_svi = mean(svi_score, na.rm = TRUE),
    median_svi = median(svi_score, na.rm = TRUE),
    sd_svi = sd(svi_score, na.rm = TRUE),
    q25_svi = quantile(svi_score, 0.25, na.rm = TRUE),
    q75_svi = quantile(svi_score, 0.75, na.rm = TRUE)
  )

svi_2014_summary <- tracts_2014_affected %>%
  st_drop_geometry() %>%
  group_by(affected_by_similar_fire) %>%
  summarize(
    mean_svi = mean(RPL_THEMES, na.rm = TRUE),
    median_svi = median(RPL_THEMES, na.rm = TRUE),
    sd_svi = sd(RPL_THEMES, na.rm = TRUE),
    q25_svi = quantile(RPL_THEMES, 0.25, na.rm = TRUE),
    q75_svi = quantile(RPL_THEMES, 0.75, na.rm = TRUE)
  )

svi_2018_summary <- tracts_2018_affected %>%
  st_drop_geometry() %>%
  group_by(affected_by_similar_fire) %>%
  summarize(
    mean_svi = mean(RPL_THEMES, na.rm = TRUE),
    median_svi = median(RPL_THEMES, na.rm = TRUE),
    sd_svi = sd(RPL_THEMES, na.rm = TRUE),
    q25_svi = quantile(RPL_THEMES, 0.25, na.rm = TRUE),
    q75_svi = quantile(RPL_THEMES, 0.75, na.rm = TRUE)
  )

# Print comparison summaries
print("Current fires - SVI comparison:")
print(svi_current_summary)
print("Historical fires (2011-2014) - SVI comparison:")
print(svi_2014_summary)
print("Historical fires (2015-2018) - SVI comparison:")
print(svi_2018_summary)

# Plot SVI distribution comparison
current_svi_data <- current_affected_tracts %>%
  st_drop_geometry() %>%
  select(affected_by_fire, svi_score) %>%
  mutate(period = "Current (2025)")

hist_2014_svi_data <- tracts_2014_affected %>%
  st_drop_geometry() %>%
  select(affected_by_similar_fire, RPL_THEMES) %>%
  rename(affected_by_fire = affected_by_similar_fire, svi_score = RPL_THEMES) %>%
  mutate(period = "Historical (2011-2014)")

hist_2018_svi_data <- tracts_2018_affected %>%
  st_drop_geometry() %>%
  select(affected_by_similar_fire, RPL_THEMES) %>%
  rename(affected_by_fire = affected_by_similar_fire, svi_score = RPL_THEMES) %>%
  mutate(period = "Historical (2015-2018)")

combined_svi_data <- bind_rows(current_svi_data, hist_2014_svi_data, hist_2018_svi_data) %>%
  filter(!is.na(svi_score))

# Create violin plot of SVI distribution
ggplot(combined_svi_data, aes(x = factor(affected_by_fire), y = svi_score, fill = period)) +
  geom_violin(position = position_dodge(width = 0.7), alpha = 0.7) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.7), alpha = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Social Vulnerability Index Distribution",
    subtitle = "Comparison between current and historical fire-affected areas",
    x = "Affected by Fire",
    y = "Social Vulnerability Index Score",
    fill = "Time Period"
  ) +
  scale_x_discrete(labels = c("0" = "Not Affected", "1" = "Fire-Affected")) +
  theme_minimal()

# Compare Health Outcomes in Similar Historical Fires --------------------------

# Extract asthma ED rates from CalEnviroScreen for historical periods
asthma_2014 <- tracts_2014_affected %>%
  st_drop_geometry() %>%
  select(GEOID, affected_by_similar_fire, Asthma) %>%
  rename(asthma_rate = Asthma) %>%
  mutate(period = "2011-2014")

asthma_2018 <- tracts_2018_affected %>%
  st_drop_geometry() %>%
  select(GEOID, affected_by_similar_fire, Asthma) %>%
  rename(asthma_rate = Asthma) %>%
  mutate(period = "2015-2018")

asthma_current <- current_affected_tracts %>%
  st_drop_geometry() %>%
  select(GEOID, affected_by_fire, asthma_rate) %>%
  mutate(period = "2025") %>%
  rename(affected_by_similar_fire = affected_by_fire)

# Combine datasets
combined_asthma <- bind_rows(asthma_2014, asthma_2018, asthma_current)

# Calculate summary statistics for asthma rates
asthma_summary <- combined_asthma %>%
  group_by(period, affected_by_similar_fire) %>%
  summarize(
    mean_asthma = mean(asthma_rate, na.rm = TRUE),
    median_asthma = median(asthma_rate, na.rm = TRUE),
    sd_asthma = sd(asthma_rate, na.rm = TRUE),
    count = n()
  )

# Create bar plot comparing asthma rates
ggplot(asthma_summary, aes(x = period, y = mean_asthma, fill = factor(affected_by_similar_fire))) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = mean_asthma - sd_asthma/sqrt(count), 
                   ymax = mean_asthma + sd_asthma/sqrt(count)),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "firebrick"),
                   labels = c("0" = "Not Fire-Affected", "1" = "Fire-Affected"),
                   name = "Fire Exposure") +
  labs(
    title = "Asthma ED Rates Comparison",
    subtitle = "Current vs. Historical Fire-Affected Areas",
    x = "Time Period",
    y = "Asthma ED Rate"
  ) +
  theme_minimal()

# Matching Analysis to Control for SVI Differences -----------------------------

# Create combined dataset for matching
matching_data_2014 <- tracts_2014_affected %>%
  st_drop_geometry() %>%
  select(
    GEOID, 
    affected_by_similar_fire,
    RPL_THEMES, # Overall SVI
    EP_POV, # % Poverty
    EP_UNEMP, # % Unemployment
    EP_MINRTY, # % Minority
    Asthma # Asthma ED rate
  ) %>%
  rename(
    fire_affected = affected_by_similar_fire,
    svi_score = RPL_THEMES,
    poverty_pct = EP_POV,
    unemp_pct = EP_UNEMP,
    minority_pct = EP_MINRTY,
    asthma_rate = Asthma
  ) %>%
  filter(!is.na(asthma_rate) & !is.na(svi_score))

# Perform propensity score matching
ps_model <- glm(fire_affected ~ svi_score + poverty_pct + minority_pct, 
               data = matching_data_2014, 
               family = binomial())

# Match fire-affected and non-affected tracts
match_result <- matchit(
  fire_affected ~ svi_score + poverty_pct + minority_pct,
  data = matching_data_2014,
  method = "nearest",
  ratio = 1
)

# Extract matched data
matched_data <- match.data(match_result)

# Compare asthma rates in matched data
t.test(asthma_rate ~ fire_affected, data = matched_data)

# Plot results
matched_summary <- matched_data %>%
  group_by(fire_affected) %>%