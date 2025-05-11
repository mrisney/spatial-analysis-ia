# Fire–Asthma Prediction Master Plan

This repository implements an end-to-end workflow to predict census-tract–level asthma ED visits following LA-area wildfires, combining spatial interpolation (Kriging), temporal forecasting (Prophet), and machine-learning (XGBoost).

---

## Repository Structure

```
├── data/
│ ├── fires/ # Raw fire perimeter files (FRAP, Eaton, Palisades)
│ ├── census_tracts/ # 2024 tract shapefile
│ ├── aqi/ # Raw AQI/PM₂.₅ CSVs
│ └── svi/ # Social Vulnerability Index data
├── results/
│ ├── svi_by_tract.csv
│ ├── fire_tracts.csv
│ ├── ces_asthma_baseline.csv
│ ├── aqi_filtered.csv
│ ├── variogram_parameters.csv
│ ├── kriged_aqi_by_tract.csv
│ ├── aqi_forecast_2025.csv
│ ├── ed_trend_feature.csv
│ ├── model_features.csv
│ ├── xgb_performance.csv
│ ├── shap_summary.csv
│ └── … PNG/HTML outputs
├── models/
│ └── xgb_best_model.json
├── code/
│ ├── 01_load_fire_and_tracts.py
│ ├── 02_list_layers.py
│ ├── 03_fetch_svi.py
│ ├── 04_extract_fire_affected_tracts.py
│ ├── 05_extract_ces_baseline_asthma.py
│ ├── 06_prepare_aqi_timeseries.py
│ ├── 07_variogram_and_kriging.R
│ ├── 08_plot_semivariogram.R
│ ├── 09_prophet_forecast_aqi.py
│ ├── 10_prophet_trend_ed.R
│ ├── 11_build_feature_matrix.py
│ ├── 12_train_xgboost.py
│ ├── 13_evaluate_and_shap.py
│ ├── 14_create_maps_and_charts.py
│ └── 15_render_report.Rmd
└── README.md # This document
```

---

## Phase 1: Data Ingestion & Cataloguing

1. **`01_load_fire_and_tracts.py`** (Python)  
   - **Purpose:** Load all fire perimeters and the 2024 census-tract shapefile into GeoPandas.  
   - **Inputs:**  
     - FRAP GDB layers + Eaton/Palisades shapefiles  
     - `data/census_tracts/tl_2024_06_tract.shp`  
   - **Outputs:**  
     - Serialized GeoDataFrames for fires and tracts.

2. **`02_list_layers.py`** (Python)  
   - **Purpose:** Inspect layers in any `.gdb` (e.g. CalEnviroScreen, SVI).  
   - **Inputs:** Path to GDB.  
   - **Outputs:** Console printout of available layers.

3. **`03_fetch_svi.py`** (Python)  
   - **Purpose:** Pull tract-level Social Vulnerability Index from CDC/Surgo API or local CSV.  
   - **Inputs:** API credentials or raw SVI file.  
   - **Outputs:** `results/svi_by_tract.csv`.

---

## Phase 2: Spatial Preprocessing

4. **`04_extract_fire_affected_tracts.py`** (Python)  
   - **Purpose:** Spatial‐join fire perimeters to tracts; classify direct vs. smoke‐only exposure (with buffer).  
   - **Inputs:** Outputs of **01**.  
   - **Outputs:**  
     - `results/fire_tracts.csv` (GEOID, fire name, exposure type)  
     - Per‐fire shapefiles of affected tracts.

5. **`05_extract_ces_baseline_asthma.py`** (Python)  
   - **Purpose:** Extract baseline annual asthma ED rates from CalEnviroScreen GDB.  
   - **Inputs:** `data/calenviroscreen40gdb_F_2021.gdb`.  
   - **Outputs:** `results/ces_asthma_baseline.csv` (GEOID, baseline rates).

---

## Phase 3: Air-Quality Filtering & Interpolation

6. **`06_prepare_aqi_timeseries.py`** (Python)  
   - **Purpose:** Load AQI CSV, parse timestamps, attach monitor coords; filter to monitors within or near each fire’s buffer.  
   - **Inputs:**  
     - `data/aqi/annual_conc_by_monitor_2016_California_sample.csv`  
     - `results/fire_tracts.csv`  
   - **Outputs:** `results/aqi_filtered.csv` (monitor, date, concentration, fire, tract).

7. **`07_variogram_and_kriging.R`** (R)  
   - **Purpose:**  
     1. Compute empirical variogram of AQI residuals.  
     2. Fit theoretical model (nugget, sill, range).  
     3. Ordinary kriging to interpolate AQI at tract centroids (block kriging).  
   - **Inputs:** `results/aqi_filtered.csv`, tract centroids.  
   - **Outputs:**  
     - `results/variogram_parameters.csv`  
     - `results/kriged_aqi_by_tract.csv`.

8. **`08_plot_semivariogram.R`** (R)  
   - **Purpose:** Plot empirical vs. fitted variogram for diagnostics.  
   - **Inputs:** `results/variogram_parameters.csv`.  
   - **Outputs:** `results/semivariogram.png`.

---

## Phase 4: Temporal Modeling (Optional)

9. **`09_prophet_forecast_aqi.py`** (Python)  
   - **Purpose:** Fit Prophet to daily kriged‐AQI per tract; forecast into 2025.  
   - **Inputs:** `results/kriged_aqi_by_tract.csv`.  
   - **Outputs:** `results/aqi_forecast_2025.csv` (predictions + intervals).

10. **`10_prophet_trend_ed.R`** (R)  
    - **Purpose:** Fit Prophet to annual ED rates per tract; extract “expected trend” feature.  
    - **Inputs:** `results/ces_asthma_baseline.csv` (multi-year ED data).  
    - **Outputs:** `results/ed_trend_feature.csv`.

---

## Phase 5: Feature Engineering & Modeling

11. **`11_build_feature_matrix.py`** (Python)  
    - **Purpose:** Merge features for each (fire, tract):  
      - Baseline asthma, SVI subindices  
      - Kriged AQI summaries (mean, exceedance days)  
      - Prophet AQI forecasts, ED trend  
      - Exposure labels  
    - **Inputs:** Outputs of **03**, **05**, **07**, **09**, **10**.  
    - **Outputs:** `results/model_features.csv`.

12. **`12_train_xgboost.py`** (Python)  
    - **Purpose:**  
      - Load features, define train/test splits (e.g. leave-one-fire-out).  
      - Hyperparameter tuning & training with XGBoost.  
      - Save best model & metrics.  
    - **Inputs:** `results/model_features.csv`.  
    - **Outputs:**  
      - `models/xgb_best_model.json`  
      - `results/xgb_performance.csv`.

13. **`13_evaluate_and_shap.py`** (Python)  
    - **Purpose:** Generate SHAP values, feature-importance plots, and spatial residual maps.  
    - **Inputs:** Best model + test dataset.  
    - **Outputs:**  
      - `results/shap_summary.csv`  
      - `results/shap_summary.png`  
      - `results/residual_map.png`.

---

## Phase 6: Visualization & Reporting

14. **`14_create_maps_and_charts.py`** (Python)  
    - **Purpose:**  
      - Produce final choropleths (predicted vs. observed ED, kriged AQI, exposure overlays).  
      - Integrate Mapbox basemaps via existing basemap scripts.  
    - **Inputs:** Results from previous phases + basemap tiles.  
    - **Outputs:** Presentation-ready PNGs.

15. **`15_render_report.Rmd`** (R Markdown)  
    - **Purpose:** Knit all tables, figures, and maps into a cohesive HTML/PDF report.  
    - **Inputs:** All CSVs & images in `results/`.  
    - **Outputs:** `report/fire_asthma_prediction_report.html`

---

> **Next steps:**  
> 1. Define CLI arguments for each script (e.g. `--input`, `--output`) so you can orchestrate the workflow via Makefile or Snakemake.  
> 2. Slot your existing code into these scripts and validate outputs incrementally.  
> 3. Iterate on feature engineering and model tuning once the data pipeline is solid.  
> 4. Use the generated plots and report for stakeholder presentations.  

---