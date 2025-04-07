# Spatial Analysis of LA Wildfires and Health Impacts

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains data, analysis code, and documentation for a research project examining the relationship between Los Angeles wildfire events (focusing on the January 2025 Eaton and Palisades fires) and respiratory health outcomes, with special attention to how social vulnerability factors modify this relationship.

The study integrates multiple spatial datasets including fire perimeters, CalEnviroScreen asthma emergency department (ED) indicators, ATSDR social vulnerability indices, air quality measurements, and building damage assessments to model the spatial association between wildfire exposure and health impacts.

## Research Questions

1. Are areas with higher social vulnerability indices experiencing greater increases in asthma-related ED visits post-fire compared to areas with lower vulnerability indices, even when controlling for similar proximity to fire perimeters?

2. How does the relationship between fire proximity and respiratory health effects vary across different levels of social vulnerability factors including socioeconomic status, housing type/quality, and race/ethnicity?

3. What spatial patterns emerge when comparing current fire impacts with historical California fires having similar perimeter characteristics?

## Repository Structure

```
├── data/
│   ├── boundaries/        # Geographic boundary files for analysis units
│   ├── fires/             # Fire perimeter data (current and historical)
│   ├── health/            # Health outcome data including asthma indicators
│   ├── social/            # Social vulnerability indices and demographic data
│   ├── air_quality/       # Air quality measurements during fire events
│   └── processed/         # Processed/derived datasets
├── code/
│   ├── 01_data_prep.R     # Data preparation and cleaning
│   ├── 02_eda.R           # Exploratory data analysis
│   ├── 03_spatial_analysis.R  # Main spatial statistical analysis
│   ├── 04_historical_comparison.R  # Historical fire comparison analysis
│   └── 05_visualization.R  # Visualization and mapping functions
├── results/
│   ├── maps/              # Generated map outputs
│   ├── figures/           # Other data visualizations
│   └── tables/            # Statistical output tables
├── docs/
│   ├── proposal.md        # Research proposal
│   └── methods.md         # Detailed methodology
└── README.md              # This readme file
```

## Data Sources

The project utilizes data from the following sources:

1. **Fire Perimeter Data**:
   - LA County January 2025 fire perimeters ([LA County GIS Data Portal](https://data.lacounty.gov/datasets/4da742952a1a4b958218f57a53d7393e/explore))
   - Historical fire perimeters from [CAL FIRE's FRAP](https://frap.fire.ca.gov/mapping/gis-data/)

2. **Health Outcome Data**:
   - [CalEnviroScreen Asthma ED indicators](https://oehha.ca.gov/calenviroscreen/indicator/asthma)
   - CDC emergency department syndromic surveillance data (January 2025)

3. **Social Vulnerability Data**:
   - [ATSDR Social Vulnerability Index](https://www.atsdr.cdc.gov/place-health)
   - CalEnviroScreen 4.0 socioeconomic indicators

4. **Environmental Exposure Data**:
   - Air quality data from PurpleAir sensors
   - Building damage assessments for the Eaton and Palisades fires
   - Post-fire soil and ash testing data (Caltech Science Exchange)

## Methodology

The analysis employs several spatial statistical methods:

1. **Buffer and Distance Analysis**: Creating multiple buffer zones around fire perimeters to establish proximity-based exposure measures

2. **Geostatistical Methods**: Spatial interpolation of air quality measurements using kriging and inverse distance weighting

3. **Spatial Autocorrelation Analysis**: Moran's I and LISA statistics to identify spatial clustering of health outcomes

4. **Spatial Regression Models**: Accounting for spatial dependencies when modeling the relationship between fire exposure, social vulnerability, and health outcomes

5. **Historical Comparative Analysis**: Matching census tracts affected by current fires with areas affected by similar historical fires to compare outcomes

## Getting Started

### Prerequisites

This project requires R (version 4.0.0 or higher) with the following packages:
- sf
- dplyr
- ggplot2
- tmap
- spdep
- geoR
- raster
- spatstat
- spatialreg

### Installation

Clone this repository:
```bash
git clone https://github.com/mrisney/spatial-analysis-ia.git
cd spatial-analysis-ia
```

### Running the Analysis

1. First, prepare the data:
```r
source("code/01_data_prep.R")
```

2. Run the exploratory data analysis:
```r
source("code/02_eda.R")
```

3. Perform the main spatial analysis:
```r
source("code/03_spatial_analysis.R")
```

## Results

[This section will be populated as analyses are completed]

## Contributing

Contributions to this project are welcome. Please feel free to open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- CalEnviroScreen for providing public health and environmental data
- ATSDR for developing the Social Vulnerability Index
- LA County GIS Data Portal for providing fire perimeter and building damage data
- Dr. Francois Tissot (Caltech) for post-fire soil and ash testing data

## Contact

For questions about this repository, please contact:
- [Marc Risney](https://github.com/mrisney)