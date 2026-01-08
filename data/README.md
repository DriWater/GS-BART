This directory contains precipitation and simulation data used in the paper.

## King County House Sales

Home sales price data for **20,149** properties in King County, Seattle. The response variable is the log-transformed house sale price. Spatial features include latitude and longitude, along with **15** additional unstructured housing-related covariates.  

**Source:** https://geodacenter.github.io/data-and-lab/KingCounty-HouseSales2015/

---

## NYC Education

New York City census data with **1,690** observations. The response variable is income per capita. Geographic location is treated as a spatial feature, while educational attainment and four additional socioeconomic variables are used as unstructured covariates.  

**Source:** https://geodacenter.github.io/data-and-lab/NYC-Census-2000

---

## U.S. Presidential Election (2016)

County-level election data for the **48 contiguous United States**. The response is the logit-transformed proportion of Republican votes. A county adjacency graph with **3,071** nodes and **8,669** edges is used, together with **55** demographic covariates.  

**Source:** https://github.com/tonmcg/US_County_Level_Election_Results_08-20

---

## U.S. Cancer Incidence

County-level cancer count data for **3,103** U.S. counties. A Poisson regression framework is used with county population as an offset. Geographic coordinates and socioeconomic variables are incorporated to construct spatial graph structures.  

**Source:** https://ghdx.healthdata.org/record/ihme-data/united-states-cancer-mortality-rates-county-1980-2014

---

## Air Pollution (2023)

County-level air quality data from **893** counties. The response is a binary indicator of whether polluted days were observed. Spatial graphs are constructed using county adjacency, latitude, longitude, and socioeconomic covariates.  

**Source:** https://www.epa.gov/air-trends/air-quality-cities-and-counties

---

## Flood Events

County-level binary indicators of flood occurrence for **3,037** counties. The same spatial adjacency structure and covariates as in the Air Pollution dataset are used.

## Simulation

Simulation data is provided in `sim_input.RData`.

## Model outputs

Pre-computed model outputs are stored in other `RData` files.
