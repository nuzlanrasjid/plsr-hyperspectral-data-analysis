# Predicting Soil Moisture & Temperature from Hyperspectral Data (PLSR)

Partial Least Squares Regression (PLSR) models in R to predict soil
moisture and temperature from hyperspectral reflectance data
(454–950 nm), comparing raw spectra against SNV-preprocessed spectra.

## Dataset

Sourced from [Kaggle: Hyperspectral Benchmark Dataset on Soil Moisture](https://www.kaggle.com/datasets/binaryjoker/hyperspectral-benchmark-dataset-on-soil-moisture).

## Method

1. Normality check on target variables (histograms, boxplots, Shapiro-Wilk)
2. Spectral filtering to 454–950 nm range
3. SNV preprocessing (compared against raw spectra)
4. PLSR with cross-validation, component count chosen by lowest RMSECV
5. Evaluation via R² (CV), RMSECV, and RPD

## Results

| Model | R² (CV) | RMSECV | RPD |
|---|---|---|---|
| Soil Moisture | 0.875331727 | 1.287714091 | 2.830872324 |
| Soil Temperature | 0.834133176 | 1.899797441 | 2,453210554 |

![Measured vs predicted soil moisture](/output/Measured%20vs%20Predicted%20soil%20moisture.png)

![Measured vs predicted soil temperature](/output/Measured%20vs%20Predicted%20soil%20temperature.png)

## How to run

```r
source("R/plsr_analysis.R")
```

Requires: `prospectr`, `pls`, `openxlsx`

## Repo structure

```
├── R/
│   └── plsr_analysis.R
├── data/
│   └── soilmoisture_dataset.csv
├── output/
│   └── (plots & exported Excel results)
└── README.md
```
