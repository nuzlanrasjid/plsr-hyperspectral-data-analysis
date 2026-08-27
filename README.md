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
| Soil Moisture | ... | ... | ... |
| Soil Temperature | ... | ... | ... |

![Measured vs predicted](output/measured-vs-predicted-moisture.png)

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
