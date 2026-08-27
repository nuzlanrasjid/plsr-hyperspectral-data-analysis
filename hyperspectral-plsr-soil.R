# Hyperspectral analysis with PLSR
# setwd("D:/Hyperspectral data")

### 1. Load libraries
library(prospectr)
library(prospectr)
library(pls)
library(stats)
library(openxlsx)


### 2. Load data
ds <- read.csv("soilmoisture_dataset.csv", check.names = FALSE)

### 3. Load Y-axes
y1 <- ds$soil_moisture
y2 <- ds$soil_temperature
y1_dens_score <- density(y1)
y2_dens_score <- density(y2)
# Normality test visual histogram 
hist(y1, freq = FALSE, col = "blue", density = c (85), 
     main = "Soil Moisture Distribution") 
hist(y2, freq = FALSE, col = "green", density = c (85), 
     main = "Soil Temperature Distribution")
# normality test visual-histogram
polygon(y1_dens_score, border = "black") 
polygon(y2_dens_score, border = "black") 
# normality test visual-boxplot
boxplot(y1, main = "Soil Moisture Distribution")
boxplot(y2, main = "Soil Temperature Distribution")
# normality test assumption_Shapiro-Wilk (normal if p value > 0.05), however
# it still can be running
shapiro.test(y1) 
shapiro.test(y2)

### 4. Load X-axes
# take all columns except the first, to check
X <- ds[, -1]
# check column names
cols <- colnames(X)
# remove columns that aren't numeric
spektrum <- suppressWarnings(as.numeric(cols))
# keep only valid spectrum columns, drop non-numeric ones
X <- X[, !is.na(spektrum)]
#update spectrum
spektrum <- spektrum[!is.na(spektrum)] # update wavelength
# filter spectrum
X <- X[, spektrum >= 454 & spektrum <= 950]
# convert to matrix
X <- as.matrix(X)


### 5. Pre-processing SNV (Optional)
#normalize spectrum
X_snv <- standardNormalVariate(as.matrix(X))
#check the visualization
matplot(spektrum[spektrum >= 454 & spektrum <= 950],
        t(X), type="l", main="Before SNV analysis")
matplot(spektrum[spektrum >= 454 & spektrum <= 950],
        t(X_snv), type="l", main="After SNV analysis")

### 6. Joining data to build the model
# joining dataframe specifically for the non-SNV model
ds_y1 <- data.frame(Target = y1, Spectra = I(X))
ds_y2 <- data.frame(Target = y2, Spectra = I(X))
# joining dataframe specifically for the SNV model
ds_y1_snv <- data.frame(Target = y1, Spectra = I(X_snv))
ds_y2_snv <- data.frame(Target = y2, Spectra = I(X_snv))


### 7. Running PLSR
# non-SNV
model_y1 <- plsr(Target ~ Spectra, 
              data = ds_y1, 
              scale = TRUE, 
              validation = "CV")
model_y2 <- plsr(Target ~ Spectra, 
                 data = ds_y2, 
                 scale = TRUE, 
                 validation = "CV")
# using SNV
model_y1_snv <- plsr(Target ~ Spectra, 
                 data = ds_y1_snv, 
                 scale = TRUE, 
                 validation = "CV")
model_y2_snv <- plsr(Target ~ Spectra, 
                 data = ds_y2_snv, 
                 scale = TRUE, 
                 validation = "CV")

### 8. summary model
# weigh based on adjCV, pick the one with the smallest adjCV
summary(model_y1) #ncomp_y1 <- 10
summary(model_y1_snv) #ncomp_y1_snv <- 7
summary(model_y2) #ncomp_y2 <- 9
summary(model_y2_snv) #ncomp_y2_snv <- 9

### 9. Choose the number of components
validationplot(model_y1, val.type = "RMSEP", main = "PLSR Soil moisture")


### 10. # evaluation 
#model_y1
for (i in 8:12) {
  pred_cv_y1 <- model_y1$validation$pred[, 1, i] # cross-validation prediction
  R2_cv_y1 <- cor(y1, pred_cv_y1)^2   # R2 CV
  rmsecv_val_y1 <- RMSEP(model_y1,
                     estimate = "CV")$val[1,1,i+1]   # R2 CV
  cat("ncomp =", i,
      "| R2 CV =", round(R2_cv_y1, 3),
      "| RMSECV =", round(rmsecv_val_y1, 3),
      "\n")
}
# model_y2
for (i in 5:10) {
  pred_cv_y2 <- model_y2$validation$pred[, 1, i] # cross-validation prediction
  R2_cv_y2 <- cor(y2, pred_cv_y2)^2   # R2 CV
  rmsecv_val_y2 <- RMSEP(model_y2,
                        estimate = "CV")$val[1,1,i+1]   # R2 CV
  cat("ncomp =", i,
      "| R2 CV =", round(R2_cv_y2, 3),
      "| RMSECV =", round(rmsecv_val_y2, 3),
      "\n")
}
# Final Plot y1
n_opt_y1 <- 9 
plot(y1,
     model_y1$validation$pred[,1,n_opt_y1], #cross-validation prediction from the PLSR model.
     main = paste("Measured vs Predicted (ncomp =", n_opt_y1, ")"), 
     xlab = "Measured Soil Moisture",
     ylab = "Predicted Soil Moisture",
     pch = 19,
     col = "darkblue")
abline(0,1,
       col = "red",
       lwd = 2)

# final plot y2
n_opt_y2 <- 9
plot(y2,
     model_y2$validation$pred[,1,n_opt_y2], #cross-validation prediction from the PLSR model.
     main = paste("Measured vs Predicted (ncomp =", n_opt_y2, ")"), 
     xlab = "Measured Soil Temperature",
     ylab = "Predicted Temperature",
     pch = 19,
     col = "darkgreen")
abline(0,1,
       col = "brown",
       lwd = 2)

# Evaluate Y1 values
pred_y1 <- model_y1$validation$pred[,1,n_opt_y1]
R2_y1 <- cor(y1, pred_y1)^2
RMSECV_y1 <- RMSEP(model_y1, estimate = "CV")$val[1,1,n_opt_y1 + 1]
RPD_y1 <- sd(y1) / RMSECV_y1
cat("=== Model Y1 (Soil Moisture) ===\n")
cat("R2 (CV) =", round(R2_y1,3), "\n")
cat("RMSECV =", round(RMSECV_y1,3), "\n")
cat("RPD =", round(RPD_y1,3), "\n\n")


#Evaluate Y2 values
pred_y2 <- model_y2$validation$pred[,1,n_opt_y2]
R2_y2 <- cor(y2, pred_y2)^2
RMSECV_y2 <- RMSEP(model_y2, estimate = "CV")$val[1,1,n_opt_y2 + 1]
RPD_y2 <- sd(y2) / RMSECV_y2
cat("=== Model Y2 (Soil Temperature) ===\n")
cat("R2 (CV) =", round(R2_y2,3), "\n")
cat("RMSECV =", round(RMSECV_y2,3), "\n")
cat("RPD =", round(RPD_y2,3), "\n\n")

### 11. Regression Plot
# model y1
coef_pls_y1 <- coef(model_y1, ncomp = n_opt_y1)

plot(spektrum,
     coef_pls_y1[,1,1],
     type = "l",
     lwd = 2,
     col = "darkblue",
     main = "PLSR Regression Coefficient - Soil Moisture",
     xlab = "Spectrum",
     ylab = "Coefficient")
abline(h = 0, col = "black", lty = 2)  # reference line at 0

#model y2
coef_pls_y2 <- coef(model_y2, ncomp = n_opt_y2)

plot(spektrum,
     coef_pls_y2[,1,1],
     type = "l",
     lwd = 2,
     col = "darkgreen",
     main = "PLSR Regression Coefficient - Soil Temperature",
     xlab = "Spectrum",
     ylab = "Coefficient")
abline(h = 0, col = "black", lty = 2)


### 13. Export workbook
## Create workbook
wb <- createWorkbook()

## Sheet 1 - Prediction
pred_df <- data.frame(
  Actual_SoilMoisture = y1,
  Predicted_SoilMoisture = pred_y1,
  Actual_SoilTemp = y2,
  Predicted_SoilTemp = pred_y2
)
addWorksheet(wb, "Prediction")
writeData(wb, "Prediction", pred_df)

## Sheet 2 - Metrics
metrics_df <- data.frame(
  Model = c("Soil_Moisture", "Soil_Temperature"),
  ncomp = c(n_opt_y1, n_opt_y2),
  R2_CV = c(R2_y1, R2_y2),
  RMSECV = c(RMSECV_y1, RMSECV_y2),
  RPD = c(RPD_y1, RPD_y2)
)
addWorksheet(wb, "Metrics")
writeData(wb, "Metrics", metrics_df)

## Sheet 3 - Regression Coefficient
coef_df <- data.frame(
  Wavelength = spektrum,
  Coefficient_SoilMoisture = coef_pls_y1[,1,1],
  Coefficient_SoilTemp = coef_pls_y2[,1,1]
)
addWorksheet(wb, "Regression_Coefficient")
writeData(wb, "Regression_Coefficient", coef_df)

### Save to Excel
saveWorkbook(wb, "PLSR_SoilMoisture_Temperature_Result.xlsx", overwrite = TRUE)
cat("Excel export completed\n")
