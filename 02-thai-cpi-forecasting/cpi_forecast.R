# getwd()
# ==============================================================================
# FINAL PROJECT: Forecasting Thai Macroeconomic Variable (CPI)
# ==============================================================================

# 1. SETUP & DATA LOADING ------------------------------------------------------
# Install packages if not already installed (uncomment to run)
# install.packages(c("ggplot2", "car", "lmtest", "forecast", "tseries"))

library(ggplot2)
library(car)      # For VIF
library(lmtest)   # For Breusch-Pagan, DW
# install.packages("forecast")  # run once if not installed
library(forecast) # For accuracy metrics

# Load Data
data <- read.csv("Data.csv")

# Convert Date
head(data$date)
data$date <- as.Date(data$date, format="%m/%d/%Y")
data$year <- as.numeric(format(data$date, "%Y"))

# Sort by year just in case
data <- data[order(data$year), ]

# Inspect Data
str(data)
summary(data)

# ==============================================================================
# (b) DATA DESCRIPTION & PLOTS
# ==============================================================================

# Time Series Plot (CPI)
plot(data$year, data$cpi, type="o", col="blue", lwd=2, 
     main="Thai CPI Over Time (2005-2024)", xlab="Year", ylab="CPI Growth (%)")
grid()

# Scatter Plots
par(mfrow=c(1,3))  # 1 row, 3 plots

plot(data$GDPTH, data$cpi,
     pch=19, col="orange",
     main="CPI vs GDPTH",
     xlab="GDPTH", ylab="CPI")

plot(data$POPTH, data$cpi,
     pch=19, col="green",
     main="CPI vs POPTH",
     xlab="Population", ylab="CPI")

plot(data$GAN, data$cpi,
     pch=19, col="red",
     main="CPI vs GAN",
     xlab="GAN", ylab="CPI")

par(mfrow=c(1,1))


# ==============================================================================
# (c) BASELINE REGRESSION MODEL (Train: 2005-2018)
# ==============================================================================
# Train/Test Split
train <- subset(data, year <= 2018)
test  <- subset(data, year >= 2019)

# Model 1: Simple OLS (CPI ~ Rice + Oil + FX)
model_base <- lm(cpi ~ GDPTH + POPTH + fx, data=train)

cat("\n--- BASELINE MODEL SUMMARY ---\n")
summary(model_base)

# Interpretation: Check coefficients for Rice, Oil, FX.

# ==============================================================================
# (d) DIAGNOSTIC TESTS
# ==============================================================================
cat("\n--- DIAGNOSTICS ---\n")

# 1. Multicollinearity (VIF)
# If VIF > 5 or 10, high multicollinearity.
vif_vals <- vif(model_base)
print(vif_vals)

# 2. Heteroscedasticity (Breusch-Pagan Test)
# H0: Homoscedasticity (Constant variance). If p < 0.05, reject H0 (Bad).
bp_test <- bptest(model_base)
print(bp_test)

# 3. Autocorrelation (Durbin-Watson)
# Range 0-4. Near 2 is good. Near 0 or 4 is bad.
dw_test <- dwtest(model_base)
print(dw_test)

# ==============================================================================
# (e) IMPROVED MODEL (Re-Specification)
# ==============================================================================
# Strategy: Use AR(1) Model (Autoregressive) because inflation is persistent.
# We add a lagged variable: CPI(t-1)

# Create Lag in full dataset
data$cpi_lag1 <- c(NA, head(data$cpi, -1))

# Re-split data to include the lag column
train_imp <- subset(data, year <= 2018)
test_imp  <- subset(data, year >= 2019)

# Model 2: AR(1) Model
model_imp <- lm(cpi ~ cpi_lag1, data=train_imp)

cat("\n--- IMPROVED MODEL (AR1) SUMMARY ---\n")
summary(model_imp)

# Justification:
# The lag term is likely significant (p < 0.05), capturing the trend better 
# than the structural variables in this small sample.

# ==============================================================================
# (f) PREDICTION & FORECAST EVALUATION
# ==============================================================================

# Predict on Test Set (2019-2024)
pred_base <- predict(model_base, newdata=test)

# For AR(1), we need dynamic prediction or static one-step.
# Here we use the actual lag values from the test set (Static Forecast).
pred_imp  <- predict(model_imp, newdata=test_imp)

# Calculate Accuracy Metrics
cat("\n--- FORECAST ACCURACY (Test Set: 2019-2024) ---\n")
acc_base <- accuracy(pred_base, test$cpi)
acc_imp  <- accuracy(pred_imp, test_imp$cpi)

print("Baseline Model Metrics:")
print(acc_base[, c("MAE", "RMSE", "MAPE")])

print("Improved Model (AR1) Metrics:")
print(acc_imp[, c("MAE", "RMSE", "MAPE")])

# ==============================================================================
# (g) CONCLUSION
# ==============================================================================
# Summarize which model had lower MAE/RMSE.
# Discuss limitations: Small sample size (N=20), annual frequency is low.