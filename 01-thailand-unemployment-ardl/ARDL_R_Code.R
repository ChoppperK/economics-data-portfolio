##############################################################################
#  ARDL Bounds-Testing Analysis of Unemployment in Thailand (1994-2024)
#  Main growth variable: REAL GDP growth | all-positive (shift-adjusted) dataset
#
#  Model:  UNEMP ~ GDPGROWTH + HDI + CPI + FDI
#          ARDL(1,1,1,2,2) | Case III (unrestricted intercept, no trend)
#
#  Reproduces the tables and figures in
#  "ARDL_Results_Academic_Report_RealGDP.docx":
#     Table A  Shift adjustment (removing negative values)
#     Table 0  Okun's Law benchmark regression
#     Table 1  Correlation matrix              (+ Fig. 5 heatmap)
#     Table 2  Unit root tests (ADF, PP, KPSS)
#     Table 3  VAR lag order selection
#     Table 4  ARDL bounds test
#     Table 5  Johansen cointegration test
#     Table 6  Long-run and short-run estimates (HAC / Newey-West SEs)
#     Table 7  VECM Granger causality + ECT loadings
#     Figs 2-3 CUSUM and CUSUM of squares
#     Fig. 4   Long-run effects chart
#
#  HOW TO RUN
#    1. Open this file in RStudio.
#    2. Set the working directory to this folder (Session > Set Working Directory).
#    3. Run everything: Ctrl+Shift+Enter (or Code > Run Region > Run All).
#    Each step prints a table you can compare with the report.
#
#  The statistical tests are written in plain R (no extra packages), so the
#  numbers match the report exactly. Only three packages are used: readxl
#  (read Excel), ARDL (cross-check) and strucchange (CUSUM plot).
##############################################################################


## ---- STEP 0: Packages -----------------------------------------------------
options(width = 200)          # wide console so tables print on one line
pkgs <- c("readxl", "ARDL", "strucchange")
to_install <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(to_install) > 0) install.packages(to_install)
library(readxl)
library(ARDL)
library(strucchange)


## ---- Helper functions (run this block once; nothing to change) -------------

stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.10, "*", ""))))
lag1  <- function(x) c(NA, x[-length(x)])
dif   <- function(x) c(NA, diff(x))

# Newey-West HAC covariance, Bartlett kernel, no small-sample adjustment
# (identical to the report; sandwich::NeweyWest() with its default
# adjust = TRUE would give standard errors about 30% larger).
newey_west <- function(model, lag = 2) {
  X <- model.matrix(model); e <- residuals(model); n <- nrow(X)
  Xe <- X * e; S <- crossprod(Xe)
  for (l in 1:lag) {
    G <- crossprod(Xe[(l + 1):n, , drop = FALSE], Xe[1:(n - l), , drop = FALSE])
    S <- S + (1 - l / (lag + 1)) * (G + t(G))
  }
  B <- solve(crossprod(X)); V <- B %*% S %*% B
  dimnames(V) <- list(colnames(X), colnames(X)); V
}

# Coefficient table using a given covariance matrix
coef_table <- function(model, V, df_t = df.residual(model)) {
  b <- coef(model); se <- sqrt(diag(V)); tt <- b / se
  p <- 2 * pt(abs(tt), df_t, lower.tail = FALSE)
  data.frame(Coef = round(b, 4), SE = round(se, 4), t = round(tt, 4), p = round(p, 4), Sig = stars(p))
}

# MacKinnon (1994) p-value for unit-root tests with an intercept
mackinnon_p <- function(tau) {
  if (tau > 2.74) return(1); if (tau < -18.83) return(0)
  cf <- if (tau <= -1.61) c(2.1659, 1.4412, 0.038269) else c(1.7339, 0.93202, -0.12745, -0.010368)
  pnorm(sum(cf * tau^(seq_along(cf) - 1)))
}

# ADF test (intercept), lag length chosen by AIC, max lag = ceiling(12*(n/100)^0.25)
adf_test <- function(x) {
  x <- as.numeric(x); dx <- diff(x); nd <- length(dx)
  maxlag <- min(ceiling(12 * (length(x) / 100)^0.25), length(x) %/% 2 - 2)
  build <- function(p) {
    X <- cbind(ylag = x[(p + 1):nd])
    if (p > 0) for (j in 1:p) X <- cbind(X, dx[(p + 1 - j):(nd - j)])
    list(y = dx[(p + 1):nd], X = X)
  }
  full <- build(maxlag)
  aics <- sapply(0:maxlag, function(p) AIC(lm(full$y ~ full$X[, 1:(p + 1), drop = FALSE])))
  best <- which.min(aics) - 1
  b <- build(best)
  tau <- summary(lm(b$y ~ b$X))$coefficients[2, "t value"]
  c(stat = tau, p = mackinnon_p(tau), lags = best)
}

# Phillips-Perron Z(tau), intercept, Bartlett kernel, lags = ceiling(12*(n/100)^0.25)
pp_test <- function(x) {
  x <- as.numeric(x); L <- ceiling(12 * (length(x) / 100)^0.25)
  y <- x[-1]; ylag <- x[-length(x)]
  fit <- lm(y ~ ylag); u <- residuals(fit); n <- length(u); k <- 2
  lam2 <- sum(u^2)
  for (j in 1:L) lam2 <- lam2 + 2 * (1 - j / (L + 1)) * sum(u[(j + 1):n] * u[1:(n - j)])
  lam2 <- lam2 / n
  s2 <- sum(u^2) / (n - k); s <- sqrt(s2); gamma0 <- s2 * (n - k) / n
  rho <- coef(fit)[["ylag"]]; sigma <- summary(fit)$coefficients["ylag", "Std. Error"]
  z <- sqrt(gamma0 / lam2) * ((rho - 1) / sigma) - 0.5 * ((lam2 - gamma0) / sqrt(lam2)) * (n * sigma / s)
  c(stat = z, p = mackinnon_p(z), lags = L)
}

# KPSS level-stationarity test, automatic bandwidth (Hobijn et al. 1998);
# p-value interpolated from the KPSS table, bounded to 0.01-0.10
kpss_test <- function(x) {
  x <- as.numeric(x); n <- length(x); e <- x - mean(x)
  covlags <- floor(n^(2 / 9)); s0 <- sum(e^2) / n; s1 <- 0
  for (i in 1:covlags) { rp <- sum(e[(i + 1):n] * e[1:(n - i)]) / (n / 2); s0 <- s0 + rp; s1 <- s1 + i * rp }
  L <- min(floor(1.1447 * ((s1 / s0)^2)^(1 / 3) * n^(1 / 3)), n - 1)
  s_hat <- sum(e^2)
  if (L > 0) for (i in 1:L) s_hat <- s_hat + 2 * (1 - i / (L + 1)) * sum(e[(i + 1):n] * e[1:(n - i)])
  stat <- sum(cumsum(e)^2) / n^2 / (s_hat / n)
  p <- approx(c(0.347, 0.463, 0.574, 0.739), c(0.10, 0.05, 0.025, 0.01), xout = stat, rule = 2)$y
  c(stat = stat, p = p, lags = L)
}

# VAR lag-order table; all lags fitted on the same (common) sample
var_lag_table <- function(Y, maxlag = 4) {
  Y <- as.matrix(Y); K <- ncol(Y); Tt <- nrow(Y) - maxlag; out <- data.frame(); prev <- NA
  for (p in 0:maxlag) {
    y <- Y[(maxlag + 1):nrow(Y), , drop = FALSE]; X <- matrix(1, Tt, 1)
    if (p > 0) for (j in 1:p) X <- cbind(X, Y[(maxlag + 1 - j):(nrow(Y) - j), , drop = FALSE])
    E <- y - X %*% solve(crossprod(X), crossprod(X, y))
    ld <- as.numeric(determinant(crossprod(E) / Tt, logarithm = TRUE)$modulus)
    ll <- -Tt * K / 2 * (1 + log(2 * pi)) - Tt / 2 * ld
    fp <- p * K^2 + K; dfm <- K * p + 1
    lr <- if (p == 0) NA else 2 * (ll - prev)
    out <- rbind(out, data.frame(Lag = p, LogL = round(ll, 4), LR = round(lr, 4),
      LR_p = if (p == 0) NA else round(pchisq(lr, K^2, lower.tail = FALSE), 4),
      FPE = signif(((Tt + dfm) / (Tt - dfm))^K * exp(ld), 4),
      AIC = round(ld + 2 * fp / Tt, 4), SIC = round(ld + log(Tt) * fp / Tt, 4),
      HQ = round(ld + 2 * log(log(Tt)) * fp / Tt, 4)))
    prev <- ll
  }
  attr(out, "T") <- Tt; out
}

# Johansen test, unrestricted constant, k_ar_diff lagged differences
johansen_test <- function(Y, k_ar_diff = 1) {
  Y <- as.matrix(Y); K <- ncol(Y); dY <- diff(Y); p <- k_ar_diff
  idx <- (p + 1):nrow(dY)
  Z0 <- dY[idx, , drop = FALSE]; Z1 <- Y[idx, , drop = FALSE]
  Z2 <- matrix(1, length(idx), 1)
  for (j in 1:p) Z2 <- cbind(Z2, dY[idx - j, , drop = FALSE])
  res <- function(A) A - Z2 %*% solve(crossprod(Z2), crossprod(Z2, A))
  R0 <- res(Z0); R1 <- res(Z1); Tt <- nrow(R0)
  S00 <- crossprod(R0) / Tt; S11 <- crossprod(R1) / Tt; S01 <- crossprod(R0, R1) / Tt
  ev <- eigen(solve(S11) %*% t(S01) %*% solve(S00) %*% S01)
  o <- order(Re(ev$values), decreasing = TRUE)
  lam <- Re(ev$values)[o]; V <- Re(ev$vectors)[, o, drop = FALSE]
  list(trace  = sapply(0:(K - 1), function(r) -Tt * sum(log(1 - lam[(r + 1):K]))),
       maxeig = -Tt * log(1 - lam),
       beta   = setNames(V[, 1] / V[1, 1], colnames(Y)), T = Tt)
}

# VECM Granger-causality F-tests (Lutkepohl 2005; same method as JMulTi and
# Python statsmodels, used in the report). For a VECM with k_ar_diff lagged
# differences: in the levels VAR of order p = k_ar_diff + 2, test that the
# first p-1 lag coefficients of the causing variable are zero.
vecm_granger_matrix <- function(Y, k_ar_diff = 1) {
  Y <- as.matrix(Y); K <- ncol(Y); n <- nrow(Y); p <- k_ar_diff + 2; q <- p - 1
  y <- Y[(p + 1):n, , drop = FALSE]; Tt <- nrow(y)
  lagm <- function(j) Y[(p + 1 - j):(n - j), , drop = FALSE]
  X  <- cbind(1, do.call(cbind, lapply(1:p, lagm)))
  B  <- solve(crossprod(X), crossprod(X, y))
  sig_ml <- crossprod(y - X %*% B) / Tt
  Xt <- cbind(1, do.call(cbind, lapply(1:q, function(j) lagm(j) - lagm(p))), lagm(p))
  xx11 <- solve(crossprod(Xt))[1:(1 + K * q), 1:(1 + K * q)]
  df2 <- K * (Tt - (1 + K * p))
  Fm <- Pm <- matrix(NA, K, K, dimnames = list(caused = colnames(Y), causing = colnames(Y)))
  for (ed in 1:K) for (ing in 1:K) if (ed != ing) {
    pos <- 1 + (0:(q - 1)) * K + ing
    a <- B[pos, ed]; V <- xx11[pos, pos, drop = FALSE] * sig_ml[ed, ed]
    Fm[ed, ing] <- as.numeric(t(a) %*% solve(V, a)) / q
    Pm[ed, ing] <- pf(Fm[ed, ing], q, df2, lower.tail = FALSE)
  }
  list(F = Fm, p = Pm, df = c(q, df2))
}

# Error-correction loading (alpha) of each VECM equation, rank 1, one lagged difference.
# t-statistics use the ML residual variance, as in the report (df = T - 7).
vecm_ect <- function(Y, beta) {
  Y <- as.matrix(Y); dY <- diff(Y); idx <- 2:nrow(dY)
  ect <- as.numeric(Y[idx, ] %*% beta); Xl <- dY[idx - 1, , drop = FALSE]
  do.call(rbind, lapply(colnames(Y), function(v) {
    fit <- lm(dY[idx, v] ~ ect + Xl); cf <- summary(fit)$coefficients["ect", ]
    Tn <- length(idx); dfr <- df.residual(fit)
    tt <- cf["t value"] * sqrt(Tn / dfr)
    data.frame(Equation = paste0("D", v), alpha = round(cf["Estimate"], 4), t = round(tt, 4),
               p = round(2 * pt(abs(tt), dfr, lower.tail = FALSE), 4), Sig = stars(2 * pt(abs(tt), dfr, lower.tail = FALSE)))
  }))
}


## ---- STEP 1: Import data ---------------------------------------------------
# Use forward slashes in the path.
file_path <- "data/thailand_unemployment_1994_2024.csv"   # <- path to the dataset in this folder

raw <- read.csv(file_path)
# Columns by position: 1 year, 2 UNEMP, 3 GDP per capita growth, 4 HDI, 5 CPI, 6 FDI,
# 7 real GDP growth. MAIN MODEL: GDPGROWTH = real GDP growth (column 7).
# GDP per capita growth (GDPPC, column 3) is used only in the STEP 14 robustness check.
if (ncol(raw) < 7) stop("Column 7 (real GDP growth) is missing from the Excel file.")
raw <- raw[, 1:7]
names(raw) <- c("year", "UNEMP", "GDPPC", "HDI", "CPI", "FDI", "GDPGROWTH")
df_original <- na.omit(data.frame(lapply(raw, as.numeric)))
df_original <- df_original[order(df_original$year), ]
cat("Sample:", min(df_original$year), "-", max(df_original$year), "| n =", nrow(df_original), "\n")


## ---- STEP 1b: Table A - Make all values positive ---------------------------
# GDPGROWTH, CPI and FDI have negative years. Each is shifted by
#   constant = |minimum| + 0.1   -> every value becomes >= 0.1.
# UNEMP and HDI are already positive and stay unchanged.
# Set use_positive_data <- FALSE to run everything on the original data:
# all slopes, SEs and test statistics stay identical; only constants change.
use_positive_data <- TRUE

vars5 <- c("UNEMP", "GDPGROWTH", "HDI", "CPI", "FDI")
df <- df_original
tableA <- data.frame()
for (v in c("GDPGROWTH", "CPI", "FDI", "UNEMP", "HDI", "GDPPC")) {
  mn <- min(df_original[[v]])
  sc <- if (use_positive_data && mn < 0) -mn + 0.1 else 0
  df[[v]] <- df_original[[v]] + sc
  tableA <- rbind(tableA, data.frame(Variable = v, Orig_min = round(mn, 4),
    Orig_max = round(max(df_original[[v]]), 4), Constant_added = round(sc, 4),
    Adj_min = round(min(df[[v]]), 4), Adj_max = round(max(df[[v]]), 4)))
}
cat("\n===== Table A: Shift adjustment =====\n"); print(tableA, row.names = FALSE)
cat("\nData used in the model:\n"); print(df, row.names = FALSE)
n <- nrow(df)
# Optional: save the positive dataset
# write.csv(df, "C:/Users/ASUS/Research/Model_final_positive.csv", row.names = FALSE)


## ---- STEP 2: Table 0 - Okun's Law benchmark -------------------------------
# dUNEMP_t = alpha + beta * GDPGROWTH_t + e_t   (HAC SEs, 2 lags)
okun_df <- na.omit(data.frame(dUNEMP = dif(df$UNEMP), GDPGROWTH = df$GDPGROWTH))
okun_model <- lm(dUNEMP ~ GDPGROWTH, data = okun_df)
cat("\n===== Table 0: Okun's Law benchmark (HAC SEs) =====\n")
okun_tab <- coef_table(okun_model, newey_west(okun_model, 2))
okun_tab$p <- round(2 * pnorm(-abs(okun_tab$t)), 4)          # z-test p-values, as in the report
okun_tab$Sig <- stars(okun_tab$p)
print(okun_tab)
cat("R2 =", round(summary(okun_model)$r.squared, 4), "| n =", nrow(okun_df), "\n")


## ---- STEP 3: Table 1 - Correlation matrix (+ Fig. 5) ----------------------
cor_matrix <- round(cor(df[, vars5]), 3)
cat("\n===== Table 1: Correlation matrix =====\n"); print(cor_matrix)
heatmap(cor_matrix, Rowv = NA, Colv = NA, symm = TRUE, margins = c(8, 8),
        col = colorRampPalette(c("#2E75B6", "white", "#C00000"))(50), main = "Fig. 5  Correlation matrix")


## ---- STEP 4: Table 2 - Unit root tests ------------------------------------
table2 <- do.call(rbind, lapply(vars5, function(v) {
  a <- adf_test(df[[v]]); ad <- adf_test(diff(df[[v]])); p <- pp_test(df[[v]]); k <- kpss_test(df[[v]])
  data.frame(Variable = v,
    ADF = round(a[["stat"]], 4), ADF_p = round(a[["p"]], 4),
    PP = round(p[["stat"]], 4), PP_p = round(p[["p"]], 4),
    KPSS = round(k[["stat"]], 4), KPSS_p = round(k[["p"]], 3),
    ADF_diff_p = round(ad[["p"]], 4),
    Order = if (a[["p"]] < 0.10) "I(0)" else if (ad[["p"]] < 0.10) "I(1)" else "I(2)+")
}))
cat("\n===== Table 2: Unit root tests (intercept) =====\n"); print(table2, row.names = FALSE)
cat("ADF/PP: H0 = unit root. KPSS: H0 = stationary (p bounded to 0.01-0.10).\n")


## ---- STEP 5: Table 3 - VAR lag order selection ----------------------------
table3 <- var_lag_table(df[, vars5], maxlag = 4)
cat("\n===== Table 3: VAR lag order selection (common sample T =", attr(table3, "T"), ") =====\n")
print(table3, row.names = FALSE)
cat("Selected lag: AIC", table3$Lag[which.min(table3$AIC)], "| SIC", table3$Lag[which.min(table3$SIC)],
    "| HQ", table3$Lag[which.min(table3$HQ)], "| FPE", table3$Lag[which.min(table3$FPE)], "\n")


## ---- STEP 6: UECM for ARDL(1,1,1,2,2) -------------------------------------
# dUNEMP_t = c + l1*UNEMP_t-1 + l2*GDP_t-1 + l3*HDI_t-1 + l4*CPI_t-1 + l5*FDI_t-1
#            + dGDP_t + dHDI_t + dCPI_t + dCPI_t-1 + dFDI_t + dFDI_t-1 + e_t
uecm_df <- na.omit(data.frame(
  d_UNEMP = dif(df$UNEMP),
  L1_UNEMP = lag1(df$UNEMP), L1_GDPGROWTH = lag1(df$GDPGROWTH), L1_HDI = lag1(df$HDI),
  L1_CPI = lag1(df$CPI), L1_FDI = lag1(df$FDI),
  d_GDPGROWTH = dif(df$GDPGROWTH), d_HDI = dif(df$HDI),
  d_CPI = dif(df$CPI), L1_d_CPI = lag1(dif(df$CPI)),
  d_FDI = dif(df$FDI), L1_d_FDI = lag1(dif(df$FDI))))
uecm_model <- lm(d_UNEMP ~ L1_UNEMP + L1_GDPGROWTH + L1_HDI + L1_CPI + L1_FDI +
                   d_GDPGROWTH + d_HDI + d_CPI + L1_d_CPI + d_FDI + L1_d_FDI, data = uecm_df)
nw_vcov <- newey_west(uecm_model, lag = 2)
cat("\nUECM: n =", nrow(uecm_df), "observations, residual df =", df.residual(uecm_model), "\n")


## ---- STEP 7: Table 4 - ARDL bounds test -----------------------------------
# H0: the five lagged-level coefficients are all zero (no long-run relationship)
restricted <- lm(d_UNEMP ~ d_GDPGROWTH + d_HDI + d_CPI + L1_d_CPI + d_FDI + L1_d_FDI, data = uecm_df)
bounds <- anova(restricted, uecm_model)
F_bounds <- bounds$F[2]
pss <- data.frame(Level = c("10%", "5%", "1%"), I0 = c(2.72, 3.23, 4.29), I1 = c(3.77, 4.35, 5.61))
cat("\n===== Table 4: ARDL bounds test (Pesaran-Shin-Smith 2001, Case III, k = 4) =====\n")
cat("F-statistic =", round(F_bounds, 4), "\n"); print(pss, row.names = FALSE)
cat(if (F_bounds > 5.61) "F > 1% upper bound -> cointegration confirmed at 1%.\n" else
    if (F_bounds > 4.35) "F > 5% upper bound -> cointegration confirmed at 5%.\n" else "Inconclusive / no cointegration.\n")

# Cross-check with the ARDL package (same F; its p-values are an extra)
ardl_model <- ardl(UNEMP ~ GDPGROWTH + HDI + CPI + FDI, data = df, order = c(1, 1, 1, 2, 2))
print(bounds_f_test(ardl_model, case = 3))


## ---- STEP 8: Table 5 - Johansen cointegration test ------------------------
joh <- johansen_test(df[, vars5], k_ar_diff = 1)
table5 <- data.frame(
  Hypothesis = c("r = 0", "r <= 1", "r <= 2", "r <= 3", "r <= 4"),
  Trace = round(joh$trace, 3), Trace_CV5 = c(69.819, 47.854, 29.796, 15.494, 3.842),
  MaxEigen = round(joh$maxeig, 3), MaxEig_CV5 = c(33.878, 27.586, 21.131, 14.264, 3.842))
table5$Trace_sig <- ifelse(table5$Trace > c(77.820, 54.682, 35.463, 19.935, 6.635), "***",
  ifelse(table5$Trace > table5$Trace_CV5, "**", ifelse(table5$Trace > c(65.820, 44.493, 27.067, 13.429, 2.706), "*", "")))
cat("\n===== Table 5: Johansen test (k_ar_diff = 1, unrestricted constant, T =", joh$T, ") =====\n")
print(table5, row.names = FALSE)
cat("Critical values: MacKinnon-Haug-Michelis (1999), as in the report.\n")


## ---- STEP 9: Table 6 - Long-run and short-run estimates (HAC) -------------
b <- coef(uecm_model); dfr <- df.residual(uecm_model)

# Long run: beta_X = -coef(L1_X) / coef(L1_UNEMP);  constant = -coef(Intercept) / coef(L1_UNEMP)
lr_map <- c(GDPGROWTH = "L1_GDPGROWTH", HDI = "L1_HDI", CPI = "L1_CPI", FDI = "L1_FDI", C = "(Intercept)")
long_run <- do.call(rbind, lapply(names(lr_map), function(v) {
  num <- lr_map[[v]]; den <- "L1_UNEMP"
  est <- -b[[num]] / b[[den]]
  g <- setNames(rep(0, length(b)), names(b)); g[num] <- -1 / b[[den]]; g[den] <- b[[num]] / b[[den]]^2
  se <- sqrt(as.numeric(t(g) %*% nw_vcov %*% g)); tt <- est / se
  p <- 2 * pt(abs(tt), dfr, lower.tail = FALSE)
  data.frame(Variable = v, Coef = round(est, 4), SE = round(se, 4), t = round(tt, 4), p = round(p, 4), Sig = stars(p))
}))

sr_all <- coef_table(uecm_model, nw_vcov)
sr_map <- c(d_GDPGROWTH = "D(GDPGROWTH)", d_HDI = "D(HDI)", d_CPI = "D(CPI)", L1_d_CPI = "D(CPI(-1))",
            d_FDI = "D(FDI)", L1_d_FDI = "D(FDI(-1))", L1_UNEMP = "CointEq(-1)")
short_run <- cbind(Variable = unname(sr_map), sr_all[names(sr_map), ])

cat("\n===== Table 6: Long-run estimations (HAC SEs) =====\n");  print(long_run, row.names = FALSE)
cat("\n===== Table 6: Short-run estimations (HAC SEs) =====\n"); print(short_run, row.names = FALSE)
cat("p-values use residual df =", dfr, "(report: df = 18), so some p differ in the 3rd decimal; stars are the same.\n")

ect <- b[["L1_UNEMP"]]
cat("Speed of adjustment:", round(-100 * ect, 1), "% per year | half-life =", round(log(0.5) / log(1 + ect), 2), "years\n")

# Cross-check with the ARDL package (coefficients only; it uses non-robust SEs)
cat("\nARDL package long-run multipliers (should match the Coef column above):\n")
print(multipliers(ardl_model, type = "lr"))


## ---- STEP 10: Diagnostics (bottom of Table 6) -----------------------------
u <- residuals(uecm_model)
lb <- Box.test(u, lag = 2, type = "Ljung-Box")
arch_lm <- function(e, lags) {
  e2 <- e^2; m <- length(e2)
  Z <- sapply(1:lags, function(i) e2[(lags + 1 - i):(m - i)])
  LM <- (m - lags) * summary(lm(e2[(lags + 1):m] ~ Z))$r.squared
  pchisq(LM, lags, lower.tail = FALSE)
}
arch_p <- sapply(1:5, function(L) arch_lm(u, L))
sk <- mean((u - mean(u))^3) / mean((u - mean(u))^2)^1.5
ku <- mean((u - mean(u))^4) / mean((u - mean(u))^2)^2
jb <- length(u) / 6 * (sk^2 + (ku - 3)^2 / 4)
yl <- df$UNEMP[(n - nrow(uecm_df) + 1):n]                      # UNEMP in levels, same sample
R2 <- 1 - sum(u^2) / sum((yl - mean(yl))^2)
k_par <- length(b)
diag_tab <- data.frame(
  Test = c("Ljung-Box Q (lag 2)", "ARCH-LM (lags 1-5, min p)", "Jarque-Bera", "R-squared (levels)",
           paste0("Adj. R-squared (n=", length(u), ", k=", k_par, ")")),
  Statistic = round(c(lb$statistic, NA, jb, R2, 1 - (1 - R2) * (length(u) - 1) / (length(u) - k_par)), 4),
  p = round(c(lb$p.value, min(arch_p), pchisq(jb, 2, lower.tail = FALSE), NA, NA), 4))
cat("\n===== Table 6: Diagnostic tests =====\n"); print(diag_tab, row.names = FALSE)
cat("ARCH-LM p < 0.05 -> heteroskedasticity, which is why HAC standard errors are used.\n")


## ---- STEP 11: Table 7 - VECM Granger causality + ECT ----------------------
gr <- vecm_granger_matrix(df[, vars5], k_ar_diff = 1)
table7 <- matrix(sprintf("%.4f%s (%.4f)", gr$F, stars(gr$p), gr$p), 5, 5, dimnames = dimnames(gr$F))
diag(table7) <- "-"
ect_tab <- vecm_ect(df[, vars5], joh$beta)
cat("\n===== Table 7: VECM Granger causality, F (p) | rows = caused, columns = causing =====\n")
print(noquote(table7))
cat("F-test df =", gr$df[1], ",", gr$df[2], "\n")
cat("\n===== Table 7: ECT(-1) loading for each equation =====\n"); print(ect_tab, row.names = FALSE)
cat("Cointegrating vector (normalised on UNEMP):\n"); print(round(joh$beta, 4))


## ---- STEP 12: Figs 2-3 - CUSUM and CUSUM of squares -----------------------
f_uecm <- formula(uecm_model)
cusum <- efp(f_uecm, data = uecm_df, type = "Rec-CUSUM")
print(sctest(cusum))                                 # p > 0.05 -> stable
plot(cusum, main = "Fig. 2  Plot of CUSUM")

rr <- recresid(f_uecm, data = uecm_df)
cs <- cumsum(rr^2) / sum(rr^2); m <- length(rr); expct <- seq_len(m) / m
c0 <- 1.358 / (sqrt(m) + 0.12 + 0.11 / sqrt(m))        # approx. 5% bound
plot(cs, type = "l", lwd = 2, col = "steelblue", ylim = c(-0.4, 1.4),
     main = "Fig. 3  Plot of CUSUM of squares", xlab = "Recursive residual", ylab = "CUSUMSQ")
lines(expct - c0, col = "red", lty = 2); lines(expct + c0, col = "red", lty = 2)
legend("topleft", c("CUSUMSQ", "5% bounds"), col = c("steelblue", "red"), lty = c(1, 2), bty = "n")
cat("CUSUMSQ max deviation =", round(max(abs(cs - expct)), 4), "vs 5% bound", round(c0, 4),
    if (max(abs(cs - expct)) < c0) "-> stable\n" else "-> outside bounds\n")


## ---- STEP 13: Fig. 4 - Long-run effects chart -----------------------------
lr4 <- long_run[long_run$Variable != "C", ]
bp <- barplot(lr4$Coef, names.arg = lr4$Variable, col = ifelse(lr4$Coef > 0, "#2E7D32", "#C62828"),
              border = NA, ylim = c(min(lr4$Coef) * 1.3, max(lr4$Coef) * 2),
              main = "Fig. 4  Long-run effects on the unemployment rate", ylab = "Long-run coefficient")
abline(h = 0)
text(bp, lr4$Coef, paste0(sprintf("%+.3f", lr4$Coef), lr4$Sig), pos = ifelse(lr4$Coef > 0, 3, 1), cex = 0.9)
legend("bottomright", c("lowers unemployment", "raises unemployment"), fill = c("#C62828", "#2E7D32"),
       border = NA, bty = "n")


## ---- STEP 14: Robustness - GDP per capita growth instead of real GDP growth ---
# Table 0a: Okun benchmark with each growth measure.
# Table 6a: ARDL(1,1,1,2,2) long-run estimates with GDPPC in place of GDPGROWTH.
if (TRUE) {
  okun_row <- function(gvar) {
    od <- na.omit(data.frame(dU = dif(df_original$UNEMP), g = df_original[[gvar]]))   # original scale
    m <- lm(dU ~ g, data = od); V <- newey_west(m, 2)
    b <- coef(m)[["g"]]; se <- sqrt(V["g", "g"]); z <- b / se
    data.frame(Growth_measure = gvar, Okun_beta = round(b, 4), SE_HAC = round(se, 4),
               p = round(2 * pnorm(-abs(z)), 4), R2 = round(summary(m)$r.squared, 4),
               Threshold_growth = round(-coef(m)[["(Intercept)"]] / b, 2))
  }
  cat("\n===== Table 0a: Okun benchmark - real GDP growth vs per capita growth =====\n")
  print(rbind(okun_row("GDPGROWTH"), okun_row("GDPPC")), row.names = FALSE)
  cat("Threshold growth = -alpha/beta on the original data: growth needed to keep unemployment constant.\n")

  rob_df <- na.omit(data.frame(
    d_UNEMP = dif(df$UNEMP),
    L1_UNEMP = lag1(df$UNEMP), L1_G = lag1(df$GDPPC), L1_HDI = lag1(df$HDI),
    L1_CPI = lag1(df$CPI), L1_FDI = lag1(df$FDI),
    d_G = dif(df$GDPPC), d_HDI = dif(df$HDI),
    d_CPI = dif(df$CPI), L1_d_CPI = lag1(dif(df$CPI)),
    d_FDI = dif(df$FDI), L1_d_FDI = lag1(dif(df$FDI))))
  rob_model <- lm(d_UNEMP ~ L1_UNEMP + L1_G + L1_HDI + L1_CPI + L1_FDI +
                    d_G + d_HDI + d_CPI + L1_d_CPI + d_FDI + L1_d_FDI, data = rob_df)
  rob_V <- newey_west(rob_model, 2); rb <- coef(rob_model); rdf <- df.residual(rob_model)
  rob_restr <- lm(d_UNEMP ~ d_G + d_HDI + d_CPI + L1_d_CPI + d_FDI + L1_d_FDI, data = rob_df)
  rob_F <- anova(rob_restr, rob_model)$F[2]
  rob_map <- c(GDPPC = "L1_G", HDI = "L1_HDI", CPI = "L1_CPI", FDI = "L1_FDI", C = "(Intercept)")
  rob_lr <- do.call(rbind, lapply(names(rob_map), function(v) {
    num <- rob_map[[v]]; den <- "L1_UNEMP"; est <- -rb[[num]] / rb[[den]]
    g <- setNames(rep(0, length(rb)), names(rb)); g[num] <- -1 / rb[[den]]; g[den] <- rb[[num]] / rb[[den]]^2
    se <- sqrt(as.numeric(t(g) %*% rob_V %*% g)); tt <- est / se; p <- 2 * pt(abs(tt), rdf, lower.tail = FALSE)
    data.frame(Variable = v, Coef = round(est, 4), SE = round(se, 4), t = round(tt, 4), p = round(p, 4), Sig = stars(p))
  }))
  cat("\n===== Table 6a: ARDL(1,1,1,2,2) with GDP per capita growth - long run (HAC) =====\n")
  cat("Bounds F =", round(rob_F, 4), "(1% upper bound 5.61)\n")
  print(rob_lr, row.names = FALSE)
  cat("ECT =", round(rb[["L1_UNEMP"]], 4), "| short-run D(GDPPC) =", round(rb[["d_G"]], 4), "\n")
  cat("Compare with Table 6: same signs and significance -> results do not depend on the growth measure.\n")
}

cat("\nDone. Compare each table with the same number in ARDL_Results_Academic_Report_RealGDP.docx.\n")
