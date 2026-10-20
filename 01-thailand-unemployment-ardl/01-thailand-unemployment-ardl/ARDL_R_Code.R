##############################################################################
#  Determinants of Unemployment in Thailand: ARDL Bounds Testing, 1994-2024
#  Final model: ARDL(2,1,0,2,0) | Case III (unrestricted intercept, no trend)
#  UNEMP ~ GDPGROWTH (real GDP growth) + HDI + CPI (inflation) + FDI (% of GDP)
#
#  Reproduces the tables and figures of the research report.
#  Data file: Thailand_Unemployment_Data_Final.xlsx (CPI and FDI from the
#  World Bank WDI, series FP.CPI.TOTL.ZG and BX.KLT.DINV.WD.GD.ZS).
#
#  HOW TO RUN
#    1. Open this file in RStudio.
#    2. In STEP 1, set file_path to where the data file is.
#    3. Run everything: Ctrl+Shift+Enter. Each step prints one table of the report.
#
#  The tests are written in plain R, so the numbers match the report. Packages:
#  readxl (read Excel), ARDL (cross-check) and strucchange (CUSUM plot).
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

# Zivot-Andrews (1992) unit-root test with one endogenous break in the intercept
# (Model A). Same algorithm as Python statsmodels zivot_andrews (Baum 2004):
# lag length chosen once by AIC on the ADF(constant + trend) regression, max 2
# lags; break dates searched over the middle 70% of the sample (trim = 15%).
# Critical values (Zivot & Andrews 1992): 1% -5.34, 5% -4.80, 10% -4.58;
# statsmodels uses -5.28 / -4.81 / -4.57 for this model.
za_test <- function(x, years, maxlag = 2, trim = 0.15) {
  x <- as.numeric(x); n <- length(x); dx <- diff(x)
  aic <- sapply(0:maxlag, function(L) {
    i <- (maxlag + 1):(n - 1); X <- cbind(1, i, x[i])
    if (L > 0) for (j in 1:L) X <- cbind(X, dx[i - j])
    m <- lm(dx[i] ~ X - 1); -2 * as.numeric(logLik(m)) + 2 * ncol(X) })
  L <- which.min(aic) - 1
  rows <- (L + 1):(n - 1); trimcnt <- floor(n * trim)
  bps <- (trimcnt + 1):(n - trimcnt)
  stats <- sapply(bps, function(bp) {
    DU <- as.numeric(seq_along(rows) > bp - (L + 1))
    X <- cbind(1, DU, seq_along(rows), x[rows])
    if (L > 0) for (j in 1:L) X <- cbind(X, dx[rows - j])
    summary(lm(dx[rows] ~ X - 1))$coefficients[4, 3] })
  c(stat = min(stats), break_year = years[bps[which.min(stats)]], lags = L)
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
file_path <- "data/Thailand_Unemployment_Data_Final.xlsx"   # <- or the full path, e.g. "C:/Users/ASUS/.../Thailand_Unemployment_Data_Final.xlsx"
if (!file.exists(file_path) && interactive()) file_path <- file.choose()   # a window opens: pick the data file

raw <- read_excel(file_path)
# Columns are found BY HEADER NAME, so their order in the Excel file does not matter.
hdr <- tolower(gsub("\\s+", " ", trimws(names(raw))))
find_col <- function(pattern, exclude = NULL, required = TRUE, label = pattern) {
  hit <- grepl(pattern, hdr, ignore.case = TRUE)
  if (!is.null(exclude)) hit <- hit & !grepl(exclude, hdr, ignore.case = TRUE)
  if (sum(hit) == 0) {
    if (required) stop("Cannot find the '", label, "' column. Headers found: ", paste(names(raw), collapse = " | "))
    return(NA)
  }
  which(hit)[1]
}
cols <- c(year      = find_col("^year",                     label = "year"),
          UNEMP     = find_col("unemploy",                  label = "Unemployment rate"),
          GDPGROWTH = find_col("gdp.*growth", exclude = "per capita", label = "Real GDP growth"),
          HDI       = find_col("hdi|human development",     label = "HDI"),
          CPI       = find_col("cpi|consumer price|inflation", label = "Inflation (CPI)"),
          FDI       = find_col("foreign direct|fdi",        label = "FDI"),
          GDPPC     = find_col("per capita", required = FALSE))
has_gdppc <- !is.na(cols[["GDPPC"]])
cols <- cols[!is.na(cols)]
cat("Columns used:\n"); print(setNames(names(raw)[cols], names(cols)))
df_original <- data.frame(lapply(raw[, cols], function(x) suppressWarnings(as.numeric(x))))
names(df_original) <- names(cols)
df_original <- na.omit(df_original[!is.na(df_original$year), ])
df_original <- df_original[order(df_original$year), ]
cat("Sample:", min(df_original$year), "-", max(df_original$year), "| n =", nrow(df_original), "\n")


## ---- STEP 2: Shift adjustment (Chapter 3.3) --------------------------------
# GDPGROWTH, CPI and FDI have negative years; each is shifted by |min| + 0.1.
# Only the constants change; every slope, SE and test statistic is identical.
use_positive_data <- TRUE
vars5 <- c("UNEMP", "GDPGROWTH", "HDI", "CPI", "FDI")
df <- df_original
shiftA <- data.frame()
for (v in c("GDPGROWTH", "CPI", "FDI", "UNEMP", "HDI", if (has_gdppc) "GDPPC")) {
  mn <- min(df_original[[v]]); sc <- if (use_positive_data && mn < 0) -mn + 0.1 else 0
  df[[v]] <- df_original[[v]] + sc
  shiftA <- rbind(shiftA, data.frame(Variable = v, Min = round(mn, 4), Constant_added = round(sc, 4)))
}
cat("\n===== Shift adjustment =====\n"); print(shiftA, row.names = FALSE)
n <- nrow(df)


## ---- STEP 3: Table 3 - Descriptive statistics (original data) --------------
desc <- do.call(rbind, lapply(c(vars5, if (has_gdppc) "GDPPC"), function(v) {
  x <- df_original[[v]]; m <- mean(x); s2 <- mean((x - m)^2)
  sk <- mean((x - m)^3) / s2^1.5; ku <- mean((x - m)^4) / s2^2; jb <- length(x) / 6 * (sk^2 + (ku - 3)^2 / 4)
  data.frame(Variable = v, Mean = round(m, 3), SD = round(sd(x), 3), Min = round(min(x), 3), Max = round(max(x), 3),
             Skewness = round(sk, 3), Kurtosis = round(ku, 3), JB = round(jb, 3), JB_p = round(pchisq(jb, 2, lower.tail = FALSE), 3))
}))
cat("\n===== Table 3: Descriptive statistics =====\n"); print(desc, row.names = FALSE)


## ---- STEP 4: Table 4 - Correlation matrix and VIF --------------------------
cat("\n===== Table 4: Correlation matrix =====\n"); print(round(cor(df[, vars5]), 3))
xv <- vars5[-1]
vif <- sapply(xv, function(v) 1 / (1 - summary(lm(reformulate(setdiff(xv, v), v), data = df))$r.squared))
cat("VIF:\n"); print(round(vif, 3))


## ---- STEP 5: Table 5 - Unit root tests -------------------------------------
table5 <- do.call(rbind, lapply(vars5, function(v) {
  a <- adf_test(df[[v]]); ad <- adf_test(diff(df[[v]])); p <- pp_test(df[[v]]); pd <- pp_test(diff(df[[v]]))
  k <- kpss_test(df[[v]]); kd <- kpss_test(diff(df[[v]]))
  data.frame(Variable = v, ADF = round(a[["stat"]], 3), ADF_p = round(a[["p"]], 3),
    ADF_d = round(ad[["stat"]], 3), ADF_d_p = round(ad[["p"]], 3), PP = round(p[["stat"]], 3), PP_p = round(p[["p"]], 3),
    PP_d = round(pd[["stat"]], 3), KPSS = round(k[["stat"]], 3), KPSS_p = round(k[["p"]], 3), KPSS_d = round(kd[["stat"]], 3),
    Order = if (a[["p"]] < 0.10) "I(0)" else if (ad[["p"]] < 0.10) "I(1)" else "I(2)+")
}))
cat("\n===== Table 5: Unit root tests (intercept) =====\n"); print(table5, row.names = FALSE)
cat("Small differences from the report can occur in the PP/KPSS lag choice; the order of integration is the same.\n")


## ---- STEP 5b: Table 6 - Zivot-Andrews unit root test with a structural break -
yrs <- df_original$year   # break dates are reported as calendar years
za_crit <- c(`1%` = -5.28, `5%` = -4.81, `10%` = -4.57)
za_star <- function(s) if (s < za_crit[1]) "***" else if (s < za_crit[2]) "**" else if (s < za_crit[3]) "*" else ""
table6 <- do.call(rbind, lapply(vars5, function(v) {
  a <- za_test(df[[v]], yrs); b <- za_test(diff(df[[v]]), yrs[-1])
  data.frame(Variable = v, ZA_level = paste0(round(a[["stat"]], 3), za_star(a[["stat"]])), Break_level = a[["break_year"]],
    ZA_diff = paste0(round(b[["stat"]], 3), za_star(b[["stat"]])), Break_diff = b[["break_year"]])
}))
cat("\n===== Table 6: Zivot-Andrews test (break in intercept) =====\n"); print(table6, row.names = FALSE)
cat("Critical values: 1% -5.28, 5% -4.81, 10% -4.57. No series is I(2), so the bounds test is valid.\n")

## ---- STEP 6: Table 7 - Okun's law benchmark --------------------------------
okun_row <- function(gvar) {
  od <- na.omit(data.frame(dU = dif(df_original$UNEMP), g = df_original[[gvar]]))     # original data
  m <- lm(dU ~ g, data = od); V <- newey_west(m, 2); b <- coef(m); se <- sqrt(V["g", "g"]); z <- b[["g"]] / se
  data.frame(Growth = gvar, Alpha = round(b[[1]], 4), Okun_beta = round(b[["g"]], 4), SE_HAC = round(se, 4),
             p = round(2 * pnorm(-abs(z)), 4), R2 = round(summary(m)$r.squared, 4), Threshold = round(-b[[1]] / b[["g"]], 2))
}
cat("\n===== Table 7: Okun's law benchmark (HAC SEs, 2 lags) =====\n")
print(rbind(okun_row("GDPGROWTH"), if (has_gdppc) okun_row("GDPPC")), row.names = FALSE)


## ---- STEP 7: Lag selection by AIC (Appendix C) -----------------------------
# All ARDL(p, q1..q4) with p = 1..2 and q = 0..2, fitted on the same sample (1996-2024).
build_ardl <- function(d, p, q, gname = "GDPGROWTH") {
  X <- data.frame(y = d$UNEMP)
  for (j in 1:p) X[[paste0("UNEMP_L", j)]] <- c(rep(NA, j), head(d$UNEMP, -j))
  for (v in names(q)) for (j in 0:q[[v]]) X[[paste0(v, "_L", j)]] <- if (j == 0) d[[if (v == "GROWTH") gname else v]] else c(rep(NA, j), head(d[[if (v == "GROWTH") gname else v]], -j))
  X
}
grid <- expand.grid(p = 1:2, g = 0:2, h = 0:2, c = 0:2, f = 0:2)
grid$AIC <- apply(grid, 1, function(r) {
  X <- build_ardl(df, r[["p"]], c(GROWTH = r[["g"]], HDI = r[["h"]], CPI = r[["c"]], FDI = r[["f"]]))[3:n, ]
  m <- lm(y ~ ., data = X); k <- length(coef(m)) + 1; nn <- nrow(X)
  -2 * as.numeric(logLik(m)) + 2 * (k - 1)       # same ranking as statsmodels ARDL AIC
})
grid <- grid[order(grid$AIC), ]
cat("\n===== Lag selection (five lowest AIC) =====\n"); print(head(grid, 5), row.names = FALSE)
cat("Selected: ARDL(", paste(unlist(grid[1, 1:5]), collapse = ","), ")  -> the report uses ARDL(2,1,0,2,0)\n")


## ---- STEP 8: UECM for ARDL(2,1,0,2,0) --------------------------------------
# dUNEMP_t = c + l1*UNEMP_t-1 + l2*GDP_t-1 + l3*HDI_t + l4*CPI_t-1 + l5*FDI_t
#            + dUNEMP_t-1 + dGDP_t + dCPI_t + dCPI_t-1 + e_t
# (HDI and FDI have zero lags, so they enter at time t; this is an exact
#  re-parameterisation of the ARDL model.)
make_uecm <- function(d, gvar = "GDPGROWTH") na.omit(data.frame(
  d_UNEMP = dif(d$UNEMP), L1_UNEMP = lag1(d$UNEMP), L1_G = lag1(d[[gvar]]), HDI = d$HDI,
  L1_CPI = lag1(d$CPI), FDI = d$FDI, L1_d_UNEMP = lag1(dif(d$UNEMP)), d_G = dif(d[[gvar]]),
  d_CPI = dif(d$CPI), L1_d_CPI = lag1(dif(d$CPI))))
uecm_df <- make_uecm(df)
uecm_model <- lm(d_UNEMP ~ L1_UNEMP + L1_G + HDI + L1_CPI + FDI + L1_d_UNEMP + d_G + d_CPI + L1_d_CPI, data = uecm_df)
nw_vcov <- newey_west(uecm_model, lag = 2)
cat("\nUECM: n =", nrow(uecm_df), "| residual df =", df.residual(uecm_model), "\n")


## ---- STEP 9: Table 8 - ARDL bounds test ------------------------------------
restricted <- lm(d_UNEMP ~ L1_d_UNEMP + d_G + d_CPI + L1_d_CPI, data = uecm_df)
F_bounds <- anova(restricted, uecm_model)$F[2]
pss <- data.frame(Level = c("10%", "5%", "1%"), I0 = c(2.45, 2.86, 3.74), I1 = c(3.52, 4.01, 5.06))
cat("\n===== Table 8: ARDL bounds test (Pesaran et al. 2001, Case III, k = 4) =====\n")
cat("F-statistic =", round(F_bounds, 3), "\n"); print(pss, row.names = FALSE)
cat(if (F_bounds > 5.06) "F > 1% upper bound -> cointegration confirmed at 1%.\n" else "Check the bounds.\n")
ardl_model <- ardl(UNEMP ~ GDPGROWTH + HDI + CPI + FDI, data = df, order = c(2, 1, 0, 2, 0))   # cross-check
print(bounds_f_test(ardl_model, case = 3))


## ---- STEP 10: Tables 9-10 - Long-run and short-run estimates (HAC) ----------
lr_table <- function(model, V, map) {
  b <- coef(model); dfr <- df.residual(model)
  do.call(rbind, lapply(names(map), function(v) {
    num <- map[[v]]; den <- "L1_UNEMP"; est <- -b[[num]] / b[[den]]
    g <- setNames(rep(0, length(b)), names(b)); g[num] <- -1 / b[[den]]; g[den] <- b[[num]] / b[[den]]^2
    se <- sqrt(as.numeric(t(g) %*% V %*% g)); tt <- est / se; p <- 2 * pt(abs(tt), dfr, lower.tail = FALSE)
    data.frame(Variable = v, Coef = round(est, 4), SE = round(se, 4), t = round(tt, 3), p = round(p, 4), Sig = stars(p))
  }))
}
lr_map <- c(GDPGROWTH = "L1_G", HDI = "HDI", CPI = "L1_CPI", FDI = "FDI", C = "(Intercept)")
long_run <- lr_table(uecm_model, nw_vcov, lr_map)
cat("\n===== Table 9: Long-run coefficients (HAC SEs) =====\n"); print(long_run, row.names = FALSE)
# Standardized beta = coefficient x SD(X) / SD(UNEMP), SDs from the original data (Table 3)
sd_y <- sd(df_original$UNEMP)
long_run$Std_beta <- c(round(long_run$Coef[1:4] * sapply(c("GDPGROWTH", "HDI", "CPI", "FDI"), function(v) sd(df_original[[v]])) / sd_y, 3), NA)
long_run$Effect_1SD_pp <- c(round(long_run$Coef[1:4] * sapply(c("GDPGROWTH", "HDI", "CPI", "FDI"), function(v) sd(df_original[[v]])), 3), NA)
cat("Standardized coefficients (ranking of effect size):\n"); print(long_run[1:4, c("Variable", "Std_beta", "Effect_1SD_pp")], row.names = FALSE)

sr_all <- coef_table(uecm_model, nw_vcov)
sr_map <- c(L1_d_UNEMP = "D(UNEMP(-1))", d_G = "D(GDPGROWTH)", HDI = "D(HDI)", d_CPI = "D(CPI)",
            L1_d_CPI = "D(CPI(-1))", FDI = "D(FDI)", L1_UNEMP = "ECT(-1)")
short_run <- cbind(Variable = unname(sr_map), sr_all[names(sr_map), ])
cat("\n===== Table 10: Short-run coefficients and ECT (HAC SEs) =====\n"); print(short_run, row.names = FALSE)
cat("Note: for HDI and FDI (zero lags) the short-run impact equals the coefficient on HDI_t and FDI_t.\n")
ect <- coef(uecm_model)[["L1_UNEMP"]]
cat("Speed of adjustment:", round(-100 * ect, 1), "% per year\n")
print(multipliers(ardl_model, type = "lr"))                       # cross-check (coefficients only)


## ---- STEP 11: Table 11 - Diagnostic tests ----------------------------------
u <- residuals(uecm_model); nU <- length(u)
# Levels form of the same ARDL(2,1,0,2,0) (identical residuals) for BG and RESET
lev <- data.frame(y = df$UNEMP, U1 = lag1(df$UNEMP), U2 = c(NA, NA, head(df$UNEMP, -2)), G0 = df$GDPGROWTH, G1 = lag1(df$GDPGROWTH),
                  H0 = df$HDI, C0 = df$CPI, C1 = lag1(df$CPI), C2 = c(NA, NA, head(df$CPI, -2)), F0 = df$FDI)
lev <- na.omit(lev); lev_m <- lm(y ~ ., data = lev); Xl <- model.matrix(lev_m); e <- residuals(lev_m)
e1 <- c(0, head(e, -1)); e2 <- c(0, 0, head(e, -2))                       # zero-padded lags (as in statsmodels)
bg_F <- anova(lm(e ~ Xl - 1), lm(e ~ Xl - 1 + e1 + e2))
bp_LM <- nU * summary(lm(u^2 ~ Xl[, -1]))$r.squared
arch_lm <- function(e, lags) { e2 <- e^2; m <- length(e2); Z <- sapply(1:lags, function(i) e2[(lags + 1 - i):(m - i)])
  LM <- (m - lags) * summary(lm(e2[(lags + 1):m] ~ Z))$r.squared; c(LM, pchisq(LM, lags, lower.tail = FALSE)) }
sk <- mean((u - mean(u))^3) / mean((u - mean(u))^2)^1.5; ku <- mean((u - mean(u))^4) / mean((u - mean(u))^2)^2
jb <- nU / 6 * (sk^2 + (ku - 3)^2 / 4)
lev$yhat2 <- fitted(lev_m)^2; reset <- anova(lev_m, lm(y ~ ., data = lev))
kp <- length(coef(uecm_model))
yl <- df$UNEMP[(n - nU + 1):n]; R2 <- 1 - sum(u^2) / sum((yl - mean(yl))^2)
diag_tab <- data.frame(
  Test = c("R2 (levels)", "Adj. R2", "Breusch-Godfrey F (2 lags)", "Breusch-Pagan LM", "ARCH LM (1 lag)", "Jarque-Bera", "Ramsey RESET F"),
  Statistic = round(c(R2, 1 - (1 - R2) * (nU - 1) / (nU - kp), bg_F$F[2], bp_LM, arch_lm(u, 1)[1], jb, reset$F[2]), 3),
  p = round(c(NA, NA, bg_F$`Pr(>F)`[2], pchisq(bp_LM, kp - 1, lower.tail = FALSE), arch_lm(u, 1)[2], pchisq(jb, 2, lower.tail = FALSE), reset$`Pr(>F)`[2]), 3))
cat("\n===== Table 11: Diagnostic tests =====\n"); print(diag_tab, row.names = FALSE)
cat("BG and RESET use the levels form of the same model (identical residuals).\n")


## ---- STEP 12: Figures 3-4 - CUSUM and CUSUM of squares ---------------------
f_uecm <- formula(uecm_model)
cusum <- efp(f_uecm, data = uecm_df, type = "Rec-CUSUM"); print(sctest(cusum)); plot(cusum, main = "Figure 3  CUSUM")
rr <- recresid(f_uecm, data = uecm_df); cs <- cumsum(rr^2) / sum(rr^2); m <- length(rr); expct <- seq_len(m) / m
c0 <- 1.358 / (sqrt(m) + 0.12 + 0.11 / sqrt(m))
plot(cs, type = "l", lwd = 2, col = "steelblue", ylim = c(-0.4, 1.4), main = "Figure 4  CUSUM of squares", xlab = "Recursive residual", ylab = "CUSUMSQ")
lines(expct - c0, col = "red", lty = 2); lines(expct + c0, col = "red", lty = 2)
cat("CUSUMSQ max deviation =", round(max(abs(cs - expct)), 4), "vs approx. 5% bound", round(c0, 4), "\n")


## ---- STEP 13: Table 12 - Johansen cointegration test -----------------------
joh <- johansen_test(df[, vars5], k_ar_diff = 1)
table11 <- data.frame(Hypothesis = c("r = 0", "r <= 1", "r <= 2", "r <= 3", "r <= 4"),
  Trace = round(joh$trace, 3), Trace_CV5 = c(69.819, 47.854, 29.796, 15.494, 3.842),
  MaxEigen = round(joh$maxeig, 3), MaxEig_CV5 = c(33.878, 27.586, 21.131, 14.264, 3.842))
cat("\n===== Table 12: Johansen test =====\n"); print(table11, row.names = FALSE)


## ---- STEP 14: Table 13 - VECM Granger causality ----------------------------
gr <- vecm_granger_matrix(df[, vars5], k_ar_diff = 1)
table12 <- matrix(sprintf("%.3f%s (%.3f)", gr$F, stars(gr$p), gr$p), 5, 5, dimnames = dimnames(gr$F)); diag(table12) <- "-"
cat("\n===== Table 13: VECM Granger causality, F (p) | rows = caused, columns = causing =====\n"); print(noquote(table12))
print(vecm_ect(df[, vars5], joh$beta), row.names = FALSE)


## ---- STEP 15: Table 14 - Robustness checks ---------------------------------
if (has_gdppc) {
  rob_df <- make_uecm(df, "GDPPC")
  rob_m <- lm(d_UNEMP ~ L1_UNEMP + L1_G + HDI + L1_CPI + FDI + L1_d_UNEMP + d_G + d_CPI + L1_d_CPI, data = rob_df)
  rob_F <- anova(lm(d_UNEMP ~ L1_d_UNEMP + d_G + d_CPI + L1_d_CPI, data = rob_df), rob_m)$F[2]
  cat("\n===== Table 14a: GDP per capita growth instead of real GDP growth =====\n")
  cat("Bounds F =", round(rob_F, 3), "| ECT =", round(coef(rob_m)[["L1_UNEMP"]], 4), "\n")
  print(lr_table(rob_m, newey_west(rob_m, 2), c(GDPPC = "L1_G", HDI = "HDI", CPI = "L1_CPI", FDI = "FDI")), row.names = FALSE)
}
old_df <- na.omit(data.frame(d_UNEMP = dif(df$UNEMP), L1_UNEMP = lag1(df$UNEMP), L1_G = lag1(df$GDPGROWTH), L1_HDI = lag1(df$HDI),
  L1_CPI = lag1(df$CPI), L1_FDI = lag1(df$FDI), d_G = dif(df$GDPGROWTH), d_HDI = dif(df$HDI), d_CPI = dif(df$CPI),
  L1_d_CPI = lag1(dif(df$CPI)), d_FDI = dif(df$FDI), L1_d_FDI = lag1(dif(df$FDI))))
old_m <- lm(d_UNEMP ~ L1_UNEMP + L1_G + L1_HDI + L1_CPI + L1_FDI + d_G + d_HDI + d_CPI + L1_d_CPI + d_FDI + L1_d_FDI, data = old_df)
old_F <- anova(lm(d_UNEMP ~ d_G + d_HDI + d_CPI + L1_d_CPI + d_FDI + L1_d_FDI, data = old_df), old_m)$F[2]
cat("\n===== Table 14b: alternative lag structure ARDL(1,1,1,2,2) =====\n")
cat("Bounds F =", round(old_F, 3), "| ECT =", round(coef(old_m)[["L1_UNEMP"]], 4), "\n")
print(lr_table(old_m, newey_west(old_m, 2), c(GDPGROWTH = "L1_G", HDI = "L1_HDI", CPI = "L1_CPI", FDI = "L1_FDI")), row.names = FALSE)

# Sensitivity: exclude 2024, whose HDI value is an estimate (UNDP data end in 2023)
sub_df <- make_uecm(df[df_original$year <= 2023, ], "GDPGROWTH")
sub_m <- lm(d_UNEMP ~ L1_UNEMP + L1_G + HDI + L1_CPI + FDI + L1_d_UNEMP + d_G + d_CPI + L1_d_CPI, data = sub_df)
sub_F <- anova(lm(d_UNEMP ~ L1_d_UNEMP + d_G + d_CPI + L1_d_CPI, data = sub_df), sub_m)$F[2]
cat("\n===== Sensitivity: sample 1994-2023 (without the estimated 2024 HDI) =====\n")
cat("N =", nobs(sub_m), "| Bounds F =", round(sub_F, 3), "| ECT =", round(coef(sub_m)[["L1_UNEMP"]], 4), "\n")
print(lr_table(sub_m, newey_west(sub_m, 2), c(GDPGROWTH = "L1_G", HDI = "HDI", CPI = "L1_CPI", FDI = "FDI")), row.names = FALSE)

cat("\nDone. Compare each table with the same table number in Final_Research_Report.docx.\n")
