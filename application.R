############################################################
# Testing Zero Restrictions in the Mixing Matrix (ICA - VAR)
# Empirical Application
############################################################

rm(list = ls())
setwd("~/Desktop/Tesi/application")

# ---------------------------------------------------------
# 1. Setup
# ---------------------------------------------------------

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(showtext)
library(tseries)
library(vars)
library(fastICA)

source("functions.R")

# ---------------------------------------------------------
# 2. Import Data
# ---------------------------------------------------------

uncertainty_data <- read_excel("Categorical_EPU_Data.xlsx") %>%
  filter(!is.na(Year), !is.na(Month)) %>%
  transmute(
    Date = as.Date(sprintf("%d-%02d-01", as.numeric(Year), as.numeric(Month))),
    EPU_General = `1. Economic Policy Uncertainty`
  )

indpro_data <- read_excel("INDPRO.xlsx") %>%
  transmute(
    Date = as.Date(observation_date),
    INDPRO
  )

cpi_data <- read_csv("CPIAUCSL.csv", show_col_types = FALSE) %>%
  transmute(
    Date = as.Date(observation_date),
    CPI = CPIAUCSL
  )

# ---------------------------------------------------------
# 3. Merge and Transformations
# ---------------------------------------------------------

svar_data <- uncertainty_data %>%
  inner_join(indpro_data, by = "Date") %>%
  inner_join(cpi_data, by = "Date") %>%
  arrange(Date) %>%
  mutate(
    l_epu     = log(EPU_General),
    ip_growth = 100 * (log(INDPRO) - log(lag(INDPRO))),
    inflation = 100 * (log(CPI) - log(lag(CPI)))
  ) %>%
  dplyr::select(Date, l_epu, ip_growth, inflation) %>%
  filter(Date >= as.Date("1985-01-01")) %>%
  drop_na()

summary(svar_data)
colSums(is.na(svar_data))

vars_to_check <- c("l_epu", "ip_growth", "inflation")

sapply(svar_data[vars_to_check], sd)

invisible(lapply(vars_to_check, function(v) {
  plot(svar_data$Date, svar_data[[v]], type = "l", main = v)
}))

# ---------------------------------------------------------
# 4. Plot
# ---------------------------------------------------------

font_add("cmu", regular = "cmunrm.ttf", italic = "cmunti.ttf")
showtext_auto()

plot_labels <- c(
  l_epu = "Economic policy uncertainty (log)",
  ip_growth = "Industrial production growth",
  inflation = "Inflation"
)

plot_data <- svar_data %>%
  pivot_longer(
    cols = -Date,
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = recode(variable, !!!plot_labels),
    variable = factor(
      variable,
      levels = c(
        "Industrial production growth",
        "Inflation",
        "Economic policy uncertainty (log)"
      )
    )
  )

ts_plot <- ggplot(plot_data, aes(x = Date, y = value)) +
  geom_line(linewidth = 0.3, alpha = 0.9) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  theme_minimal(base_family = "cmu") +
  labs(x = NULL, y = NULL) +
  theme(
    strip.text = element_text(size = 11),
    axis.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "grey85"),
    panel.background = element_blank(),
    panel.spacing = unit(1, "lines")
  )

ts_plot

ggsave("timeseries_plot.png", plot = ts_plot, width = 7, height = 6, dpi = 300)

# ---------------------------------------------------------
# 5. Stationarity tests
# ---------------------------------------------------------

adf_results <- lapply(svar_data[vars_to_check], adf.test)
adf_results

# ---------------------------------------------------------
# 6. VAR Estimation and Diagnostics
# ---------------------------------------------------------

var_data <- svar_data %>%
  dplyr::select(ip_growth, inflation, l_epu)

lag_selection <- VARselect(var_data, lag.max = 12, type = "const")
print(lag_selection$selection)

model_var <- VAR(var_data, p = 2, type = "const")
serial.test(model_var, lags.pt = 16, type = "PT.asymptotic")

# Residual autocorrelation -> try p = 3
model_var3 <- VAR(var_data, p = 3, type = "const")
serial.test(model_var3, lags.pt = 16, type = "PT.asymptotic")

# Residual normality test
normality.test(model_var, multivariate.only = FALSE)

# ---------------------------------------------------------
# 7. ICA
# ---------------------------------------------------------

ures <- resid(model_var)
Y <- as.matrix(var_data)

T <- nrow(ures)
K <- ncol(ures)
p <- 2

AA <- Acoef(model_var)

const <- sapply(seq_len(K), function(j) {
  coef(model_var)[[j]]["const", 1]
})

A <- fAp_fastICA(ures, sseed = 123)
round(A, 3)

# ---------------------------------------------------------
# 8. Bootstrap
# ---------------------------------------------------------

N <- 1000

Aboot <- array(NA, dim = c(K, K, N))
vecAboot <- matrix(NA, nrow = K^2, ncol = N)

h <- 13
IRFboot <- array(NA, dim = c(K, K, h + 1, N))

for (b in 1:N) {
  
  unew <- ures[sample(T, replace = TRUE), ]
  
  Ynew <- matrix(0, nrow = (T + p), ncol = K)
  Ynew[1:p, ] <- Y[1:p, ]
  
  for (i in (p + 1):(T + p)) {
    for (j in 1:p) {
      Ynew[i, ] <- Ynew[i, ] + AA[[j]] %*% Ynew[i - j, ]
    }
    Ynew[i, ] <- const + Ynew[i, ] + unew[i - p, ]
  }
  
  Ynew <- as.data.frame(Ynew)
  
  varest_new <- vars::VAR(Ynew, p = p, type = "const")
  ures_star <- resid(varest_new)
  
  icares <- fastICA(ures_star, K, tol = 1e-14, maxit = 3000, verbose = FALSE)
  W <- t((icares$K) %*% (icares$W))
  A_star <- solve(W)
  
  P <- myfrob(A_star, A, K)
  A_star <- A_star %*% P
  
  IRFboot[, , , b] <- compute_struct_irf(varest_new, A_star, h = h)
  
  Aboot[, , b] <- A_star
  vecAboot[, b] <- as.vector(A_star)
}

vecA <- as.vector(A)

#----------------------------------------------------------
# Bootstrap confidence bands for IRFs
#----------------------------------------------------------

lower_irf <- apply(IRFboot, c(1, 2, 3), quantile, probs = 0.05, na.rm = TRUE)
upper_irf <- apply(IRFboot, c(1, 2, 3), quantile, probs = 0.95, na.rm = TRUE)

#----------------------------------------------------------
# 8. Bootstrap Inference
#----------------------------------------------------------

alpha <- 0.1
ncoef <- K^2

SD    <- apply(vecAboot, 1, sd)
tstat <- vecA / SD

df <- K^2 * (p + 1) + K

tvalue <- qt(1 - alpha / 2, (T + p - df - 1))

# --------------- P-VALUES ----------------
pvalues <- 2 * (1 - pt(abs(tstat), df = (T + p - df - 1)))
pval_matrix <- matrix(pvalues, nrow = K, ncol = K)
# -----------------------------------------

mcoef <- matrix("*", K, K)

for (i in 1:ncoef) {
  if (abs(tstat[i]) <= tvalue) {
    mcoef[i] <- "0"
  }
}

#----------------------------------------------------------
# 9. Output
#----------------------------------------------------------

round(A, 3)
print(mcoef)
round(pval_matrix, 3)

############### IMPULSE RESPONSE ANALYSIS #################

#----------------------------------------------------------
# Structural IRFs from the estimated VAR and ICA matrix A
#----------------------------------------------------------

h <- 13

IRF <- compute_struct_irf(model_var, A, h = h)

shock_labels <- c(
  s1 = "Industrial production shock",
  s2 = "Inflation shock",
  s3 = "Uncertainty shock"
)

response_labels <- c(
  ip_growth = "Industrial production growth",
  inflation = "Inflation",
  l_epu = "Economic policy uncertainty"
)

#----------------------------------------------------------
# IRF data for plotting
#----------------------------------------------------------

irf_band_df <- as.data.frame.table(IRF, responseName = "irf")
names(irf_band_df) <- c("response", "shock", "horizon", "irf")

lower_df <- as.data.frame.table(lower_irf, responseName = "lower")
upper_df <- as.data.frame.table(upper_irf, responseName = "upper")

irf_band_df$lower <- lower_df$lower
irf_band_df$upper <- upper_df$upper

irf_band_df$horizon <- as.numeric(as.character(irf_band_df$horizon))

irf_band_df$response <- factor(
  response_labels[irf_band_df$response],
  levels = response_labels
)

irf_band_df$shock <- factor(
  shock_labels[irf_band_df$shock],
  levels = shock_labels
)

#----------------------------------------------------------
# Plot of structural IRFs
#----------------------------------------------------------

# Fixed scale

irf_band_plot <- ggplot(irf_band_df, aes(x = horizon, y = irf)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 0.4) +
  facet_grid(response ~ shock, scales = "free_y") +
  labs(
    x = "Months",
    y = "Response"
  ) +
  theme_minimal(base_family = "cmu") +
  theme(
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

irf_band_plot

ggsave("irf_structural_ica_bands.png", plot = irf_band_plot, width = 8, height = 6, dpi = 300)

# Free scale

irf_band_plot_free <- ggplot(irf_band_df, aes(x = horizon, y = irf)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 0.4) +
  facet_grid(response ~ shock, scales = "free_y") +
  labs(
    x = "Months",
    y = "Response"
  ) +
  theme_minimal(base_family = "cmu") +
  theme(
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

irf_band_plot_free

ggsave("irf_structural_ica_bands_free.png", plot = irf_band_plot_free, width = 8, height = 6, dpi = 300)

#----------------------------------------------------------
# Bootstrap distribution of each coefficient
#----------------------------------------------------------

library(ggtext)

font_add(
  "cmu",
  regular = "cmunrm.ttf",
  italic  = "cmunti.ttf"
)

showtext_auto()

boot_long <- data.frame(
  value = as.vector(vecAboot),
  coef_id = rep(seq_len(K^2), times = N)
) %>%
  mutate(
    row = ((coef_id - 1) %% K) + 1,
    col = ((coef_id - 1) %/% K) + 1,
    label = paste0("<i>a</i><sub style='font-size:6pt'>", row, col, "</sub>"),
    Ahat = as.vector(A)[coef_id]
  )

boot_dist_plot <- ggplot(boot_long, aes(x = value)) +
  geom_density(fill = "grey80", alpha = 0.8, linewidth = 0.3) +
  geom_vline(aes(xintercept = Ahat), linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ label, scales = "free", ncol = K) +
  theme_minimal(base_family = "cmu") +
  labs(
    x = "Bootstrap coefficient value",
    y = "Density"
  ) +
  theme(
    strip.text = element_markdown(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    panel.grid.major = element_line(linewidth = 0.15, color = "grey85"),
    panel.grid.minor = element_line(linewidth = 0.08, color = "grey90")
  )

boot_dist_plot

ggsave(
  "bootstrap_distributions_A_coefficients.png",
  plot = boot_dist_plot,
  width = 8,
  height = 6,
  dpi = 300
)


########## Financial Market Uncetainty ############

# Stock market volatility index (VIX

vix_data <- read_csv("VIXCLS.csv", show_col_types = FALSE) %>%
  transmute(
    Date = as.Date(observation_date),
    VIX = VIXCLS
  ) %>%
  mutate(
    Month = format(Date, "%Y-%m")
  ) %>%
  group_by(Month) %>%
  summarise(
    VIX = mean(VIX, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Date = as.Date(paste0(Month, "-01"))
  ) %>%
  dplyr::select(Date, VIX)

# New dataset

financial_data <- indpro_data %>%
  inner_join(cpi_data, by = "Date") %>%
  inner_join(vix_data, by = "Date") %>%
  arrange(Date) %>%
  mutate(
    ip_growth = 100 * (log(INDPRO) - log(lag(INDPRO))),
    inflation = 100 * (log(CPI) - log(lag(CPI))),
    l_vix = log(VIX)
  ) %>%
  dplyr::select(Date, ip_growth, inflation, l_vix) %>%
  filter(
    Date >= as.Date("1990-01-01"),
    Date <= as.Date("2025-12-01")
  ) %>%
  drop_na()

# VAR

var_data_fin <- financial_data %>%
  dplyr::select(ip_growth, inflation, l_vix)

lag_selection_fin <- VARselect(
  var_data_fin,
  lag.max = 12,
  type = "const"
)

print(lag_selection_fin$selection)

# p = 2

model_var_fin <- VAR(
  var_data_fin,
  p = 2,
  type = "const"
)

serial.test(
  model_var_fin,
  lags.pt = 16,
  type = "PT.asymptotic"
)

normality.test(
  model_var_fin,
  multivariate.only = FALSE
)

# ICA

ures_fin <- resid(model_var_fin)

A_fin <- fAp_fastICA(
  ures_fin,
  sseed = 123
)

round(A_fin, 3)

# 8. Bootstrap

T_fin <- nrow(ures_fin)
K_fin <- ncol(ures_fin)
p_fin <- 2

Y_fin <- as.matrix(var_data_fin)

AA_fin <- Acoef(model_var_fin)

const_fin <- sapply(seq_len(K_fin), function(j) {
  coef(model_var_fin)[[j]]["const", 1]
})

N <- 1000

Aboot_fin <- array(
  NA,
  dim = c(K_fin, K_fin, N)
)

vecAboot_fin <- matrix(
  NA,
  nrow = K_fin^2,
  ncol = N
)

h <- 13

IRFboot_fin <- array(
  NA,
  dim = c(K_fin, K_fin, h + 1, N)
)

for (b in 1:N) {
  
  unew_fin <- ures_fin[sample(T_fin, replace = TRUE), ]
  
  Ynew_fin <- matrix(
    0,
    nrow = (T_fin + p_fin),
    ncol = K_fin
  )
  
  Ynew_fin[1:p_fin, ] <- Y_fin[1:p_fin, ]
  
  for (i in (p_fin + 1):(T_fin + p_fin)) {
    for (j in 1:p_fin) {
      Ynew_fin[i, ] <- Ynew_fin[i, ] +
        AA_fin[[j]] %*% Ynew_fin[i - j, ]
    }
    
    Ynew_fin[i, ] <- const_fin +
      Ynew_fin[i, ] +
      unew_fin[i - p_fin, ]
  }
  
  Ynew_fin <- as.data.frame(Ynew_fin)
  
  varest_new_fin <- vars::VAR(
    Ynew_fin,
    p = p_fin,
    type = "const"
  )
  
  ures_star_fin <- resid(varest_new_fin)
  
  icares_fin <- fastICA(
    ures_star_fin,
    K_fin,
    tol = 1e-14,
    maxit = 3000,
    verbose = FALSE
  )
  
  W_fin <- t(
    (icares_fin$K) %*% (icares_fin$W)
  )
  
  A_star_fin <- solve(W_fin)
  
  P_fin <- myfrob(
    A_star_fin,
    A_fin,
    K_fin
  )
  
  A_star_fin <- A_star_fin %*% P_fin
  
  IRFboot_fin[, , , b] <- compute_struct_irf(
    varest_new_fin,
    A_star_fin,
    h = h
  )
  
  Aboot_fin[, , b] <- A_star_fin
  
  vecAboot_fin[, b] <- as.vector(A_star_fin)
}

vecA_fin <- as.vector(A_fin)


# Bootstrap confidence bands for IRFs

lower_irf_fin <- apply(
  IRFboot_fin,
  c(1, 2, 3),
  quantile,
  probs = 0.05,
  na.rm = TRUE
)

upper_irf_fin <- apply(
  IRFboot_fin,
  c(1, 2, 3),
  quantile,
  probs = 0.95,
  na.rm = TRUE
)


# 8. Bootstrap Inference

alpha <- 0.1

ncoef_fin <- K_fin^2

SD_fin <- apply(
  vecAboot_fin,
  1,
  sd
)

tstat_fin <- vecA_fin / SD_fin

df_fin <- K_fin^2 * (p_fin + 1) + K_fin

tvalue_fin <- qt(
  1 - alpha / 2,
  (T_fin + p_fin - df_fin - 1)
)


# --------------- P-VALUES ----------------

pvalues_fin <- 2 * (
  1 - pt(
    abs(tstat_fin),
    df = (T_fin + p_fin - df_fin - 1)
  )
)

pval_matrix_fin <- matrix(
  pvalues_fin,
  nrow = K_fin,
  ncol = K_fin
)

# -----------------------------------------

mcoef_fin <- matrix(
  "*",
  K_fin,
  K_fin
)

for (i in 1:ncoef_fin) {
  if (abs(tstat_fin[i]) <= tvalue_fin) {
    mcoef_fin[i] <- "0"
  }
}


# 9. Output

round(A_fin, 3)

print(mcoef_fin)

round(pval_matrix_fin, 3)


############### IMPULSE RESPONSE ANALYSIS #################

# Structural IRFs from the estimated VAR and ICA matrix A

h <- 13

IRF_fin <- compute_struct_irf(
  model_var_fin,
  A_fin,
  h = h
)

shock_labels_fin <- c(
  s1 = "Industrial production shock",
  s2 = "Inflation shock",
  s3 = "Financial uncertainty shock"
)

response_labels_fin <- c(
  ip_growth = "Industrial production growth",
  inflation = "Inflation",
  l_vix = "Stock market volatility (log VIX)"
)


#----------------------------------------------------------
# IRF data for plotting
#----------------------------------------------------------

irf_band_df_fin <- as.data.frame.table(
  IRF_fin,
  responseName = "irf"
)

names(irf_band_df_fin) <- c(
  "response",
  "shock",
  "horizon",
  "irf"
)

lower_df_fin <- as.data.frame.table(
  lower_irf_fin,
  responseName = "lower"
)

upper_df_fin <- as.data.frame.table(
  upper_irf_fin,
  responseName = "upper"
)

irf_band_df_fin$lower <- lower_df_fin$lower

irf_band_df_fin$upper <- upper_df_fin$upper

irf_band_df_fin$horizon <- as.numeric(
  as.character(irf_band_df_fin$horizon)
)

irf_band_df_fin$response <- factor(
  response_labels_fin[irf_band_df_fin$response],
  levels = response_labels_fin
)

irf_band_df_fin$shock <- factor(
  shock_labels_fin[irf_band_df_fin$shock],
  levels = shock_labels_fin
)


# Plot of structural IRFs

# Fixed scale

irf_band_plot_fin <- ggplot(
  irf_band_df_fin,
  aes(x = horizon, y = irf)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2
  ) +
  geom_line(
    linewidth = 0.4
  ) +
  facet_grid(
    response ~ shock,
    scales = "free_y"
  ) +
  labs(
    x = "Months",
    y = "Response"
  ) +
  theme_minimal(base_family = "cmu") +
  theme(
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

irf_band_plot_fin

ggsave(
  "irf_structural_ica_bands_financial.png",
  plot = irf_band_plot_fin,
  width = 8,
  height = 6,
  dpi = 300
)


# Free scale

irf_band_plot_free_fin <- ggplot(
  irf_band_df_fin,
  aes(x = horizon, y = irf)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2
  ) +
  geom_line(
    linewidth = 0.4
  ) +
  facet_grid(
    response ~ shock,
    scales = "free_y"
  ) +
  labs(
    x = "Months",
    y = "Response"
  ) +
  theme_minimal(base_family = "cmu") +
  theme(
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

irf_band_plot_free_fin

ggsave(
  "irf_structural_ica_bands_financial_free.png",
  plot = irf_band_plot_free_fin,
  width = 8,
  height = 6,
  dpi = 300
)


# Bootstrap distribution of each coefficient

library(ggtext)

font_add(
  "cmu",
  regular = "cmunrm.ttf",
  italic = "cmunti.ttf"
)

showtext_auto()

boot_long_fin <- data.frame(
  value = as.vector(vecAboot_fin),
  coef_id = rep(
    seq_len(K_fin^2),
    times = N
  )
) %>%
  mutate(
    row = ((coef_id - 1) %% K_fin) + 1,
    col = ((coef_id - 1) %/% K_fin) + 1,
    label = paste0(
      "<i>a</i><sub style='font-size:6pt'>",
      row,
      col,
      "</sub>"
    ),
    Ahat = as.vector(A_fin)[coef_id]
  )

boot_dist_plot_fin <- ggplot(
  boot_long_fin,
  aes(x = value)
) +
  geom_density(
    fill = "grey80",
    alpha = 0.8,
    linewidth = 0.3
  ) +
  geom_vline(
    aes(xintercept = Ahat),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  facet_wrap(
    ~ label,
    scales = "free",
    ncol = K_fin
  ) +
  theme_minimal(base_family = "cmu") +
  labs(
    x = "Bootstrap coefficient value",
    y = "Density"
  ) +
  theme(
    strip.text = element_markdown(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    panel.grid.major = element_line(
      linewidth = 0.15,
      color = "grey85"
    ),
    panel.grid.minor = element_line(
      linewidth = 0.08,
      color = "grey90"
    )
  )

boot_dist_plot_fin

ggsave(
  "bootstrap_distributions_A_coefficients_financial.png",
  plot = boot_dist_plot_fin,
  width = 8,
  height = 6,
  dpi = 300
)
