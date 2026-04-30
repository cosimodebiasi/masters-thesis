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