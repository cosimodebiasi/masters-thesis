############################################################
# Testing Zero Restrictions in the Mixing Matrix (ICA - VAR)
# Simulation with K = 2
############################################################

rm(list = ls())

#----------------------------------------------------------
# 1. Setup
#----------------------------------------------------------

setwd("~/Desktop/Tesi/simulation")

# Required packages
library(vars)
library(dse)
library(rmutil)
library(fastICA)
library(gtools)
library(matrixcalc)

# Load custom functions
source("functions.R")

#----------------------------------------------------------
# 2. Simulate VAR(2) Process
#----------------------------------------------------------

p <- 2   # number of lags
K <- 2   # number of variables

# Lag polynomial A(L)
AL <- array(
  c(1, -0.5, 0.3,
    0,  0.2, 0.1,
    0, -0.2, 0.7,
    1,  0.5, -0.3),
  c(3, 2, 2)
)

# Mixing matrix (recursive structure)
B0inv <- diag(2)
B0inv[1, 2] <- 0.5

# Constant term
c <- c(1, 2)

# Define VAR(2)
var2 <- ARMA(A = AL, B = B0inv, TREND = c)

# Simulation settings
h <- 10
T  <- 1000 + h

set.seed(123)

# Generate shocks
w <- matrix(rlaplace(T * K, s = sqrt(1/2)), nrow = T, ncol = K)

# Simulate data
varsim <- dse::simulate(
  var2,
  sampleT = T,
  noise = list(w = w),
  rng = list(seed = c(46))
)

vardat <- matrix(varsim$output, nrow = T, ncol = K)
colnames(vardat) <- c("y1", "y2")

# Remove burn-in
ds <- vardat[(h + 1):T, ]

save(ds, file = "ds.Rdata")

# Plot simulated series
plot.ts(ds, main = "Simulated series", xlab = "Time")

#----------------------------------------------------------
# 3. Estimate Reduced-Form VAR
#----------------------------------------------------------

Y <- as.matrix(ds)

varest <- vars::VAR(Y, p = p, type = "const")

AA <- Acoef(varest)
ures <- resid(varest)

T <- nrow(ures)
K <- ncol(ures)

# Estimated constant
const <- sapply(1:K, function(eq) coef(varest)[[eq]]["const", 1])

# Diagnostics
plot.ts(ures)
normality.test(varest, multivariate.only = FALSE)

#----------------------------------------------------------
# 4. Apply ICA to Residuals
#----------------------------------------------------------

A <- fAp_fastICA(ures, sseed = 123)

#----------------------------------------------------------
# 5. Bootstrap
#----------------------------------------------------------

N <- 1000

Aboot <- array(NA, dim = c(K, K, N))
vecAboot <- matrix(NA, nrow = K^2, ncol = N)

for (b in 1:N) {
  
  # Resample residuals
  unew <- ures[sample(T, replace = TRUE), ]
  
  # Reconstruct time series
  Ynew <- matrix(0, nrow = (T + p), ncol = K)
  Ynew[1:p, ] <- Y[1:p, ]
  
  for (i in (p + 1):(T + p)) {
    for (j in 1:p) {
      Ynew[i, ] <- Ynew[i, ] + AA[[j]] %*% Ynew[i - j, ]
    }
    Ynew[i, ] <- const + Ynew[i, ] + unew[i - p, ]
  }
  
  Ynew <- as.data.frame(Ynew)
  
  # Re-estimate VAR
  varest_new <- vars::VAR(Ynew, p = p, type = "const")
  ures_star <- resid(varest_new)
  
  # ICA estimation
  icares <- fastICA(ures_star, K, tol = 1e-14, maxit = 3000, verbose = FALSE)
  W <- t((icares$K) %*% (icares$W))
  A_star <- solve(W)
  
  # Align permutation and signs
  P <- myfrob(A_star, A, K)
  A_star <- A_star %*% P
  
  # Store results
  Aboot[, , b] <- A_star
  vecAboot[, b] <- as.vector(A_star)
}

vecA <- as.vector(A)

#----------------------------------------------------------
# 6. Bootstrap Inference
#----------------------------------------------------------

alpha <- 0.01
ncoef <- K^2

SD    <- apply(vecAboot, 1, sd)
tstat <- vecA / SD

# Degrees of freedom
df <- K^2 * (p + 1) + K

tvalue <- qt(1 - alpha / 2, (T + p - df - 1))

# Significance matrix
mcoef <- matrix("*", K, K)

for (i in 1:ncoef) {
  if (abs(tstat[i]) <= tvalue) {
    mcoef[i] <- "0"
  }
}

# Output results
print(B0inv)
round(A, 3)
print(mcoef)

############################################################
# Testing Zero Restrictions in the Mixing Matrix (ICA - VAR)
# Simulation with K = 3
############################################################

rm(list = ls())

#----------------------------------------------------------
# 1. Setup
#----------------------------------------------------------

setwd("~/Desktop/Tesi/simulation")

library(vars)
library(dse)
library(rmutil)
library(fastICA)
library(gtools)
library(matrixcalc)

source("functions.R")

#----------------------------------------------------------
# 2. Simulate VAR(2) Process
#----------------------------------------------------------

p <- 2
K <- 3

# Lag polynomial
AL <- array(0, dim = c(p + 1, K, K))

# Lag 0 (identity)
AL[1, , ] <- diag(K)

# Lag 1
AL[2, , ] <- matrix(c(
  0.5, 0.1, 0.0,
  0.0, 0.3, 0.1,
  0.0, 0.0, 0.4
), K, K, byrow = TRUE)

# Lag 2
AL[3, , ] <- matrix(c(
  -0.2, 0.0, 0.0,
  0.0,-0.1, 0.0,
  0.0, 0.0,-0.2
), K, K, byrow = TRUE)

# Mixing matrix (recursive structure)
B0inv <- diag(K)
B0inv[1, 2] <- 0.5
B0inv[1, 3] <- 0.3
B0inv[2, 3] <- 0.4

# Constant term
c <- c(1, 2, 1.5)

# VAR model
var2 <- ARMA(A = AL, B = B0inv, TREND = c)

# Simulation settings
h <- 10
T  <- 1000 + h

set.seed(123)

w <- matrix(rlaplace(T * K, s = sqrt(1/2)), nrow = T, ncol = K)

varsim <- dse::simulate(
  var2,
  sampleT = T,
  noise = list(w = w),
  rng = list(seed = c(46))
)

vardat <- matrix(varsim$output, nrow = T, ncol = K)
colnames(vardat) <- c("y1", "y2", "y3")

# Remove burn-in
ds <- vardat[(h + 1):T, ]

#----------------------------------------------------------
# 3. Estimate Reduced-Form VAR
#----------------------------------------------------------

Y <- as.matrix(ds)

varest <- vars::VAR(Y, p = p, type = "const")

AA <- Acoef(varest)
ures <- resid(varest)

T <- nrow(ures)

# Constant
const <- sapply(1:K, function(eq) coef(varest)[[eq]]["const", 1])

# Diagnostics
normality.test(varest, multivariate.only = FALSE)

#----------------------------------------------------------
# 4. ICA
#----------------------------------------------------------

A <- fAp_fastICA(ures, sseed = 123)

#----------------------------------------------------------
# 5. Bootstrap
#----------------------------------------------------------

N <- 1000

Aboot <- array(NA, dim = c(K, K, N))
vecAboot <- matrix(NA, nrow = K^2, ncol = N)

for (b in 1:N) {
  
  # Resample residuals
  unew <- ures[sample(T, replace = TRUE), ]
  
  # Reconstruct series
  Ynew <- matrix(0, nrow = (T + p), ncol = K)
  Ynew[1:p, ] <- Y[1:p, ]
  
  for (i in (p + 1):(T + p)) {
    for (j in 1:p) {
      Ynew[i, ] <- Ynew[i, ] + AA[[j]] %*% Ynew[i - j, ]
    }
    Ynew[i, ] <- const + Ynew[i, ] + unew[i - p, ]
  }
  
  Ynew <- as.data.frame(Ynew)
  
  # Re-estimate VAR
  varest_new <- vars::VAR(Ynew, p = p, type = "const")
  ures_star <- resid(varest_new)
  
  # ICA
  icares <- fastICA(ures_star, K, tol = 1e-14, maxit = 3000, verbose = FALSE)
  W <- t((icares$K) %*% (icares$W))
  A_star <- solve(W)
  
  # Align (permutation + sign)
  P <- myfrob(A_star, A, K)
  A_star <- A_star %*% P
  
  # Store
  Aboot[, , b] <- A_star
  vecAboot[, b] <- as.vector(A_star)
}

vecA <- as.vector(A)

#----------------------------------------------------------
# 6. Inference
#----------------------------------------------------------

alpha <- 0.01
ncoef <- K^2

SD    <- apply(vecAboot, 1, sd)
tstat <- vecA / SD

# Degrees of freedom
df <- K^2 * (p + 1) + K

tvalue <- qt(1 - alpha / 2, (T + p - df - 1))

# Significance matrix
mcoef <- matrix("*", K, K)

for (i in 1:ncoef) {
  if (abs(tstat[i]) <= tvalue) {
    mcoef[i] <- "0"
  }
}

# Output results
print(B0inv)
print(round(A, 3))
print(mcoef)


########################################################################
# Recursiveness Test via Cholesky Orthogonalisation and Independence
# (Applied to the Trivariate Simulation)
########################################################################

library(gtools)
library(dHSIC)

alpha_rec <- 0.05
B_dhsic <- 499

# Test function
test_recursive_order <- function(ures, ord, alpha = 0.05, B = 499) {
  
  # 1. Reorder reduced-form residuals
  Uord <- as.matrix(ures[, ord, drop = FALSE])
  
  # 2. Covariance matrix of reordered residuals
  Sigma_u <- cov(Uord)
  
  # 3. Cholesky factorisation
  # In R: chol(Sigma_u) = upper triangular matrix R
  # such that Sigma_u = t(R) %*% R
  Rchol <- chol(Sigma_u)
  Pchol <- t(Rchol)
  
  # 4. Orthogonalised shocks
  # ehat_t = P^{-1} u_t
  Ehat <- t(solve(Pchol, t(Uord)))
  
  # 5. dHSIC test of mutual independence
  dh <- dHSIC::dhsic.test(
    X = Ehat,
    alpha = alpha,
    method = "permutation",
    B = B,
    matrix.input = TRUE
  )
  
  return(list(
    order = ord,
    Sigma_u = Sigma_u,
    Pchol = Pchol,
    shocks = Ehat,
    statistic = dh$statistic,
    p.value = dh$p.value,
    reject = (dh$p.value < alpha)
  ))
}


# Loop over all permutations of the variables
perm_mat <- gtools::permutations(n = K, r = K, v = 1:K)
nperm <- nrow(perm_mat)

rec_results <- vector("list", nperm)

for (i in 1:nperm) {
  
  rec_results[[i]] <- test_recursive_order(
    ures = ures,
    ord = perm_mat[i, ],
    alpha = alpha_rec,
    B = B_dhsic
  )
}

# Summary table
rec_summary <- data.frame(
  ordering = apply(perm_mat, 1, function(x) paste(colnames(Y)[x], collapse = " -> ")),
  statistic = sapply(rec_results, function(x) x$statistic),
  p.value = sapply(rec_results, function(x) x$p.value),
  reject_H0_independence = sapply(rec_results, function(x) x$reject)
)

# Multiple testing correction with Bonferroni
alpha_bonf <- alpha_rec / nperm
rec_summary$reject_bonf <- rec_summary$p.value < alpha_bonf

print(rec_summary)

cat("\nBonferroni corrected alpha:", alpha_bonf, "\n")

recursive_supported_raw <- any(!rec_summary$reject_H0_independence)
recursive_supported_bonf <- any(!rec_summary$reject_bonf)

# Decision without multiple testing restriction
if (recursive_supported_raw) {
  cat("Model compatible with recursive ICA\n")
  cat("Orderings where independence is NOT rejected:\n")
  print(rec_summary$ordering[!rec_summary$reject_H0_independence])
} else {
  cat("Model NOT compatible with recursive ICA\n")
}

# Decision with Bonferroni correction
if (recursive_supported_bonf) {
  cat("Model compatible with recursive ICA\n")
  cat("Orderings where independence is NOT rejected after Bonferroni correction:\n")
  print(rec_summary$ordering[!rec_summary$reject_bonf])
} else {
  cat("Model NOT compatible with recursive ICA\n")
}

# Best ordering and p-value
best_idx <- which.max(rec_summary$p.value)
best_order <- perm_mat[best_idx, ]

cat("\nBest ordering:", paste(colnames(Y)[best_order], collapse = " -> "), "\n")
cat("Best p-value:", rec_summary$p.value[best_idx], "\n")

# Plot orthogonalised shocks for the best ordering

Ubest <- as.matrix(ures[, best_order, drop = FALSE])

Sigma_best <- cov(Ubest)
P_best <- t(chol(Sigma_best))

Ehat_best <- t(solve(P_best, t(Ubest)))
colnames(Ehat_best) <- paste0("e", 1:K)

plot.ts(Ehat_best, main = "Orthogonalised shocks (best ordering)", xlab = "Time")

# Scatterplot matrix
pairs(Ehat_best, main = "Pairwise plots of orthogonalised shocks")

############################################################
# Monte Carlo - Zero Restrictions (K = 3)
############################################################

rm(list = ls())

library(vars)
library(dse)
library(rmutil)
library(fastICA)
library(gtools)
library(matrixcalc)

source("functions.R")

#----------------------------------------------------------
# Settings
#----------------------------------------------------------

M <- 100 # Monte Carlo replications
N <- 500 # lower number of bootstrap iterations to speed up the simulation

p <- 2
K <- 3
h <- 10
Tsim <- 1000 + h

# True structure
B0inv <- diag(K)
B0inv[1, 2] <- 0.5
B0inv[1, 3] <- 0.3
B0inv[2, 3] <- 0.4

# True zero positions
true_zero <- matrix(0, K, K)
true_zero[2,1] <- 1
true_zero[3,1] <- 1
true_zero[3,2] <- 1

# Storage
zero_detect <- array(0, dim = c(K, K, M))

#----------------------------------------------------------
# Monte Carlo loop
#----------------------------------------------------------

for (m in 1:M) {
  
  cat("Replication:", m, "\n")
  
  # Lag polynomial
  AL <- array(0, dim = c(p + 1, K, K))
  AL[1,,] <- diag(K)
  
  AL[2,,] <- matrix(c(
    0.5, 0.1, 0.0,
    0.0, 0.3, 0.1,
    0.0, 0.0, 0.4
  ), K, K, byrow = TRUE)
  
  AL[3,,] <- matrix(c(
    -0.2, 0.0, 0.0,
    0.0,-0.1, 0.0,
    0.0, 0.0,-0.2
  ), K, K, byrow = TRUE)
  
  cvec <- c(1,2,1.5)
  
  var2 <- ARMA(A = AL, B = B0inv, TREND = cvec)
  
  # Simulate
  w <- matrix(rlaplace(Tsim * K, s = sqrt(1/2)), nrow = Tsim)
  
  varsim <- simulate(var2, sampleT = Tsim, noise = list(w = w))
  Y <- varsim$output[(h+1):Tsim, ]
  
  # Estimate VAR
  varest <- VAR(Y, p = p, type = "const")
  AA <- Acoef(varest)
  ures <- resid(varest)
  T <- nrow(ures)
  
  const <- sapply(1:K, function(eq) coef(varest)[[eq]]["const",1])
  
  # ICA
  A <- fAp_fastICA(ures)
  
  # Bootstrap
  Aboot <- matrix(NA, nrow = K^2, ncol = N)
  
  for (b in 1:N) {
    
    unew <- ures[sample(T, replace = TRUE), ]
    
    Ynew <- matrix(0, nrow = T + p, ncol = K)
    Ynew[1:p, ] <- Y[1:p, ]
    
    for (i in (p+1):(T+p)) {
      for (j in 1:p) {
        Ynew[i, ] <- Ynew[i, ] + AA[[j]] %*% Ynew[i-j, ]
      }
      Ynew[i, ] <- const + Ynew[i, ] + unew[i-p, ]
    }
    
    varest_new <- VAR(Ynew, p = p, type = "const")
    ures_star <- resid(varest_new)
    
    icares <- fastICA(ures_star, K, tol=1e-14, maxit=3000)
    W <- t((icares$K) %*% (icares$W))
    A_star <- solve(W)
    
    P <- myfrob(A_star, A, K)
    A_star <- A_star %*% P
    
    Aboot[,b] <- as.vector(A_star)
  }
  
  vecA <- as.vector(A)
  SD <- apply(Aboot, 1, sd)
  tstat <- vecA / SD
  
  df <- K^2 * (p + 1) + K
  tvalue <- qt(0.995, (T + p - df - 1))
  
  # Detect zeros
  mcoef <- matrix(0, K, K)
  
  for (i in 1:(K^2)) {
    if (abs(tstat[i]) <= tvalue) {
      mcoef[i] <- 1
    }
  }
  
  zero_detect[,,m] <- mcoef
}

#----------------------------------------------------------
# Results
#----------------------------------------------------------

freq_zero <- apply(zero_detect, c(1,2), mean)

print("Frequency of zero detection:")
round(freq_zero, 3)

barplot(c(0.98, 0.99, 1.00),
        names.arg = c("(2,1)", "(3,1)", "(3,2)"),
        ylab = "Frequency of zero classification")

############################################################
# Strength of Coefficients (K = 3)
############################################################

rm(list = ls())

library(vars)
library(dse)
library(rmutil)
library(fastICA)
library(gtools)
library(matrixcalc)

source("functions.R")

#----------------------------------------------------------
# Settings
#----------------------------------------------------------

gamma_vals <- c(0.5, 0.3, 0.2, 0.1)

p <- 2
K <- 3
h <- 10
Tsim <- 1000 + h
N <- 1000      
alpha <- 0.01

#----------------------------------------------------------
# Loop over gamma
#----------------------------------------------------------

for (gamma in gamma_vals) {
  
  cat("\n=============================\n")
  cat("Gamma =", gamma, "\n")
  cat("=============================\n")
  
  #-------------------------------
  # 1. DGP
  #-------------------------------
  
  # Lag polynomial
  AL <- array(0, dim = c(p + 1, K, K))
  AL[1,,] <- diag(K)
  
  AL[2,,] <- matrix(c(
    0.5, 0.1, 0.0,
    0.0, 0.3, 0.1,
    0.0, 0.0, 0.4
  ), K, K, byrow = TRUE)
  
  AL[3,,] <- matrix(c(
    -0.2, 0.0, 0.0,
    0.0,-0.1, 0.0,
    0.0, 0.0,-0.2
  ), K, K, byrow = TRUE)
  
  # Mixing matrix
  B0inv <- diag(K)
  B0inv[1,2] <- gamma
  B0inv[1,3] <- gamma
  B0inv[2,3] <- gamma
  
  cvec <- c(1,2,1.5)
  
  var2 <- ARMA(A = AL, B = B0inv, TREND = cvec)
  
  #-------------------------------
  # 2. Simulation
  #-------------------------------
  
  w <- matrix(rlaplace(Tsim * K, s = sqrt(1/2)), nrow = Tsim)
  
  varsim <- simulate(var2, sampleT = Tsim, noise = list(w = w))
  Y <- varsim$output[(h+1):Tsim, ]
  
  #-------------------------------
  # 3. VAR estimation
  #-------------------------------
  
  varest <- VAR(Y, p = p, type = "const")
  AA <- Acoef(varest)
  ures <- resid(varest)
  T <- nrow(ures)
  
  const <- sapply(1:K, function(eq) coef(varest)[[eq]]["const",1])
  
  #-------------------------------
  # 4. ICA
  #-------------------------------
  
  A <- fAp_fastICA(ures)
  
  #-------------------------------
  # 5. Bootstrap
  #-------------------------------
  
  Aboot <- matrix(NA, nrow = K^2, ncol = N)
  
  for (b in 1:N) {
    
    # Resample residuals
    unew <- ures[sample(T, replace = TRUE), ]
    
    # Reconstruct series
    Ynew <- matrix(0, nrow = T + p, ncol = K)
    Ynew[1:p, ] <- Y[1:p, ]
    
    for (i in (p+1):(T+p)) {
      for (j in 1:p) {
        Ynew[i, ] <- Ynew[i, ] + AA[[j]] %*% Ynew[i-j, ]
      }
      Ynew[i, ] <- const + Ynew[i, ] + unew[i-p, ]
    }
    
    # Re-estimate VAR
    varest_new <- VAR(Ynew, p = p, type = "const")
    ures_star <- resid(varest_new)
    
    # ICA
    icares <- fastICA(ures_star, K, tol=1e-14, maxit=3000)
    W <- t((icares$K) %*% (icares$W))
    A_star <- solve(W)
    
    # Align (permutation + sign)
    P <- myfrob(A_star, A, K)
    A_star <- A_star %*% P
    
    Aboot[,b] <- as.vector(A_star)
  }
  
  #-------------------------------
  # 6. Inference
  #-------------------------------
  
  vecA <- as.vector(A)
  SD <- apply(Aboot, 1, sd)
  tstat <- vecA / SD
  
  df <- K^2 * (p + 1) + K
  tvalue <- qt(1 - alpha/2, (T + p - df - 1))
  
  mcoef <- matrix("*", K, K)
  
  for (i in 1:(K^2)) {
    if (abs(tstat[i]) <= tvalue) {
      mcoef[i] <- "0"
    }
  }
  
  #-------------------------------
  # 7. Output
  #-------------------------------
  
  cat("Estimated A:\n")
  print(round(A, 3))
  
  cat("\nSignificance matrix S:\n")
  print(mcoef)
}



