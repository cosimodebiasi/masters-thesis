############################################################
# Function: ICA with Column Reordering
############################################################

fAp_fastICA <- function(ures, sseed = 46) {
  
  set.seed(sseed)
  
  n <- ncol(ures)
  X <- t(ures)
  
  icares <- fastICA(
    t(X),
    n,
    tol = 1e-14,
    maxit = 3000,
    verbose = FALSE
  )
  
  set.seed(NULL)
  
  W <- t((icares$K) %*% (icares$W))
  A <- solve(W)
  
  # Reorder columns based on largest absolute values
  aba <- abs(A)
  CC  <- matrix(0, n, 2)
  
  for (i in 1:n) {
    cc <- which(aba == max(aba), arr.ind = TRUE)
    aba[cc[1], ] <- 0
    aba[, cc[2]] <- 0
    CC[i, ] <- cc
  }
  
  sn <- sapply(1:n, function(i) CC[CC[,1] == i, 2])
  
  Ap <- A[, sn]
  colnames(Ap) <- paste0("s", 1:n)
  
  # Enforce positive diagonal
  for (i in 1:n) {
    if (Ap[i, i] < 0) {
      Ap[, i] <- -Ap[, i]
    }
  }
  
  return(Ap)
}

############################################################
# Function: Sign and Permutation (Frobenius minimisation)
############################################################

myfrob <- function(A.hat, A, K) {
  
  perm <- permutations(ncol(A), ncol(A))
  sign <- permutations(2, ncol(A), repeats.allowed = TRUE)
  sign[sign == 2] <- -1
  
  pr <- diag(ncol(A))
  
  best_norm <- Inf
  best_P    <- NULL
  
  for (j in 1:nrow(perm)) {
    for (i in 1:nrow(sign)) {
      
      P <- pr[, perm[j, ]] * sign[i, ]
      A_tilde <- A.hat %*% P
      
      frob <- (1 / sqrt(K - 1)) * frobenius.norm(A_tilde - A)
      
      if (frob < best_norm) {
        best_norm <- frob
        best_P    <- P
      }
    }
  }
  
  return(best_P)
}


#############################################################
# Function: Structural Impulse Response Functions via ICA
#############################################################

compute_struct_irf <- function(varest_obj, A_mat, h = 13) {
  
  Phi <- Acoef(varest_obj)
  K <- ncol(A_mat)
  p <- length(Phi)
  
  # Companion matrix
  Fmat <- matrix(0, nrow = K * p, ncol = K * p)
  Fmat[1:K, 1:(K * p)] <- do.call(cbind, Phi)
  
  if (p > 1) {
    Fmat[(K + 1):(K * p), 1:(K * (p - 1))] <- diag(K * (p - 1))
  }
  
  # Selector matrix
  J <- cbind(diag(K), matrix(0, nrow = K, ncol = K * (p - 1)))
  
  # IRFs: response x shock x horizon
  IRF <- array(0, dim = c(K, K, h + 1))
  IRF[, , 1] <- A_mat
  
  Fpower <- diag(K * p)
  
  for (s in 1:h) {
    Fpower <- Fpower %*% Fmat
    IRF[, , s + 1] <- J %*% Fpower %*% t(J) %*% A_mat
  }
  
  dimnames(IRF) <- list(
    response = colnames(A_mat),
    shock = colnames(A_mat),
    horizon = 0:h
  )
  
  return(IRF)
}
