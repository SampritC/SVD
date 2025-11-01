##Algorithm 1

# Householder QR decomposition
QR = function(A) {
  m <- nrow(A)
  n <- ncol(A)
  
  unit <- function(v) v / sqrt(sum(v * v))
  hmult <- function(u, x) x - 2 * sum(u * x) * u
  shaver <- function(x) {
    x[1] <- x[1] - sqrt(sum(x * x))
    unit(x)
  }
  
  Q <- diag(1, m)
  R <- A
  
  for (i in 1:min(m, n)) {
    # Skip if the below-diagonal part of the i-th column is (nearly) zero
    if (i < m && sum(R[(i + 1):m, i]^2) > 1e-7) {
      u <- shaver(R[i:m, i])
      
      # Build the Householder reflector H_i of size m x m
      H <- diag(1, m)
      uuT <- matrix(u, ncol = 1) %*% t(matrix(u, ncol = 1))
      H_sub <- diag(1, m - i + 1) - 2 * uuT
      H[i:m, i:m] <- H_sub
      
      # Apply H on the left: R <- H R, Q <- Q H
      R <- H %*% R
      Q <- Q %*% H
    }
  }
  
  # Round R just like your original code
  R <- round(R, digits = 6)
  list(Q = Q, R = R)
}

# QR-iteration eigen-decomposition (orthogonal iteration)
eig = function(M, max_iter = 1000) {
  X <- M
  P <- diag(1, nrow(M))
  for (i in 1:max_iter) {
    decomp <- QR(X)
    Q <- decomp$Q
    R <- decomp$R
    P <- P %*% Q
    X <- R %*% Q
  }
  # X converges to (approximately) upper-triangular with eigenvalues on the diagonal
  vals <- diag(round(X, 5))
  vecs <- round(P, 5)
  list(values = vals, vectors = vecs)
}

# SVD via eigen-decomposition of M^T M and M M^T
SVDcomp = function(M) {
  D <- t(M) %*% M
  E <- M %*% t(M)
  
  eD <- eig(D)  # V from right-eigenvectors of M^T M
  eE <- eig(E)  # U from eigenvectors of M M^T
  
  # singular values are sqrt of eigenvalues of M^T M
  svals <- sqrt(as.numeric(eD$values))
  
  # Optional: sort in decreasing order (keep U,V in the same order if needed)
  ord <- order(svals, decreasing = TRUE)
  svals <- svals[ord]
  V <- eD$vectors[, ord, drop = FALSE]
  U <- eE$vectors[, ord, drop = FALSE]
  
  # Print to match your original intent
  print(U)
  print(svals)
  print(V)
  
  invisible(list(U = U, d = svals, V = V))
}

# Example call:
M <- matrix(c(3, 1, 0,
               1, 3, 1,
               0, 1, 3), nrow = 3, byrow = TRUE)
SVDcomp(M)

##Algorithm 2

SVD = function(B) {
  # A function to compute the SVD of a matrix with nrow >= 2, ncol >= 2
  
  # If B is "wide", we will process t(B) and swap U,V later; code below pads by +1
  if (nrow(B) < ncol(B)) {
    A <- matrix(0, ncol(B) + 1, nrow(B) + 1)
    A[1:ncol(B), 1:nrow(B)] <- t(B)
  } else {
    A <- matrix(0, nrow(B) + 1, ncol(B) + 1)
    A[1:nrow(B), 1:ncol(B)] <- B
  }
  
  m <- nrow(A) - 1
  n <- ncol(A) - 1
  
  unit <- function(v) v / sqrt(sum(v * v))                          # unit vector
  hmult <- function(u, x) x - 2 * sum(u * x) * u                     # (I - 2uu^T) x
  shaver <- function(x) { x[1] <- x[1] - sqrt(sum(x * x)); unit(x) } # Householder "shave"
  
  # Golub–Kahan bidiagonalization storing Householder vectors in the extra row/col of A
  for (i in 1:n) {
    # Left reflector on column i (below diagonal)
    if (i < m && sum((A[(i + 1):m, i])^2) > 10^(-12)) {
      if (i != m) {
        temp <- shaver(A[i:m, i])
      } else {
        temp <- 0
      }
      for (j in i:n) {
        A[i:m, j] <- hmult(temp, A[i:m, j])  # apply Householder on the left
      }
      A[(i + 1):(m + 1), i] <- temp          # store vector
    } else {
      if (i < m) A[(i + 1):(m + 1), i] <- 0
    }
    
    # Right reflector on row i (to zero out superdiagonal entries beyond i+1)
    if (i <= (n - 2)) {
      if (sum((A[i, (i + 2):n])^2) > 10^(-12)) {
        temp <- shaver(A[i, (i + 1):n])
        for (k in i:m) {
          A[k, (i + 1):n] <- hmult(temp, A[k, (i + 1):n])  # apply Householder on the right
        }
        A[i, (i + 2):(n + 1)] <- temp                      # store vector
      } else {
        A[i, (i + 2):(n + 1)] <- 0
      }
    }
  }
  
  # Build U_A from stored left Householder vectors
  U_A <- diag(1, m)
  for (i in 1:n) {
    newmat1 <- diag(1, m)
    uvec <- A[(i + 1):(m + 1), i]
    uuT <- matrix(uvec, ncol = 1) %*% t(matrix(uvec, ncol = 1))
    temp <- diag(1, (m + 1 - i)) - 2 * uuT
    newmat1[i:m, i:m] <- temp
    U_A <- U_A %*% newmat1
  }
  
  # Build V_A from stored right Householder vectors
  V_A <- diag(1, n)
  if (n >= 3) {
    for (i in 1:(n - 2)) {
      newmat1 <- diag(1, n)
      vvec <- A[i, (i + 2):(n + 1)]
      vvT <- matrix(vvec, ncol = 1) %*% t(matrix(vvec, ncol = 1))
      temp <- diag(1, (n - i)) - 2 * vvT
      newmat1[(i + 1):n, (i + 1):n] <- temp
      V_A <- V_A %*% newmat1
    }
  }
  
  # Bidiagonal block to iterate on
  if (nrow(B) >= ncol(B)) {
    bid <- t(U_A) %*% B %*% V_A
  } else {
    bid <- t(U_A) %*% t(B) %*% V_A
  }
  bid <- round(bid, digits = 6)  # bidiagonal (approximately)
  
  # Truncate U_A to m x n so later multiplications align with the n-by-n core
  U_A <- U_A[1:m, 1:n, drop = FALSE]
  
  # Diagonalize the n x n bidiagonal core using alternating QR/LQ steps
  U3 <- diag(1, n)
  V3 <- diag(1, n)
  
  epsilon <- max(abs(bid[1:n, 1:n] - diag(diag(bid[1:n, 1:n]))))
  while (epsilon > 1e-7) {
    Q <- qr.Q(qr(bid[1:n, 1:n]))
    R <- qr.R(qr(bid[1:n, 1:n]))
    L <- t(qr.R(qr(t(R))))
    P <- qr.Q(qr(t(R)))
    
    U3 <- U3 %*% Q
    V3 <- t(P) %*% V3
    bid[1:n, 1:n] <- L
    
    epsilon <- max(abs(bid[1:n, 1:n] - diag(diag(bid[1:n, 1:n]))))
  }
  
  V3 <- t(V3)
  
  U <- U_A %*% U3
  V <- V_A %*% V3
  sig <- round(bid[1:n, 1:n], digits = 6)
  
  # If we flipped to t(B) earlier, undo by swapping U and V and transposing reconstruction
  if (nrow(B) < ncol(B)) {
    mat1 <- t(U %*% sig %*% t(V))
    tmp <- U; U <- V; V <- tmp
  } else {
    mat1 <- U %*% sig %*% t(V)
  }
  
  # Make singular values nonnegative (flip corresponding U columns if needed)
  for (i in 1:n) {
    if (sig[i, i] < 0) {
      U[, i] <- -U[, i]
      sig[i, i] <- -sig[i, i]
    }
  }
  
  cat("\nU is:\n")
  print(round(U, digits = 6))   # U (may be tall); extend if you need full square orthonormal basis
  
  cat("\nSigma is:\n")
  print(sig)                    # diagonal matrix of singular values (n x n)
  
  cat("\nV is:\n")
  print(round(V, digits = 6))   # V (n x n block; may be non-square vs original dims)
  
  cat("\nThe product U %*% Sigma %*% t(V) is:\n")
  print(round(mat1, digits = 5))
}

# Examples
mat <- matrix(c(2, 6, 1,
                8, 5, 0,
                7, 4, 3), 3, 3, byrow = TRUE)
SVD(mat)

mat <- matrix(1:15, 5, 3)
SVD(mat)

mat <- matrix(1:40, 5, 8)
SVD(mat)
