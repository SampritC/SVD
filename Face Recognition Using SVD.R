library(jpeg)

# Parameters
H <- 200
W <- 180
n_i <- 10
n_j <- 5

# Allocate data matrix
M <- matrix(0, nrow = n_i * n_j, ncol = H * W)

# Subject names (length 10; index with i)
names_vec <- c(
  "male.pacole",
  "male.pspliu",
  "male.sjbeck",
  "male.skumar",
  "male.rsanti",
  "female.anpage",
  "female.asamma",
  "female.klclar",
  "female.ekavaz",
  "female.drbost"
)

k <- 0
for (i in 1:n_i) {
  for (j in 1:n_j) {
    k <- k + 1
    # Build filename like ".../grayfaces/male.pacole.1.jpg"
    filename <- file.path(
      "C:/Users/samah/OneDrive/Desktop/grayfaces",
      paste0(names_vec[i], ".", j, ".jpg")
    )
    if (!file.exists(filename)) stop("File not found: ", filename)
    
    img <- readJPEG(filename)               # matrix (grayscale) or 3D array (RGB)
    if (length(dim(img)) == 3L) {
      # Convert RGB → grayscale (simple average); adjust if you prefer luma
      img <- apply(img, c(1, 2), mean)
    }
    if (!all(dim(img) == c(H, W))) {
      stop("Image has dims ", paste(dim(img), collapse = "x"),
           "; expected ", H, "x", W, " for file: ", filename)
    }
    M[k, ] <- as.vector(img)
  }
}

# Column means and mean face
means <- apply(M, 2, mean)

# Show mean face
plot(NA, xlim = c(1, 2), ylim = c(1, 2), type = "n", xlab = "", ylab = "", axes = FALSE)
y <- means
dim(y) <- c(H, W)
rasterImage(as.raster(y), 1, 1, 2, 2)

# Center data, build S = N N^T, eigendecompose
N <- scale(M, scale = FALSE)
S <- N %*% t(N)
eig <- eigen(S)

# Scree
plot(eig$values, type = "b", xlab = "Component", ylab = "Eigenvalue")

# Project to top-8 components
P <- t(eig$vectors[, 1:8]) %*% N   # 8 x (H*W)

# Normalize each projected vector to unit length
shape <- function(x) x / sqrt(sum(x * x))
Q <- apply(P, 2, shape)

# Reconstruct centered images from 8 components
L <- eig$vectors[, 1:8] %*% Q      # 50 x (H*W)

# Add back mean to each reconstructed row
for (j in 1:(n_i * n_j)) {
  L[j, ] <- L[j, ] + means
}

# Display one reconstructed face (e.g., #23), rescaled to [0,1]
plot(NA, xlim = c(1, 2), ylim = c(1, 2), type = "n", xlab = "", ylab = "", axes = FALSE)
y <- abs(L[23, ])
ext <- range(y)
y <- (y - ext[1]) / (ext[2] - ext[1])
dim(y) <- c(H, W)
rasterImage(as.raster(y), 1, 1, 2, 2)

##SVD Version

# Reading images and storing them in the matrix named " images "
library(jpeg)

# Set working directory
setwd("C:/Users/User/Desktop/facesbwandcol/grayfaces")

# Create a 50 x 36000 matrix to hold image data
images <- matrix(0, 50, 36000)

# Vector of subject names
v <- c(
  "female.anpage",
  "female.asamma",
  "female.drbost",
  "female.ekavaz",
  "female.klclar",
  "male.pacole",
  "male.pspliu",
  "male.rsanti",
  "male.sjbeck",
  "male.skumar"
)
# Build file matrix -----------------------------------------------------------
for (i in 1:10) {
  for (j in 1:5) {
    filename <- paste0(v[i], ".", j, ".jpg")
    images[(5 * i - 5 + j), ] <- as.vector(readJPEG(filename))
  }
}

# Mean face -------------------------------------------------------------------
meanvec <- colMeans(images)             # 1 x 36000 numeric vector

# Center images ---------------------------------------------------------------
matrix1 <- images
for (i in 1:50) {
  matrix1[i, ] <- matrix1[i, ] - meanvec
}

# SVD (eigenfaces live in V) --------------------------------------------------
newmat <- svd(matrix1)
V <- newmat$v   # 36000 x 50 (right singular vectors)

# Display helper --------------------------------------------------------------
display <- function(num1 = 0, num2 = 0, meanface = 0, vec = matrix(0, 1, 36000)) {
  # num1: subject index in v (1..10); num2: photo index (1..5)
  # meanface = 1 to show average face
  # vec: optionally provide a 1 x 36000 vector to display directly
  
  if (meanface == 1) {
    img <- matrix(meanvec, 200, 180)
  } else {
    k <- 5 * num1 - 5 + num2
    img <- matrix(images[k, ], 200, 180)
  }
  
  if (!identical(vec, matrix(0, 1, 36000))) {
    img <- matrix(vec, 200, 180)
  }
  
  # Flip to match expected orientation (transpose + horizontal flip)
  img1 <- t(img)
  img2 <- matrix(0, 180, 200)
  for (i in 1:200) img2[, 201 - i] <- img1[, i]
  
  image(img2, useRaster = TRUE, col = grey(seq(0, 1, length = 256)), axes = FALSE, xlab = "", ylab = "")
}

# Linear combination / reconstruction ----------------------------------------
lin_comb <- function(vect, n) {
  # Reconstruct using the first n eigenfaces (columns of V)
  p <- numeric(n)
  face <- matrix(vect, 1, 36000)
  for (i in 1:n) {
    p[i] <- sum((vect - meanvec) * V[, i])
    face <- face + p[i] * matrix(V[, i], 1, 36000)
  }
  display(1, 1, 0, (face + meanvec))
}

# Example: reconstructed image using 6 leading eigenfaces
lin_comb(images[25, ], 6)
