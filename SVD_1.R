#Task: You may be knowing about a typesetting language called LaTeX. By default LaTeX
#places the things in the regular horizontal way, it can also rotate and scale (possibly
#different scaling factors in x- and y-directions): While writing a math book, we
#needed a way to place a writing on a slant surface drawn on the page:
#LaTeX has no direct way to achieve this. But we could get the effect by first doing
#a little SVD computation. How?



SVD = function(theta, phi, a, b) {
  x = c(seq(0, 1, by = 0.001), seq(0.25, 0.75, by = 0.001))
  y = c(2 * x[1:500], 2 - 2 * x[501:1001], rep(0.5, 501))
  plot(x, y, xlim = c(0, 4), ylim = c(-1, 1))
  
  A = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2)
  v = as.matrix(rbind(x, y))
  z = A %*% v
  
  B = matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), 2)
  R = diag(c(a, b))
  Q = B %*% R %*% t(B)
  
  w = Q %*% z
  plot(w[1, ], w[2, ], xlim = c(0, 4), ylim = c(-1, 1))
}

SVD(-0.9, -0.5, 1.2, 0.95)
SVD(0, 0, 1, 1)

#Reshaping a square to a parallelogram using R:
f = function(x1, y1, x2, y2) {
  a1 = (x2 - x1) / 2
  a2 = -(x1 + x2) / 2
  a3 = (y2 - y1) / 2
  a4 = -(y1 + y2) / 2
  
  A = matrix(c(a1, a2, a3, a4), 2, 2, byrow = TRUE)
  b = svd(A)
  U = b$u
  V = b$v
  
  if (abs(U[1, 1] + U[2, 2]) < 10^(-10)) {
    U[, 1] = U[, 1] * (-1)
  }
  if (abs(V[1, 1] + V[2, 2]) < 10^(-10)) {
    V[, 1] = V[, 1] * (-1)
  }
  
  sig = diag(b$d, 2)
  temp1 = acos(V[1, 1])
  temp2 = acos(U[1, 1])
  
  c = 180 / pi
  theta1 = (-c) * temp1
  theta2 = c * temp2
  
  if (abs(sin(temp1) - V[1, 2]) < 10^(-10)) {
    theta1 = -theta1
  }
  if (abs(sin(temp2) - U[1, 2]) < 10^(-10)) {
    theta2 = -theta2
  }
  
  x = sig[1, 1]
  y = sig[2, 2]
  
  cat("\n\nFirst angle of rotation (in degrees) is", theta1, "\n")
  cat("\n\nScale factor along x-axis is", x, "\n")
  cat("\n\nScale factor along y-axis is", y, "\n")
  cat("\n\nSecond angle of rotation (in degrees) is", theta2, "\n")
}

f(-6, -1, 1, -3)  # For modifying the letter alpha
f(-4, -3, 2, -3)  # For Italics Font

