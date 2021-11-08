alfa.contour <- function(m, s, a, n = 100, x = NULL, cont.line = FALSE) {
  ## the type parameter determines whether the additive or
  ## the isometric log-ratio transformation will be used.
  ## If type='alr' (the default) the additive
  ## log-ratio transformation is used.
  ## If type='ilr', the isometric log-ratio is used
  ## n is the number of points of each axis used
  nam <- c("X1", "X2", "X3")
  x1 <- seq(0.001, 0.999, length = n)
  sqrt3 <- sqrt(3)
  x2 <- seq(0.001, sqrt3/2 - 0.001, length = n)
  mat <- matrix(nrow = n, ncol = n)
  ha <- t( helm(3) )
  down <- det(2 * pi * s)^(-0.5)
  st <- solve(s)

  for ( i in 1:c(n/2) ) {
    for ( j in 1:n ) {
      if ( x2[j] < sqrt3 * x1[i] ) {
        ## This checks if the point will lie inside the triangle
        ## The next 4 lines calculate the composition
        w3 <- 2 * x2[j] / sqrt3
        w2 <- x1[i] - x2[j] / sqrt3
        w1 <- 1 - w2 - w3
        w <- c(w1, w2, w3)
        if ( abs(a) > 1e-9 ) {
          z <- w^a
          ta <- sum(z)
          z <- 3 / a * z / ta - 1/a
        } else  z <- log(w) - mean( log(w) )
        y <- as.vector( z %*% ha ) ## multiply by the Helmert sub-matrix
        can <- down * exp( -0.5 * ( ( y - m ) %*% st %*% ( y - m ) ) )
        if ( abs(can) < Inf )  mat[i, j] <- can
      }
    }
  }

  for ( i in c(n/2 + 1):n ) {
    for ( j in 1:n ) {
      ## This checks if the point will lie inside the triangle
      if ( x2[j] < sqrt3 - sqrt3 * x1[i] ) {
        ## The next 4 lines calculate the composition
        w3 <- 2 * x2[j] / sqrt3
        w2 <- x1[i] - x2[j] / sqrt3
        w1 <- 1 - w2 - w3
        w <- c(w1, w2, w3)
        if ( abs(a) > 1e-9 ) {
          z <- w^a
          ta <- sum(z)
          z <- 3 / a * z / ta - 1/a
        } else  z <- log(w) - mean( log(w) )
        y <- as.vector( z %*% ha ) ## multiply by the Helmert sub-matrix
        can <- down * exp( -0.5 * ( ( y - m ) %*% st %*% ( y - m ) ) )
        if ( abs(can) < Inf )   mat[i, j] <- can
      }
    }
  }


  # Create triangle corners
  b1 <- c(0.5, 0, 1, 0.5)

  b2 <- c(sqrt3/2, 0, 0, sqrt3/2)
  b <- cbind(b1, b2)

  # Axes
  b_x1 <- seq(from = 0, to = 1, length.out = 11)
  b_y1 <- rep(0, times = 11)

  b_x2 <- seq(from = 0.5, to = 0, length.out = 11)
  b_y2 <- seq(from = sqrt3/2, to = 0, length.out = 11)

  b_x4 <- seq(from = 1, to = 0.5, length.out = 11)
  b_y4 <- seq(from = 0, to = sqrt3/2, length.out = 11)

  # Plot window

  # Continuous color legend
  # Note that it disappears EVERY BLACK LINE!!!!!!
  # So, for the ones you want, you must do col = "black"
  # For more, see here
  # https://stackoverflow.com/questions/8068366/removing-lines-within-filled-contour-legend

  par(fg = NA)


  # Filled contoure plot in base R
  filled.contour(x1, x2, mat,

                 # Number of levels
                 # the greater the more interpolate
                 nlevels = 200,

                 # Colors with base R

                 # Select color function
                 color.palette = colorRampPalette( c( "blue",
                                                      "cyan",
                                                      "yellow",
                                                      "red") ),


                 # Adjust axes to points
                 plot.axes = {
                   ## Manual axes
                   # Axis 1
                   text(b_x1, b_y1,
                        c("","0.1", "0.2", "0.3", "0.4", "0.5",
                          "0.6", "0.7", "0.8", "0.9", ""),
                        adj = c(0.5, 1.5),
                        col = "black",
                        cex = 1);


                   # Axis 2
                   text(b_x2, b_y2,
                        c("","0.1", "0.2", "0.3", "0.4", "0.5",
                          "0.6", "0.7", "0.8", "0.9", ""),
                        adj = c(1.25, -0.15),
                        col = "black",
                        cex = 1);


                   # Axis 4
                   text(b_x4, b_y4,
                        c("","0.1", "0.2", "0.3", "0.4", "0.5",
                          "0.6", "0.7", "0.8", "0.9", ""),
                        adj = c(-0.25, -0.15),
                        col = "black",
                        cex = 1);

                   # Add compositional data
                   if ( !is.null(x) ) {
                     proj <- matrix(c(0, 1, 0.5, 0, 0, sqrt3/2), ncol = 2)
                     xa <- x %*% proj
                     points(xa[, 1], xa[, 2], col = "black")
                     nam2 <- colnames(x)
                     if ( !is.null(nam2) )  nam <- nam2
                   };

                   # Show corner titles
                   text( b[1, 1], b[1, 2] + 0.07, nam[3], cex = 1, col = "black", font = 2 );
                   text( b[2:3, 1], b[2:3, 2] - 0.07, nam[1:2], cex = 1, col = "black", font = 2 );

                   # Draw triangle in two dimensions
                   points(b[, 1], b[, 2],
                          type = "l",
                          lwd = 4.5,
                          col = "black");

                   # Add contour lines
                   if ( cont.line ) {
                     contour(x1, x2, mat,
                             pt = "s",
                             # Color of contour lines
                             # Otherwise par(fg = NA) will
                             # disappear them...

                             col="black",


                             # Number of levels
                             nlevels = 10,

                             # Size of contour numbers
                             labcex = 0.8,

                             # Width of contour lines
                             lwd = 1,

                             add = TRUE) }
                 },

                 # Legend tick lines
                 key.axes = {axis(4, col = "black")},

                 # Axes labs
                 xlab = "",
                 ylab = "",

                 # Plot limits
                 xlim = c(-0.1, 1.1),
                 ylim = c(-0.1, 1.1) )

  # Plot window
}
