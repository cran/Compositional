 ################################
#### Contour plot of the Dirichlet distribution in S^2
#### Tsagris Michail 1/2013
#### mtsagris@yahoo.gr
################################
diri.contour <- function(a, n = 100, x = NULL, cont.line = FALSE) {
  ## a are the estimated Dirichlet parameters
  ## n shows the number of points at which the density is calculated
  ## so, n^2 points are used.
  ## x should be a 3-part compositional data or NULL for no data
  nam <- c("X1", "X2", "X3")
  x1 <- seq(0.001, 0.999, length = n)  ## coordinates of x
  sqrt3 <- sqrt(3)
  x2 <- seq(0.001, sqrt3/2 - 1e-03, length = n)  ## coordinates of y
  # old code
  # mat <- matrix(nrow = n, ncol = n)
  # be <- prod( gamma(a)) / gamma(sum(a) )  ## beta function
  #
  # for ( i in 1:c(n/2) ) {
  #   for (j in 1:n) {
  #     if ( x2[j] < sqrt3 * x1[i] ) {
  #       ## This checks if the point will lie inside the triangle
  #       ## the next three lines invert the points which lie inside
  #       ## the triangle back into the composition in S^2
  #       w3 <- 2 * x2[j] / sqrt3
  #       w2 <- x1[i] - x2[j]/sqrt3
  #       w1 <- 1 - w2 - w3
  #       w <- c(w1, w2, w3)
  #       can <- prod( w^(a - 1) ) / be
  #       if (abs(can) < Inf)  mat[i, j] <- can
  #     }
  #   }
  # }
  # for (i in c(n/2 + 1):n) {
  #   for (j in 1:n) {
  #     ## This checks if the point will lie inside the triangle
  #     if ( x2[j] < sqrt3 - sqrt3 * x1[i] ) {
  #       ## the next three lines invert the points which lie inside
  #       ## the triangle back into the composition in S^2
  #       w3 <- 2 * x2[j] / sqrt3
  #       w2 <- x1[i] - x2[j]/sqrt3
  #       w1 <- 1 - w2 - w3
  #       w <- c(w1, w2, w3)
  #       can <- prod( w^(a - 1) ) / be
  #       if (abs(can) < Inf)  mat[i, j] <- can
  #     }
  #   }
  # }
  wa <- NULL
  for ( i in 1:n ) {
    w3 <- 2 * x2 / sqrt3
    w2 <- x1[i] - x2/sqrt3
    w1 <- 1 - w2 - w3
    wa <- rbind(wa, cbind(w1, w2, w3) )
  }

  can <- Compositional::ddiri(wa, a, logged = FALSE)
  mat <- matrix(can, byrow = TRUE, nrow = n, ncol = n)

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
##### C A U T I O N !!!
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
                   color.palette =  colorRampPalette( c( "blue",
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
                               # Draw triangle in two dimensions
                                points(b[, 1], b[, 2],
                                       type = "l",
                                       lwd = 2,
                                       col = "black");
                       # Add 3-part compositional data
                       if ( !is.null(x) ) {
                          proj <- matrix(c(0, 1, 0.5, 0, 0, sqrt3/2), ncol = 2)
                          xa <- x %*% proj
                          points(xa[, 1], xa[, 2], col = "black")
                          nam2 <- colnames(x)
                          if ( !is.null(nam2) )  nam2 <- nam
                       };
                          # Show corner titles
                           text( b[1, 1], b[1, 2] + 0.07, nam[3], cex = 1, col = "black", font = 2 );
                           text( b[2:3, 1], b[2:3, 2] - 0.07, nam[1:2], cex = 1, col = "black", font = 2 );
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
                                         lwd = 1.5,
                                         add = TRUE)
                               }
                  },
                  # Legend tick lines
                   key.axes = {axis(4, col = "black")},
                  # Axes labs
                   xlab = "",
                   ylab = "",
                  # Plot limits
                   xlim = c(-0.1, 1.1),
                   ylim = c(-0.1, 1.1) )
# Reset Plot window
}












# kent.contour <- function(G, param, n = 100, x = NULL, cont.line = FALSE) {
#
#   nam <- c("X1", "X2", "X3")
#   x1 <- seq(0.001, 0.999, length = n)  ## coordinates of x
#   sqrt3 <- sqrt(3)
#   x2 <- seq(0.001, sqrt3/2 - 1e-03, length = n)  ## coordinates of y
#
#   wa <- NULL
#   for ( i in 1:n ) {
#     w3 <- 2 * x2 / sqrt3
#     w2 <- x1[i] - x2/sqrt3
#     w1 <- 1 - w2 - w3
#     wa <- rbind(wa, cbind(w1, w2, w3) )
#   }
#
#   y <- sqrt(wa)
#   can <- Directional::dkent(y, G, param)
#   mat <- matrix(can, byrow = TRUE, nrow = n, ncol = n)
#
#   # Create triangle corners
#   b1 <- c(0.5, 0, 1, 0.5)
#
#   b2 <- c(sqrt3/2, 0, 0, sqrt3/2)
#   b <- cbind(b1, b2)
#
#   # Axes
#   b_x1 <- seq(from = 0, to = 1, length.out = 11)
#   b_y1 <- rep(0, times = 11)
#
#   b_x2 <- seq(from = 0.5, to = 0, length.out = 11)
#   b_y2 <- seq(from = sqrt3/2, to = 0, length.out = 11)
#
#   b_x4 <- seq(from = 1, to = 0.5, length.out = 11)
#   b_y4 <- seq(from = 0, to = sqrt3/2, length.out = 11)
#
#
#   # Plot window
#   ##### C A U T I O N !!!
#   # Continuous color legend
#   # Note that it disappears EVERY BLACK LINE!!!!!!
#   # So, for the ones you want, you must do col = "black"
#   # For more, see here
#   # https://stackoverflow.com/questions/8068366/removing-lines-within-filled-contour-legend
#
#   par(fg = NA)
#
#
#   # Filled contoure plot in base R
#   filled.contour(x1, x2, mat,
#
#                  # Number of levels
#                  # the greater the more interpolate
#                  nlevels = 200,
#
#                  # Colors with base R
#
#                  # Select color function
#                  color.palette =  colorRampPalette( c( "blue",
#                                                        "cyan",
#                                                        "yellow",
#                                                        "red") ),
#
#                  # Adjust axes to points
#                  plot.axes = {
#                    ## Manual axes
#                    # Axis 1
#                    text(b_x1, b_y1,
#                         c("","0.1", "0.2", "0.3", "0.4", "0.5",
#                           "0.6", "0.7", "0.8", "0.9", ""),
#                         adj = c(0.5, 1.5),
#                         col = "black",
#                         cex = 1);
#
#
#                    # Axis 2
#                    text(b_x2, b_y2,
#                         c("","0.1", "0.2", "0.3", "0.4", "0.5",
#                           "0.6", "0.7", "0.8", "0.9", ""),
#                         adj = c(1.25, -0.15),
#                         col = "black",
#                         cex = 1);
#
#
#                    # Axis 4
#                    text(b_x4, b_y4,
#                         c("","0.1", "0.2", "0.3", "0.4", "0.5",
#                           "0.6", "0.7", "0.8", "0.9", ""),
#                         adj = c(-0.25, -0.15),
#                         col = "black",
#                         cex = 1);
#
#
#                    # Draw triangle in two dimensions
#                    points(b[, 1], b[, 2],
#                           type = "l",
#                           lwd = 2,
#                           col = "black");
#
#
#                    # Add 3-part compositional data
#                    if ( !is.null(x) ) {
#                      proj <- matrix(c(0, 1, 0.5, 0, 0, sqrt3/2), ncol = 2)
#                      xa <- x %*% proj
#                      points(xa[, 1], xa[, 2], col = "black")
#                      nam2 <- colnames(x)
#                      if ( !is.null(nam2) )  nam <- nam2
#                    };
#
#
#                    # Show corner titles
#                    text( b[1, 1], b[1, 2] + 0.07, nam[3], cex = 1, col = "black", font = 2 );
#                    text( b[2:3, 1], b[2:3, 2] - 0.07, nam[1:2], cex = 1, col = "black", font = 2 );
#
#                    # Add contour lines
#
#                    if ( cont.line ) {
#                      contour(x1, x2, mat,
#                              pt = "s",
#                              # Color of contour lines
#                              # Otherwise par(fg = NA) will
#                              # disappear them...
#
#                              col="black",
#
#                              # Number of levels
#                              nlevels = 10,
#
#                              # Size of contour numbers
#                              labcex = 0.8,
#
#                              # Width of contour lines
#                              lwd = 1.5,
#
#                              add = TRUE)
#                    }
#                  },
#
#                  # Legend tick lines
#                  key.axes = {axis(4, col = "black")},
#
#                  # Axes labs
#                  xlab = "",
#                  ylab = "",
#
#                  # Plot limits
#                  xlim = c(-0.1, 1.1),
#                  ylim = c(-0.1, 1.1) )
#
#   # Reset Plot window
# }

