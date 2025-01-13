acor <- function(y, x, a, type = "dcor") {

  Dx <- dim(x)        ;    Dy <- dim(y)
  ncx <- Dx[2] - 1    ;    ncy <- Dy[2] - 1   ;   nr <- Dx[1]

  if ( length(a) == 1 )  {

    if ( type == "cancor" ) {
      res <- 0
      y <- Compositional::alfa(y, a)$aff
      qy <- qr(y)
      dy <- qy$rank
      if ( dy > 0 ) {
        x <- Compositional::alfa(x, a)$aff
        qx <- qr(x)
        dx <- qx$rank
        if ( dx > 0 ) {
          res <- svd(qr.qty(qx, qr.qy( qy, diag(1, nr, dy)) )[1L:dx, ,drop = FALSE], dx, dy)$d
        }
      }  ##  end  if ( dy > 0 ) {

    } else if (type == "dcor") {
      y <- Compositional::alfa(y, a)$aff
      x <- Compositional::alfa(x, a)$aff
      res <- Rfast::dcor(y, x)$dcor

    } else {
      res <- numeric(2)
      y <- Compositional::alfa(y, a)$aff
      x <- Compositional::alfa(x, a)$aff
      res[1] <- Rfast::dcor(y, x)$dcor
      names(res) <- c("dcor", "cancor")

      qy <- qr(y)
      dy <- qy$rank
      if ( dy > 0 ) {
        qx <- qr(x)
        dx <- qx$rank
        if ( dx > 0 ) {
          res[2] <- svd(qr.qty(qx, qr.qy( qy, diag(1, nr, dy)) )[1L:dx, ,drop = FALSE], dx, dy)$d[1]
        }
      }  ##  end  if ( dy > 0 ) {
    }  ## end  if ( type == "cancor" ) {

  } else {

    if ( min(x) == 0  |  min(y) == 0 )  a <- a[ a > 0 ]
    nr <- dim(x)[1]
    res <- numeric( length(a) )
    names(res) <- paste("alpha=", a, sep = "")

    if ( type == "cancor" ) {

      for ( i in 1:length(a) ) {
        ep <- NULL
        y1 <- Compositional::alfa(y, a[i])$aff
        qy1 <- qr(y1)
        dy1 <- qy1$rank
        if ( dy1 > 0 )  {
          x1 <- Compositional::alfa(x, a[i])$aff
          qx1 <- qr(x1)
          dx1 <- qx1$rank
          if ( dx1 > 0 )  {
            ep <- svd(qr.qty(qx1, qr.qy( qy1, diag(1, nr, dy)) )[1L:dx, ,drop = FALSE], dx1, dy1)$d[1]
          }
        }  ##  end  if ( dy1 > 0 ) {
          if ( !is.null(ep) )  res[i] <- ep
      }  ##  end  for ( i in 1:length(a) ) {

    } else if ( type =="dcor" ) {
      res <- numeric( length(a) )
      names(res) <- paste("alpha=", a, sep = "")
      for ( i in 1:length(a) ) {
        res[i] <- Rfast::dcor( alfa(y, a[i])$aff, alfa(x, a[i])$aff )$dcor
      }

    } else {
      res <- matrix(NA, nrow = 2, ncol = length(a) )
      colnames(res) <- paste("alpha=", a, sep = "")
      rownames(res) <- c("dcor", "cancor")

      for ( i in 1:length(a) ) {
        y1 <- Compositional::alfa(y, a[i])$aff
        x1 <- Compositional::alfa(x, a[i])$aff
        res[1, i] <- Rfast::dcor(y1, x1)$dcor
        ep <- NULL
        qy1 <- qr(y1)
        dy1 <- qy1$rank
        if ( dy1 > 0 )  {
          qx1 <- qr(x1)
          dx1 <- qx1$rank
          if ( dx1 > 0 )  {
            ep <- svd(qr.qty(qx1, qr.qy( qy1, diag(1, nr, dy)) )[1L:dx, ,drop = FALSE], dx1, dy1)$d[1]
          }
        }  ##  end  if ( dy1 > 0 ) {
          if ( !is.null(ep) )  res[2, i] <- ep
        }
    }  ##  end  if ( type == "cancor" ) {

  }  ##  end   if ( length(a) == 1 )  {

  res
}

