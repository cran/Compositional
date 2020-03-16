################################
#### Classification for compositional data using the alpha-transformation
#### The k-NN algorithm
#### Tsagris Michail 8/2015
#### References: Tsagris, M., Preston S. and Wood A.T.A. (2016).
#### Improved classication for compositional data using the alpha-transformation
#### Journal of Classification (To appear)
#### http://arxiv.org/pdf/1506.04976v2.pdf
#### mtsagris@yahoo.gr
################################
alfa.knn <- function(xnew, x, ina, a = 1, k = 5, type = "S", mesos = TRUE, apostasi = "euclidean", rann = FALSE) {
  ## x is the matrix containing the data
  ## ina indicates the groups
  ## a is the value of the power parameter
  ## k in the number of nearest neighbours
  ## apostasi is the type of metric used, either "euclidean" or "manhattan"
  ## type is either S or NS. Should the standard k-NN be use or not
  ## if mesos is TRUE, then the arithmetic mean distance of the
  ## k nearest points will be used
  ## Both of these apply for the non-standard,
  ## algorithm, that is when type=NS
  ## xnew is the new dataset. It can be a single vector or a matrix
  p <- dim(x)[2]
  xnew <- as.matrix(xnew)
  xnew <- matrix(xnew, ncol = p)  ## makes sure xnew is a matrix
  ina <- as.numeric(ina)
  nc <- max(ina) ## The number of groups
  nu <- dim(xnew)[1]
  if ( !is.null(a) ) {
    znew <- Compositional::alfa(xnew, a, h = FALSE)$aff
    z <- Compositional::alfa(x, a, h = FALSE)$aff
  } else {
    znew <- xnew
	  z <- x
  }
  if (type == "NS") {
    ## Non Standard algorithm
    klen <- length(k)
    g <- matrix(0, nu, klen)
    ta <- matrix(nrow = nu, ncol = nc)
    apo <- list()
    for (m in 1:nc) {
      disa <- Rfast::dista(znew, z[ina == m,], type = apostasi, trans = FALSE)
      apo[[ m ]] <- Rfast::colSort(disa)[1:max(k), ]
    }
    if ( mesos ) {
      for (j in 1:klen) {
        for (m in 1:nc) {
          ta[, m] <- Rfast::colmeans( apo[[ m ]][1:k[j], , drop = FALSE] )
        }
      }
      g[, j] <- Rfast::rowMins(ta)
    } else {  ## mesos = FALSE
      for (j in 1:klen) {
        for (m in 1:nc) {
          ta[, m] <- Rfast::colhameans( apo[[ m ]][1:k[j], , drop = FALSE] )
		}
      }
      g[, j] <- Rfast::rowMins(ta)
    }  ## end if (mesos)

  } else if (type == "S") {
    ## Standard algorithm
    if ( rann ) {
     klen <- length(k)
     di <- RANN::nn2( data = x, query = xnew, k = max(k) )$nn.idx
     g <- matrix(nrow = nu, ncol = klen)
     m1 <- matrix(nrow = max(k), ncol = nu)
     for ( i in 1:nu )  m1[, i] <- ina[ di[i, ] ]
     for ( j in 1:klen ) g[, j] <- Rfast::colMaxs( Rfast::colTabulate( m1[1:k[j], ] ) )
    } else  g <- Rfast::knn(xnew = znew, y = ina, x = z, k = k, dist.type = apostasi, type = "C", freq.option = 1)
  }  ## end if (type == "S")
  colnames(g) <- paste("k=", k, sep = "")
  g
}



