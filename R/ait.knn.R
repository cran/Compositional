ait.knn <- function (xnew, x, ina, a = 1, k = 5, type = "S", mesos = TRUE, 
    apostasi = "euclidean", rann = FALSE) {

    p <- dim(x)[2]
    xnew <- as.matrix(xnew)
    xnew <- matrix(xnew, ncol = p)
    ina <- as.numeric(ina)
    nc <- max(ina)
    nu <- dim(xnew)[1]
	
    if (!is.null(a)) {
        znew <- Compositional::ait(xnew, a, h = FALSE)
        z <- Compositional::ait(x, a, h = FALSE)
    }   else {
        znew <- xnew
        z <- x
    }
	
    if ( type == "NS" ) {
        klen <- length(k)
        g <- matrix(0, nu, klen)
        ta <- matrix(nrow = nu, ncol = nc)
        apo <- list()
        for (m in 1:nc) {
            disa <- Rfast::dista(znew, z[ina == m, ], type = apostasi, 
                trans = FALSE)
            apo[[m]] <- Rfast::colSort(disa)[1:max(k), ]
        }
        if (mesos) {
            for (j in 1:klen) {
                for (m in 1:nc) {
                  ta[, m] <- Rfast::colmeans(apo[[m]][1:k[j], , drop = FALSE])
                }
            }
            g[, j] <- Rfast::rowMins(ta)
        }   else {
            for (j in 1:klen) {
                for (m in 1:nc) {
                  ta[, m] <- Rfast::colhameans(apo[[m]][1:k[j], , drop = FALSE])
                }
            }
            g[, j] <- Rfast::rowMins(ta)
        }
    }   else if ( type == "S" ) {
        if (rann) {
            klen <- length(k)
            di <- RANN::nn2(data = x, query = xnew, k = max(k))$nn.idx
            g <- matrix(nrow = nu, ncol = klen)
            m1 <- matrix(nrow = max(k), ncol = nu)
            for (i in 1:nu) m1[, i] <- ina[di[i, ]]
            for (j in 1:klen) g[, j] <- Rfast::colMaxs( Rfast::colTabulate(m1[1:k[j], ]) )
        }   else g <- Rfast::knn(xnew = znew, y = ina, x = z, k = k, 
            dist.type = apostasi, type = "C", freq.option = 1)
    }
    colnames(g) <- paste("k=", k, sep = "")
    g
}
