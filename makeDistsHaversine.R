
makeDistsHaversine<-function (datacoords, knotcoords, knotmat = TRUE, polys = NULL, 
          type = "A", plot.transition = FALSE, grid.dim = c(100, 100)) 
{
  if (is.null(polys)) {
    if (length(which(is.na(knotcoords[, 1]))) > 0) 
      stop("remove NAs from knotcoords")
    d2k <- matrix(0, ncol = dim((knotcoords))[1], nrow = length(datacoords[, 
                                                                           1]))
    for (i in 1:dim(knotcoords)[1]) {
      # d2k[, i] <- sqrt((datacoords[, 1] - knotcoords[i,  1])^2 + 
      #                    (datacoords[, 2] - knotcoords[i, 2])^2)
      
      d2k[, i] <- distm(datacoords, knotcoords[i, ], fun = distHaversine)
    }
    if (knotmat == T) {
      knotDist = as.matrix(dist(na.omit(knotcoords), method = "euclidean", 
                                diag = TRUE, upper = TRUE))
      return(list(dataDist = d2k, knotDist = knotDist))
    }
    else {
      return(list(dataDist = d2k))
    }
  }
  if (!is.null(polys)) {
    Nx = seq(min(c(datacoords[, 1], knotcoords[, 1])), 
             max(c(datacoords[,1], knotcoords[, 1])), length = grid.dim[1])
    Ny = seq(min(c(datacoords[, 2], knotcoords[, 2])), 
             max(c(datacoords[,2], knotcoords[, 2])), length = grid.dim[2])
    xygrid <- expand.grid(x = Nx, y = Ny)
    if (type == "B") {
      if (!is.null(names(datacoords))) {
        knotcoords <- data.frame(knotcoords)
        names(knotcoords) <- names(datacoords)
      }
      datacoords2 = rbind(datacoords, knotcoords)
      geodistsoutput <- getGeoDist(xygrid = xygrid, polys = polys, 
                                   datalocations = datacoords2, plot.transition = plot.transition)
      d2k <- geodistsoutput$distance[(1:nrow(datacoords)), 
                                     ((nrow(datacoords) + 1):nrow(datacoords2))]
      if (knotmat == T) {
        knotDist = geodistsoutput$distance[(nrow(datacoords) + 
                                              1):nrow(datacoords2), (nrow(datacoords) + 1):nrow(datacoords2)]
        return(list(dataDist = d2k, knotDist = knotDist))
      }
      else {
        return(list(dataDist = d2k))
      }
    }
    if (type == "A") {
      geodistsoutput <- getGeoDist(xygrid = xygrid, polys = polys, 
                                   datalocations = datacoords, plot.transition = plot.transition)
      d2k <- geodistsoutput$distance[, attr(knotcoords, 
                                            "points.selected")]
      if (knotmat == T) {
        knotDist = geodistsoutput$distance[attr(knotcoords, 
                                                "points.selected"), attr(knotcoords, "points.selected")]
        return(list(dataDist = d2k, knotDist = knotDist))
      }
      else {
        return(list(dataDist = d2k))
      }
    }
  }
}


