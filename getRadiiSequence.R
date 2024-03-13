
getRadiiSequence<-function(method=NULL, numberofradii=10, 
                           distMatrix, basis,
                           rvin=NULL,
                           xydata, response, 
                           alpha = 0, vgmmodel = "Sph", 
                           showplots = FALSE, 
                           ...){
  
  if(is.null(distMatrix)) stop("**** No distance matrix provided.****")
  
  ifelse(method %in% c("original", "variogram"),1 , stop("*** method not one of original or variogram"))
  
  distMatrix[which(is.infinite(distMatrix), arr.ind = T)]<-NA
  
  # ~~~~~~~ variogram method ~~~~~~~~~ # 
  
  if(method == "variogram"){
    data <- data.frame(xydata, response)
    names(data) <- c("x", "y", "response")
    
    sp::coordinates(data) = ~x + y
    suppressWarnings({
      vg <- gstat::variogram(response ~ 1, data, alpha= alpha, ...)
      fit.vg <- try(gstat::fit.variogram(vg, model=gstat::vgm(vgmmodel)), silent = TRUE)
      
      if(class(fit.vg)[1]=="try-error"){
        fit.vg <- gstat::fit.variogram(vg, model=gstat::vgm(vgmmodel), fit.method = 6)
      }
    })
    #print(fit.vg)
    if(showplots == TRUE){
      par(mfrow=c(1,2))
      print(plot(vg))
      print(plot(fit.vg, cutoff=max(vg$dist)))
    }
    best.s <- fit.vg$range[2]
    
    
    if(max(distMatrix, na.rm = TRUE)<best.s){
      method = "original"
      print(paste0("Range (", 
                   round(best.s, 2), 
                   ") is greater than maximum distance in distMatrix (",
                   round(max(distMatrix, na.rm = TRUE), 3), 
                   ") so r_seq created using method = 'original'"
      ))
    }else{
      vg.gap <- vg$dist[2] - vg$dist[1]
      nr <- floor(numberofradii/2)
      l <- seq(best.s - (vg.gap * nr), best.s, length=nr)
      
      if(min(l)<0){
        l <- seq(min(abs(l)), best.s, length=nr)
      }
      
      u <- seq(best.s, best.s + (vg.gap * nr), length=nr)
      s_seq <- unique(c(l, u))
      
      
      if(basis=='gaussian'){
        
        r_seq <- 1/((sqrt(2) * s_seq))
        attr(r_seq, "Method") <- "Variogram"
        attr(r_seq, "vg.fit") <- fit.vg
      }
      
      if(basis=='exponential'){
        r_seq <- sqrt(s_seq)
        attr(r_seq, "Method") <- "Variogram"
        attr(r_seq, "vg.fit") <- fit.vg
      }
    } 
  }
  
  # ~~~~~~~~~ Original method ~~~~~~~~~~~~~~~~~
  # from CRESS paper
  
  if(method == "original"){
    if(basis=='gaussian'){
      minDist <- mean(apply(distMatrix,2,min, na.rm=TRUE))
      meanDist <- mean(apply(distMatrix,2,mean, na.rm=TRUE))
      if (is.null(rvin)) {
        rval_max<- sqrt(-log(0.7)/meanDist**2)
        rval_min<- sqrt(-log(0.001)/meanDist**2)
      } else {
        rval_max = rvin[1]
        rval_min <- rvin[2]
      }
      r_seq<- exp(seq(log(rval_min), log(rval_max), length=numberofradii))
    }
    
    if(basis=='exponential'){
      numberofradii = numberofradii+2
      # establish smallest observation-knot distance
      rmin<- sqrt(max(distMatrix, na.rm=TRUE)/21)
      rmax<- sqrt(max(distMatrix, na.rm=TRUE)/3e-7)
      r_seq <- exp(seq(log(rmin), log(rmax), length=numberofradii))[-c(1,numberofradii)]
    }
    attr(r_seq, "Method") <- "Original"
    
  }
  return(r_seq)
  
}


getRadiiChoices.vario<-function(numberofradii=10, xydata, response, 
                                basis, alpha = 0, vgmmodel = "Sph", 
                                showplots = FALSE, distMatrix = NULL, ...){
  
  if(is.null(distMatrix)) stop("**** No distance matrix provided.****")
  
  
  
  return(r_seq) 
}