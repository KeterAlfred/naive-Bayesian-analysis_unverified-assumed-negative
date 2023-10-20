

# setwd("C:\\Users\\Keter\\OneDrive - ITG\\EDCTP\\PhD\\lca\\Vukuzazi\\revised analysis\\vb_sim_analysis\\naive")

#===========================================================================
#
#
# Analysis of simulated data - assuming unverified are negative
#
#===========================================================================






#Questions
#What do you think is the prevalence of PTB in those not tested with Xpert Ultra & culture?
#What do you think is the sensitivity of Xpert Ultra and culture among those not tested using Xpert Ultra and culture?
#Should we assume perfect specificity for Xpert Ulta and culture? among the unverified?


# Create the design matrix
K=6
DH = matrix( NA, nrow = 2^K, ncol = K  )  #

DH[,1] <- rep(c(1, 0),2^K/2)
DH[,2] <- rep(c(rep(1,2), rep(0,2)),2^K/4)
DH[,3] <- rep(c(rep(1,4), rep(0,4)),2^K/8)
DH[,4] <- rep(c(rep(1,8), rep(0,8)),2^K/16)
DH[,5] <- rep(c(rep(1,16), rep(0,16)),2^K/32)
DH[,6] <- rep(c(rep(1,32), rep(0,32)),2^K/64)




naive_dat <- foreach(m=c(1:M)+j, .inorder = FALSE) %do% {
  
  
  
  
  #Data for the model
  
  dat.tested = as.data.frame(sim_dat_Final1[,,m])
  dat.tested[,"Y50"] <- ifelse( dat.tested[,"Y50"] %in% 9, 0, dat.tested[,"Y50"])
  dat.tested[,"Y60"] <- ifelse( dat.tested[,"Y60"] %in% 9, 0, dat.tested[,"Y60"])
  dat.tested$crs = ifelse( dat.tested$Y50 %in% 1 | dat.tested$Y60 %in% 1, 1, 0 )
  
  
  #First create a variable with 1s 
  ones = rep(1, length.out = dim(dat.tested)[1])
  
  var.nm  <- c("Y1", "Y2", "Y4",  "Y50", "Y60","crs")
  

  
  pattern1 <- aggregate( ones ~   Y1 + Y2  + Y4 + Y50 + Y60 + crs, data = dat.tested,function(x)NROW(x))
  
  pattern1.s = pattern1[ order(-pattern1[,7],-pattern1[,6],-pattern1[,5],-pattern1[,4],-pattern1[,3],-pattern1[,2],-pattern1[,1]), ]
  D.s = DH[ order(-DH[,6],-DH[,5],-DH[,4],-DH[,3],-DH[,2],-DH[,1]), ]
  
  
  colnames(D.s) <- var.nm
  colnames(pattern1.s) <- c(var.nm,"Freqobs")
  
  pred.data = merge(as.data.frame(D.s), as.data.frame(pattern1.s),by = var.nm, all = TRUE)
  pred.data[,"Freqobs"] <- ifelse( is.na( pred.data[,"Freqobs"]),0,pred.data[,"Freqobs"])

  
  jags.datax = pred.data[order(-pred.data[,6],-pred.data[,5],
                               -pred.data[,4],-pred.data[,3],
                               -pred.data[,2],-pred.data[,1]),c(var.nm,"Freqobs")]
  
  jags.data = jags.datax[,c(var.nm,"Freqobs")]
  
  
  

  
  
  #Prepare the collapsed data as a list; this is important for jags to run
  jags.sim.data_cc1 = list( N = sum(jags.data[,"Freqobs"]),
                            Freqobs = jags.data[,"Freqobs"],
                            y = jags.data[,c("Y1", "Y2", "Y4",  "Y50", "Y60","crs")],
                            J = dim(jags.data)[1],
                            K = 5,
                            d = jags.data[,"crs"])
  #Assign names to the objects in the list
  names(jags.sim.data_cc1) <- c("N","Freqobs","y","J","K","d")
  

  
  
  return( jags.sim.data_cc1 )
}

