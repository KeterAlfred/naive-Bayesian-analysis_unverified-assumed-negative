

#
#
# Model for performing naive evaluation of diagnostic tests
#
#=========================================================




#Load required packages
require(R2jags)

# Create the design matrix
K=6   #Actual number of tests: here any TB symp., radiologist concl., CAD4TB v7, 
#     Culture, Xpert, CRS 
DH = matrix( NA, nrow = 2^K, ncol = K  )  #

DH[,1] <- rep(c(1, 0),2^K/2)
DH[,2] <- rep(c(rep(1,2), rep(0,2)),2^K/4)
DH[,3] <- rep(c(rep(1,4), rep(0,4)),2^K/8)
DH[,4] <- rep(c(rep(1,8), rep(0,8)),2^K/16)
DH[,5] <- rep(c(rep(1,16), rep(0,16)),2^K/32)
DH[,6] <- rep(c(rep(1,32), rep(0,32)),2^K/64)
# DH[,7] <- rep(c(rep(1,64), rep(0,64)),2^K/128)
# DH[,8] <- rep(c(rep(1,128), rep(0,128)),2^K/256)


model.naive  <- function(){
  
  #Likelihood
  Freqobs[1:J] ~ dmulti(pp[1:J],N)
  
  for(j in 1 : J){
    # d[j] ~ dbin(prev,N)
    pp[j] = (pi*prod(p1[j,1:K]^y[j,1:K]*(1-p1[j,1:K])^(1-y[j,1:K])))^d[j] * 
      ((1-pi)*prod(p0[j,1:K]^y[j,1:K]*(1-p0[j,1:K])^(1-y[j,1:K])))^(1-d[j])
    
    for(k in 1 : K){
      p1[j,k] = phi( a[k] + b[k]*d[j])
      p0[j,k] = phi( a[k] )
    }
  }
  for(k in 1 : K){
    a[k] ~ dnorm(0,1)
    b[k] ~ dnorm(0,1)
    se[k] = phi( a[k] + b[k])
    sp[k] = 1 - phi( a[k])
  }
  theta ~ dnorm(0,1)
  pi = phi(theta)
}




# Determine prevalence in different groups
#==============================================================
model.naive.groups  <- function(){
  
  #Likelihood
  Freqobs[1:J] ~ dmulti(pp[1:J],N)
  nv[1:G] ~ dmulti(pv[1:G],N)
  
  for(j in 1 : J){

    pp[j] =      (pij[j]*prod(p1[j,1:K]^y[j,1:K]*(1-p1[j,1:K])^(1-y[j,1:K])))^d[j] * 
             ((1-pij[j])*prod(p0[j,1:K]^y[j,1:K]*(1-p0[j,1:K])^(1-y[j,1:K])))^(1-d[j])
      

    
    for(k in 1 : K){
      p1[j,k] = phi( a[k,1] + a[k,4]*d[j] + a[k,2]*(y[j,7]==2) + 
                                            a[k,3]*(y[j,7]==3) ) # + 
                                            # a[k,4]*(y[j,7]==4) )
      
      p0[j,k] = phi( a[k,1] +               a[k,2]*(y[j,7]==2) + 
                                            a[k,3]*(y[j,7]==3) ) # + 
                                            # a[k,4]*(y[j,7]==4)  )
    }
    
    pij[j] = phi( theta[1] + theta[2]*(y[j,7]==2) + 
                             theta[3]*(y[j,7]==3) ) # + 
                             # theta[4]*(y[j,7]==4) )
  }
# for( g in 1 : 5 ){
  for(k in 1 : K){
    a[k,1] ~ dnorm(0,1)
    a[k,2] ~ dnorm(0,1)
    a[k,3] ~ dnorm(0,1)
    a[k,4] ~ dnorm(0,1)
    # a[k,5] ~ dnorm(0,1)
    se[k,1] =     phi( a[k,1] + a[k,4])
    sp[k,1] = 1 - phi( a[k,1])
    
    se[k,2] =     phi( a[k,2] + a[k,4] + a[k,2])
    sp[k,2] = 1 - phi( a[k,1] + a[k,2])
    
    se[k,3] =     phi( a[k,2] + a[k,4] + a[k,3])
    sp[k,3] = 1 - phi( a[k,1] + a[k,3])
    
    # se[k,4] =     phi( a[k,2] + a[k,5] + a[k,4])
    # sp[k,4] = 1 - phi( a[k,1] + a[k,4])
    
    se[k,4] = se[k,1]*pv[1] + se[k,2]*pv[2] + se[k,3]*pv[3]
    sp[k,4] = sp[k,1]*pv[1] + sp[k,2]*pv[2] + sp[k,3]*pv[3]
  }
# }
  


  for(g in 1:G){
    theta[g] ~ dnorm(0,1)
  }
  pi[1] = phi(theta[1])
  pi[2] = phi(theta[1] + theta[2])
  pi[3] = phi(theta[1] + theta[3])
  # pi[4] = phi(theta[1] + theta[4])
  # pi[5] = pi[1]*pv[1] + pi[2]*pv[2] + pi[3]*pv[3] + pi[4]*pv[4]
  pi[4] = pi[1]*pv[1] + pi[2]*pv[2] + pi[3]*pv[3]
  
  pv[1:G] ~ ddirich(c(1,1,1))
}


#===============================================================




model.naive.v  <- function(){
  
  #Likelihood
  Freqobs[1:J] ~ dmulti(pp[1:J],N)
  
  for(j in 1 : J){
    # d[j] ~ dbin(prev,N)
    pp[j] = (pv*(prevj[j]*prod(p11[j,1:K]^y[j,1:K]*(1-p11[j,1:K])^(1-y[j,1:K])))^d[j] * 
               ((1-prevj[j])*prod(p10[j,1:K]^y[j,1:K]*(1-p10[j,1:K])^(1-y[j,1:K])))^(1-d[j]))^y[j,9]*
      ((1-pv)*(prevj[j]*prod(p01[j,1:K]^y[j,1:K]*(1-p01[j,1:K])^(1-y[j,1:K])))^d[j] * 
         ((1-prevj[j])*prod(p00[j,1:K]^y[j,1:K]*(1-p00[j,1:K])^(1-y[j,1:K])))^(1-d[j]))^(1-y[j,9])
    
    for(k in 1 : K){
      p11[j,k] = phi( a[1,k] + a[2,k]*d[j] + a[3,k]*y[j,9])
      p10[j,k] = phi( a[1,k] + a[3,k]*y[j,9])
      p01[j,k] = phi( a[1,k] + a[2,k]*d[j])
      p00[j,k] = phi( a[1,k] )
    }
    prevj[j] = phi(theta[1] + theta[2]*y[j,9])
  }
  for(k in 1 : K){
    a[1,k] ~ dnorm(0,1)
    a[2,k] ~ dnorm(0,1)
    a[3,k] ~ dnorm(0,1)
    se[k,1] = phi( a[1,k] + a[2,k] + a[3,k])
    sp[k,1] = 1 - phi( a[1,k] + a[3,k])
    se[k,2] = phi( a[1,k] + a[2,k])
    sp[k,2] = 1 - phi( a[1,k])
    se[k,3] = (se[k,1]*prev[1]*pv + se[k,2]*prev[2]*(1-pv))/(prev[1]*pv + prev[2]*(1-pv))
    sp[k,3] = 1 - ( (1-sp[k,1])*prev[1]*pv + (1-sp[k,2])*prev[2]*(1-pv) )/(prev[1]*pv + prev[2]*(1-pv))
  }
  theta[1] ~ dnorm(0,1)
  theta[2] ~ dnorm(0,1)

  mu ~ dnorm(0,1)
  pv= phi( mu )
  
  prev[1] = phi(theta[1] + theta[2])
  prev[2] = phi( theta[1] )
  prev[3] = prev[2]*(1-pv) + prev[1]*pv
}





#=================================================================

K=6   #Actual number of tests: here any TB symp., radiologist concl., CAD4TB v7, Culture, Xpert, CRS 
DH0 = matrix( NA, nrow = 2^K, ncol = K  )  #

DH0[,1] <- rep(c(1, 0),2^K/2)
DH0[,2] <- rep(c(rep(1,2), rep(0,2)),2^K/4)
DH0[,3] <- rep(c(rep(1,4), rep(0,4)),2^K/8)
DH0[,4] <- rep(c(rep(1,8), rep(0,8)),2^K/16)
DH0[,5] <- rep(c(rep(1,16), rep(0,16)),2^K/32)
DH0[,6] <- rep(c(rep(1,32), rep(0,32)),2^K/64)
# DH0[,7] <- rep(c(rep(1,64), rep(0,64)),2^K/128)
# DH0[,8] <- rep(c(rep(1,128), rep(0,128)),2^K/256)



model.naive.imp  <- function(){
  
  #Likelihood
  Freqobs[1:J] ~ dmulti(pp[1:J],N)
  
  for(j in 1 : J){
    
    pp[j] = (pred_prev*prod(p1[j,1:K]^y[j,1:K]*(1-p1[j,1:K])^(1-y[j,1:K])))^d[j] * 
      ((1-pred_prev)*prod(p0[j,1:K]^y[j,1:K]*(1-p0[j,1:K])^(1-y[j,1:K])))^(1-d[j])
    
    for(k in 1 : 5){
      p1[j,k] = phi( a[k] + b[k]*d[j])
      p0[j,k] = phi( a[k] )
    }
  
  for(k in 1 : 2){
    p1[j,k+5] = phi( c[1,k] + c[2,k]*y[j,1] + c[3,k]*y[j,3] + c[4,k]*y[j,4] + c[5,k]*y[j,5]  + c[6,k]*d[j] )
    p0[j,k+5] = phi( c[1,k] ) # + c[2,k]*y[j,1] + c[3,k]*y[j,3] + c[4,k]*y[j,4] + c[5,k]*y[j,5] )
  }
  
  # y[j,5] = round( pp1[j,1]*Freqobs )
  # y[j,6] = round( pp1[j,2]*Freqobs )
  
  }
  
  for(j in 1 : J0){
  y0[j,5] ~ dbin(q1[j,1],N0)
  y0[j,6] ~ dbin(q1[j,2],N0)
  for(k in 1 : 2){
    q1[j,k] = phi( c[1,k] + c[2,k]*y0[j,1] + c[3,k]*y0[j,2] + c[4,k]*y0[j,3] + c[5,k]*y0[j,4] + c[6,k]*y0[j,7] )
    q0[j,k] = phi( c[1,k] ) #+ c[2,k]*y0[j,1] + c[3,k]*y0[j,2] + c[4,k]*y0[j,3] + c[5,k]*y0[j,4] )
    
  }
}
  for(k in 1 : 6){
    c[k,1] ~ dnorm(0,1)
    c[k,2] ~ dnorm(0,1)
  }
  for(k in 1 : 7){
    a[k] ~ dnorm(0,1)
    b[k] ~ dnorm(0,1)
    se[k] = phi( a[k] + b[k])
    sp[k] = 1 - phi( a[k])
  }
  
  # theta ~ dnorm(0,1)
  # prev = phi(theta)
  
  pred_prev = mean( p1[,6] + (1-p1[,6])*p1[,7] )
}
