M<- 1000
# mom
# This is to get the mean and variance under Kimura
Kimura.model <- function(x0,t,Q){
  e <- rep(1,4)
  E <- e%*%t(e)
  alpha <- Q[1,2]
  beta <- Q[1,3]
  a1 <- 1/4
  a2 <- 1/4
  a3 <- 1/2
  a <- c(a1,a2,a3)
  b1 <- 0
  b2 <- 4*beta
  b3 <- 2*(alpha+beta)
  b <- c(b1,b2,b3)
  A1 <- E
  A2 <- matrix(c(1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1),nrow=4)
  A3 <- matrix(c(1,-1,0,0,-1,1,0,0,0,0,1,-1,0,0,-1,1),nrow=4)
  A <- list(A1,A2,A3)
  eQt <- a1*exp(-b1*t)*A1 + a2*exp(-b2*t)*A2 + a3*exp(-b3*t)*A3
  Ext <- t(x0)%*%eQt
  gijk <- function(i,j,k){
    if(1+b[i]+b[k] == b[j]){
      return(t*exp(-b[j]*t))
    }else{
      mygijk <- (exp(-b[j]*t)-exp(-(1+b[i]+b[k])*t))/(1+b[i]+b[k]-b[j])
      return(mygijk)
    }
  }
  sum <- 0
  for (i in 1:3){
    for(j in 1:3){
      for (k in 1:3){
        sum <- sum + a[i]*a[j]*a[k]*gijk(i,j,k)*A[[i]]%*%diag(c(t(x0)%*%A[[j]]))%*%A[[k]]
      }
    }
  }
  Varxt <- sum - t(eQt)%*%x0%*%t(x0)%*%eQt*(1-exp(-t))
  mylist <- NULL
  mylist$E <- Ext
  mylist$V <- Varxt
  return(mylist)
}

# This is to estimate the parameters of hierarchical beta
Hierarch.Beta.approx <- function(Ext,Varxt){
  muomega <- Ext[1]+Ext[2]
  mueta1 <- Ext[1]/muomega
  mueta2 <- Ext[3]/(1-muomega)
  Varomega <- Varxt[1,1]+Varxt[2,2]+2*Varxt[1,2]
  
  #Vareta1 <- (Varxt[1,1]-Varomega*mueta1^2)/(Varomega+muomega^2)
  #Vareta2 <- (Varxt[3,3]-Varomega*mueta2^2)/(Varomega+(1-muomega)^2)
  Vareta1 <- (Varxt[1,1] +Varxt[2,2] - (mueta1*mueta1+(1-mueta1)*(1-mueta1))*Varomega)/(2*(muomega*muomega+Varomega))
  Vareta2 <- (Varxt[3,3] +Varxt[4,4] - (mueta2*mueta2+(1-mueta2)*(1-mueta2))*Varomega)/(2*((1-muomega)*(1-muomega)+Varomega))
  phiomega <- muomega*(1-muomega)/(Varomega)-1
  phieta1 <- mueta1*(1-mueta1)/Vareta1-1
  phieta2 <- mueta2*(1-mueta2)/Vareta2-1
  return(list(c(muomega,Varomega),c(mueta1,Vareta1),c(mueta2,Vareta2)))
  
}


integration_simulationT <- function(parameter.list,nbins,M){
  muomega <- parameter.list[[1]][1]
  Varomega <- parameter.list[[1]][2]
  mueta1 <- parameter.list[[2]][1]
  Vareta1 <- parameter.list[[2]][2]
  mueta2 <- parameter.list[[3]][1]
  Vareta2 <- parameter.list[[3]][2]
  phiomega <- abs(muomega*(1-muomega)/(Varomega)-1)
  phieta1 <- abs(mueta1*(1-mueta1)/Vareta1-1)
  phieta2 <- abs(mueta2*(1-mueta2)/Vareta2-1)
  
  alphaAG <- muomega*phiomega
  betaAG <- phiomega*(1-muomega)
  alphaA <- mueta1*phieta1
  betaA <- phieta1*(1-mueta1)
  alphaC <- mueta2*phieta2
  betaC <- phieta2*(1-mueta2)
  
  # probability of spike = 0 
  
  prob_spike0 <- function(alpha,beta,nConstant = 1){
    
    return((1/nConstant) * pbeta(1/(M),alpha,beta))
    
  }
  
  # probability of spike = 1
  
  prob_spike1 <- function(alpha,beta,nConstant = 1){
    
    return((1/nConstant) * (1-pbeta(1-1/(M),alpha,beta)))
    
  }
  
  
  ibeta <- function(alpha,beta,left,right){
    return(pbeta(right,alpha,beta) - pbeta(left,alpha,beta))
  }
  
  # Given A + G = 0
  hierarchical_A_Plus_G_0T <- function(alleleCi){
    left <- alleleCi[1]
    right <- alleleCi[2]
    
    return(ibeta(alphaC,betaC,left,right) * prob_spike0(alphaAG,betaAG))
    
  }
  
  # Given A + G = 1
  hierarchical_C_Plus_T_0T <- function(alleleAi){
    left <- alleleAi[1]
    right <- alleleAi[2]
    
    return(ibeta(alphaA,betaA,left,right) * prob_spike1(alphaAG,betaAG))
    
  }
  
  # Given A + C = 0
  hierarchical_A_Plus_C_0T <- function(alleleGi){
    left <- alleleGi[1]
    right <- alleleGi[2]
    
    alleleG <- (left+right)/2
    
    return(ibeta(alphaAG,betaAG,left,right) * prob_spike0(alphaA,betaA,alleleG) * prob_spike0(alphaC,betaC,(1-alleleG)))
    
  }
  
  # Given A + T = 0
  hierarchical_A_Plus_T_0T <- function(alleleGi){
    left <- alleleGi[1]
    right <- alleleGi[2]
    
    alleleG <- (left+right)/2
    
    return(ibeta(alphaAG,betaAG,left,right) * prob_spike0(alphaA,betaA,alleleG) * prob_spike1(alphaC,betaC,(1-alleleG)))
    
  }
  
  # Given G + C = 0
  hierarchical_G_Plus_C_0T <- function(alleleAi){
    left <- alleleAi[1]
    right <- alleleAi[2]
    
    alleleA <- (left+right)/2
    
    return(ibeta(alphaAG,betaAG,left,right) * prob_spike1(alphaA,betaA,alleleA) * prob_spike0(alphaC,betaC,1-alleleA))
    
    
  }
  
  # Given G + T = 0T
  hierarchical_G_Plus_T_0T <- function(alleleAi){
    left <- alleleAi[1]
    right <- alleleAi[2]
    
    alleleA <- (left+right)/2
    
    return(ibeta(alphaAG,betaAG,left,right) * prob_spike1(alphaA,betaA,alleleA) * prob_spike1(alphaC,betaC,1-alleleA))
    
  }
  
  prob_spike_vector <- c()
  
  # For one allele equals to 1
  # A = 1, G = 0, C = 0, T = 0
  spike_1 <- function(alphaAG,betaAG,alphaA,betaA){
    return(prob_spike1(alphaAG,betaAG) * prob_spike1(alphaA,betaA))
  }
  
  # A = 0, G = 1, C = 0, T = 0
  spike_2 <- function(alphaAG,betaAG,alphaA,betaA){
    return(prob_spike1(alphaAG,betaAG) * prob_spike0(alphaA,betaA))
  }
  
  # A = 0, G = 0, C = 1, T = 0
  spike_3 <- function(alphaAG,betaAG,alphaC,betaC){
    return(prob_spike0(alphaAG,betaAG) * prob_spike1(alphaC,betaC))
  }
  
  # A = 0, G = 0, C = 0, T = 1
  spike_4 <- function(alphaAG,betaAG,alphaC,betaC){
    return(prob_spike0(alphaAG,betaAG) * prob_spike0(alphaC,betaC))
  }
  
  
  # for four spikes
  # A = 1
  prob_spike_vector <- c(prob_spike_vector,spike_1(alphaAG,betaAG,alphaA,betaA))
  
  # G = 1
  prob_spike_vector <- c(prob_spike_vector,spike_2(alphaAG,betaAG,alphaA,betaA))
  
  # C = 1
  prob_spike_vector <- c(prob_spike_vector,spike_3(alphaAG,betaAG,alphaC,betaC))
  
  # T = 1
  prob_spike_vector <- c(prob_spike_vector,spike_4(alphaAG,betaAG,alphaC,betaC))
  
  # Then we calculate the probability for two allele cases
  CL_bin <- CL_bins(nbins,M)
  
  bin_points <- c()
  
  for (i in 1: nrow(CL_bin)){
    bin_points <- c(bin_points,(CL_bin[i,1]+CL_bin[i,2])/2)
  }
  
  interval_length <- CL_bin[,3]
  
  bin_vector <- c()
  # Given A + G = 0
  count = 1
  for (i in bin_points){
    alleleCi <- CL_bin[count,c(1,2)]
    bin_vector <- c(bin_vector,hierarchical_A_Plus_G_0T(alleleCi))
    count = count + 1
  }
  
  # Given A + G = 1
  count = 1
  for (i in bin_points){
    alleles <- CL_bin[count,c(1,2)]
    bin_vector <- c(bin_vector,hierarchical_C_Plus_T_0T(alleles))
    count = count + 1
  }
  
  # Given A + C = 0
  count = 1
  for (i in bin_points){
    alleles <- CL_bin[count,c(1,2)]
    bin_vector <- c(bin_vector,hierarchical_A_Plus_C_0T(alleles))
    count = count + 1
  }
  
  # Given A + T = 0
  count = 1
  for (i in bin_points){
    alleles <- CL_bin[count,c(1,2)]
    bin_vector <- c(bin_vector,hierarchical_A_Plus_T_0T(alleles))
    count = count + 1
  }
  
  # Given G + C = 0
  count = 1
  for (i in bin_points){
    alleles <- CL_bin[count,c(1,2)]
    bin_vector <- c(bin_vector,hierarchical_G_Plus_C_0T(alleles))
    count = count + 1
  }
  
  # Given G + T = 0
  count = 1
  for (i in bin_points){
    alleles <- CL_bin[count,c(1,2)]
    bin_vector <- c(bin_vector,hierarchical_G_Plus_T_0T(alleles))
    count = count + 1
  }
  
  mylist <- list()
  
  mylist[[1]] <- c(normalize(bin_vector))
  
  mylist[[2]] <- c(prob_spike_vector,normalize(bin_vector) *(1-sum(prob_spike_vector)))
  
  names(mylist) <- c("pure_compare_bins","consider_spike")
  
  return(mylist)
  
}

wonderful_simulation <- function(parameter.list,nbins,M){
  muomega <- parameter.list[[1]][1]
  Varomega <- parameter.list[[1]][2]
  mueta1 <- parameter.list[[2]][1]
  Vareta1 <- parameter.list[[2]][2]
  mueta2 <- parameter.list[[3]][1]
  Vareta2 <- parameter.list[[3]][2]
  phiomega <- abs(muomega*(1-muomega)/(Varomega)-1)
  phieta1 <- abs(mueta1*(1-mueta1)/Vareta1-1)
  phieta2 <- abs(mueta2*(1-mueta2)/Vareta2-1)
  
  alphaAG <- muomega*phiomega
  betaAG <- phiomega*(1-muomega)
  alphaA <- mueta1*phieta1
  betaA <- phieta1*(1-mueta1)
  alphaC <- mueta2*phieta2
  betaC <- phieta2*(1-mueta2)
  
  # probability of spike = 0 
  
  prob_spike0 <- function(alpha,beta,nConstant = 1){
    
    if (alpha > 10^3 && beta > 10^3){
      return(0)
    } else{
      return((1/nConstant)*M^(-alpha)/(alpha*beta(alpha,beta)))
    }
  }
  
  # probability of spike = 1
  
  prob_spike1 <- function(alpha,beta,nConstant = 1){
    
    if (alpha > 10^3 && beta > 10^3){
      return(0)
    } else{
      return((1/nConstant)*M^(-beta)/(beta*beta(alpha,beta)))
    }
  }
  
  
  # For one allele equals to 1
  # A = 1, G = 0, C = 0, T = 0
  spike_1 <- function(alphaAG,betaAG,alphaA,betaA){
    return(prob_spike1(alphaAG,betaAG) * prob_spike1(alphaA,betaA))
  }
  
  # A = 0, G = 1, C = 0, T = 0
  spike_2 <- function(alphaAG,betaAG,alphaA,betaA){
    return(prob_spike1(alphaAG,betaAG) * prob_spike0(alphaA,betaA))
  }
  
  # A = 0, G = 0, C = 1, T = 0
  spike_3 <- function(alphaAG,betaAG,alphaC,betaC){
    return(prob_spike0(alphaAG,betaAG) * prob_spike1(alphaC,betaC))
  }
  
  # A = 0, G = 0, C = 0, T = 1
  spike_4 <- function(alphaAG,betaAG,alphaC,betaC){
    return(prob_spike0(alphaAG,betaAG) * prob_spike0(alphaC,betaC))
  }
  
  prob_spike_vector <- c()
  # for four spikes
  # A = 1
  prob_spike_vector <- c(prob_spike_vector,spike_1(alphaAG,betaAG,alphaA,betaA))
  
  # G = 1
  prob_spike_vector <- c(prob_spike_vector,spike_2(alphaAG,betaAG,alphaA,betaA))
  
  # C = 1
  prob_spike_vector <- c(prob_spike_vector,spike_3(alphaAG,betaAG,alphaC,betaC))
  
  # T = 1
  prob_spike_vector <- c(prob_spike_vector,spike_4(alphaAG,betaAG,alphaC,betaC))
  
  
  
  # Given A + G = 0
  hierarchical_A_Plus_G_0 <- function(alleles,alphaC,betaC,alphaAG,betaAG){
    alleleC <- alleles[3]
    if (alleleC != 0 && alleleC != 1){
      return(dbeta(alleleC,alphaC,betaC) * prob_spike0(alphaAG,betaAG))
    }
  }
  
  # Given A + G = 1
  hierarchical_C_Plus_T_0 <- function(alleles,alphaA,betaA,alphaAG,betaAG){
    alleleA <- alleles[1]
    if (alleleA != 0 && alleleA != 1){
      return(dbeta(alleleA,alphaA,betaA) * prob_spike1(alphaAG,betaAG))
    } 
  }
  
  # Given A + C = 0
  hierarchical_A_Plus_C_0 <- function(alleles,alphaA,betaA,alphaC,betaC,alphaAG,betaAG){
    alleleG <- alleles[2]
    if (alleleG != 0 && alleleG != 1){
      return(dbeta(alleleG,alphaAG,betaAG) * prob_spike0(alphaA,betaA,alleleG) * prob_spike0(alphaC,betaC,(1-alleleG)))
    } 
  }
  
  # Given A + T = 0
  hierarchical_A_Plus_T_0 <- function(alleles,alphaA,betaA,alphaC,betaC,alphaAG,betaAG){
    alleleG <- alleles[2]
    if (alleleG != 0 && alleleG != 1){
      return(dbeta(alleleG,alphaAG,betaAG) * prob_spike0(alphaA,betaA,alleleG) * prob_spike1(alphaC,betaC,(1-alleleG)))
    } 
  }
  
  # Given G + C = 0
  hierarchical_G_Plus_C_0 <- function(alleles,alphaA,betaA,alphaC,betaC,alphaAG,betaAG){
    alleleA <- alleles[1]
    if (alleleA != 0 && alleleA != 1){
      return(dbeta(alleleA,alphaAG,betaAG) * prob_spike1(alphaA,betaA,alleleA) * prob_spike0(alphaC,betaC,1-alleleA))
    } 
    
  }
  
  # Given G + T = 0
  hierarchical_G_Plus_T_0 <- function(alleles,alphaA,betaA,alphaC,betaC,alphaAG,betaAG){
    alleleA <- alleles[1]
    if (alleleA != 0 && alleleA != 1){
      return(dbeta(alleleA,alphaAG,betaAG) * prob_spike1(alphaA,betaA,alleleA) * prob_spike1(alphaC,betaC,1-alleleA))
    } 
  }
  
  # Then we calculate the probability for two allele cases
  
  CL_bin <- CL_bins(nbins,M)
  
  bin_points <- c()
  
  for (i in 1: nrow(CL_bin)){
    bin_points <- c(bin_points,(CL_bin[i,1]+CL_bin[i,2])/2)
  }
  
  interval_length <- CL_bin[,3]
  
  bin_vector <- c()
  # Given A + G = 0
  count = 1
  for (i in bin_points){
    alleles <- c(0,0,i,1-i)
    bin_vector <- c(bin_vector,hierarchical_A_Plus_G_0(alleles,alphaC,betaC,alphaAG,betaAG)*interval_length[count])
    count = count + 1
  }
  
  # Given A + G = 1
  count = 1
  for (i in bin_points){
    alleles <- c(i,1-i,0,0)
    bin_vector <- c(bin_vector,hierarchical_C_Plus_T_0(alleles,alphaA,betaA,alphaAG,betaAG)*interval_length[count])
    count = count + 1
  }
  
  # Given A + C = 0
  count = 1
  for (i in bin_points){
    alleles <- c(0,i,0,1-i)
    bin_vector <- c(bin_vector,hierarchical_A_Plus_C_0(alleles,alphaA,betaA,alphaC,betaC,alphaAG,betaAG)*interval_length[count])
    count = count + 1
  }
  
  # Given A + T = 0
  count = 1
  for (i in bin_points){
    alleles <- c(0,i,1-i,0)
    bin_vector <- c(bin_vector,hierarchical_A_Plus_T_0(alleles,alphaA,betaA,alphaC,betaC,alphaAG,betaAG)*interval_length[count])
    count = count + 1
  }
  
  # Given G + C = 0
  count = 1
  for (i in bin_points){
    alleles <- c(i,0,0,1-i)
    bin_vector <- c(bin_vector,hierarchical_G_Plus_C_0(alleles,alphaA,betaA,alphaC,betaC,alphaAG,betaAG)*interval_length[count])
    count = count + 1
  }
  
  # Given G + T = 0
  count = 1
  for (i in bin_points){
    alleles <- c(i,0,1-i,0)
    bin_vector <- c(bin_vector,hierarchical_G_Plus_T_0(alleles,alphaA,betaA,alphaC,betaC,alphaAG,betaAG)*interval_length[count])
    count = count + 1
  }
  
  mylist <- list()
  
  mylist[[1]] <- c(normalize(bin_vector))
  
  mylist[[2]] <- c(prob_spike_vector,normalize(bin_vector) *(1-sum(prob_spike_vector)))
  
  names(mylist) <- c("pure_compare_bins","consider_spike")
  
  return(mylist)
  
}


normalize <- function(v){
  return(v/sum(v))
}


filter_and_normalize <- function(v,N){
  for (i in 1:length(v)){
    if (v[i] < 1/N){
      v[i] <- 0
    }
  }
  return(normalize(v))
}


hdd <- function(x,y){
  sqrt(sum((sqrt(x) -sqrt(y))^2)) / (2^-0.5)
}



all_simulation <- function(x0,t,Q,nbins,n.samples){
  N <- M
  EV <- Kimura.model(x0,t,Q*N)
  
  parameter.list <- Hierarch.Beta.approx(EV$E,EV$V)
  
  theoretical <- wonderful_simulation(parameter.list,nbins,N)
  
  integration <- integration_simulationT(parameter.list,nbins,N)
  
  return(abs(theoretical$consider_spike-integration$consider_spike)/integration$consider_spike)
}

#### Chebyshev bins
CL_bins <- function(nsize,M){
  mymatrix <- matrix(nrow = nsize+1,ncol=3)
  myc = c()
  for (n in 1:nsize){
    myc <- c(myc,cos((2*n-1)*pi / (2*nsize)))
  }
  b = 1 - 1/M
  a = 1/M
  myc <- rev((b-a)/2 * myc + (a+b)/2)
  
  i = 1/M
  for (n in 1:nsize){
    mymatrix[n,1] <- i
    mymatrix[n,2] <- myc[n] 
    i = myc[n]
  }
  mymatrix[n+1,1] <- i
  mymatrix[n+1,2] <- 1-1/M
  mymatrix[,3] <- mymatrix[,2]-mymatrix[,1]
  return(mymatrix)
}


# equal size bins
CL_bins <- function(nsize,M){
  mymatrix <- matrix(nrow = nsize+1,ncol=3)
  
  myc <- seq(1/M,1-1/M,(1-1/M-1/M)/(nsize+1))
  myc <- myc[-c(1,nsize+2)]
  i = 1/M
  for (n in 1:nsize){
    mymatrix[n,1] <- i
    mymatrix[n,2] <- myc[n] 
    i = myc[n]
  }
  mymatrix[n+1,1] <- i
  mymatrix[n+1,2] <- 1-1/M
  mymatrix[,3] <- mymatrix[,2]-mymatrix[,1]
  return(mymatrix)
}





###################################
our_sim <- function(nbins){
  x0 <- c(0.5,0.5,0,0)
  ts <- seq(1,100,1)
  n.samples <- 100000
  heat.map1 <- matrix(nrow=6*(nbins+1)+4,ncol=length(ts))
  k = 2
  u = 10^-8
  Q <- matrix(c(-(k+2)*u, k * u, u, u,
                k * u, -(k+2)*u, u, u,
                u, u, -(k+2)*u, k * u,
                u, u, k * u, -(k+2)*u),nrow=4)
  Q <- Q*M
  nn = 1
  for (t in ts){
    heat.map1[,nn] <- all_simulation(x0,t,Q,nbins,n.samples)
    nn = nn + 1
  }
  return(heat.map1)
}

our_sim_multiple <- function(nbins){
  x0 <- c(0.5,0.5,0,0)
  ts <- seq(0.1,10,0.1)
  n.samples <- 100000
  
  k = 2
  u = 10^-8
  Q <- matrix(c(-(k+2)*u, k * u, u, u,
                k * u, -(k+2)*u, u, u,
                u, u, -(k+2)*u, k * u,
                u, u, k * u, -(k+2)*u),nrow=4)
  heat.map2 <- matrix(nrow=length(nbins),ncol=length(ts))
  bb <- 1
  for (bin in nbins){
    nn = 1
    heat.map1 <- matrix(nrow=6*(bin+1)+4,ncol=length(ts))
    for (t in ts){
      heat.map1[,nn] <- all_simulation(x0,t,Q,bin,n.samples)
      nn = nn + 1
    }
    heat.map2[bb,] <- colSums(heat.map1)/dim(heat.map1)[1]
    bb <- bb + 1
  }
  
  return(heat.map2)
}


# equal size bins
CL_bins <- function(nsize,M){
  mymatrix <- matrix(nrow = nsize+1,ncol=3)
  
  myc <- seq(1/M,1-1/M,(1-1/M-1/M)/(nsize+1))
  myc <- myc[-c(1,nsize+2)]
  i = 1/M
  for (n in 1:nsize){
    mymatrix[n,1] <- i
    mymatrix[n,2] <- myc[n] 
    i = myc[n]
  }
  mymatrix[n+1,1] <- i
  mymatrix[n+1,2] <- 1-1/M
  mymatrix[,3] <- mymatrix[,2]-mymatrix[,1]
  return(mymatrix)
}

heat.map1<-our_sim_multiple(seq(2,100,1))
#### Chebyshev bins
CL_bins <- function(nsize,M){
  mymatrix <- matrix(nrow = nsize+1,ncol=3)
  myc = c()
  for (n in 1:nsize){
    myc <- c(myc,cos((2*n-1)*pi / (2*nsize)))
  }
  b = 1 - 1/M
  a = 1/M
  myc <- rev((b-a)/2 * myc + (a+b)/2)
  
  i = 1/M
  for (n in 1:nsize){
    mymatrix[n,1] <- i
    mymatrix[n,2] <- myc[n] 
    i = myc[n]
  }
  mymatrix[n+1,1] <- i
  mymatrix[n+1,2] <- 1-1/M
  mymatrix[,3] <- mymatrix[,2]-mymatrix[,1]
  return(mymatrix)
}
heat.map2<-our_sim_multiple(seq(2,100,1))




##### equal size 
d = data.frame(x=rep(seq(0.1,10,0.1),99), 
               y=rep(seq(2,100,1),each=100), 
               z=c(t(heat.map1)))
levelplot(z~x*y,data=d,col.regions=colorRampPalette(c("yellow","orange","red", "black")),at=c(seq(0,0.8,length=30)),xlab="Time",ylab = "Number of Bins",main="Discretization by Equally Spaced Bin Points")
levelplot(z~x*y,data=d,col.regions=colorRampPalette(c("blue","light blue","light green", "yellow","orange","red","black")),at=c(seq(0,0.8,length=400)),xlab="Time",ylab = "Number of Bins",main="Errors Introduced by Using Uniformly Spaced Bin Points (Alternative Formulas)")


##### chebyshev
d = data.frame(x=rep(seq(0.1,10,0.1),99), 
               y=rep(seq(2,100,1),each=100), 
               z=c(t(heat.map2)))
levelplot(z~x*y,data=d,col.regions=colorRampPalette(c("yellow","orange","red", "black")),at=c(seq(0,0.8,length=20)),xlab="Time",ylab = "Number of Bins",main="Discretization by Chebyshev-Lobatto Points")

levelplot(z~x*y,data=d,col.regions=colorRampPalette(c("blue","light blue","light green", "yellow","orange","red","black")),at=c(seq(0,0.8,length=400)),xlab="Time",ylab = "Number of Bins",main="Errors Introduced by Using Chebyshev-Lobatto Points (Alternative Formulas)")

########### chebyshev one

second <- our_sim(15)
levelplot(t(second),col.regions=colorRampPalette(c("yellow","orange","red", "black")),at=c(seq(0,0.8,length=100)),xlab="Time",ylab="Bin Index",main = "Chebyshev-Lobatto Points Discretization for Bin Size = 15")



