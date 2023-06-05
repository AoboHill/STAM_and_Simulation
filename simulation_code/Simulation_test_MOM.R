# The simulation function, time scale: t = r/N
library(FNN)
K <- 4
set.seed(1)
true.simulation <- function(x0, U, K, t, N){
  for (i in 1 : (t * N)){
    x0 <- c(rmultinom(1,N,c(x0%*%U)))/N
  }
  return(x0)
}
h.d <- function(X,Y,nn) {
  d <- ncol(X) # number of dimensions, must be the same in X and Y
  n <- nrow(X) # number of samples in X
  m <- nrow(Y) # number of samples in Y
  Ni <- c()
  Mi <- c()
  gx <- c()
  eta <- n/m
  Z <- rbind(X,Y)
  knns <- get.knnx(Z,Y,k=nn)[[1]]
  
  find.freq <- function(data){
    data=as.data.frame(data)
    data$index=c(1:nrow(data))
    d=unique(data[,-ncol(data)])
    d$index1=c(1:nrow(d))
    data1=merge(data,d)
    data1$length=1
    data2=aggregate(.~index1,data1[,c("index1","length")],length)
    result=merge(data1[,-c(ncol(data1))],data2,by="index1")
    result=result[order(result$index),-1]
    return(result)
  }
  
  X.freq <- find.freq(X)
  Y.freq <- find.freq(Y)
  Z.freq <- find.freq(Z)
  
  from.Y <- Z.freq[1:n,6] - X.freq[,6]
  X.freq<-cbind(X.freq,from.Y)
  
  from.X <- Z.freq[seq(n+1,n+m),6] - Y.freq[,6]
  Y.freq <- cbind(Y.freq,from.X)
  
  Z.freq <- cbind(Z.freq,c(X.freq[,6],from.X))
  Z.freq <- cbind(Z.freq,c(from.Y,Y.freq[,6]))
  #find_ni <- function(data){
  #  ni <- sum(Z.freq[c(data),7])
  #  mi <- sum(Z.freq[c(data),8])
  #  return(ni)
  #}
  #find_mi <- function(data){
  #  ni <- sum(Z.freq[c(data),7])
  #  mi <- sum(Z.freq[c(data),8])
  #  return(mi)
  #}
  find_ni <- function(data){
    ni <- sum(Z.freq[c(data),7]/(Z.freq[c(data),7]+Z.freq[c(data),8]))
    return(ni)
  }
  
  Ni <- apply(knns,1,find_ni)
  Mi <- rep(nn,10000) - Ni
  
  
  gx <- (sqrt((eta)^2 * Ni/(Mi+1))-1)^2/2
  gx[which(is.na(gx))] <- 0
  gx[which(gx<0)] <- 0
  Dg <- max(1/m*sum(gx),0)
  return(sqrt(Dg))
  
  
}

##### Below functions are used to sample the allele frequency from mom method
# This is to sample for the hierarchical beta distribution
simulate.Hierarch.Beta <- function(parameter.list){
  muomega <- parameter.list[[1]][1]
  Varomega <- parameter.list[[1]][2]
  mueta1 <- parameter.list[[2]][1]
  Vareta1 <- parameter.list[[2]][2]
  mueta2 <- parameter.list[[3]][1]
  Vareta2 <- parameter.list[[3]][2]
  phiomega <- abs(muomega*(1-muomega)/(Varomega)-1)
  phieta1 <- abs(mueta1*(1-mueta1)/Vareta1-1)
  phieta2 <- abs(mueta2*(1-mueta2)/Vareta2-1)
  
  w <- rbeta(1,muomega*phiomega, phiomega*(1-muomega))
  eta1 <- rbeta(1,mueta1*phieta1, phieta1*(1-mueta1))
  eta2 <- rbeta(1,mueta2*phieta2, phieta2*(1-mueta2))
  
  x1 <- w*eta1
  x2 <- w*(1-eta1)
  x3 <- (1-w)*eta2
  x4 <- 1 - x1 - x2 - x3
  
  return(c(x1,x2,x3,x4))
}

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

Hierarch.Beta.approx_before<-function(Ext,Varxt){
  muomega <- Ext[1]+Ext[2]
  mueta1 <- Ext[1]/muomega
  mueta2 <- Ext[3]/(1-muomega)
  Varomega <- Varxt[1,1]+Varxt[2,2]+2*Varxt[1,2]
  Vareta1 <- (Varxt[1,1]-Varomega*mueta1^2)/(Varomega+muomega^2)
  Vareta2 <- (Varxt[3,3]-Varomega*mueta2^2)/(Varomega+(1-muomega)^2)
  phiomega <- muomega*(1-muomega)/(Varomega)-1
  phieta1 <- mueta1*(1-mueta1)/Vareta1-1
  phieta2 <- mueta2*(1-mueta2)/Vareta2-1
  return(list(c(muomega,Varomega),c(mueta1,Vareta1),c(mueta2,Vareta2)))
  
}
######

##### The main function for sampling from mom and true
simulation.Kimura <- function(x0,t,N){
  k = 2
  u = 10^-8
  Q <- matrix(c(-(k+2)*u, k * u, u, u,
                k * u, -(k+2)*u, u, u,
                u, u, -(k+2)*u, k * u,
                u, u, k * u, -(k+2)*u),nrow=4)
  Q1 <- Q
  diag(Q1) <- diag(Q1) + 1
  U <- Q1
  Q <- Q*N
  
  
  EV <- Kimura.model(x0,t,Q)
  parameter.list <- Hierarch.Beta.approx(EV$E,EV$V)
  n.samples <- 10000
  mom <- matrix(rep(0,n.samples * K),nrow = n.samples,ncol = K)
  true <- matrix(rep(0,n.samples * K),nrow = n.samples,ncol = K)
  for (i in 1:n.samples){
    mom[i,] <- simulate.Hierarch.Beta(parameter.list)
    true[i,] <- true.simulation(x0,U,K,t,N)
  }
  return(list(mom,true))
  
}



## this section is for playing with the bins and calculating the H distance
discretize.H.distance <- function(alleles,alleles1){
  c.matrix <- matrix(nrow=4,ncol=1)
  c.matrix[,1] <- c(1,0.1,0.1,0.1)
  count = c(0)
  for (x in 1:nrow(alleles)){
    alleles[x,][which(alleles[x,]==0)] <- alleles[x,][which(alleles[x,]==0)] + 10e-8
    indictor <- FALSE
    for (i in 1:ncol(c.matrix)){
      if (alleles[x,][1]<=c.matrix[1,i] && alleles[x,][1]>c.matrix[1,i]-0.1 
          && alleles[x,][2]<=c.matrix[2,i] && alleles[x,][2]>c.matrix[2,i]-0.1 
          && alleles[x,][3]<=c.matrix[3,i] && alleles[x,][3]>c.matrix[3,i]-0.1
          && alleles[x,][4]<=c.matrix[4,i] && alleles[x,][4]>c.matrix[4,i]-0.1){
        count[i] = count[i] + 1
        indictor <- TRUE
        break
      } 
    }
    if (!indictor){
      c.matrix <- cbind(c.matrix,findBin(alleles[x,]))
      count <- c(count,1)
    }
  }
  count1 <- rep(0,length(count))
  for (x in 1:nrow(alleles1)){
    alleles1[x,][which(alleles1[x,]==0)] <- alleles1[x,][which(alleles1[x,]==0)] + 10e-8
    indictor <- FALSE
    for (i in 1:ncol(c.matrix)){
      if (alleles1[x,][1]<=c.matrix[1,i] && alleles1[x,][1]>c.matrix[1,i]-0.1 
          && alleles1[x,][2]<=c.matrix[2,i] && alleles1[x,][2]>c.matrix[2,i]-0.1 
          && alleles1[x,][3]<=c.matrix[3,i] && alleles1[x,][3]>c.matrix[3,i]-0.1
          && alleles1[x,][4]<=c.matrix[4,i] && alleles1[x,][4]>c.matrix[4,i]-0.1){
        count1[i] = count1[i] + 1
        indictor <- TRUE
        break
      }
    }
    if (!indictor){
      c.matrix <- cbind(c.matrix,findBin(alleles1[x,]))
      count1 <- c(count1,1)
      count <- c(count,0)
    }
  }
  
  prob.mom <- count/sum(count)
  prob.true <- count1/sum(count1)
  Hellinger.distance <- sqrt(sum((sqrt(prob.mom)-sqrt(prob.true))^2))/sqrt(2)
  return(Hellinger.distance)
  
}


### run the simulation

## Find the bin that the value should belong to
findBin <- function(xs){
  mVector <- c(0,0,0,0)
  for (i in 1:4){
    mVector[i] <- ceiling(xs[i]*10)*0.1
  }
  return(mVector)
}


# a little bit computational trick
computation.trick <- function(mom,true,N){
  for (i in 1:nrow(mom)){
    for(j in 1:ncol(mom)){
      if (mom[i,j] < 1/N){
        mom[i,j] <- 0
      }
    }
    if (sum(mom[i,]) != 1){
      mom[i,] <- mom[i,]/sum(mom[i,])
    }
  }
  return(list(mom,true))
}






N <- 1000
ts <- c(10)
A.G <- c(0.5)
A.O.G <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
C.O.T <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
A.G.list <- list()
time.list <- list()
time = 1
for (t in ts){
  a.plus.g = 1
  for (c in A.G){
    heat.map.matrix <- matrix(rep(0,81),nrow = 9)
    a.of.g = 1
    for (x in A.O.G){
      c.of.t = 1
      for(y in C.O.T){
        x0 <- c(c*x,(1-x)*c,(1-c)*(y),(1-c)*(1-y))
        
        simus<-simulation.Kimura(x0,t,N)
        mom <-simus[[1]]
        true <-simus[[2]]
        distribution.list <- computation.trick(mom,true,N)
        continuous <- max(h.d(distribution.list[[1]],distribution.list[[2]],nn=100))
        print(continuous)
        heat.map.matrix[a.of.g,c.of.t] <- continuous
        print(c(time,a.plus.g,a.of.g,c.of.t))
        c.of.t = c.of.t + 1
        print(heat.map.matrix)
      }
      a.of.g = a.of.g + 1
    }
    A.G.list[[a.plus.g]] <- heat.map.matrix
    a.plus.g = a.plus.g + 1
  }
  time.list[[time]] <- A.G.list
  time = time + 1
}



for (i in 1:21){

  d = data.frame(x=rep(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),9), 
                 y=rep(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),each=9), 
                 z=c(time.list[[i]][[1]]))
  levelplot(z~x*y,data=d,col.regions=colorRampPalette(c("yellow","orange","red", "black")),at=c(seq(0,0.4,length=400)),xlab="C/(C+T)",ylab = "A/A+G",main=paste("t = ",t[i]))
  
}

t <- c(5,10)
x = 2
d = data.frame(x=rep(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),9), 
               y=rep(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),each=9), 
               z=c(time.list[[x]][[1]]))
levelplot(z~x*y,data=d,col.regions=colorRampPalette(c("yellow","red","pink", "blue","black")),at=c(seq(0,0.4,length=400)),xlab="C/(C+T)",ylab = "A/A+G",main=paste("t = ",t[x]))


