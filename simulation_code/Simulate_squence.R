# simulation of tree sampling
library(TreeSim)
library(ape)
library(expm) 



# The death rate is 0 here
birth_rate <- 0.00001 # lambda
death_rate <- 0 # mu




# We generate a pure brith tree (Yule Model)
yule_tree <- sim.bdtree(birth_rate,0,n = 4)
yule_tree<-ladderize(yule_tree)
plot(yule_tree,show.tip.label = TRUE,show.node.label = TRUE, edge.width = 2, cex = 2,use.edge.length = TRUE,root.edge =T,main="Simulated Pure Birth Species Tree")
x <- yule_tree
x$edge.length


# normalize the Q matrix (Normalized)
k <- 2
u <- 0.00000001 
Q <- matrix(c(-(k+2)*u, k * u, u, u,
              k * u, -(k+2)*u, u, u,
              u, u, -(k+2)*u, k * u,
              u, u, k * u, -(k+2)*u),nrow = 4)
#stationary <- c(0.25,0.25,0.25,0.25)
#scale <- sum(diag(Q) * stationary)
#Q <- Q/(-scale)
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

#U <- Q
#diag(U) <- c(0,0,0,0)


true.simulation <- function(x0, Q, N, t){
  EV <- Kimura.model(x0,t,Q*N)
  parameters.list <-Hierarch.Beta.approx(EV$E,EV$V)
  xt <- simulate.Hierarch.Beta(parameters.list)
  return(filter_and_normalize(xt))
}



# This function is for us to sample for one site

simSeq.phylo <- function(x,Q,root.freq,N,sampled_individuals=10) {
  x <- reorder(x)
  edge <- x$edge
  nNodes <- max(edge)
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parent[!match(parent, child, 0)][1])
  tl <- x$edge.length
  x0 = root.freq
  res_freq <- matrix(NA,nNodes,4)
  res_freq[root,]<-root.freq
  for (i in seq_along(tl)) {
    from <- parent[i]
    to <- child[i]
    t <- tl[i]
    x0 <- res_freq[from,]
    res_freq[to,] <- true.simulation(x0,Q,N,t) 
  }
  k <- length(x$tip.label)
  label <- c(x$tip.label, as.character( (k + 1):nNodes))
  rownames(res_freq) <- label
  n.tip <- k
  res_freq <- res_freq[1:n.tip,]
  res1 <- matrix(NA, nNodes, sampled_individuals)
  label <- c(x$tip.label, as.character( (k + 1):nNodes))
  rownames(res1) <- label
  res1 <- res1[1:n.tip,]
  for (i in 1:n.tip){
    res1[i,] <- sample(c("A","G","C","T"),sampled_individuals,T,res_freq[i,])
  }
  return(res1)
}

normalize <- function(v){
  return(v/sum(v))
}

filter_and_normalize <- function(v,N){
  for (i in 1:length(v)){
    if (v[i] < 10^-16){
      v[i] <- 0
    }
  }
  return(normalize(v))
}

# This function is for us to get the final simulation
sim.sites <- function(x,Q,N,n.sites,sampled_individual=10){
  
  
  datalist <- list()
  
  for (i in 1:n.sites){
    root.freq <- true.simulation(c(0.25,0.25,0.25,0.25),Q,N,100000000000)
    datalist[[i]] <- simSeq.phylo(x,Q,root.freq,N)
  }
  
  n.tip <- length(x$tip.label)
  
  allAlignment <- list()
  
  for (i in 1:n.tip){
    seq.list <- list()
    for (j in 1:sampled_individual){
      one.seq <- ''
      for (k in 1:n.sites){
        one.seq <- paste(one.seq,datalist[[k]][i,j],sep = '')
      }
      seq.list[[j]] <- one.seq
    }
    allAlignment[[i]] <- seq.list
  }
  
  names(allAlignment) <- rownames(datalist[[1]])[1:n.tip]
  
  return(allAlignment)
}

result<-sim.sites(x,Q,1000,1000)

x = c(1,0,0,0)
for (i in 1:1000){
  print(true.simulation(c(1,0,0,0),Q,1000,40000))
  x = x + true.simulation(c(1,0,0,0),Q,1000,40000)
}
print(x)
