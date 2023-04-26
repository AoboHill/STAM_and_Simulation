# simulation of tree sampling
library(TreeSim)
library(ape)
library(expm) 


set.seed(1)

# The death rate is 0 here
birth_rate <- 0.01 # lambda
death_rate <- 0 # mu




# We generate a pure brith tree (Yule Model)
yule_tree <- sim.bdtree(birth_rate,0,n = 4)
plot(yule_tree,show.tip.label = TRUE,show.node.label = TRUE, edge.width = 2, cex = 2)
x <- yule_tree



# normalize the Q matrix (Normalized)
k <- 2
u <- 0.05
Q <- matrix(c(-(k+2)*u, k * u, u, u,
              k * u, -(k+2)*u, u, u,
              u, u, -(k+2)*u, k * u,
              u, u, k * u, -(k+2)*u),nrow = 4)
#stationary <- c(0.25,0.25,0.25,0.25)
#scale <- sum(diag(Q) * stationary)
#Q <- Q/(-scale)


#U <- Q
#diag(U) <- c(0,0,0,0)


# Wright-Fisher simulation for one branch
#true.simulation <- function(x0, Q, N, t){
#  return(c(rmultinom(1,N,c(x0 %*% expm(N*Q*t)))))
#}


true.simulation <- function(x0, Q, N, t){
  Q <- Q/N 
  diag(Q) <- diag(Q) + 1
  U <- Q
  for (i in 1 : round(t * N)){
    x0 <- c(rmultinom(1,N,c(x0%*%U)))/N
  }
  return(x0)
}



# This function is to convert the number of the allele to the allele frequency given
# the vector of characters
countToFreq <- function(x){
  return(c(sum(x == "A"),sum(x == "G"),sum(x == "C"),sum(x == "T"))/length(x))
}

# This function is to help convert the number of alleles to characters
# which can help us identify the meaning of the numbers
FreqToCount <- function(x){
  return(c(rep("A",x[1]),rep("G",,x[2]),rep("C",x[3]),rep("T",x[4])))
}



# This function is for us to sample for one site

simSeq.phylo <- function(x,Q,rootsite,N,sampled_individuals=10) {
  x <- reorder(x)
  edge <- x$edge
  nNodes <- max(edge)
  res <- matrix(NA, nNodes, N)
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parent[!match(parent, child, 0)][1])
  res[root, ] <- rootsite
  tl <- x$edge.length
  x0 = rootsite
  for (i in seq_along(tl)) {
    from <- parent[i]
    to <- child[i]
    t <- tl[i]
    x0 <- countToFreq(res[from,])
    alleleCount <- true.simulation(x0,Q,N,t) * N
    res[to,] = FreqToCount(alleleCount)
  }
  k <- length(x$tip.label)
  label <- c(x$tip.label, as.character( (k + 1):nNodes))
  rownames(res) <- label
  res <- res[1:n.tip,]
  res1 <- matrix(NA, nNodes, sampled_individuals)
  label <- c(x$tip.label, as.character( (k + 1):nNodes))
  rownames(res1) <- label
  res1 <- res1[1:n.tip,]
  for (i in 1:n.tip){
    res1[i,] <- sample(c("A","G","C","T"),sampled_individuals,T,countToFreq(res[i,]))
  }
  return(res1)
}





# This function is for us to get the final simulation
sim.sites <- function(x,Q,N,n.sites,sampled_individual=10){
  
  
  datalist <- list()
  
  for (i in 1:n.sites){
    rootsite <- rep(sample(c("A","G","C","T"),1,T,rep(0.25,4)),N)
    datalist[[i]] <- simSeq.phylo(x,Q,rootsite,N)
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
