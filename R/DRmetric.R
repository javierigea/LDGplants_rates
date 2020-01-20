library(ape)
library(geiger)
library(plyr)
library(taxonlookup)
library(Taxonstand)
library(BAMMtools)
library(phangorn)
library(phytools)
library(caper)
library(picante)
library(diversitree)
library(RColorBrewer)

DR_statistic <- function(x, return.mean = FALSE){
  
  rootnode <- length(x$tip.label) + 1
  
  sprates <- numeric(length(x$tip.label))
  for (i in 1:length(sprates)){
    node <- i
    index <- 1
    qx <- 0
    while (node != rootnode){
      el <- x$edge.length[x$edge[,2] == node]
      node <- x$edge[,1][x$edge[,2] == node]
      
      qx <- qx + el* (1 / 2^(index-1))
      
      index <- index + 1
    }
    sprates[i] <- 1/qx
  }
  
  if (return.mean){
    return(mean(sprates))		
  }else{
    names(sprates) <- x$tip.label
    return(sprates)
  }
  
}