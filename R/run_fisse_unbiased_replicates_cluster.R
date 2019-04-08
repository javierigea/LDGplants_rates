.libPaths('/home_old/ji247/R/x86_64-pc-linux-gnu-library/3.5')
setwd('/home_old/ji247/ldg_plants/')

library(ape)
library(geiger)
library(plyr)
library(BAMMtools)
library(caper)
library(picante)
library(diversitree)
library(phangorn)

####function by L Revell to fix errors in is.ultrametric
force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

args<-commandArgs(trailingOnly = TRUE)
print(args)
interval.min<-as.numeric(args[1])
interval.max<-as.numeric(args[2])


source('./R/run_fisse/traitDependent_functions.R')
unbiased.tables.files<-list.files('./results/subsampled_unbiased/input_tables/',pattern='.txt')
unbiased.tables<-lapply(unbiased.tables.files,function(x)read.table(paste('./results/subsampled_unbiased/input_tables/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
names(unbiased.tables)<-gsub(unbiased.tables.files,pattern='GBIFdata.BAMM.subsampled_',replacement='')
for (i in c(interval.min:interval.max)){
  cat(gsub(names(unbiased.tables)[i],pattern='_table.txt',replacement=''),'\n')
  tree<-read.tree("./raw_data/Qian_GBIFdata.tree")
  #this is here to prevent is.ultrametric type errors with ape >3.5
  tree<-check_and_fix_ultrametric(tree)
  traits.subsampled<-unbiased.tables[[i]]$tropical
  names(traits.subsampled)<-unbiased.tables[[i]]$binomial
  tree.subsampled<-drop.tip(tree,setdiff(tree$tip.label,names(traits.subsampled)))
  cat('running fisse on ',gsub(names(unbiased.tables)[i],pattern='_table.txt',replacement=''),'\n')
  system.time(fisse.binary.subsampled<- FISSE.binary.median(phy=tree.subsampled, states = traits.subsampled,reps=100,tol=0.1))
  saveRDS(fisse.binary.subsampled,file=paste('./results/subsampled_unbiased/fisse/fisse.binary.subsampled_',gsub(names(unbiased.tables)[i],pattern='_table.txt',replacement=''),'_median.RDS',sep=''))
  
}
