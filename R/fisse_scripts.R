library(ape)
library(diversitree)
source('./R/run_fisse/traitDependent_functions.R')
run_fisse_unbiased_replicates<-function(GBIFdata.BAMM,tree,sampling.proportions,replicates){
  cat('subsampling datasets','\n')
  GBIFdata.BAMM.subsampling<-lapply(c(1:replicates),function(x)GBIFdata.BAMM[c(sample(which(GBIFdata.BAMM$tropical==0),size = round(sampling.proportions[1]*size.dataset)),sample(which(GBIFdata.BAMM$tropical==1),size = round(sampling.proportions[2]*size.dataset))),])
  cat('subsampling trees','\n')
  trees.subsampling<-lapply(GBIFdata.BAMM.subsampling,function(x){dropped.tree<-drop.tip(tree,setdiff(tree$tip.label,x$binomial));return(dropped.tree)})
  cat('subsampling traits','\n')
  traits.subsampling<-lapply(GBIFdata.BAMM.subsampling,function(x){traits<-x$tropical;names(traits)<-x$binomial;return(traits)})
  cat('running FISSE','\n')
  fisse.binary.tropical.subsamplings<- lapply(c(1:length(traits.subsampling)),function(x){cat('running FISSE for ',x,' replicate','\n');FISSE.binary(trees.subsampling[[x]], traits.subsampling[[x]],reps=100)})
  saveRDS(fisse.binary.tropical.subsamplings,file=paste('./results/fisse.binary.tropical.subsamplings.',replicates,'.replicates.RDS',sep=''))
}
  