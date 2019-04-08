#################################################################################
# (c) Javier Igea
# igea.javier@gmail.com
# Last modified: 18/06/2017

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#################################################################################
library(plyr)
library(picante)
library(BAMMtools)
library(RPANDA)
library(geiger)
library(taxonlookup)
library(ape)
library(parallel)
library(phytools)
library(PBD)
library(caper)
library(MuMIn)
library(phylolm)
###################################################################
###################################################################
###################################################################
###################################################################
#from Rabosky 2016 "Challenges in estimation" in the supplementary info
###   This equation estimates the probability that a sample of k taxa from a clade of n total taxa
###       includes the root node, under a Yule process.
###
###   The equation is taken from:
###
###   Sanderson, M. J. 1996. How many taxa must be sampled to identify the root node of a large clade?
###                    Systematic Biology 45: 168 - 173
crownCaptureProbability <- function(n, k){
  p1 <- 2*(n-k)/((n-1)*(k+1));
  return(1-p1);
}

######get phylogroups function by Victor Soria-Carrasco
######modified from extract_phylogroups at https://bitbucket.org/visoca/compdivrate/
get.phylogroups<-function(tree, minage, maxage, mincladesize, ncores){
  root.age<-max(node.depth.edgelength(tree))
  nodes<-seq(1,tree$Nnode)+length(tree$tip.label)
  nodes.ages<-as.numeric(unlist(mclapply(nodes, 
                                         function(x) root.age-nodeheight(tree,x), mc.cores=ncores)))
  nodes.desc<-as.numeric(unlist(mclapply(nodes, 
                                         function(x) length(getDescendants(tree,x)), mc.cores=ncores)))
  nodes<-data.frame(node=nodes,age=nodes.ages,ndesc=nodes.desc)
  chosen.nodes<-nodes[(nodes$age>=minage & nodes$age<=maxage &
                         nodes$ndesc>mincladesize),]
  
  # sort nodes by clade size and remove the descendendants
  # in order to keep the most inclusive clades only
  chosen.nodes.sort<-chosen.nodes[order(-chosen.nodes[,3]),]
  for(i in 1:nrow(chosen.nodes.sort)){
    if (!is.na(chosen.nodes.sort[i,1])){
      desc<-getDescendants(tree, chosen.nodes.sort[i,1])
      chosen.nodes.sort<-chosen.nodes.sort[!(chosen.nodes.sort[,1] %in% desc),]
    }
  }
  
  # extract clades
  phylogroups<-list()
  for (i in 1:nrow(chosen.nodes.sort)){
    n<-chosen.nodes.sort$node[i]
    phylogroups[[i]]<-extract.clade(tree,n)
  }
  class(phylogroups)<-"multiPhylo"
  
  return(list(nodes=chosen.nodes.sort, phylogroups=phylogroups))
  #return(phylogroups)
}


######this runs RPANDA across phylogroups
####run_RPANDA_BM_phylogroup_select_size<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  f.lamb.cst<-function(t,y){y[1]}
####  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  f.mu.cst.0<-function(t,y){0}
####  f.mu.cst<-function(t,y){y[1]}
####  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t);y[1]*exp(y[2]*t)}
####  lamb_par_init.exp<-c(0.5,0.01)
####  mu_par_init.exp<-c(0.05,-0.01)
####  lamb_par_init.cst<-c(0.5)
####  mu_par_init.cst<-c(0.05)
####  mu_par_init.0<-c()
####  res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.cst.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par = c(5),mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  }
####  res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par_init.exp,mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.exp.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par=c(5,0.05),mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  }
####  res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par_init.cst,mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  if(res.lambda.cst.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par=c(5),mu_par=c(0.01),f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  }  
####  #res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par_init.cst,mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  #if(res.lambda.cst.mu.exp$mu_par[1]<0){
####  #  cat('negative','\n')
####  #  res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = c(5),mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  #}  
####  #res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par_init.exp,mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  #if(res.lambda.exp.mu.exp$mu_par[1]<0){
####  #  res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par=c(5,0.05),mu_par=c(0.05,-0.01),f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  ##  cat('negative','\n')
####  #}  
####  res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par_init.exp,mu_par_init.cst,f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  if(res.lambda.exp.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par=c(5,0.05),mu_par=c(0.01),f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  }  
####  #aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
####  #names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
####  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.cst$aicc)
####  names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.cst','lambda.cst.mu.cst')
####  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
####  if (names(aicc.weights)[1]=='lambda.cst.mu0'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.0'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####  }
####  #else if (names(aicc.weights)[1]=='lambda.exp.mu.exp'){
####  #  genus.table$lambda.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]
####  #  genus.table$mu.rpanda1<-res.lambda.exp.mu.exp$mu_par[1]
####  #  genus.table$ndr.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]-res.lambda.exp.mu.exp$mu_par[1]
####  #}
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.exp.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]-res.lambda.exp.mu.cst$mu_par[1]
####  }
####  #else if (names(aicc.weights)[1]=='lambda.cst.mu.exp'){
####  #  genus.table$lambda.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]
####  #  genus.table$mu.rpanda1<-res.lambda.cst.mu.exp$mu_par[1]
####  #  genus.table$ndr.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]-res.lambda.cst.mu.exp$mu_par[1]
####  #}
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
####  }
####  if (genus.table$mu.rpanda1[1]<0){
####    cat('*****NEGATIVE EXTINCTION RATE','\n')
####  }
####  genus.table$aicc<-names(aicc.weights)[1]
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if(length((trait))>1){
####    res.trait1<-fitContinuous(phylogroup,trait,model="BM")
####    genus.table$sigsq<-res.trait1$opt$sigsq
####    genus.table$z0<-res.trait1$opt$z0
####  }else{
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####  }
####  
####  return(genus.table)
####}
##########this runs RPANDA across phylogroups
####run_RPANDA_EB_phylogroup_select_size<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  f.lamb.cst<-function(t,y){y[1]}
####  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  f.mu.cst.0<-function(t,y){0}
####  f.mu.cst<-function(t,y){y[1]}
####  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t);y[1]*exp(y[2]*t)}
####  lamb_par_init.exp<-c(0.5,0.01)
####  mu_par_init.exp<-c(0.05,-0.01)
####  lamb_par_init.cst<-c(0.5)
####  mu_par_init.cst<-c(0.05)
####  mu_par_init.0<-c()
####  res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.cst.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par = c(5),mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  }
####  res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par_init.exp,mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.exp.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par=c(5,0.05),mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  }
####  res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par_init.cst,mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  if(res.lambda.cst.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par=c(5),mu_par=c(0.01),f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  }  
####  #res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par_init.cst,mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  #if(res.lambda.cst.mu.exp$mu_par[1]<0){
####  #  cat('negative','\n')
####  #  res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = c(5),mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  #}  
####  #res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par_init.exp,mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  #if(res.lambda.exp.mu.exp$mu_par[1]<0){
####  #  res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par=c(5,0.05),mu_par=c(0.05,-0.01),f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  ##  cat('negative','\n')
####  #}  
####  res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par_init.exp,mu_par_init.cst,f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  if(res.lambda.exp.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par=c(5,0.05),mu_par=c(0.01),f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  }  
####  #aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
####  #names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
####  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.cst$aicc)
####  names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.cst','lambda.cst.mu.cst')
####  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
####  if (names(aicc.weights)[1]=='lambda.cst.mu0'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.0'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####  }
####  #else if (names(aicc.weights)[1]=='lambda.exp.mu.exp'){
####  #  genus.table$lambda.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]
####  #  genus.table$mu.rpanda1<-res.lambda.exp.mu.exp$mu_par[1]
####  #  genus.table$ndr.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]-res.lambda.exp.mu.exp$mu_par[1]
####  #}
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.exp.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]-res.lambda.exp.mu.cst$mu_par[1]
####  }
####  #else if (names(aicc.weights)[1]=='lambda.cst.mu.exp'){
####  #  genus.table$lambda.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]
####  #  genus.table$mu.rpanda1<-res.lambda.cst.mu.exp$mu_par[1]
####  #  genus.table$ndr.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]-res.lambda.cst.mu.exp$mu_par[1]
####  #}
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
####  }
####  if (genus.table$mu.rpanda1[1]<0){
####    cat('*****NEGATIVE EXTINCTION RATE','\n')
####  }
####  genus.table$aicc<-names(aicc.weights)[1]
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if (length(trait)<2){
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####  }else{
####    name.check<-name.check(phylogroup,trait)
####    if(name.check!='OK'){
####      tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
####    }else{
####      tips.with.data<-phylogroup
####    }
####    if(length(tips.with.data$tip.label)>mincladesize){
####      res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
####      res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
####      genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
####      genus.table$aicc.BM<-res.trait1.BM$opt$aicc
####      genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
####      genus.table$aicc.EB<-res.trait1.EB$opt$aicc
####      genus.table$a.EB<-res.trait1.EB$opt$a
####    }else{
####      genus.table$sigsq.BM<-NA
####      genus.table$aicc.BM<-NA
####      genus.table$sigsq.EB<-NA
####      genus.table$aicc.EB<-NA
####      genus.table$a.EB<-NA
####    }
####    
####  }
####  return(genus.table)
####}
####
####run_RPANDA_EB_phylogroup_select_size_6models<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  f.lamb.cst<-function(t,y){y[1]}
####  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  f.mu.cst.0<-function(t,y){0}
####  f.mu.cst<-function(t,y){y[1]}
####  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  lamb_par_init.exp<-c(0.5,0.01)
####  mu_par_init.exp<-c(0.05,-0.01)
####  lamb_par_init.cst<-c(0.5)
####  mu_par_init.cst<-c(0.05)
####  mu_par_init.0<-c()
####  res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.cst.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par = c(5),mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  }
####  if(res.lambda.cst.mu.0$lamb_par[1]>3){
####    lamb_par_init.exp<-c(5,0.01)
####    mu_par_init.exp<-c(5,0.01)
####    lamb_par_init.cst<-c(0.5)
####    mu_par_init.cst<-c(0.2)
####    mu_par_init.0<-c()
####    
####  }
####  cat('lambdaexpmu0','\n')
####  res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par_init.exp,mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.exp.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par=c(5,0.05),mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE)
####  }
####  cat('lambdacstmucst','\n')
####  res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par_init.cst,mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  if(res.lambda.cst.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par=c(5),mu_par=c(0.01),f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE)
####  }
####  cat('lambdacstmuexp','\n')
####  res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par_init.cst,mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  if(res.lambda.cst.mu.exp$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = c(5),mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE)
####  }  
####  cat('lambdaexpmuexp','\n')
####  res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par_init.exp,mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  if(res.lambda.exp.mu.exp$mu_par[1]<0){
####    res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par=c(5,0.05),mu_par=c(0.05,-0.01),f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####  #  cat('negative','\n')
####  }
####  cat('lambdaexpmucst','\n')
####  res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par_init.exp,mu_par_init.cst,f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  if(res.lambda.exp.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par=c(5,0.05),mu_par=c(0.01),f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE)
####  }  
####  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
####  names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
####  #aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.cst$aicc)
####  #names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.cst','lambda.cst.mu.cst')
####  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
####  if (names(aicc.weights)[1]=='lambda.cst.mu0'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.0'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.exp'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.exp.mu.exp$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]-res.lambda.exp.mu.exp$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.exp.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]-res.lambda.exp.mu.cst$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.exp'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.exp$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]-res.lambda.cst.mu.exp$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
####  }
####  if (genus.table$mu.rpanda1[1]<0){
####    cat('*****NEGATIVE EXTINCTION RATE','\n')
####  }
####  genus.table$aicc<-names(aicc.weights)[1]
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if (length(trait)<2){
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####  }else{
####    name.check<-name.check(phylogroup,trait)
####    if(name.check!='OK'){
####      tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
####    }else{
####      tips.with.data<-phylogroup
####    }
####    if(length(tips.with.data$tip.label)>mincladesize){
####      res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
####      res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
####      genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
####      genus.table$aicc.BM<-res.trait1.BM$opt$aicc
####      genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
####      genus.table$aicc.EB<-res.trait1.EB$opt$aicc
####      genus.table$a.EB<-res.trait1.EB$opt$a
####    }else{
####      genus.table$sigsq.BM<-NA
####      genus.table$aicc.BM<-NA
####      genus.table$sigsq.EB<-NA
####      genus.table$aicc.EB<-NA
####      genus.table$a.EB<-NA
####    }
####    
####  }
####  return(genus.table)
####}
####

####run_RPANDA_EB_z0_phylogroup_select_lat<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$aicc<-NA
####   #genus.table$sigsq.BM<-NA
####   #genus.table$aicc.BM<-NA
####   #genus.table$sigsq.EB<-NA
####   #genus.table$aicc.EB<-NA
####   #genus.table$a.EB<-NA
####   
####   #genus.table$z0.EB<-NA
####    #genus.table$z0.BM<-NA
####    
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$aicc<-NA
####    #genus.table$sigsq.BM<-NA
####    #genus.table$aicc.BM<-NA
####    #genus.table$sigsq.EB<-NA
####    #genus.table$aicc.EB<-NA
####    #genus.table$a.EB<-NA
####    #genus.table$aicc<-NA
####    #genus.table$z0.EB<-NA
####    #genus.table$z0.BM<-NA
####    return(genus.table)
####  }
####  f.lamb.cst<-function(t,y){y[1]}
####  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  f.mu.cst.0<-function(t,y){0}
####  f.mu.cst<-function(t,y){y[1]}
####  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  lamb_par_init.exp<-c(0.5,0.01)
####  mu_par_init.exp<-c(0.05,-0.01)
####  lamb_par_init.cst<-c(0.5)
####  mu_par_init.cst<-c(0.05)
####  mu_par_init.0<-c()
####  res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.cst.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par = c(5),mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  }
####  if(res.lambda.cst.mu.0$lamb_par[1]>0.75){
####    lamb_par_init.exp<-c(5,0.01)
####    mu_par_init.exp<-c(5,0.01)
####    lamb_par_init.cst<-c(5)
####    mu_par_init.cst<-c(0.5)
####    mu_par_init.0<-c()
####    
####  }
####  cat('lambdaexpmu0','\n')
####  res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par_init.exp,mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE,dt=1e-3)
####  if(res.lambda.exp.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par=c(5,0.05),mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE,dt=1e-3)
####  }
####  cat('lambdacstmucst','\n')
####  res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par_init.cst,mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE,dt=1e-3)
####  if(res.lambda.cst.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par=c(5),mu_par=c(0.01),f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE,dt=1e-3)
####  }
####  cat('lambdacstmuexp','\n')
####  res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par_init.cst,mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE,dt=1e-3)
####  if(res.lambda.cst.mu.exp$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = c(5),mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE,dt=1e-3)
####  }  
####  cat('lambdaexpmuexp','\n')
####  res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par_init.exp,mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE,dt=1e-3)
####  if(res.lambda.exp.mu.exp$mu_par[1]<0){
####    res.lambda.exp.mu.exp<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par=c(5,0.05),mu_par=c(0.05,-0.01),f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE)
####    #  cat('negative','\n')
####  }
####  cat('lambdaexpmucst','\n')
####  res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par_init.exp,mu_par_init.cst,f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE,dt=1e-3)
####  if(res.lambda.exp.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.exp.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par=c(5,0.05),mu_par=c(0.01),f=sampling.fraction,expo.lamb=TRUE,cst.mu=TRUE,dt=1e-3)
####  }  
####  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
####  names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
####  #aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.cst$aicc)
####  #names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.cst','lambda.cst.mu.cst')
####  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
####  if (names(aicc.weights)[1]=='lambda.cst.mu0'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.0'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.exp'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.exp.mu.exp$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]-res.lambda.exp.mu.exp$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.exp.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.exp.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]-res.lambda.exp.mu.cst$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.exp'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.exp$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]-res.lambda.cst.mu.exp$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
####  }
####  if (genus.table$mu.rpanda1[1]<0){
####    cat('*****NEGATIVE EXTINCTION RATE','\n')
####  }
####  genus.table$aicc<-names(aicc.weights)[1]
####  #trait<-as.numeric(genus.table$Median.Latitude)
####  #names(trait)<-genus.table$tip
####  #trait<-trait[!is.na(trait)]
####  #if (length(trait)<2){
####  #  genus.table$sigsq.BM<-NA
####  #  genus.table$aicc.BM<-NA
####  #  genus.table$sigsq.EB<-NA
####  #  genus.table$aicc.EB<-NA
####  #  genus.table$a.EB<-NA
####  #  genus.table$z0.EB<-NA
####  #  genus.table$z0.BM<-NA
####  #}else{
####  #  name.check<-name.check(phylogroup,trait)
####  #  if(name.check!='OK'){
####  #    tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
####  #  }else{
####  #    tips.with.data<-phylogroup
####  #  }
####  #  if(length(tips.with.data$tip.label)>mincladesize){
####  #    res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
####  #    res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
####  #    genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
####  #    genus.table$aicc.BM<-res.trait1.BM$opt$aicc
####  #    genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
####  #    genus.table$aicc.EB<-res.trait1.EB$opt$aicc
####  #    genus.table$a.EB<-res.trait1.EB$opt$a
####  #    genus.table$z0.EB<-res.trait1.EB$opt$z0
####  #    genus.table$z0.BM<-res.trait1.BM$opt$z0
####  #  }else{
####  #    genus.table$sigsq.BM<-NA
####  #    genus.table$aicc.BM<-NA
####  #    genus.table$sigsq.EB<-NA
####  #    genus.table$aicc.EB<-NA
####  #    genus.table$a.EB<-NA
####  #    genus.table$z0.EB<-NA
####  #    genus.table$z0.BM<-NA
####  #  }
####  #  
####  #}
####  return(genus.table)
####}


####run_RPANDA_phylogroup<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$aicc<-NA
####    return(genus.table)
####  }
####  
####  f.lamb.cst<-function(t,y){y[1]}
####  f.mu.cst<-function(t,y){y>0;y[1]}
####  f.mu.cst.0<-function(t,y){0}
####  mu_par_init.0<-c()
####  lamb_par_init.cst<-0.5
####  mu_par_init.cst<-0.1
####  cat('lambda.cst.mu0','\n')
####  res.lambda.cst.mu.0<-try(fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE))
####  cat('lambda.cst.mu.cst','\n')
####  res.lambda.cst.mu.cst<-try(fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par = c(res.lambda.cst.mu.0$lamb_par[1]),mu_par = mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE))
####  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  
####  lamb_par_init.exp<-c(res.lambda.cst.mu.cst$lamb_par[1],0.05)
####  mu_par_init.exp<-c(res.lambda.cst.mu.cst$mu_par[1],0.1)
####  if(res.lambda.cst.mu.cst$lamb_par[1]<0){
####    lamb_par_init.exp<-c(0,0.1)
####  }else{
####    lamb_par_init.exp<-c(res.lambda.cst.mu.cst$lamb_par[1],0.05)
####  }
####  if(res.lambda.cst.mu.cst$mu_par[1]<0){
####    mu_par_init.exp<-c(0,0.1)
####  }else{
####    mu_par_init.exp<-c(res.lambda.cst.mu.cst$mu_par[1],0.1)
####  }
####  
####  if(f.lamb.exp(tot_time,lamb_par_init.exp)<f.mu.exp(tot_time,mu_par_init.exp)){
####    lamb_par_init.exp<-c(mu_par_init.exp[1]*1.0005,mu_par_init.exp[2]*0.999)
####  }
####  
####  
####  
####  
####  cat('lambda.exp.mu.0','\n')
####  res.lambda.exp.mu.0<-try(fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par = lamb_par_init.exp,mu_par = mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE,dt=1e-3))
####  cat('lambda.exp.mu.cst','\n')
####  res.lambda.exp.mu.cst<-try(fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par = lamb_par_init.exp,mu_par = mu_par_init.cst,f=sampling.fraction,expo.lamb =TRUE,cst.mu=TRUE,dt=1e-3))
####  cat('lambda.cst.mu.exp','\n')
####  res.lambda.cst.mu.exp<-try(fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = lamb_par_init.cst,mu_par = mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE,dt=1e-3))
####  cat('lambda.exp.mu.exp','\n')
####  res.lambda.exp.mu.exp<-try(fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par = lamb_par_init.exp,mu_par = mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE,dt=1e-3))
####  cat('aicvector parameters','\n')
####  if(class(res.lambda.exp.mu.0)=='try-error'){
####    res.lambda.exp.mu.0<-try(fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par = c(res.lambda.cst.mu.cst$lamb_par[1],0.01),mu_par = mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE,dt=1e-3))
####  }
####  if(class(res.lambda.exp.mu.cst)=='try-error'){
####    res.lambda.exp.mu.cst<-try(fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par = c(res.lambda.cst.mu.cst$lamb_par[1],0.01),mu_par = mu_par_init.cst,f=sampling.fraction,expo.lamb =TRUE,cst.mu=TRUE,dt=1e-3))
####  }
####  if(class(res.lambda.cst.mu.exp)=='try-error'){
####    res.lambda.cst.mu.exp<-try(fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = 0.5,mu_par =c(0.1,-0.1) ,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE,dt=1e-3))
####  }
####  if(class(res.lambda.cst.mu.0)=='try-error'){
####    class(res.lambda.cst.mu.0)<-'fit.bd'
####    res.lambda.cst.mu.0$aicc<-1000
####  }
####  if(is.infinite(res.lambda.cst.mu.0$aicc)){
####    res.lambda.cst.mu.0$aicc<-1000
####  }
####  if(class(res.lambda.exp.mu.0)=='try-error'){
####    class(res.lambda.exp.mu.0)<-'fit.bd'
####    res.lambda.exp.mu.0$aicc<-1000
####  }
####  if(is.infinite(res.lambda.exp.mu.0$aicc)||class(res.lambda.exp.mu.0)=='try-error'){
####    res.lambda.exp.mu.0$aicc<-1000
####  }
####  if(class(res.lambda.exp.mu.exp)=='try-error'){
####    class(res.lambda.exp.mu.exp)<-'fit.bd'
####    res.lambda.exp.mu.exp$aicc<-1000
####  }
####  if(is.infinite(res.lambda.exp.mu.exp$aicc)||class(res.lambda.exp.mu.exp)=='try-error'){
####    res.lambda.exp.mu.exp$aicc<-1000
####  }
####  if(class(res.lambda.exp.mu.cst)=='try-error'){
####    class(res.lambda.exp.mu.cst)<-'fit.bd'
####    res.lambda.exp.mu.cst$aicc<-1000
####  }
####  if(is.infinite(res.lambda.exp.mu.cst$aicc)||class(res.lambda.exp.mu.cst)=='try-error'){
####    res.lambda.exp.mu.cst$aicc<-1000
####  }
####  if(class(res.lambda.cst.mu.exp)=='try-error'){
####    class(res.lambda.cst.mu.exp)<-'fit.bd'
####    res.lambda.cst.mu.exp$aicc<-1000
####  }
####  if(is.infinite(res.lambda.cst.mu.exp$aicc)||class(res.lambda.cst.mu.exp)=='try-error'){
####    res.lambda.cst.mu.exp$aicc<-1000
####  }
####  if(class(res.lambda.cst.mu.cst)=='try-error'){
####    class(res.lambda.cst.mu.cst)<-'fit.bd'
####    res.lambda.cst.mu.cst$aicc<-1000
####  }
####  if(is.infinite(res.lambda.cst.mu.cst$aicc)||class(res.lambda.cst.mu.cst)=='try-error'){
####    res.lambda.cst.mu.cst$aicc<-1000
####  }
####  
####  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
####  names(aicc.vector)<-c('lambda.cst.mu0','lambda.exp.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
####  if(length(which(aicc.vector==1000))>0){
####    aicc.vector[which(aicc.vector==1000)]<-max(aicc.vector[-which(aicc.vector==1000)])
####  }
####  #aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.cst.mu.cst$aicc)
####  #names(aicc.vector)<-c('lambda.cst.mu0','lambda.cst.mu.cst')
####  cat('getting parameters','\n')
####  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
####  if (names(aicc.weights)[1]=='lambda.cst.mu0'){
####    lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####    mu.rpanda1<-0
####    ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####  }else if (names(aicc.weights)[1]=='lambda.exp.mu.0'){
####    lambda.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####    mu.rpanda1<-0
####    ndr.rpanda1<-res.lambda.exp.mu.0$lamb_par[1]
####  }else if (names(aicc.weights)[1]=='lambda.exp.mu.exp'){
####    lambda.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]
####    mu.rpanda1<-res.lambda.exp.mu.exp$mu_par[1]
####    ndr.rpanda1<-res.lambda.exp.mu.exp$lamb_par[1]-res.lambda.exp.mu.exp$mu_par[1]
####  }else if (names(aicc.weights)[1]=='lambda.exp.mu.cst'){
####    lambda.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]
####    mu.rpanda1<-res.lambda.exp.mu.cst$mu_par[1]
####    ndr.rpanda1<-res.lambda.exp.mu.cst$lamb_par[1]-res.lambda.exp.mu.cst$mu_par[1]
####  }else if (names(aicc.weights)[1]=='lambda.cst.mu.exp'){
####    lambda.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]
####    mu.rpanda1<-res.lambda.cst.mu.exp$mu_par[1]
####    ndr.rpanda1<-res.lambda.cst.mu.exp$lamb_par[1]-res.lambda.cst.mu.exp$mu_par[1]
####  }else if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
####    lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
####    mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
####    ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
####  }
####  if (mu.rpanda1<0){
####    cat('*****NEGATIVE EXTINCTION RATE','\n')
####  }
####  aicc<-names(aicc.weights)[1]
####  results.vector<-c(lambda.rpanda1,mu.rpanda1,ndr.rpanda1,aicc)
####  names(results.vector)<-c('lambda.rpanda1','mu.rpanda1','ndr.rpanda1','aicc')
####  
####  cat('outputting','\n')
####  genus.table$lambda.rpanda1<-as.numeric(results.vector['lambda.rpanda1'])
####  genus.table$mu.rpanda1<-as.numeric(results.vector['mu.rpanda1'])
####  genus.table$ndr.rpanda1<-as.numeric(results.vector['ndr.rpanda1'])
####  genus.table$aicc<-results.vector['aicc']
####  return(genus.table)
####  
####}

run_RPANDA_phylogroup_modelave<-function(phylogroup,table,sampling,mincladesize){
  genus.table<-table[match(phylogroup$tip.label,table$tip),]
  cat(genus.table$New_Species[1],'\n')
  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
  genus.table.genera<-unique(genus.table.genera)
  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
  colnames(genus.counts)<-c('genus','n.of.tips')
  genus.table.genera<-merge(genus.table.genera,genus.counts)
  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
  genus.table$phylogroup.sampling.fraction<-sampling.fraction
  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
  tot_time<-max(node.age(phylogroup)$ages)
  genus.table$phylogroup.size<-nrow(genus.table)
  genus.table$phylogroup.name<-phylogroup$tip.label[1]
  #skip for small or undersampled clades
  if(genus.table$phylogroup.size[1]<mincladesize){
    genus.table$lambda.rpanda1<-NA
    genus.table$mu.rpanda1<-NA
    genus.table$ndr.rpanda1<-NA
    genus.table$aicc<-NA
    return(genus.table)
  }
  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
    genus.table$lambda.rpanda1<-NA
    genus.table$mu.rpanda1<-NA
    genus.table$ndr.rpanda1<-NA
    genus.table$aicc<-NA
    return(genus.table)
  }
  
  f.lamb.cst<-function(t,y){y[1]}
  f.mu.cst<-function(t,y){y>0;y[1]}
  f.mu.cst.0<-function(t,y){0}
  mu_par_init.0<-c()
  lamb_par_init.cst<-0.5
  mu_par_init.cst<-0.1
  cat('lambda.cst.mu0','\n')
  res.lambda.cst.mu.0<-try(fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE))
  cat('lambda.cst.mu.cst','\n')
  res.lambda.cst.mu.cst<-try(fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par = c(res.lambda.cst.mu.0$lamb_par[1]),mu_par = mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE))
  f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
  f.mu.exp<-function(t,y){y[1]*exp(y[2]*t)}
  
  lamb_par_init.exp<-c(res.lambda.cst.mu.cst$lamb_par[1],0.05)
  mu_par_init.exp<-c(res.lambda.cst.mu.cst$mu_par[1],0.1)
  if(res.lambda.cst.mu.cst$lamb_par[1]<0){
    lamb_par_init.exp<-c(0,0.1)
  }else{
    lamb_par_init.exp<-c(res.lambda.cst.mu.cst$lamb_par[1],0.05)
  }
  if(res.lambda.cst.mu.cst$mu_par[1]<0){
    mu_par_init.exp<-c(0,0.1)
  }else{
    mu_par_init.exp<-c(res.lambda.cst.mu.cst$mu_par[1],0.1)
  }
  
  if(f.lamb.exp(tot_time,lamb_par_init.exp)<f.mu.exp(tot_time,mu_par_init.exp)){
    lamb_par_init.exp<-c(mu_par_init.exp[1]*1.0005,mu_par_init.exp[2]*0.999)
  }
  
  
  
  
  cat('lambda.exp.mu.0','\n')
  res.lambda.exp.mu.0<-try(fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par = lamb_par_init.exp,mu_par = mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE,dt=1e-3))
  cat('lambda.exp.mu.cst','\n')
  res.lambda.exp.mu.cst<-try(fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par = lamb_par_init.exp,mu_par = mu_par_init.cst,f=sampling.fraction,expo.lamb =TRUE,cst.mu=TRUE,dt=1e-3))
  cat('lambda.cst.mu.exp','\n')
  res.lambda.cst.mu.exp<-try(fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = lamb_par_init.cst,mu_par = mu_par_init.exp,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE,dt=1e-3))
  cat('lambda.exp.mu.exp','\n')
  res.lambda.exp.mu.exp<-try(fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.exp,lamb_par = lamb_par_init.exp,mu_par = mu_par_init.exp,f=sampling.fraction,expo.lamb=TRUE,expo.mu=TRUE,dt=1e-3))
  cat('aicvector parameters','\n')
  if(class(res.lambda.exp.mu.0)=='try-error'){
    res.lambda.exp.mu.0<-try(fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst.0,lamb_par = c(res.lambda.cst.mu.cst$lamb_par[1],0.01),mu_par = mu_par_init.0,f=sampling.fraction,expo.lamb =TRUE,fix.mu=TRUE,dt=1e-3))
  }
  if(class(res.lambda.exp.mu.cst)=='try-error'){
    res.lambda.exp.mu.cst<-try(fit_bd(phylogroup,tot_time,f.lamb.exp,f.mu.cst,lamb_par = c(res.lambda.cst.mu.cst$lamb_par[1],0.01),mu_par = mu_par_init.cst,f=sampling.fraction,expo.lamb =TRUE,cst.mu=TRUE,dt=1e-3))
  }
  if(class(res.lambda.cst.mu.exp)=='try-error'){
    res.lambda.cst.mu.exp<-try(fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.exp,lamb_par = 0.5,mu_par =c(0.1,-0.1) ,f=sampling.fraction,cst.lamb =TRUE,expo.mu=TRUE,dt=1e-3))
  }
  #setting AIC to 1000 (arbitrarily) for the models that failed to be evaluated, then I will set the AIC of these models to the maxAIC
  if(class(res.lambda.cst.mu.0)=='try-error'){
    class(res.lambda.cst.mu.0)<-'fit.bd'
    res.lambda.cst.mu.0$aicc<-1000
  }
  if(is.infinite(res.lambda.cst.mu.0$aicc)){
    res.lambda.cst.mu.0$aicc<-1000
  }
  if(class(res.lambda.exp.mu.0)=='try-error'){
    class(res.lambda.exp.mu.0)<-'fit.bd'
    res.lambda.exp.mu.0$aicc<-1000
  }
  if(is.infinite(res.lambda.exp.mu.0$aicc)||class(res.lambda.exp.mu.0)=='try-error'){
    res.lambda.exp.mu.0$aicc<-1000
  }
  if(class(res.lambda.exp.mu.exp)=='try-error'){
    class(res.lambda.exp.mu.exp)<-'fit.bd'
    res.lambda.exp.mu.exp$aicc<-1000
  }
  if(is.infinite(res.lambda.exp.mu.exp$aicc)||class(res.lambda.exp.mu.exp)=='try-error'){
    res.lambda.exp.mu.exp$aicc<-1000
  }
  if(class(res.lambda.exp.mu.cst)=='try-error'){
    class(res.lambda.exp.mu.cst)<-'fit.bd'
    res.lambda.exp.mu.cst$aicc<-1000
  }
  if(is.infinite(res.lambda.exp.mu.cst$aicc)||class(res.lambda.exp.mu.cst)=='try-error'){
    res.lambda.exp.mu.cst$aicc<-1000
  }
  if(class(res.lambda.cst.mu.exp)=='try-error'){
    class(res.lambda.cst.mu.exp)<-'fit.bd'
    res.lambda.cst.mu.exp$aicc<-1000
  }
  if(is.infinite(res.lambda.cst.mu.exp$aicc)||class(res.lambda.cst.mu.exp)=='try-error'){
    res.lambda.cst.mu.exp$aicc<-1000
  }
  if(class(res.lambda.cst.mu.cst)=='try-error'){
    class(res.lambda.cst.mu.cst)<-'fit.bd'
    res.lambda.cst.mu.cst$aicc<-1000
  }
  if(is.infinite(res.lambda.cst.mu.cst$aicc)||class(res.lambda.cst.mu.cst)=='try-error'){
    res.lambda.cst.mu.cst$aicc<-1000
  }
  
  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.exp.mu.0$aicc,res.lambda.exp.mu.exp$aicc,res.lambda.exp.mu.cst$aicc,res.lambda.cst.mu.exp$aicc,res.lambda.cst.mu.cst$aicc)
  names(aicc.vector)<-c('lambda.cst.mu.0','lambda.exp.mu.0','lambda.exp.mu.exp','lambda.exp.mu.cst','lambda.cst.mu.exp','lambda.cst.mu.cst')
  if(length(which(aicc.vector==1000))>0){
    aicc.vector[which(aicc.vector==1000)]<-max(aicc.vector[-which(aicc.vector==1000)])
  }
  #aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.cst.mu.cst$aicc)
  #names(aicc.vector)<-c('lambda.cst.mu0','lambda.cst.mu.cst')
  cat('getting parameters','\n')
  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
  #this does the model averaging. Take each weight and do the weighted.mean
  aicc.weighted.lambda<-vector()
  aicc.weighted.mu<-vector()
  for(e in 1:length(aicc.weights)){
    lambda.aux<-0
    mu.aux<-0
    
    if(class(get(paste('res.',names(aicc.weights[e]),sep='')))=='fit.bd'){
      if((length(get(paste('res.',names(aicc.weights[e]),sep=''))$lamb_par[1])>0)&(is.na(get(paste('res.',names(aicc.weights[e]),sep=''))$lamb_par[1])==F)){
        lambda.aux<-get(paste('res.',names(aicc.weights[e]),sep=''))$lamb_par[1]
        
      }
      aicc.weighted.lambda[e]<-lambda.aux
      if(length(grep('mu.0',names(aicc.weights[e])))>0){
        mu.aux<-0
        aicc.weighted.mu[e]<-mu.aux
      }else{
        if((length(get(paste('res.',names(aicc.weights[e]),sep=''))$mu_par[1])>0)&(is.na(get(paste('res.',names(aicc.weights[e]),sep=''))$mu_par[1])==F)){
          mu.aux<-get(paste('res.',names(aicc.weights[e]),sep=''))$mu_par[1]
          aicc.weighted.mu[e]<-mu.aux
        }
      }
    }else{
      aicc.weighted.lambda[e]<-0
      aicc.weighted.mu[e]<-0
    }
    
  }
  aicc.weighted.lambda[aicc.weighted.lambda<0]<-0
  aicc.weighted.mu[aicc.weighted.mu<0]<-0
  lambda.rpanda1<-weighted.mean(x=aicc.weighted.lambda,w=aicc.weights)
  mu.rpanda1<-weighted.mean(x=aicc.weighted.mu,w=aicc.weights)
  ndr.rpanda1<-lambda.rpanda1-mu.rpanda1
  
  aicc<-names(aicc.weights)[1]
  results.vector<-c(lambda.rpanda1,mu.rpanda1,ndr.rpanda1,aicc)
  names(results.vector)<-c('lambda.rpanda1','mu.rpanda1','ndr.rpanda1','aicc')
  
  cat('outputting','\n')
  genus.table$lambda.rpanda1<-as.numeric(results.vector['lambda.rpanda1'])
  genus.table$mu.rpanda1<-as.numeric(results.vector['mu.rpanda1'])
  genus.table$ndr.rpanda1<-as.numeric(results.vector['ndr.rpanda1'])
  genus.table$aicc<-results.vector['aicc']
  return(genus.table)
  
}
##########this runs MS across phylogroups
####run_MS_BM_phylogroup_select_size<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####    return(genus.table)
####  }
####  genus.table$lambda.ms0<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0)
####  genus.table$lambda.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)/(1-0.5)
####  genus.table$lambda.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)/(1-0.9)
####  genus.table$ndr.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)
####  genus.table$ndr.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if(length((trait))>1){
####    res.trait1<-fitContinuous(phylogroup,trait,model="BM")
####    genus.table$sigsq<-res.trait1$opt$sigsq
####    genus.table$z0<-res.trait1$opt$z0
####  }else{
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####  }
####  return(genus.table)
####}
####
##########this runs MS + trait with BM & EB across phylogroups
####run_MS_EB_phylogroup_select_size<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    return(genus.table)
####  }
####  genus.table$lambda.ms0<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0)
####  genus.table$lambda.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)/(1-0.5)
####  genus.table$lambda.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)/(1-0.9)
####  genus.table$ndr.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)
####  genus.table$ndr.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if (length(trait)<2){
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####  }else{
####    name.check<-name.check(phylogroup,trait)
####    if(name.check!='OK'){
####      tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
####    }else{
####      tips.with.data<-phylogroup
####    }
####    if(length(tips.with.data$tip.label)>mincladesize){
####      res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
####      res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
####      genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
####      genus.table$aicc.BM<-res.trait1.BM$opt$aicc
####      genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
####      genus.table$aicc.EB<-res.trait1.EB$opt$aicc
####      genus.table$a.EB<-res.trait1.EB$opt$a
####    }else{
####      genus.table$sigsq.BM<-NA
####      genus.table$aicc.BM<-NA
####      genus.table$sigsq.EB<-NA
####      genus.table$aicc.EB<-NA
####      genus.table$a.EB<-NA
####    }
####      
####  }
####  return(genus.table)
####}

####run_RPANDA_linear_phylogroup_select_lat<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$aicc<-NA
####    #genus.table$sigsq.BM<-NA
####    #genus.table$aicc.BM<-NA
####    #genus.table$sigsq.EB<-NA
####    #genus.table$aicc.EB<-NA
####    #genus.table$a.EB<-NA
####    
####    #genus.table$z0.EB<-NA
####    #genus.table$z0.BM<-NA
####    
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.rpanda1<-NA
####    genus.table$mu.rpanda1<-NA
####    genus.table$ndr.rpanda1<-NA
####    genus.table$aicc<-NA
####    #genus.table$sigsq.BM<-NA
####    #genus.table$aicc.BM<-NA
####    #genus.table$sigsq.EB<-NA
####    #genus.table$aicc.EB<-NA
####    #genus.table$a.EB<-NA
####    #genus.table$aicc<-NA
####    #genus.table$z0.EB<-NA
####    #genus.table$z0.BM<-NA
####    return(genus.table)
####  }
####  f.lamb.cst<-function(t,y){y[1]}
####  #f.lamb.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  f.lamb.lin<-function(t,y){y[1]+y[2]*t}
####  f.mu.cst.0<-function(t,y){0}
####  f.mu.cst<-function(t,y){y[1]}
####  #f.mu.exp<-function(t,y){y[1]*exp(y[2]*t)}
####  f.mu.lin<-function(t,y){y[1]+y[2]*t}
####  lamb_par_init.lin<-c(0.5,0.01)
####  mu_par_init.lin<-c(0.05,-0.01)
####  lamb_par_init.cst<-c(0.5)
####  mu_par_init.cst<-c(0.05)
####  mu_par_init.0<-c()
####  res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par_init.cst,mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  if(res.lambda.cst.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst.0,lamb_par = c(5),mu_par_init.0,f=sampling.fraction,cst.lamb =TRUE,fix.mu=TRUE)
####  }
####  if(res.lambda.cst.mu.0$lamb_par[1]>3){
####    lamb_par_init.lin<-c(5,0.01)
####    mu_par_init.lin<-c(5,0.01)
####    lamb_par_init.cst<-c(0.5)
####    mu_par_init.cst<-c(0.2)
####    mu_par_init.0<-c()
####    
####  }
####  cat('lambdalinmu0','\n')
####  res.lambda.lin.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.lin,f.mu.cst.0,lamb_par_init.lin,mu_par_init.0,f=sampling.fraction,fix.mu=TRUE,dt=1e-3)
####  if(res.lambda.lin.mu.0$lamb_par[1]<0){
####    cat('negative','\n')
####    res.lambda.lin.mu.0<-fit_bd(phylogroup,tot_time,f.lamb.lin,f.mu.cst.0,lamb_par=c(5,0.05),mu_par_init.0,f=sampling.fraction,fix.mu=TRUE,dt=1e-3)
####  }
####  cat('lambdacstmucst','\n')
####  res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par_init.cst,mu_par_init.cst,f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE,dt=1e-3)
####  if(res.lambda.cst.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.cst,lamb_par=c(5),mu_par=c(0.01),f=sampling.fraction,cst.lamb =TRUE,cst.mu=TRUE,dt=1e-3)
####  }
####  cat('lambdacstmulin','\n')
####  res.lambda.cst.mu.lin<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.lin,lamb_par_init.cst,mu_par_init.lin,f=sampling.fraction,cst.lamb =TRUE,dt=1e-3)
####  if(res.lambda.cst.mu.lin$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.cst.mu.lin<-fit_bd(phylogroup,tot_time,f.lamb.cst,f.mu.lin,lamb_par = c(5),mu_par_init.lin,f=sampling.fraction,cst.lamb =TRUE,dt=1e-3)
####  }  
####  cat('lambdalinmulin','\n')
####  res.lambda.lin.mu.lin<-fit_bd(phylogroup,tot_time,f.lamb.lin,f.mu.lin,lamb_par_init.lin,mu_par_init.lin,f=sampling.fraction,dt=1e-3)
####  if(res.lambda.lin.mu.lin$mu_par[1]<0){
####    res.lambda.lin.mu.lin<-fit_bd(phylogroup,tot_time,f.lamb.lin,f.mu.lin,lamb_par=c(5,0.05),mu_par=c(0.05,-0.01),f=sampling.fraction,dt=1e-3)
####    #  cat('negative','\n')
####  }
####  cat('lambdalinmucst','\n')
####  res.lambda.lin.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.lin,f.mu.cst,lamb_par_init.lin,mu_par_init.cst,f=sampling.fraction,cst.mu=TRUE,dt=1e-3)
####  if(res.lambda.lin.mu.cst$mu_par[1]<0){
####    cat('negative','\n')
####    res.lambda.lin.mu.cst<-fit_bd(phylogroup,tot_time,f.lamb.lin,f.mu.cst,lamb_par=c(5,0.05),mu_par=c(0.01),f=sampling.fraction,cst.mu=TRUE,dt=1e-3)
####  }  
####  aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.lin.mu.0$aicc,res.lambda.lin.mu.lin$aicc,res.lambda.lin.mu.cst$aicc,res.lambda.cst.mu.lin$aicc,res.lambda.cst.mu.cst$aicc)
####  names(aicc.vector)<-c('lambda.cst.mu0','lambda.lin.mu.0','lambda.lin.mu.lin','lambda.lin.mu.cst','lambda.cst.mu.lin','lambda.cst.mu.cst')
####  #aicc.vector<-c(res.lambda.cst.mu.0$aicc,res.lambda.lin.mu.0$aicc,res.lambda.lin.mu.cst$aicc,res.lambda.cst.mu.cst$aicc)
####  #names(aicc.vector)<-c('lambda.cst.mu0','lambda.lin.mu.0','lambda.lin.mu.cst','lambda.cst.mu.cst')
####  aicc.weights<-sort(Weights(aicc.vector),decreasing=T)
####  if (names(aicc.weights)[1]=='lambda.cst.mu0'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.lin.mu.0'){
####    genus.table$lambda.rpanda1<-res.lambda.lin.mu.0$lamb_par[1]
####    genus.table$mu.rpanda1<-0
####    genus.table$ndr.rpanda1<-res.lambda.lin.mu.0$lamb_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.lin.mu.lin'){
####    genus.table$lambda.rpanda1<-res.lambda.lin.mu.lin$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.lin.mu.lin$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.lin.mu.lin$lamb_par[1]-res.lambda.lin.mu.lin$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.lin.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.lin.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.lin.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.lin.mu.cst$lamb_par[1]-res.lambda.lin.mu.cst$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.lin'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.lin$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.lin$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.lin$lamb_par[1]-res.lambda.cst.mu.lin$mu_par[1]
####  }
####  else if (names(aicc.weights)[1]=='lambda.cst.mu.cst'){
####    genus.table$lambda.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]
####    genus.table$mu.rpanda1<-res.lambda.cst.mu.cst$mu_par[1]
####    genus.table$ndr.rpanda1<-res.lambda.cst.mu.cst$lamb_par[1]-res.lambda.cst.mu.cst$mu_par[1]
####  }
####  if (genus.table$mu.rpanda1[1]<0){
####    cat('*****NEGATIVE EXTINCTION RATE','\n')
####  }
####  genus.table$aicc<-names(aicc.weights)[1]
####  #trait<-as.numeric(genus.table$Median.Latitude)
####  #names(trait)<-genus.table$tip
####  #trait<-trait[!is.na(trait)]
####  #if (length(trait)<2){
####  #  genus.table$sigsq.BM<-NA
####  #  genus.table$aicc.BM<-NA
####  #  genus.table$sigsq.EB<-NA
####  #  genus.table$aicc.EB<-NA
####  #  genus.table$a.EB<-NA
####  #  genus.table$z0.EB<-NA
####  #  genus.table$z0.BM<-NA
####  #}else{
####  #  name.check<-name.check(phylogroup,trait)
####  #  if(name.check!='OK'){
####  #    tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
####  #  }else{
####  #    tips.with.data<-phylogroup
####  #  }
####  #  if(length(tips.with.data$tip.label)>mincladesize){
####  #    res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
####  #    res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
####  #    genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
####  #    genus.table$aicc.BM<-res.trait1.BM$opt$aicc
####  #    genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
####  #    genus.table$aicc.EB<-res.trait1.EB$opt$aicc
####  #    genus.table$a.EB<-res.trait1.EB$opt$a
####  #    genus.table$z0.EB<-res.trait1.EB$opt$z0
####  #    genus.table$z0.BM<-res.trait1.BM$opt$z0
####  #  }else{
####  #    genus.table$sigsq.BM<-NA
####  #    genus.table$aicc.BM<-NA
####  #    genus.table$sigsq.EB<-NA
####  #    genus.table$aicc.EB<-NA
####  #    genus.table$a.EB<-NA
####  #    genus.table$z0.EB<-NA
####  #    genus.table$z0.BM<-NA
####  #  }
####  #  
####  #}
####  return(genus.table)
####}

##########this runs MS across phylogroups
####run_MS_BM_phylogroup_select_size<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####    return(genus.table)
####  }
####  genus.table$lambda.ms0<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0)
####  genus.table$lambda.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)/(1-0.5)
####  genus.table$lambda.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)/(1-0.9)
####  genus.table$ndr.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)
####  genus.table$ndr.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if(length((trait))>1){
####    res.trait1<-fitContinuous(phylogroup,trait,model="BM")
####    genus.table$sigsq<-res.trait1$opt$sigsq
####    genus.table$z0<-res.trait1$opt$z0
####  }else{
####    genus.table$sigsq<-NA
####    genus.table$z0<-NA
####  }
####  return(genus.table)
####}
####
##########this runs MS + trait with BM & EB across phylogroups
####run_MS_EB_phylogroup_select_size<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####    return(genus.table)
####  }
####  genus.table$lambda.ms0<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0)
####  genus.table$lambda.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)/(1-0.5)
####  genus.table$lambda.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)/(1-0.9)
####  genus.table$ndr.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)
####  genus.table$ndr.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)
####  trait<-as.numeric(genus.table$Log10.Seed.Weight)
####  names(trait)<-genus.table$tip
####  trait<-trait[!is.na(trait)]
####  if (length(trait)<2){
####    genus.table$sigsq.BM<-NA
####    genus.table$aicc.BM<-NA
####    genus.table$sigsq.EB<-NA
####    genus.table$aicc.EB<-NA
####    genus.table$a.EB<-NA
####  }else{
####    name.check<-name.check(phylogroup,trait)
####    if(name.check!='OK'){
####      tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
####    }else{
####      tips.with.data<-phylogroup
####    }
####    if(length(tips.with.data$tip.label)>mincladesize){
####      res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
####      res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
####      genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
####      genus.table$aicc.BM<-res.trait1.BM$opt$aicc
####      genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
####      genus.table$aicc.EB<-res.trait1.EB$opt$aicc
####      genus.table$a.EB<-res.trait1.EB$opt$a
####    }else{
####      genus.table$sigsq.BM<-NA
####      genus.table$aicc.BM<-NA
####      genus.table$sigsq.EB<-NA
####      genus.table$aicc.EB<-NA
####      genus.table$a.EB<-NA
####    }
####      
####  }
####  return(genus.table)
####}
######this runs MS + trait with BM & EB across phylogroups
####run_MS_EB_z0_phylogroup_select_lat<-function(phylogroup,table,sampling,mincladesize){
####  genus.table<-table[match(phylogroup$tip.label,table$tip),]
####  cat(genus.table$New_Species[1],'\n')
####  #to get the sampling fraction, I'll use a weighted mean of the genus.species.sampling.fraction
####  genus.table.genera<-genus.table[,c('genus','genus.species.sampling.fraction')]
####  genus.table.genera<-unique(genus.table.genera)
####  genus.counts<-as.data.frame(table(as.character(genus.table$genus)))
####  colnames(genus.counts)<-c('genus','n.of.tips')
####  genus.table.genera<-merge(genus.table.genera,genus.counts)
####  sampling.fraction<-weighted.mean(genus.table.genera$genus.species.sampling.fraction,w=genus.table.genera$n.of.tips)
####  genus.table$phylogroup.sampling.fraction<-sampling.fraction
####  genus.table$phylogroup.crowncapture<-crownCaptureProbability(sum(genus.table.genera$n.of.tips)/sampling.fraction,sum(genus.table.genera$n.of.tips))
####  tot_time<-max(node.age(phylogroup)$ages)
####  genus.table$phylogroup.size<-nrow(genus.table)
####  genus.table$phylogroup.name<-phylogroup$tip.label[1]
####  #skip for small or undersampled clades
####  if(genus.table$phylogroup.size[1]<mincladesize){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    #genus.table$sigsq.BM<-NA
####    #genus.table$aicc.BM<-NA
####    #genus.table$sigsq.EB<-NA
####    #genus.table$aicc.EB<-NA
####    #genus.table$a.EB<-NA
####    #genus.table$z0.EB<-NA
####    #genus.table$z0.BM<-NA
####    return(genus.table)
####  }
####  if(genus.table$phylogroup.sampling.fraction[1]<sampling){
####    genus.table$lambda.ms0<-NA
####    genus.table$lambda.ms05<-NA
####    genus.table$lambda.ms09<-NA
####    genus.table$ndr.ms05<-NA
####    genus.table$ndr.ms09<-NA
####    #genus.table$sigsq.BM<-NA
####    #genus.table$aicc.BM<-NA
####    #genus.table$sigsq.EB<-NA
####    #genus.table$aicc.EB<-NA
####    #genus.table$a.EB<-NA
####    #genus.table$z0.EB<-NA
####    #genus.table$z0.BM<-NA
####    return(genus.table)
####  }
####  genus.table$lambda.ms0<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0)
####  genus.table$lambda.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)/(1-0.5)
####  genus.table$lambda.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)/(1-0.9)
####  genus.table$ndr.ms05<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.5)
####  genus.table$ndr.ms09<-bd.ms(time=tot_time,n=(nrow(genus.table))/sampling.fraction,epsilon=0.9)
####  #trait<-as.numeric(genus.table$Median.Latitude)
####  #names(trait)<-genus.table$tip
####  #trait<-trait[!is.na(trait)]
####  #if (length(trait)<2){
####  #  genus.table$sigsq.BM<-NA
####  #  genus.table$aicc.BM<-NA
####  #  genus.table$sigsq.EB<-NA
####  #  genus.table$aicc.EB<-NA
####  #  genus.table$a.EB<-NA
####  #  genus.table$z0.EB<-NA
####  #  genus.table$z0.BM<-NA
####  #}else{
####  #  name.check<-name.check(phylogroup,trait)
####  #  if(name.check!='OK'){
####  #    tips.with.data<-drop.tip(phylogroup,name.check$tree_not_data)
####  #  }else{
####  #    tips.with.data<-phylogroup
####  #  }
####  #  if(length(tips.with.data$tip.label)>mincladesize){
####  #    res.trait1.BM<-fitContinuous(phylogroup,trait,model="BM")
####  #    res.trait1.EB<-fitContinuous(phylogroup,trait,model="EB")
####  #    genus.table$sigsq.BM<-res.trait1.BM$opt$sigsq
####  #    genus.table$aicc.BM<-res.trait1.BM$opt$aicc
####  #    genus.table$sigsq.EB<-res.trait1.EB$opt$sigsq
####  #    genus.table$aicc.EB<-res.trait1.EB$opt$aicc
####  #    genus.table$a.EB<-res.trait1.EB$opt$a
####  #    genus.table$z0.EB<-res.trait1.EB$opt$z0
####  #    genus.table$z0.BM<-res.trait1.BM$opt$z0
####  #    
####  #  }else{
####  #    genus.table$sigsq.BM<-NA
####  #    genus.table$aicc.BM<-NA
####  #    genus.table$sigsq.EB<-NA
####  #    genus.table$aicc.EB<-NA
####  #    genus.table$a.EB<-NA
####  #    genus.table$z0.EB<-NA
####  #    genus.table$z0.BM<-NA
####  #    
####  #  }
####  #  
####  #}
####  return(genus.table)
####}
####
###wrapper to run clade analysis (with Magallon & Sanderson estimator) + pgls
#this creates the phylogroups (there's another script to run if phyogroups have previously been saved)
####run_clades_MS_pgls<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  #save phylogroup object to Rsave file
####  save(phylogroups,file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_BM_phylogroup_select_size(x,table,mincladesize,sampling))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############this is for lambda
####  pdf(paste('./output/plots/clade_analyses_new/MSlambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSlambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  ###############this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/MSndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.ms0,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms0','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####}
#this will load presaved phylogroups
####run_clades_MS_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  #phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_BM_phylogroup_select_size(x,table,mincladesize=mincladesize,sampling=sampling))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############this is for lambda
####  pdf(paste('./output/plots/clade_analyses_new/MSlambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSlambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  ###############this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/MSndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####}
####
#####this will load presaved phylogroups
####run_clades_MS_EB_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  #phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_EB_phylogroup_select_size(x,table,mincladesize=mincladesize,sampling=sampling))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #AIC model selection,AICcBM-AICcEB
####  df.phylogroups$diff.AICc.BM.EB<-df.phylogroups$aicc.BM-df.phylogroups$aicc.EB
####  df.phylogroups$best.model.sigsq<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'sigsq.BM']
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'sigsq.EB']
####  df.phylogroups$output.best.model<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'output.best.model']<-'BM'
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'output.best.model']<-'EB'
####  output.best.model<-as.data.frame(table(df.phylogroups$output.best.model))
####  write.table(output.best.model,file=paste('./output/clade_analyses/MS_EB/MS_EB_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'best_fit_BMEB_model.txt',sep=''),quote=F,sep='\t',row.names=F)
####  
####  ###############this is for lambda
####  pdf(paste('./output/clade_analyses/MS_EB/MS_EB_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$best.model.sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/MS_EB/MS_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  ###############this is for ndr
####  pdf(paste('./output/clade_analyses/MS_EB/MS_EB_ndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$best.model.sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/MS_EB/MS_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_ndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####}



#this will load presaved phylogroups
#####run_clades_MS_EB_z0_pgls_savedphylogroups_abslat<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  #phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  phylogroups<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'.RDS',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_EB_z0_phylogroup_select_lat(x,table,mincladesize=mincladesize,sampling=sampling))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of lat data points per phylogroup
####  phylogroup.species.lat<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.lat[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude)])))
####  }
####  phylogroup.species.lat<-as.data.frame(phylogroup.species.lat,stringsAsFactors = F)
####  colnames(phylogroup.species.lat)<-c('phylogroup.name','phylogroup.species.with.lat.data')
####  phylogroup.species.lat$phylogroup.species.with.lat.data<-as.numeric(phylogroup.species.lat$phylogroup.species.with.lat.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.lat,all.x=TRUE)
####  
####  #get the mean lat in each phylogroup
####  mean.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  df.phylogroups<-merge(df.phylogroups,mean.lat.phylogroups)
####  phylogroups.counts.stricts.df<-matrix(NA,ncol=4,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  colnames(phylogroups.counts.stricts.df)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup','phylogroup.name')
####  #get the counts of strict tropical
####  for(i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.table<-df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]
####    phylogroup.table<-phylogroup.table[complete.cases(phylogroup.table),]
####    strict.counts<-count(phylogroup.table, c("phylogroup.name","strict.tropical"))
####    strict.counts<-strict.counts[,c('strict.tropical','freq')]
####    if((0%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(0,0))
####    }
####    if((1%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(1,0))
####    }
####    if((2%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(2,0))
####    }
####    phylogroup.table.strict.counts<-as.data.frame(matrix(NA,nrow=1,ncol=3),stringsAsFactors = F)
####    colnames(phylogroup.table.strict.counts)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup')
####    phylogroup.table.strict.counts$n.strict.temperate.species.phylogroup<-strict.counts[strict.counts$strict.tropical==0,'freq']
####    phylogroup.table.strict.counts$n.strict.tropical.species.phylogroup<-strict.counts[strict.counts$strict.tropical==1,'freq']
####    phylogroup.table.strict.counts$n.strict.widespread.species.phylogroup<-strict.counts[strict.counts$strict.tropical==2,'freq']
####    phylogroup.table.strict.counts$phylogroup.name<-phylogroup.table$phylogroup.name[1]
####    phylogroups.counts.stricts.df<-rbind(phylogroups.counts.stricts.df,phylogroup.table.strict.counts)
####  }
####  count.strict.tropical.phylogroups<-aggregate(strict.tropical~phylogroup.name,data=df.phylogroups,count)
####  phylogroups.counts.stricts.df<-phylogroups.counts.stricts.df[complete.cases(phylogroups.counts.stricts.df),]
####  df.phylogroups<-merge(df.phylogroups,phylogroups.counts.stricts.df,all.x=T)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  
####  
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.lat.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #AIC model selection,AICcBM-AICcEB
####  
####  #strict tropical phylogroups =  phylogroups with > .8 of species strict.tropical
####  #strict temperate phylogroups =  phylogroups with > .8 of species strict.tropical
####  df.phylogroups$strict.tropical.character<-NA
####  df.phylogroups[(df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>.7,'strict.tropical.character']<-1
####  df.phylogroups[(df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>.7,'strict.tropical.character']<-0
####  
####  ###############this is for lambda
####  pdf(paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/MS_EB_z0/MS_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Median.Latitude,xlab='phylogroup.Median.Latitude')
####  wilcox.test.trop.vs.temp<-wilcox.test(df.phylogroups[df.phylogroups$strict.tropical.character==0,'lambda.ms05'],df.phylogroups[df.phylogroups$strict.tropical.character==1,'lambda.ms05'])
####  boxplot(df.phylogroups[df.phylogroups$strict.tropical.character==0,'lambda.ms05'],df.phylogroups[df.phylogroups$strict.tropical.character==1,'lambda.ms05'],main=paste('strict.temp vs strict.trop (wilcox test pval = ',round(wilcox.test.trop.vs.temp$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups[df.phylogroups$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups[df.phylogroups$strict.tropical.character==1,]),')',sep='')),notch=T)
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09')
####  
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.Median.Latitude.lambda<-try(pgls(log10(lambda.ms05)~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  pgls.Median.Latitude.lambda<-try(pgls(log10(lambda.ms05)~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  #pgls.Median.Latitude.sigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.Median.Latitude.lambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-90,60),ylim=c(-2,1))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.lambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.Median.Latitude.lambda<-gls(log10(lambda.ms05) ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-90,90),ylim=c(-2,1))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.lambda',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],8),round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],8))
####  }
####  
####  ##if(class(pgls.seedsigsq)!='try-error'){
####  ##  plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####  ##  abline(pgls.seedsigsq)
####  ##  pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  ##}
####  ##if(class(pgls.seedsigsq)=='try-error'){
####  ##  row.names(phylogroups.data)<-phylogroups.data$Tip
####  ##  pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####  ##  plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####  ##  abline(pgls.seedsigsq)
####  ##  pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  ##}
####  
####  results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/MS_EB_z0/MS_EB_z0_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results_medianlatitude.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  ##this is for strict trop vs strict temp
####  #delete na (phylogroups that are not strict trop or strict temp)
####  df.phylogroups<-df.phylogroups[complete.cases(df.phylogroups),]
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$strict.tropical.character,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09)
####  colnames(phylogroups.data)<-c('Tip','strict.tropical.character','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09')
####  
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.strict.tropical.character.lambda<-try(pgls(log10(lambda.ms05)~as.character(strict.tropical.character), data = phylogroups.tree.object, lambda='ML'))
####  pgls.strict.tropical.characterlambda<-try(pgls(log10(lambda.ms05)~as.character(strict.tropical.character), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  #pgls.Median.Latitude.sigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.strict.tropical.character.lambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$strict.tropical.character,xlab='strict.tropical',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.strict.tropical.character.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.strict.tropical.character.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-1,2),ylim=c(-2,1))
####    abline(pgls.strict.tropical.character.lambda)
####    pgls.strict.tropical.character.lambda.lambdaresults<-c('pgls.strict.tropical.character',round(summary(pgls.strict.tropical.character.lambda)$coefficients[2,1],8),round(summary(pgls.strict.tropical.character.lambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.lambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.strict.tropical.character.lambda<-gls(log10(lambda.ms05) ~ as.character(strict.tropical.character), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$strict.tropical.character,xlab='strict.tropical',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.strict.tropical.character.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.strict.tropical.character.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-1,2),ylim=c(-2,1))
####    abline(pgls.strict.tropical.character.lambda)
####    pgls.strict.tropical.character.lambda.lambdaresults<-c('pgls.strict.tropical.character.lambda',round(summary(pgls.strict.tropical.character.lambda)$tTable[2,1],8),round(summary(pgls.strict.tropical.character.lambda)$tTable[2,4],8))
####  }
####  
####  results.df<-as.data.frame(rbind(pgls.strict.tropical.character.lambda.lambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/MS_EB_z0/MS_EB_z0_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results_stricttropical.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  dev.off()
####  
####  ###############this is for ndr
####  ##pdf(paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/MS_EB_z0/MS_EB_z0_ndr_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  ##par(mfrow=c(2,2))
####  ##hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  ##hist(df.phylogroups$phylogroup.Median.Latitude,xlab='phylogroup.Median.Latitude')
####  ##phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  ##phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$z0,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$best.model.sigsq)
####  ##colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  ##phylogroups.tree$node.label<-NULL
####  ###this is for caper
####  ##phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  ###for mean family seed weight
####  ##cat('running pgls','\n')
####  ##pgls.Median.Latitudendr<-try(pgls(log10(ndr.ms05)~phylogroup.Median.Latitude, data = phylogroups.tree.object, lambda='ML'))
####  ###summary(pgls.seedndr)
####  ###pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  ###summary(pgls.seedsigsq)
####  ######try both optimisations (caper + nlme if caper gives an error)
####  ##if(class(pgls.Median.Latitudendr)!='try-error'){
####  ##  plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Median.Latitude,xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitudendr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitudendr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-90,60),ylim=c(-2,1))
####  ##  abline(pgls.Median.Latitudendr)
####  ##  pgls.Median.Latitudendrresults<-c('pgls.Median.Latitudendr',round(summary(pgls.Median.Latitudendr)$coefficients[2,1],8),round(summary(pgls.Median.Latitudendr)$coefficients[2,4],8))
####  ##}
####  ##if(class(pgls.Median.Latitudendr)=='try-error'){
####  ##  row.names(phylogroups.data)<-phylogroups.data$Tip
####  ##  pgls.Median.Latitudendr<-gls(log10(ndr.ms05) ~ phylogroup.Median.Latitude, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####  ##  plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Median.Latitude,xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitudendr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitudendr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-90,60),ylim=c(-2,1))
####  ##  abline(pgls.Median.Latitudendr)
####  ##  pgls.Median.Latitudendrresults<-c('pgls.Median.Latitudendr',round(summary(pgls.Median.Latitudendr)$tTable[2,1],8),round(summary(pgls.Median.Latitudendr)$tTable[2,4],8))
####  ##}
####  ###if(class(pgls.seedsigsq)!='try-error'){
####  ###  plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####  ###  abline(pgls.seedsigsq)
####  ###  pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  ###}
####  ###if(class(pgls.seedsigsq)=='try-error'){
####  ###  row.names(phylogroups.data)<-phylogroups.data$Tip
####  ###  pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####  ###  plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####  ###  abline(pgls.seedsigsq)
####  ###  pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  ###}
####  ##
####  ##dev.off()
####  ##results.df<-as.data.frame(rbind(pgls.Median.Latitudendrresults))
####  ##colnames(results.df)<-c('analysis','slope','pvalue')
####  ##rownames(results.df)<-NULL
####  ##write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/MS_EB_z0/MS_EB_z0_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_ndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####}
####
####


####run_clades_RPANDA_EB_pgls_savedphylogroups_6models<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_EB_phylogroup_select_size_6models(x,table,sampling,mincladesize))
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #AIC model selection,AICcBM-AICcEB
####  df.phylogroups$diff.AICc.BM.EB<-df.phylogroups$aicc.BM-df.phylogroups$aicc.EB
####  df.phylogroups$best.model.sigsq<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'sigsq.BM']
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'sigsq.EB']
####  
####  df.phylogroups$output.best.model<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'output.best.model']<-'BM'
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'output.best.model']<-'EB'
####  output.best.model<-as.data.frame(table(df.phylogroups$output.best.model))
####  write.table(output.best.model,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_lambda_6models_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'best_fit_BMEB_model.txt',sep=''),quote=F,sep='\t',row.names=F)
####  
####  
####  ###############correlation with lambda
####  pdf(paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_lambda_6models_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$best.model.sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #output a table with the best fitting models for RPANDA for each clade
####  best.fitting.RPANDA<-as.data.frame(table(df.phylogroups$aicc),stringsAsFactors = F)
####  write.table(best.fitting.RPANDA,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_6models_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_bestfitmodel_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_6models_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  ###########this is for ndr
####  pdf(paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_ndr_analysis_clades_6models_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$best.model.sigsq)
####  
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #####this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_6models_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_ndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####}
####run_clades_RPANDA_EB_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_EB_phylogroup_select_size(x,table,sampling,mincladesize))
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #AIC model selection,AICcBM-AICcEB
####  df.phylogroups$diff.AICc.BM.EB<-df.phylogroups$aicc.BM-df.phylogroups$aicc.EB
####  df.phylogroups$best.model.sigsq<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'sigsq.BM']
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'sigsq.EB']
####  
####  df.phylogroups$output.best.model<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'output.best.model']<-'BM'
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'output.best.model']<-'EB'
####  output.best.model<-as.data.frame(table(df.phylogroups$output.best.model))
####  write.table(output.best.model,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'best_fit_BMEB_model.txt',sep=''),quote=F,sep='\t',row.names=F)
####  
####  
####  ###############correlation with lambda
####  pdf(paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_lambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$best.model.sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #output a table with the best fitting models for RPANDA for each clade
####  best.fitting.RPANDA<-as.data.frame(table(df.phylogroups$aicc),stringsAsFactors = F)
####  write.table(best.fitting.RPANDA,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_bestfitmodel_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  ###########this is for ndr
####  pdf(paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_ndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$best.model.sigsq)
####  
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #####this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/RPANDA_EB/RPANDA_EB_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_ndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####}

#####run_clades_RPANDA_pgls_savedphylogroups_abslat<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  phylogroups<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'.RDS',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  #phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_EB_z0_phylogroup_select_lat(x,table,sampling,mincladesize))
####  #phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_EB_z0_phylogroup_select_lat(x,table,sampling,mincladesize))
####  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_phylogroup(x,table,sampling,mincladesize))
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of lat data points per phylogroup
####  phylogroup.species.lat<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.lat[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude)])))
####  }
####  phylogroup.species.lat<-as.data.frame(phylogroup.species.lat,stringsAsFactors = F)
####  colnames(phylogroup.species.lat)<-c('phylogroup.name','phylogroup.species.with.lat.data')
####  phylogroup.species.lat$phylogroup.species.with.lat.data<-as.numeric(phylogroup.species.lat$phylogroup.species.with.lat.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.lat,all.x=TRUE)
####  
####  #get the mean lat in each phylogroup
####  mean.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  df.phylogroups<-merge(df.phylogroups,mean.lat.phylogroups)
####  phylogroups.counts.stricts.df<-matrix(NA,ncol=4,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  colnames(phylogroups.counts.stricts.df)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup','phylogroup.name')
####  #get the counts of strict tropical
####  for(i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.table<-df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]
####    phylogroup.table<-phylogroup.table[complete.cases(phylogroup.table),]
####    strict.counts<-count(phylogroup.table, c("phylogroup.name","strict.tropical"))
####    strict.counts<-strict.counts[,c('strict.tropical','freq')]
####    if((0%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(0,0))
####    }
####    if((1%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(1,0))
####    }
####    if((2%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(2,0))
####    }
####    phylogroup.table.strict.counts<-as.data.frame(matrix(NA,nrow=1,ncol=3),stringsAsFactors = F)
####    colnames(phylogroup.table.strict.counts)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup')
####    phylogroup.table.strict.counts$n.strict.temperate.species.phylogroup<-strict.counts[strict.counts$strict.tropical==0,'freq']
####    phylogroup.table.strict.counts$n.strict.tropical.species.phylogroup<-strict.counts[strict.counts$strict.tropical==1,'freq']
####    phylogroup.table.strict.counts$n.strict.widespread.species.phylogroup<-strict.counts[strict.counts$strict.tropical==2,'freq']
####    phylogroup.table.strict.counts$phylogroup.name<-phylogroup.table$phylogroup.name[1]
####    phylogroups.counts.stricts.df<-rbind(phylogroups.counts.stricts.df,phylogroup.table.strict.counts)
####  }
####  count.strict.tropical.phylogroups<-aggregate(strict.tropical~phylogroup.name,data=df.phylogroups,count)
####  phylogroups.counts.stricts.df<-phylogroups.counts.stricts.df[complete.cases(phylogroups.counts.stricts.df),]
####  df.phylogroups<-merge(df.phylogroups,phylogroups.counts.stricts.df,all.x=T)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  
####  
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.lat.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #AIC model selection,AICcBM-AICcEB
####  
####  #strict tropical phylogroups =  phylogroups with > .8 of species strict.tropical
####  #strict temperate phylogroups =  phylogroups with > .8 of species strict.tropical
####  df.phylogroups$strict.tropical.character<-NA
####  df.phylogroups[(df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>.7,'strict.tropical.character']<-1
####  df.phylogroups[(df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>.7,'strict.tropical.character']<-0
####  
####  ###############this is for lambda
####  pdf(paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Median.Latitude,xlab='phylogroup.Median.Latitude')
####  wilcox.test.trop.vs.temp<-wilcox.test(df.phylogroups[df.phylogroups$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups[df.phylogroups$strict.tropical.character==1,'lambda.rpanda1'])
####  boxplot(df.phylogroups[df.phylogroups$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups[df.phylogroups$strict.tropical.character==1,'lambda.rpanda1'],main=paste('strict.temp vs strict.trop (wilcox test pval = ',round(wilcox.test.trop.vs.temp$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups[df.phylogroups$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups[df.phylogroups$strict.tropical.character==1,]),')',sep='')),notch=T)
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$lambda.rpanda1)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','lambda.rpanda1')
####  
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.Median.Latitude.lambda<-try(pgls(log10(lambda.rpanda1)~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  pgls.Median.Latitude.lambda<-try(pgls(log10(lambda.rpanda1)~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  #pgls.Median.Latitude.sigsq<-try(pgls(log10(lambda.rpanda1)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.Median.Latitude.lambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda1)~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-90,60),ylim=c(-2,1))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.lambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.Median.Latitude.lambda<-gls(log10(lambda.rpanda1) ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda1)~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-90,90),ylim=c(-2,1))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.lambda',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],8),round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],8))
####  }
####  
####  ##if(class(pgls.seedsigsq)!='try-error'){
####  ##  plot(log10(phylogroups.tree.object$data$lambda.rpanda1)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####  ##  abline(pgls.seedsigsq)
####  ##  pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  ##}
####  ##if(class(pgls.seedsigsq)=='try-error'){
####  ##  row.names(phylogroups.data)<-phylogroups.data$Tip
####  ##  pgls.seedsigsq<-gls(log10(lambda.rpanda1) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####  ##  plot(log10(phylogroups.data$lambda.rpanda1)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####  ##  abline(pgls.seedsigsq)
####  ##  pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  ##}
####  
####  results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results_medianlatitude.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  ##this is for strict trop vs strict temp
####  #delete na (phylogroups that are not strict trop or strict temp)
####  df.phylogroups<-df.phylogroups[complete.cases(df.phylogroups),]
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$strict.tropical.character,df.phylogroups$lambda.rpanda1,df.phylogroups$ndr.rpanda1)
####  colnames(phylogroups.data)<-c('Tip','strict.tropical.character','lambda.rpanda1','ndr.rpanda1')
####  
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.strict.tropical.character.lambda<-try(pgls(log10(lambda.rpanda1)~as.character(strict.tropical.character), data = phylogroups.tree.object, lambda='ML'))
####  pgls.strict.tropical.characterlambda<-try(pgls(log10(lambda.rpanda1)~as.character(strict.tropical.character), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  #pgls.Median.Latitude.sigsq<-try(pgls(log10(lambda.rpanda1)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.strict.tropical.character.lambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda1)~phylogroups.tree.object$data$strict.tropical.character,xlab='strict.tropical',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.strict.tropical.character.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.strict.tropical.character.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-1,2),ylim=c(-2,1))
####    abline(pgls.strict.tropical.character.lambda)
####    pgls.strict.tropical.character.lambda.lambdaresults<-c('pgls.strict.tropical.character',round(summary(pgls.strict.tropical.character.lambda)$coefficients[2,1],8),round(summary(pgls.strict.tropical.character.lambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.lambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.strict.tropical.character.lambda<-gls(log10(lambda.rpanda1) ~ as.character(strict.tropical.character), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda1)~phylogroups.data$strict.tropical.character,xlab='strict.tropical',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.strict.tropical.character.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.strict.tropical.character.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-1,2),ylim=c(-2,1))
####    abline(pgls.strict.tropical.character.lambda)
####    pgls.strict.tropical.character.lambda.lambdaresults<-c('pgls.strict.tropical.character.lambda',round(summary(pgls.strict.tropical.character.lambda)$tTable[2,1],8),round(summary(pgls.strict.tropical.character.lambda)$tTable[2,4],8))
####  }
####  
####  results.df<-as.data.frame(rbind(pgls.strict.tropical.character.lambda.lambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_lambda_results_stricttropical.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  dev.off()
####  
####}  
####
####run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold<-function(tree,minage,maxage,mincladesize,ncores,sampling,table,strict.tropical.threshold,GBIF.sampling,loadRPANDA){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  if(loadRPANDA==F){
####    phylogroups<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'.RDS',sep=''))
####    phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####    length(unique(phylogroups.species))
####    cat('measuring diversification/trait evolution','\n')
####    phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_phylogroup(x,table,sampling,mincladesize=mincladesize))
####    
####  }else if (loadRPANDA==T){
####    cat('loading RPANDA phylogroups object','\n')
####    phylogroups.RPANDA<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'_RPANDA.RDS',sep=''))
####  }
####  # phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_EB_z0_phylogroup_select_lat(x,table,sampling,mincladesize))
####  #phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_linear_phylogroup_select_lat(x,table,sampling,mincladesize))
####  
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of lat data points per phylogroup
####  phylogroup.species.lat<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.lat[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude)])))
####  }
####  phylogroup.species.lat<-as.data.frame(phylogroup.species.lat,stringsAsFactors = F)
####  colnames(phylogroup.species.lat)<-c('phylogroup.name','phylogroup.species.with.lat.data')
####  phylogroup.species.lat$phylogroup.species.with.lat.data<-as.numeric(phylogroup.species.lat$phylogroup.species.with.lat.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.lat,all.x=TRUE)
####  
####  #get the mean lat in each phylogroup
####  #mean.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,mean)
####  #colnames(mean.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  #df.phylogroups<-merge(df.phylogroups,mean.lat.phylogroups)
####
####  #median is very correlated with mean
####  median.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,median)
####  colnames(median.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  df.phylogroups<-merge(df.phylogroups,median.lat.phylogroups)
####  
####  #sd.Median.Latitude is quite big in phylogroups
####  sd.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,function(x) sd(abs(x)))
####  colnames(sd.lat.phylogroups)[2]<-'phylogroup.sd.Median.Latitude'
####  sd.lat.phylogroups$phylogroup.sd.Median.Latitude[is.na(sd.lat.phylogroups$phylogroup.sd.Median.Latitude)]<-0
####  df.phylogroups<-merge(df.phylogroups,sd.lat.phylogroups)
####  
####  phylogroups.counts.stricts.df<-matrix(NA,ncol=4,nrow=0)
####  colnames(phylogroups.counts.stricts.df)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup','phylogroup.name')
####  #get the counts of strict tropical
####  for(d in 1:length(unique(df.phylogroups$phylogroup.name))){
####    cat(d,'\n')
####    phylogroup.table<-df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[d],]
####    phylogroup.table<-phylogroup.table[complete.cases(phylogroup.table),]
####    strict.counts<-count(phylogroup.table, c("phylogroup.name","strict.tropical"))
####    strict.counts<-strict.counts[,c('strict.tropical','freq')]
####    #add missing 0,1,2 counts (e.g if there are not strict tropical -  1s -  add a row with strict tropical =1 and freq =0)
####    if((0%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(0,0))
####    }
####    if((1%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(1,0))
####    }
####    if((2%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(2,0))
####    }
####    phylogroup.table.strict.counts<-as.data.frame(matrix(NA,nrow=1,ncol=3),stringsAsFactors = F)
####    colnames(phylogroup.table.strict.counts)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup')
####    phylogroup.table.strict.counts$n.strict.temperate.species.phylogroup<-strict.counts[strict.counts$strict.tropical==0,'freq']
####    phylogroup.table.strict.counts$n.strict.tropical.species.phylogroup<-strict.counts[strict.counts$strict.tropical==1,'freq']
####    phylogroup.table.strict.counts$n.strict.widespread.species.phylogroup<-strict.counts[strict.counts$strict.tropical==2,'freq']
####    phylogroup.table.strict.counts$phylogroup.name<-phylogroup.table$phylogroup.name[1]
####    phylogroups.counts.stricts.df<-rbind(phylogroups.counts.stricts.df,phylogroup.table.strict.counts)
####  }
####  count.strict.tropical.phylogroups<-aggregate(strict.tropical~phylogroup.name,data=df.phylogroups,count)
##### phylogroups.counts.stricts.df<-phylogroups.counts.stricts.df[complete.cases(phylogroups.counts.stricts.df),]
####  df.phylogroups<-merge(df.phylogroups,phylogroups.counts.stricts.df,all.x=T)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  
####  df.phylogroups$turnover.rpanda1<-df.phylogroups$lambda.rpanda1+df.phylogroups$mu.rpanda1
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #filter by GBIF sampling
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>GBIF.sampling,]
####  
####  #dataset with strict tropical and strict temperate species
####  #strict tropical phylogroups =  phylogroups with > .7 of species strict.tropical
####  #strict temperate phylogroups =  phylogroups with > .7 of species strict.tropical
####  df.phylogroups$strict.tropical.character<-NA
####  df.phylogroups[(df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-1
####  df.phylogroups[(df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-0
####  
####  df.phylogroups.strict<-df.phylogroups[!is.na(df.phylogroups$strict.tropical.character),]
####  
####  df.phylogroups$proportion.strict.tropical<-df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  df.phylogroups$proportion.strict.temperate<-df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  
####
####  ###############this is for lambda
####  pdf(paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_new_GBIFsampling.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Median.Latitude,xlab='phylogroup.Median.Latitude')
####  if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####    wilcox.test.trop.vs.temp.lambda<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
####    t.test.trop.vs.temp.lambda<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.mu<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
####    t.test.trop.vs.temp.mu<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.ndr<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
####    t.test.trop.vs.temp.ndr<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.turnover<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    t.test.trop.vs.temp.turnover<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    
####    wilcox.t.results<-as.data.frame(rbind(c('lambda',wilcox.test.trop.vs.temp.lambda$p.value,t.test.trop.vs.temp.lambda$p.value),c('mu',wilcox.test.trop.vs.temp.mu$p.value,t.test.trop.vs.temp.mu$p.value),c('ndr',wilcox.test.trop.vs.temp.ndr$p.value,t.test.trop.vs.temp.ndr$p.value),c('turnover',wilcox.test.trop.vs.temp.turnover$p.value,t.test.trop.vs.temp.turnover$p.value)))
####    colnames(wilcox.t.results)<-c('parameter','pvalue.wilcox','pvalue.ttest')
####    write.table(wilcox.t.results,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_wilcox.ttest_table_new_GBIFsampling.txt',sep=''),sep='\t',quote=F,row.names=F)
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],main=paste('lambda strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.lambda$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T)
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],main=paste('mu strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.mu$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T)
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],main=paste('ndr strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.ndr$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T)
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],main=paste('turnover strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.turnover$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T)
####    
####  }
####  
####  
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  #phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$strict.tropical.character,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1)
####  #colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1,df.phylogroups$proportion.strict.tropical,df.phylogroups$proportion.strict.temperate)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1','proportion.strict.tropical','proportion.strict.temperate')
####  
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.Median.Latitude.lambda<-try(pgls(lambda.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  pgls.Median.Latitude.mu<-try(pgls(mu.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  pgls.Median.Latitude.ndr<-try(pgls(ndr.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  pgls.Median.Latitude.turnover<-try(pgls(turnover.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  
####  pgls.prop.trop.lambda<-try(pgls(lambda.rpanda1~proportion.strict.tropical, data = phylogroups.tree.object, lambda='ML'))
####  pgls.prop.trop.mu<-try(pgls(mu.rpanda1~proportion.strict.tropical, data = phylogroups.tree.object, lambda='ML'))
####  pgls.prop.trop.ndr<-try(pgls(ndr.rpanda1~proportion.strict.tropical, data = phylogroups.tree.object, lambda='ML'))
####  pgls.prop.trop.turnover<-try(pgls(turnover.rpanda1~proportion.strict.tropical, data = phylogroups.tree.object, lambda='ML'))
####  
####  pgls.prop.temp.lambda<-try(pgls(lambda.rpanda1~proportion.strict.temperate, data = phylogroups.tree.object, lambda='ML'))
####  pgls.prop.temp.mu<-try(pgls(mu.rpanda1~proportion.strict.temperate, data = phylogroups.tree.object, lambda='ML'))
####  pgls.prop.temp.ndr<-try(pgls(ndr.rpanda1~proportion.strict.temperate, data = phylogroups.tree.object, lambda='ML'))
####  pgls.prop.temp.turnover<-try(pgls(turnover.rpanda1~proportion.strict.temperate, data = phylogroups.tree.object, lambda='ML'))
####  
####  if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####    
####    phylogroups.tree.strict<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups.strict$Tip))
####    phylogroups.data.strict<-data.frame(df.phylogroups.strict$Tip,df.phylogroups.strict$phylogroup.Median.Latitude,df.phylogroups.strict$strict.tropical.character,df.phylogroups.strict$lambda.rpanda1,df.phylogroups.strict$mu.rpanda1,df.phylogroups.strict$ndr.rpanda1,df.phylogroups.strict$turnover.rpanda1)
####    colnames(phylogroups.data.strict)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####    phylogroups.tree.strict.object<-comparative.data(phy = phylogroups.tree.strict,data=phylogroups.data.strict,names.col='Tip',vcv=TRUE)
####    
####    pgls.strict.tropical.lambda<-try(pgls(lambda.rpanda1~strict.tropical.character, data = phylogroups.tree.strict.object, lambda='ML'))
####    pgls.strict.tropical.mu<-try(pgls(mu.rpanda1~strict.tropical.character, data = phylogroups.tree.strict.object, lambda='ML'))
####    pgls.strict.tropical.ndr<-try(pgls(ndr.rpanda1~strict.tropical.character, data = phylogroups.tree.strict.object, lambda='ML'))
####    pgls.strict.tropical.turnover<-try(pgls(turnover.rpanda1~strict.tropical.character, data = phylogroups.tree.strict.object, lambda='ML'))
####    
####    strict.tropical.phylanova<-phylogroups.tree.strict.object$data$strict.tropical.character
####    names(strict.tropical.phylanova)<-rownames(phylogroups.tree.strict.object$data)
####    
####    strict.tropical.phylanova.lambda<-phylogroups.tree.strict.object$data$lambda.rpanda1
####    names(strict.tropical.phylanova.lambda)<-rownames(phylogroups.tree.strict.object$data)
####    
####    strict.tropical.phylanova.mu<-phylogroups.tree.strict.object$data$mu.rpanda1
####    names(strict.tropical.phylanova.mu)<-rownames(phylogroups.tree.strict.object$data)
####    
####    strict.tropical.phylanova.ndr<-phylogroups.tree.strict.object$data$ndr.rpanda1
####    names(strict.tropical.phylanova.ndr)<-rownames(phylogroups.tree.strict.object$data)
####    
####    strict.tropical.phylanova.turnover<-phylogroups.tree.strict.object$data$turnover.rpanda1
####    names(strict.tropical.phylanova.turnover)<-rownames(phylogroups.tree.strict.object$data)
####    
####    phylanova.strict.tropical.lambda<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.lambda,nsim = 100,posthoc = T)
####    phylanova.strict.tropical.mu<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.mu,nsim = 100,posthoc = T)
####    phylanova.strict.tropical.ndr<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.ndr,nsim = 100,posthoc = T)
####    phylanova.strict.tropical.turnover<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.turnover,nsim = 100,posthoc = T)
####  }
####  #summary(pgls.seedlambda)
####  #pgls.Median.Latitude.sigsq<-try(pgls(log10(lambda.rpanda1)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if((class(pgls.Median.Latitude.lambda)!='try-error')&(class(pgls.Median.Latitude.mu)!='try-error')&(class(pgls.Median.Latitude.ndr)!='try-error')&(class(pgls.Median.Latitude.turnover)!='try-error')&(class(pgls.prop.trop.lambda)!='try-error')&(class(pgls.prop.trop.mu)!='try-error')&(class(pgls.prop.trop.ndr)!='try-error')&(class(pgls.prop.trop.turnover)!='try-error')&(class(pgls.prop.temp.lambda)!='try-error')&(class(pgls.prop.temp.mu)!='try-error')&(class(pgls.prop.temp.ndr)!='try-error')&(class(pgls.prop.temp.turnover)!='try-error')){
####    plot(phylogroups.tree.object$data$lambda.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    plot(phylogroups.tree.object$data$mu.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))  
####    abline(pgls.Median.Latitude.mu)
####    plot(phylogroups.tree.object$data$ndr.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    plot(phylogroups.tree.object$data$turnover.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))  
####    abline(pgls.Median.Latitude.turnover)
####    
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.lambda',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],8))
####    pgls.Median.Latitude.muresults<-c('pgls.Median.Latitude.mu',round(summary(pgls.Median.Latitude.mu)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.mu)$coefficients[2,4],8))
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],8))
####    pgls.Median.Latitude.turnoverresults<-c('pgls.Median.Latitude.turnover',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.turnover)$coefficients[2,4],8))
####    
####    plot(phylogroups.tree.object$data$lambda.rpanda1~phylogroups.tree.object$data$proportion.strict.tropical,xlab='prop.trop',ylab='lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    plot(phylogroups.tree.object$data$mu.rpanda1~phylogroups.tree.object$data$proportion.strict.tropical,xlab='prop.trop',ylab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)  
####    abline(pgls.prop.trop.mu)
####    plot(phylogroups.tree.object$data$ndr.rpanda1~phylogroups.tree.object$data$proportion.strict.tropical,xlab='prop.trop',ylab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    plot(phylogroups.tree.object$data$turnover.rpanda1~phylogroups.tree.object$data$proportion.strict.tropical,xlab='prop.trop',ylab='turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)  
####    abline(pgls.prop.trop.turnover)
####    
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.lambda',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],8),round(summary(pgls.prop.trop.lambda)$coefficients[2,4],8))
####    pgls.prop.trop.muresults<-c('pgls.prop.trop.mu',round(summary(pgls.prop.trop.mu)$coefficients[2,1],8),round(summary(pgls.prop.trop.mu)$coefficients[2,4],8))
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],8),round(summary(pgls.prop.trop.ndr)$coefficients[2,4],8))
####    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.turnover',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],8),round(summary(pgls.prop.trop.turnover)$coefficients[2,4],8))
####    
####    plot(phylogroups.tree.object$data$lambda.rpanda1~phylogroups.tree.object$data$proportion.strict.temperate,xlab='prop.temp',ylab='lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    plot(phylogroups.tree.object$data$mu.rpanda1~phylogroups.tree.object$data$proportion.strict.temperate,xlab='prop.temp',ylab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)  
####    abline(pgls.prop.temp.mu)
####    plot(phylogroups.tree.object$data$ndr.rpanda1~phylogroups.tree.object$data$proportion.strict.temperate,xlab='prop.temp',ylab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    plot(phylogroups.tree.object$data$turnover.rpanda1~phylogroups.tree.object$data$proportion.strict.temperate,xlab='prop.temp',ylab='turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)  
####    abline(pgls.prop.temp.turnover)
####    
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.lambda',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],8),round(summary(pgls.prop.temp.lambda)$coefficients[2,4],8))
####    pgls.prop.temp.muresults<-c('pgls.prop.temp.mu',round(summary(pgls.prop.temp.mu)$coefficients[2,1],8),round(summary(pgls.prop.temp.mu)$coefficients[2,4],8))
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],8),round(summary(pgls.prop.temp.ndr)$coefficients[2,4],8))
####    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.turnover',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],8),round(summary(pgls.prop.temp.turnover)$coefficients[2,4],8))
####    if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####      
####      pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.lambda',round(summary(pgls.strict.tropical.lambda)$coefficients[2,1],8),round(summary(pgls.strict.tropical.lambda)$coefficients[2,4],8))
####      pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$coefficients[2,1],8),round(summary(pgls.strict.tropical.mu)$coefficients[2,4],8))
####      pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$coefficients[2,1],8),round(summary(pgls.strict.tropical.ndr)$coefficients[2,4],8))
####      pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.turnover',round(summary(pgls.strict.tropical.turnover)$coefficients[2,1],8),round(summary(pgls.strict.tropical.turnover)$coefficients[2,4],8))
####      
####      phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####      phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####      phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####      phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####      
####    }
####    
####  }else{
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    
####    pgls.Median.Latitude.lambda<-try(gls(lambda.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    pgls.Median.Latitude.mu<-try(gls(mu.rpanda1~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    pgls.Median.Latitude.ndr<-try(gls(ndr.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    pgls.Median.Latitude.turnover<-try(gls(turnover.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    
####    pgls.prop.trop.lambda<-try(gls(lambda.rpanda1 ~ proportion.strict.tropical, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    pgls.prop.trop.mu<-try(gls(mu.rpanda1~ proportion.strict.tropical, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    pgls.prop.trop.ndr<-try(gls(ndr.rpanda1 ~ proportion.strict.tropical, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    pgls.prop.trop.turnover<-try(gls(turnover.rpanda1 ~ proportion.strict.tropical, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    
####    pgls.prop.temp.lambda<-try(gls(lambda.rpanda1 ~ proportion.strict.temperate, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    pgls.prop.temp.mu<-try(gls(mu.rpanda1~ proportion.strict.temperate, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    pgls.prop.temp.ndr<-try(gls(ndr.rpanda1 ~ proportion.strict.temperate, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    pgls.prop.temp.turnover<-try(gls(turnover.rpanda1 ~ proportion.strict.temperate, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    
####    if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####      
####      pgls.strict.tropical.lambda<-try(gls(lambda.rpanda1 ~ strict.tropical.character, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####      pgls.strict.tropical.mu<-try(gls(mu.rpanda1 ~ strict.tropical.character, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####      pgls.strict.tropical.ndr<-try(gls(ndr.rpanda1 ~ strict.tropical.character, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####      pgls.strict.tropical.turnover<-try(gls(turnover.rpanda1 ~ strict.tropical.character, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####      
####    }
####    plot(phylogroups.data$lambda.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='lambda^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    plot(phylogroups.data$mu.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='mu^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.mu)
####    plot(phylogroups.data$lndr.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='ndr^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    plot(phylogroups.data$turnover.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='turnover^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.turnover)
####    
####    plot(phylogroups.data$lambda.rpanda1~phylogroups.data$proportion.strict.tropical,xlab='prop.trop',ylab='lambda^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    plot(phylogroups.data$mu.rpanda1~phylogroups.data$proportion.strict.tropical,xlab='prop.trop',ylab='mu^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.mu)
####    plot(phylogroups.data$lndr.rpanda1~phylogroups.data$proportion.strict.tropical,xlab='prop.trop',ylab='ndr^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    plot(phylogroups.data$turnover.rpanda1~phylogroups.data$proportion.strict.tropical,xlab='prop.trop',ylab='turnover^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.turnover)
####    
####    plot(phylogroups.data$lambda.rpanda1~phylogroups.data$proportion.strict.temperate,xlab='prop.temp',ylab='lambda^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    plot(phylogroups.data$mu.rpanda1~phylogroups.data$proportion.strict.temperate,xlab='prop.temp',ylab='mu^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.mu)
####    plot(phylogroups.data$lndr.rpanda1~phylogroups.data$proportion.strict.temperate,xlab='prop.temp',ylab='ndr^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    plot(phylogroups.data$turnover.rpanda1~phylogroups.data$proportion.strict.temperate,xlab='prop.temp',ylab='turnover^(1/3)',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.turnover)
####    
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.lambda',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],8),round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],8))
####    pgls.Median.Latitude.muresults<-c('pgls.Median.Latitude.mu',round(summary(pgls.Median.Latitude.mu)$tTable[2,1],8),round(summary(pgls.Median.Latitude.mu)$tTable[2,4],8))
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],8),round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],8))
####    pgls.Median.Latitude.turnoverresults<-c('pgls.Median.Latitude.turnover',round(summary(pgls.Median.Latitude.turnover)$tTable[2,1],8),round(summary(pgls.Median.Latitude.turnover)$tTable[2,4],8))
####    
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.lambda',round(summary(pgls.prop.trop.lambda)$tTable[2,1],8),round(summary(pgls.prop.trop.lambda)$tTable[2,4],8))
####    pgls.prop.trop.muresults<-c('pgls.prop.trop.mu',round(summary(pgls.prop.trop.mu)$tTable[2,1],8),round(summary(pgls.prop.trop.mu)$tTable[2,4],8))
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$tTable[2,1],8),round(summary(pgls.prop.trop.ndr)$tTable[2,4],8))
####    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.turnover',round(summary(pgls.prop.trop.turnover)$tTable[2,1],8),round(summary(pgls.prop.trop.turnover)$tTable[2,4],8))
####    
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.lambda',round(summary(pgls.prop.temp.lambda)$tTable[2,1],8),round(summary(pgls.prop.temp.lambda)$tTable[2,4],8))
####    pgls.prop.temp.muresults<-c('pgls.prop.temp.mu',round(summary(pgls.prop.temp.mu)$tTable[2,1],8),round(summary(pgls.prop.temp.mu)$tTable[2,4],8))
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$tTable[2,1],8),round(summary(pgls.prop.temp.ndr)$tTable[2,4],8))
####    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.turnover',round(summary(pgls.prop.temp.turnover)$tTable[2,1],8),round(summary(pgls.prop.temp.turnover)$tTable[2,4],8))
####    
####    if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####      
####      pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.lambda',round(summary(pgls.strict.tropical.lambda)$tTable[2,1],8),round(summary(pgls.strict.tropical.lambda)$tTable[2,4],8))
####      pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$tTable[2,1],8),round(summary(pgls.strict.tropical.mu)$tTable[2,4],8))
####      pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$tTable[2,1],8),round(summary(pgls.strict.tropical.ndr)$tTable[2,4],8))
####      pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.turnover',round(summary(pgls.strict.tropical.turnover)$tTable[2,1],8),round(summary(pgls.strict.tropical.turnover)$tTable[2,4],8))
####      
####      phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####      phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####      phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####      phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####    }
####  }
####  
####  if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####    results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.muresults,pgls.Median.Latitude.ndrresults,pgls.Median.Latitude.turnoverresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.muresults,pgls.prop.trop.ndrresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.muresults,pgls.prop.temp.ndrresults,pgls.prop.temp.turnoverresults,pgls.strict.tropical.lambdaresults,pgls.strict.tropical.muresults,pgls.strict.tropical.ndrresults,pgls.strict.tropical.turnoverresults,phylanova.strict.tropical.lambda.results,phylanova.strict.tropical.mu.results,phylanova.strict.tropical.ndr.results,phylanova.strict.tropical.turnover.results))
####  }else{
####    results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.muresults,pgls.Median.Latitude.ndrresults,pgls.Median.Latitude.turnoverresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.muresults,pgls.prop.trop.ndrresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.muresults,pgls.prop.temp.ndrresults,pgls.prop.temp.turnoverresults))
####  }
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_',strict.tropical.threshold,'.tropthreshold_pgls_table_new_GBIFsampling.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####    
####    plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='lambda')
####    points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####    
####    plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='mu')
####    points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####    
####    plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='ndr')
####    points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####    
####    plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='turnover')
####    points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  }
####  dev.off()
####}

####run_clades_RPANDA_pgls_savedphylogroups_and_save<-function(tree,minage,maxage,mincladesize,ncores,sampling,table,strict.tropical.threshold,GBIF.sampling){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  phylogroups<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'.RDS',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_phylogroup(x,table,sampling,mincladesize=mincladesize))
####  saveRDS(phylogroups.RPANDA,file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'_RPANDA.RDS',sep=''))
####}
####

run_clades_RPANDAmodelave_pgls_savedphylogroups_and_save<-function(tree,minage,maxage,mincladesize,ncores,sampling,table,strict.tropical.threshold,GBIF.sampling){
  colnames(table)[4]<-'tip'
  colnames(table)[9]<-'genus.species.sampling.fraction'
  cat('getting clades','\n')
  phylogroups<-readRDS(file=paste('./results/clade_analyses//phylogroups_',minage,'_',maxage,'_size',mincladesize,'.RDS',sep=''))
  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
  length(unique(phylogroups.species))
  cat('measuring diversification/trait evolution','\n')
  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_phylogroup_modelave(x,table,sampling,mincladesize=mincladesize))
  saveRDS(phylogroups.RPANDA,file=paste('./results/clade_analyses/phylogroups_',minage,'_',maxage,'_size',mincladesize,'_RPANDAmodelave.RDS',sep=''))
}
####run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log<-function(tree,minage,maxage,mincladesize,ncores,sampling,table,strict.tropical.threshold,GBIF.sampling,loadRPANDA){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  if(loadRPANDA==F){
####    phylogroups<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'.RDS',sep=''))
####    phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####    length(unique(phylogroups.species))
####    cat('measuring diversification/trait evolution','\n')
####    phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_phylogroup(x,table,sampling,mincladesize=mincladesize))
####    
####  }else if(loadRPANDA==T){
####    cat('loading RPANDA phylogroups object','\n')
####    phylogroups.RPANDA<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'_RPANDA.RDS',sep=''))
####  }
####  
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of lat data points per phylogroup
####  phylogroup.species.lat<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.lat[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude)])))
####  }
####  phylogroup.species.lat<-as.data.frame(phylogroup.species.lat,stringsAsFactors = F)
####  colnames(phylogroup.species.lat)<-c('phylogroup.name','phylogroup.species.with.lat.data')
####  phylogroup.species.lat$phylogroup.species.with.lat.data<-as.numeric(phylogroup.species.lat$phylogroup.species.with.lat.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.lat,all.x=TRUE)
####  
####  #get the mean lat in each phylogroup
####  #mean.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,mean)
####  #colnames(mean.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  #df.phylogroups<-merge(df.phylogroups,mean.lat.phylogroups)
####  
####  #median is very correlated with mean
####  median.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,median)
####  colnames(median.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  df.phylogroups<-merge(df.phylogroups,median.lat.phylogroups)
####  
####  #sd.Median.Latitude is quite big in phylogroups
####  sd.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,function(x) sd(abs(x)))
####  colnames(sd.lat.phylogroups)[2]<-'phylogroup.sd.Median.Latitude'
####  sd.lat.phylogroups$phylogroup.sd.Median.Latitude[is.na(sd.lat.phylogroups$phylogroup.sd.Median.Latitude)]<-0
####  df.phylogroups<-merge(df.phylogroups,sd.lat.phylogroups)
####  
####  phylogroups.counts.stricts.df<-matrix(NA,ncol=4,nrow=0)
####  colnames(phylogroups.counts.stricts.df)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup','phylogroup.name')
####  #get the counts of strict tropical
####  for(d in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.table<-df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[d],]
####    phylogroup.table<-phylogroup.table[complete.cases(phylogroup.table),]
####    strict.counts<-count(phylogroup.table, c("phylogroup.name","strict.tropical"))
####    strict.counts<-strict.counts[,c('strict.tropical','freq')]
####    #add missing 0,1,2 counts (e.g if there are not strict tropical -  1s -  add a row with strict tropical =1 and freq =0)
####    if((0%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(0,0))
####    }
####    if((1%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(1,0))
####    }
####    if((2%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(2,0))
####    }
####    phylogroup.table.strict.counts<-as.data.frame(matrix(NA,nrow=1,ncol=3),stringsAsFactors = F)
####    colnames(phylogroup.table.strict.counts)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup')
####    phylogroup.table.strict.counts$n.strict.temperate.species.phylogroup<-strict.counts[strict.counts$strict.tropical==0,'freq']
####    phylogroup.table.strict.counts$n.strict.tropical.species.phylogroup<-strict.counts[strict.counts$strict.tropical==1,'freq']
####    phylogroup.table.strict.counts$n.strict.widespread.species.phylogroup<-strict.counts[strict.counts$strict.tropical==2,'freq']
####    phylogroup.table.strict.counts$phylogroup.name<-phylogroup.table$phylogroup.name[1]
####    phylogroups.counts.stricts.df<-rbind(phylogroups.counts.stricts.df,phylogroup.table.strict.counts)
####  }
####  count.strict.tropical.phylogroups<-aggregate(strict.tropical~phylogroup.name,data=df.phylogroups,count)
####  # phylogroups.counts.stricts.df<-phylogroups.counts.stricts.df[complete.cases(phylogroups.counts.stricts.df),]
####  df.phylogroups<-merge(df.phylogroups,phylogroups.counts.stricts.df,all.x=T)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  
####  df.phylogroups$turnover.rpanda1<-df.phylogroups$lambda.rpanda1+df.phylogroups$mu.rpanda1
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #filter by GBIF sampling
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>GBIF.sampling,]
####  
####  #dataset with strict tropical and strict temperate species
####  #strict tropical phylogroups =  phylogroups with > .7 of species strict.tropical
####  #strict temperate phylogroups =  phylogroups with > .7 of species strict.tropical
####  df.phylogroups$strict.tropical.character<-NA
####  df.phylogroups[(df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-1
####  df.phylogroups[(df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-0
####  
####  
####  
####  df.phylogroups$proportion.strict.tropical<-df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  df.phylogroups$proportion.strict.temperate<-df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  
####  df.phylogroups$lambda.rpanda1<-log10(df.phylogroups$lambda.rpanda1)
####  df.phylogroups$turnover.rpanda1<-log10(df.phylogroups$turnover.rpanda1)
####  
####  df.phylogroups.strict<-df.phylogroups[!is.na(df.phylogroups$strict.tropical.character),]
####  
####  
####  ###############this is for lambda
####  pdf(paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_new_GBIFsampling_log.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Median.Latitude,xlab='phylogroup.Median.Latitude')
####  
####  if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####    wilcox.test.trop.vs.temp.lambda<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
####    t.test.trop.vs.temp.lambda<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.mu<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
####    t.test.trop.vs.temp.mu<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.ndr<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
####    t.test.trop.vs.temp.ndr<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.turnover<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    t.test.trop.vs.temp.turnover<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    
####    wilcox.t.results<-as.data.frame(rbind(c('log10.lambda',wilcox.test.trop.vs.temp.lambda$p.value,t.test.trop.vs.temp.lambda$p.value),c('mu',wilcox.test.trop.vs.temp.mu$p.value,t.test.trop.vs.temp.mu$p.value),c('ndr',wilcox.test.trop.vs.temp.ndr$p.value,t.test.trop.vs.temp.ndr$p.value),c('log10.turnover',wilcox.test.trop.vs.temp.turnover$p.value,t.test.trop.vs.temp.turnover$p.value)))
####    colnames(wilcox.t.results)<-c('parameter','pvalue.wilcox','pvalue.ttest')
####    write.table(wilcox.t.results,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_wilcox.ttest_table_new_GBIFsampling_log.txt',sep=''),sep='\t',quote=F,row.names=F)
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],main=paste('lambda strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.lambda$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.lambda.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],main=paste('mu strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.mu$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='mu.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],main=paste('ndr strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.ndr$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='ndr.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],main=paste('turnover strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.turnover$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.turnover.rpanda1')
####    
####  }
####  
####  
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  #phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$strict.tropical.character,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1)
####  #colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$phylogroup.sd.Median.Latitude,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1,df.phylogroups$proportion.strict.tropical,df.phylogroups$proportion.strict.temperate)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','phylogroup.sd.Median.Latitude','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1','proportion.strict.tropical','proportion.strict.temperate')
####  
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  row.names(phylogroups.data)<-phylogroups.data$Tip
####  
####  #create binary trait (trop = 1, temp =0) for phyloglm
####  phylogroups.data$binary.tropical<-NA
####  
####  phylogroups.data[phylogroups.data$proportion.strict.temperate>strict.tropical.threshold&phylogroups.data$proportion.strict.tropical<strict.tropical.threshold,'binary.tropical']<-0
####  phylogroups.data[phylogroups.data$proportion.strict.temperate<strict.tropical.threshold&phylogroups.data$proportion.strict.tropical>strict.tropical.threshold,'binary.tropical']<-1
####  phylogroups.data2<-phylogroups.data[complete.cases(phylogroups.data),]
####  phylogroups.tree2<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data2$Tip))
####  phylogroups.tree.object2<-comparative.data(phy = phylogroups.tree2,data=phylogroups.data2,names.col='Tip',vcv=TRUE)
####  
####  
####  
####  pgls.Median.Latitude.lambda<-try(pgls(lambda.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.lambda)=='try-error'){
####    pgls.Median.Latitude.lambda<-try(gls(lambda.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.Median.Latitude.mu<-try(pgls(mu.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.mu)=='try-error'){
####    pgls.Median.Latitude.mu<-try(gls(mu.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.Median.Latitude.ndr<-try(pgls(ndr.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.ndr)=='try-error'){
####    pgls.Median.Latitude.ndr<-try(gls(ndr.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.Median.Latitude.turnover<-try(pgls(turnover.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.turnover)=='try-error'){
####    pgls.Median.Latitude.turnover<-try(gls(turnover.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  pgls.prop.trop.lambda<-try(pgls(proportion.strict.tropical~lambda.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.lambda)=='try-error'){
####    pgls.prop.trop.lambda<-try(gls(proportion.strict.tropical~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.trop.mu<-try(pgls(proportion.strict.tropical~mu.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.mu)=='try-error'){
####    pgls.prop.trop.mu<-try(gls(proportion.strict.tropical~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.trop.ndr<-try(pgls(proportion.strict.tropical~ndr.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.ndr)=='try-error'){
####    pgls.prop.trop.ndr<-try(gls(proportion.strict.tropical~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.trop.turnover<-try(pgls(proportion.strict.tropical~turnover.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.turnover)=='try-error'){
####    pgls.prop.trop.turnover<-try(gls(proportion.strict.tropical~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  pgls.prop.temp.lambda<-try(pgls(proportion.strict.temperate~lambda.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.lambda)=='try-error'){
####    pgls.prop.temp.lambda<-try(gls(proportion.strict.temperate~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.mu<-try(pgls(proportion.strict.temperate~mu.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.mu)=='try-error'){
####    pgls.prop.temp.mu<-try(gls(proportion.strict.temperate~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.ndr<-try(pgls(proportion.strict.temperate~ndr.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.ndr)=='try-error'){
####    pgls.prop.temp.ndr<-try(gls(proportion.strict.temperate~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.turnover<-try(pgls(proportion.strict.temperate~turnover.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.turnover)=='try-error'){
####    pgls.prop.temp.turnover<-try(gls(proportion.strict.temperate~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  
####  #  phylogroups.tree.strict<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups.strict$Tip))
####  #  phylogroups.data.strict<-data.frame(df.phylogroups.strict$Tip,df.phylogroups.strict$phylogroup.Median.Latitude,df.phylogroups.strict$strict.tropical.character,df.phylogroups.strict$lambda.rpanda1,df.phylogroups.strict$mu.rpanda1,df.phylogroups.strict$ndr.rpanda1,df.phylogroups.strict$turnover.rpanda1)
####  #  colnames(phylogroups.data.strict)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####  #  phylogroups.tree.strict.object<-comparative.data(phy = phylogroups.tree.strict,data=phylogroups.data.strict,names.col='Tip',vcv=TRUE)
####  #  
####  #  pgls.strict.tropical.lambda<-try(pgls(strict.tropical.character~lambda.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.mu<-try(pgls(strict.tropical.character~mu.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.ndr<-try(pgls(strict.tropical.character~ndr.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.turnover<-try(pgls(strict.tropical.character~turnover.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  
####  #  
####  #  pgls.strict.tropical.lambda<-try(gls(strict.tropical.character~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.mu<-try(gls(strict.tropical.character~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.ndr<-try(gls(strict.tropical.character~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.turnover<-try(gls(strict.tropical.character~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  
####  #  strict.tropical.phylanova<-phylogroups.tree.strict.object$data$strict.tropical.character
####  #  names(strict.tropical.phylanova)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.lambda<-phylogroups.tree.strict.object$data$lambda.rpanda1
####  #  names(strict.tropical.phylanova.lambda)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.mu<-phylogroups.tree.strict.object$data$mu.rpanda1
####  #  names(strict.tropical.phylanova.mu)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.ndr<-phylogroups.tree.strict.object$data$ndr.rpanda1
####  #  names(strict.tropical.phylanova.ndr)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.turnover<-phylogroups.tree.strict.object$data$turnover.rpanda1
####  #  names(strict.tropical.phylanova.turnover)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  phylanova.strict.tropical.lambda<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.lambda,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.mu<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.mu,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.ndr<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.ndr,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.turnover<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.turnover,nsim = 100,posthoc = T)
####  #}
####  #summary(pgls.seedlambda)
####  #pgls.Median.Latitude.sigsq<-try(pgls(log10(lambda.rpanda1)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.Median.Latitude.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.log10.lambda',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.Median.Latitude.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.log10.lambda',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],8),round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.mu)
####    pgls.Median.Latitude.muresults<-c('pgls.Median.Latitude.mu',round(summary(pgls.Median.Latitude.mu)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.mu)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.mu)
####    pgls.Median.Latitude.muresults<-c('pgls.Median.Latitude.mu',round(summary(pgls.Median.Latitude.mu)$tTable[2,1],8),round(summary(pgls.Median.Latitude.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],8),round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.turnover)
####    pgls.Median.Latitude.turnoverresults<-c('pgls.Median.Latitude.log10.turnover',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.turnover)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.turnover)
####    pgls.Median.Latitude.turnoverresults<-c('pgls.Median.Latitude.log10.turnover',round(summary(pgls.Median.Latitude.turnover)$tTable[2,1],8),round(summary(pgls.Median.Latitude.turnover)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],8),round(summary(pgls.prop.trop.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$tTable[2,1],8),round(summary(pgls.prop.trop.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.mu)
####    pgls.prop.trop.muresults<-c('pgls.prop.trop.mu',round(summary(pgls.prop.trop.mu)$coefficients[2,1],8),round(summary(pgls.prop.trop.mu)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.mu)
####    pgls.prop.trop.muresults<-c('pgls.prop.trop.mu',round(summary(pgls.prop.trop.mu)$tTable[2,1],8),round(summary(pgls.prop.trop.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],8),round(summary(pgls.prop.trop.ndr)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$tTable[2,1],8),round(summary(pgls.prop.trop.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.turnover)
####    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.log10.turnover',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],8),round(summary(pgls.prop.trop.turnover)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.turnover)
####    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.log10.turnover',round(summary(pgls.prop.trop.turnover)$tTable[2,1],8),round(summary(pgls.prop.trop.turnover)$tTable[2,4],8))
####  }
####  
####  if(class(pgls.prop.temp.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],8),round(summary(pgls.prop.temp.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$tTable[2,1],8),round(summary(pgls.prop.temp.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.mu)
####    pgls.prop.temp.muresults<-c('pgls.prop.temp.mu',round(summary(pgls.prop.temp.mu)$coefficients[2,1],8),round(summary(pgls.prop.temp.mu)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.mu)
####    pgls.prop.temp.muresults<-c('pgls.prop.temp.mu',round(summary(pgls.prop.temp.mu)$tTable[2,1],8),round(summary(pgls.prop.temp.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],8),round(summary(pgls.prop.temp.ndr)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$tTable[2,1],8),round(summary(pgls.prop.temp.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.turnover)
####    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.log10.turnover',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],8),round(summary(pgls.prop.temp.turnover)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.turnover)
####    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.log10.turnover',round(summary(pgls.prop.temp.turnover)$tTable[2,1],8),round(summary(pgls.prop.temp.turnover)$tTable[2,4],8))
####  }
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #    
####  #   pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.log10.lambda',round(summary(pgls.strict.tropical.lambda)$coefficients[2,1],8),round(summary(pgls.strict.tropical.lambda)$coefficients[2,4],8))
####  #   pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$coefficients[2,1],8),round(summary(pgls.strict.tropical.mu)$coefficients[2,4],8))
####  #   pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$coefficients[2,1],8),round(summary(pgls.strict.tropical.ndr)$coefficients[2,4],8))
####  #   pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.log10.turnover',round(summary(pgls.strict.tropical.turnover)$coefficients[2,1],8),round(summary(pgls.strict.tropical.turnover)$coefficients[2,4],8))
####  #   
####  #   phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.log10.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####  #   phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####  #   phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####  #   phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.log10.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####  #   
####  # }
####  # 
####  
####  #   pgls.strict.tropical.lambda<-try(gls(strict.tropical.character~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.mu<-try(gls(strict.tropical.character~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.ndr<-try(gls(strict.tropical.character~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.turnover<-try(gls(strict.tropical.character~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   
####  # }
####  # 
####  # if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #   
####  #   pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.log10.lambda',round(summary(pgls.strict.tropical.lambda)$tTable[2,1],8),round(summary(pgls.strict.tropical.lambda)$tTable[2,4],8))
####  #   pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$tTable[2,1],8),round(summary(pgls.strict.tropical.mu)$tTable[2,4],8))
####  #   pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$tTable[2,1],8),round(summary(pgls.strict.tropical.ndr)$tTable[2,4],8))
####  #   pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.log10.turnover',round(summary(pgls.strict.tropical.turnover)$tTable[2,1],8),round(summary(pgls.strict.tropical.turnover)$tTable[2,4],8))
####  #   
####  #   phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.log10.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####  #   phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####  #   phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####  #   phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.log10.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####  #}
####  
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.muresults,pgls.Median.Latitude.ndrresults,pgls.Median.Latitude.turnoverresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.muresults,pgls.prop.trop.ndrresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.muresults,pgls.prop.temp.ndrresults,pgls.prop.temp.turnoverresults,pgls.strict.tropical.lambdaresults,pgls.strict.tropical.muresults,pgls.strict.tropical.ndrresults,pgls.strict.tropical.turnoverresults,phylanova.strict.tropical.lambda.results,phylanova.strict.tropical.mu.results,phylanova.strict.tropical.ndr.results,phylanova.strict.tropical.turnover.results))
####  #}else{
####    results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.muresults,pgls.Median.Latitude.ndrresults,pgls.Median.Latitude.turnoverresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.muresults,pgls.prop.trop.ndrresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.muresults,pgls.prop.temp.ndrresults,pgls.prop.temp.turnoverresults))
####    colnames(results.df)<-c('analysis','slope','pvalue')
####    if((min(table(phylogroups.tree.object2$data$binary.tropical))>1)&&(length(table(phylogroups.tree.object2$data$binary.tropical))>1)){
####      pgls.binary.trop.lambda<-phyloglm(binary.tropical~lambda.rpanda1, data = phylogroups.tree.object2$data, phy=phylogroups.tree2)
####      pgls.binary.trop.turnover<-phyloglm(binary.tropical~turnover.rpanda1, data = phylogroups.tree.object2$data, phy=phylogroups.tree2)
####      
####      pgls.binary.trop.lambda.results<-c('pgls.binary.trop.log10.lambda',round(summary(pgls.binary.trop.lambda)$coefficients[2,1],8),round(summary(pgls.binary.trop.lambda)$coefficients[2,4],8))
####      pgls.binary.trop.turnover.results<-c('pgls.binary.trop.log10.turnover',round(summary(pgls.binary.trop.turnover)$coefficients[2,1],8),round(summary(pgls.binary.trop.turnover)$coefficients[2,4],8))
####      results.binary.trop<-as.data.frame(rbind(pgls.binary.trop.lambda.results,pgls.binary.trop.turnover.results))
####      colnames(results.binary.trop)<-c('analysis','slope','pvalue')
####      results.df<-rbind(results.df,results.binary.trop)
####    }
####    
####  
####  
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_',strict.tropical.threshold,'.tropthreshold_pgls_table_new_GBIFsampling_log.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='log10.lambda')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='mu')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='ndr')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='log10.turnover')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #}
####  dev.off()
####}

####run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave<-function(tree,minage,maxage,mincladesize,ncores,sampling,table,strict.tropical.threshold,GBIF.sampling,loadRPANDA){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  if(loadRPANDA==F){
####    phylogroups<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'.RDS',sep=''))
####    phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####    length(unique(phylogroups.species))
####    cat('measuring diversification/trait evolution','\n')
####    phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_phylogroup_modelave(x,table,sampling,mincladesize=mincladesize))
####    
####  }else if(loadRPANDA==T){
####    cat('loading RPANDA phylogroups object','\n')
####    phylogroups.RPANDA<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'_RPANDAmodelave.RDS',sep=''))
####  }
####  
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of lat data points per phylogroup
####  phylogroup.species.lat<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.lat[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude)])))
####  }
####  phylogroup.species.lat<-as.data.frame(phylogroup.species.lat,stringsAsFactors = F)
####  colnames(phylogroup.species.lat)<-c('phylogroup.name','phylogroup.species.with.lat.data')
####  phylogroup.species.lat$phylogroup.species.with.lat.data<-as.numeric(phylogroup.species.lat$phylogroup.species.with.lat.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.lat,all.x=TRUE)
####  
####  #get the mean lat in each phylogroup
####  #mean.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,mean)
####  #colnames(mean.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  #df.phylogroups<-merge(df.phylogroups,mean.lat.phylogroups)
####  
####  #median is very correlated with mean
####  median.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,median)
####  colnames(median.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  df.phylogroups<-merge(df.phylogroups,median.lat.phylogroups)
####  
####  #sd.Median.Latitude is quite big in phylogroups
####  sd.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,function(x) sd(abs(x)))
####  colnames(sd.lat.phylogroups)[2]<-'phylogroup.sd.Median.Latitude'
####  sd.lat.phylogroups$phylogroup.sd.Median.Latitude[is.na(sd.lat.phylogroups$phylogroup.sd.Median.Latitude)]<-0
####  df.phylogroups<-merge(df.phylogroups,sd.lat.phylogroups)
####  
####  phylogroups.counts.stricts.df<-matrix(NA,ncol=4,nrow=0)
####  colnames(phylogroups.counts.stricts.df)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup','phylogroup.name')
####  #get the counts of strict tropical
####  for(d in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.table<-df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[d],]
####    phylogroup.table<-phylogroup.table[complete.cases(phylogroup.table),]
####    strict.counts<-count(phylogroup.table, c("phylogroup.name","strict.tropical"))
####    strict.counts<-strict.counts[,c('strict.tropical','freq')]
####    #add missing 0,1,2 counts (e.g if there are not strict tropical -  1s -  add a row with strict tropical =1 and freq =0)
####    if((0%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(0,0))
####    }
####    if((1%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(1,0))
####    }
####    if((2%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(2,0))
####    }
####    phylogroup.table.strict.counts<-as.data.frame(matrix(NA,nrow=1,ncol=3),stringsAsFactors = F)
####    colnames(phylogroup.table.strict.counts)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup')
####    phylogroup.table.strict.counts$n.strict.temperate.species.phylogroup<-strict.counts[strict.counts$strict.tropical==0,'freq']
####    phylogroup.table.strict.counts$n.strict.tropical.species.phylogroup<-strict.counts[strict.counts$strict.tropical==1,'freq']
####    phylogroup.table.strict.counts$n.strict.widespread.species.phylogroup<-strict.counts[strict.counts$strict.tropical==2,'freq']
####    phylogroup.table.strict.counts$phylogroup.name<-phylogroup.table$phylogroup.name[1]
####    phylogroups.counts.stricts.df<-rbind(phylogroups.counts.stricts.df,phylogroup.table.strict.counts)
####  }
####  count.strict.tropical.phylogroups<-aggregate(strict.tropical~phylogroup.name,data=df.phylogroups,count)
####  # phylogroups.counts.stricts.df<-phylogroups.counts.stricts.df[complete.cases(phylogroups.counts.stricts.df),]
####  df.phylogroups<-merge(df.phylogroups,phylogroups.counts.stricts.df,all.x=T)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  
####  df.phylogroups$turnover.rpanda1<-df.phylogroups$lambda.rpanda1+df.phylogroups$mu.rpanda1
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #filter by GBIF sampling
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>GBIF.sampling,]
####  
####  #dataset with strict tropical and strict temperate species
####  #strict tropical phylogroups =  phylogroups with > .7 of species strict.tropical
####  #strict temperate phylogroups =  phylogroups with > .7 of species strict.tropical
####  df.phylogroups$strict.tropical.character<-NA
####  df.phylogroups[(df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-1
####  df.phylogroups[(df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-0
####  
####  
####  
####  df.phylogroups$proportion.strict.tropical<-df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  df.phylogroups$proportion.strict.temperate<-df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  
####  df.phylogroups$lambda.rpanda1<-log10(df.phylogroups$lambda.rpanda1)
####  df.phylogroups$turnover.rpanda1<-log10(df.phylogroups$turnover.rpanda1)
####  
####  df.phylogroups.strict<-df.phylogroups[!is.na(df.phylogroups$strict.tropical.character),]
####  
####  
####  ###############this is for lambda
####  pdf(paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_new_GBIFsampling_log_modelave.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Median.Latitude,xlab='phylogroup.Median.Latitude')
####  
####  if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####    wilcox.test.trop.vs.temp.lambda<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
####    t.test.trop.vs.temp.lambda<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.mu<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
####    t.test.trop.vs.temp.mu<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.ndr<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
####    t.test.trop.vs.temp.ndr<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.turnover<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    t.test.trop.vs.temp.turnover<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    
####    wilcox.t.results<-as.data.frame(rbind(c('log10.lambda',wilcox.test.trop.vs.temp.lambda$p.value,t.test.trop.vs.temp.lambda$p.value),c('mu',wilcox.test.trop.vs.temp.mu$p.value,t.test.trop.vs.temp.mu$p.value),c('ndr',wilcox.test.trop.vs.temp.ndr$p.value,t.test.trop.vs.temp.ndr$p.value),c('log10.turnover',wilcox.test.trop.vs.temp.turnover$p.value,t.test.trop.vs.temp.turnover$p.value)))
####    colnames(wilcox.t.results)<-c('parameter','pvalue.wilcox','pvalue.ttest')
####    write.table(wilcox.t.results,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_wilcox.ttest_table_new_GBIFsampling_log_modelave.txt',sep=''),sep='\t',quote=F,row.names=F)
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],main=paste('lambda strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.lambda$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.lambda.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],main=paste('mu strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.mu$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='mu.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],main=paste('ndr strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.ndr$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='ndr.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],main=paste('turnover strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.turnover$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.turnover.rpanda1')
####    
####  }
####  
####  
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  #phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$strict.tropical.character,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1)
####  #colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$phylogroup.sd.Median.Latitude,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1,df.phylogroups$proportion.strict.tropical,df.phylogroups$proportion.strict.temperate)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','phylogroup.sd.Median.Latitude','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1','proportion.strict.tropical','proportion.strict.temperate')
####  
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  row.names(phylogroups.data)<-phylogroups.data$Tip
####  
####  #create binary trait (trop = 1, temp =0) for phyloglm
####  phylogroups.data$binary.tropical<-NA
####  
####  phylogroups.data[phylogroups.data$proportion.strict.temperate>strict.tropical.threshold&phylogroups.data$proportion.strict.tropical<strict.tropical.threshold,'binary.tropical']<-0
####  phylogroups.data[phylogroups.data$proportion.strict.temperate<strict.tropical.threshold&phylogroups.data$proportion.strict.tropical>strict.tropical.threshold,'binary.tropical']<-1
####  phylogroups.data2<-phylogroups.data[complete.cases(phylogroups.data),]
####  phylogroups.tree2<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data2$Tip))
####  phylogroups.tree.object2<-comparative.data(phy = phylogroups.tree2,data=phylogroups.data2,names.col='Tip',vcv=TRUE)
####  
####  
####  
####  pgls.Median.Latitude.lambda<-try(pgls(lambda.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.lambda)=='try-error'){
####    pgls.Median.Latitude.lambda<-try(gls(lambda.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.Median.Latitude.mu<-try(pgls(mu.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.mu)=='try-error'){
####    pgls.Median.Latitude.mu<-try(gls(mu.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.Median.Latitude.ndr<-try(pgls(ndr.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.ndr)=='try-error'){
####    pgls.Median.Latitude.ndr<-try(gls(ndr.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.Median.Latitude.turnover<-try(pgls(turnover.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.turnover)=='try-error'){
####    pgls.Median.Latitude.turnover<-try(gls(turnover.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  pgls.prop.trop.lambda<-try(pgls(proportion.strict.tropical~lambda.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.lambda)=='try-error'){
####    pgls.prop.trop.lambda<-try(gls(proportion.strict.tropical~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    if (class(pgls.prop.trop.lambda)=='try-error'){
####      lam<-seq(0,1,0.1)
####      fit.prop.trop.lambda<-list()
####      form.prop.trop.lambda<-proportion.strict.tropical~lambda.rpanda1
####      phylogroups.data<-phylogroups.data[complete.cases(phylogroups.data[,c('lambda.rpanda1')]),]
####      for (i in seq_along(lam)){
####        cor <- corPagel(lam[i], phy = drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data$Tip)), fixed = TRUE)
####        fit.prop.trop.lambda[[i]] <- gls(form.prop.trop.lambda, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML")
####      }
####      pgls.prop.trop.lambda<-fit.prop.trop.lambda[[which.min(sapply(fit.prop.trop.lambda, logLik))]]
####    }
####  }
####  pgls.prop.trop.mu<-try(pgls(proportion.strict.tropical~mu.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.mu)=='try-error'){
####    pgls.prop.trop.mu<-try(gls(proportion.strict.tropical~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.trop.ndr<-try(pgls(proportion.strict.tropical~ndr.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.ndr)=='try-error'){
####    pgls.prop.trop.ndr<-try(gls(proportion.strict.tropical~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.trop.turnover<-try(pgls(proportion.strict.tropical~turnover.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.turnover)=='try-error'){
####    pgls.prop.trop.turnover<-try(gls(proportion.strict.tropical~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    if (class(pgls.prop.trop.turnover)=='try-error'){
####      lam<-seq(0,1,0.1)
####      fit.prop.trop.turnover<-list()
####      form.prop.trop.turnover<-proportion.strict.tropical~turnover.rpanda1
####      phylogroups.data<-phylogroups.data[complete.cases(phylogroups.data[,c('turnover.rpanda1')]),]
####      for (i in seq_along(lam)){
####        cor <- corPagel(lam[i], phy = drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data$Tip)), fixed = TRUE)
####        fit.prop.trop.turnover[[i]] <- gls(form.prop.trop.turnover, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML")
####      }
####      pgls.prop.trop.turnover<-fit.prop.trop.turnover[[which.min(sapply(fit.prop.trop.turnover, logLik))]]
####    }
####  }
####  
####  pgls.prop.temp.lambda<-try(pgls(proportion.strict.temperate~lambda.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.lambda)=='try-error'){
####    pgls.prop.temp.lambda<-try(gls(proportion.strict.temperate~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.mu<-try(pgls(proportion.strict.temperate~mu.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.mu)=='try-error'){
####    pgls.prop.temp.mu<-try(gls(proportion.strict.temperate~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.ndr<-try(pgls(proportion.strict.temperate~ndr.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.ndr)=='try-error'){
####    pgls.prop.temp.ndr<-try(gls(proportion.strict.temperate~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    if(class(pgls.prop.temp.ndr)=='try-error'){
####      pgls.prop.temp.ndr<-try(gls(proportion.strict.tropical~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####      if (class(pgls.prop.temp.ndr)=='try-error'){
####        lam<-seq(0,1,0.1)
####        fit.prop.temp.ndr<-list()
####        form.prop.temp.ndr<-proportion.strict.tropical~ndr.rpanda1
####        phylogroups.data<-phylogroups.data[complete.cases(phylogroups.data[,c('ndr.rpanda1')]),]
####        for (i in seq_along(lam)){
####          cor <- corPagel(lam[i], phy = drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data$Tip)), fixed = TRUE)
####          fit.prop.temp.ndr[[i]] <- gls(form.prop.temp.ndr, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML")
####        }
####        pgls.prop.temp.ndr<-fit.prop.temp.ndr[[which.min(sapply(fit.prop.temp.ndr, logLik))]]
####      }
####    }
####  }
####  pgls.prop.temp.turnover<-try(pgls(proportion.strict.temperate~turnover.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.turnover)=='try-error'){
####    pgls.prop.temp.turnover<-try(gls(proportion.strict.temperate~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  
####  #  phylogroups.tree.strict<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups.strict$Tip))
####  #  phylogroups.data.strict<-data.frame(df.phylogroups.strict$Tip,df.phylogroups.strict$phylogroup.Median.Latitude,df.phylogroups.strict$strict.tropical.character,df.phylogroups.strict$lambda.rpanda1,df.phylogroups.strict$mu.rpanda1,df.phylogroups.strict$ndr.rpanda1,df.phylogroups.strict$turnover.rpanda1)
####  #  colnames(phylogroups.data.strict)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####  #  phylogroups.tree.strict.object<-comparative.data(phy = phylogroups.tree.strict,data=phylogroups.data.strict,names.col='Tip',vcv=TRUE)
####  #  
####  #  pgls.strict.tropical.lambda<-try(pgls(strict.tropical.character~lambda.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.mu<-try(pgls(strict.tropical.character~mu.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.ndr<-try(pgls(strict.tropical.character~ndr.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.turnover<-try(pgls(strict.tropical.character~turnover.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  
####  #  
####  #  pgls.strict.tropical.lambda<-try(gls(strict.tropical.character~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.mu<-try(gls(strict.tropical.character~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.ndr<-try(gls(strict.tropical.character~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.turnover<-try(gls(strict.tropical.character~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  
####  #  strict.tropical.phylanova<-phylogroups.tree.strict.object$data$strict.tropical.character
####  #  names(strict.tropical.phylanova)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.lambda<-phylogroups.tree.strict.object$data$lambda.rpanda1
####  #  names(strict.tropical.phylanova.lambda)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.mu<-phylogroups.tree.strict.object$data$mu.rpanda1
####  #  names(strict.tropical.phylanova.mu)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.ndr<-phylogroups.tree.strict.object$data$ndr.rpanda1
####  #  names(strict.tropical.phylanova.ndr)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.turnover<-phylogroups.tree.strict.object$data$turnover.rpanda1
####  #  names(strict.tropical.phylanova.turnover)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  phylanova.strict.tropical.lambda<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.lambda,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.mu<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.mu,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.ndr<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.ndr,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.turnover<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.turnover,nsim = 100,posthoc = T)
####  #}
####  #summary(pgls.seedlambda)
####  #pgls.Median.Latitude.sigsq<-try(pgls(log10(lambda.rpanda1)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.Median.Latitude.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.log10.lambda',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.Median.Latitude.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.log10.lambda',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],8),round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.mu)
####    pgls.Median.Latitude.muresults<-c('pgls.Median.Latitude.mu',round(summary(pgls.Median.Latitude.mu)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.mu)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.mu)
####    pgls.Median.Latitude.muresults<-c('pgls.Median.Latitude.mu',round(summary(pgls.Median.Latitude.mu)$tTable[2,1],8),round(summary(pgls.Median.Latitude.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],8),round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.turnover)
####    pgls.Median.Latitude.turnoverresults<-c('pgls.Median.Latitude.log10.turnover',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.turnover)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.turnover)
####    pgls.Median.Latitude.turnoverresults<-c('pgls.Median.Latitude.log10.turnover',round(summary(pgls.Median.Latitude.turnover)$tTable[2,1],8),round(summary(pgls.Median.Latitude.turnover)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],8),round(summary(pgls.prop.trop.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$tTable[2,1],8),round(summary(pgls.prop.trop.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.mu)
####    pgls.prop.trop.muresults<-c('pgls.prop.trop.mu',round(summary(pgls.prop.trop.mu)$coefficients[2,1],8),round(summary(pgls.prop.trop.mu)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.mu)
####    pgls.prop.trop.muresults<-c('pgls.prop.trop.mu',round(summary(pgls.prop.trop.mu)$tTable[2,1],8),round(summary(pgls.prop.trop.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],8),round(summary(pgls.prop.trop.ndr)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$tTable[2,1],8),round(summary(pgls.prop.trop.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.turnover)
####    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.log10.turnover',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],8),round(summary(pgls.prop.trop.turnover)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.turnover)
####    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.log10.turnover',round(summary(pgls.prop.trop.turnover)$tTable[2,1],8),round(summary(pgls.prop.trop.turnover)$tTable[2,4],8))
####  }
####  
####  if(class(pgls.prop.temp.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],8),round(summary(pgls.prop.temp.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$tTable[2,1],8),round(summary(pgls.prop.temp.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.mu)
####    pgls.prop.temp.muresults<-c('pgls.prop.temp.mu',round(summary(pgls.prop.temp.mu)$coefficients[2,1],8),round(summary(pgls.prop.temp.mu)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.mu)
####    pgls.prop.temp.muresults<-c('pgls.prop.temp.mu',round(summary(pgls.prop.temp.mu)$tTable[2,1],8),round(summary(pgls.prop.temp.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],8),round(summary(pgls.prop.temp.ndr)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$tTable[2,1],8),round(summary(pgls.prop.temp.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.turnover)
####    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.log10.turnover',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],8),round(summary(pgls.prop.temp.turnover)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.turnover)
####    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.log10.turnover',round(summary(pgls.prop.temp.turnover)$tTable[2,1],8),round(summary(pgls.prop.temp.turnover)$tTable[2,4],8))
####  }
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #    
####  #   pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.log10.lambda',round(summary(pgls.strict.tropical.lambda)$coefficients[2,1],8),round(summary(pgls.strict.tropical.lambda)$coefficients[2,4],8))
####  #   pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$coefficients[2,1],8),round(summary(pgls.strict.tropical.mu)$coefficients[2,4],8))
####  #   pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$coefficients[2,1],8),round(summary(pgls.strict.tropical.ndr)$coefficients[2,4],8))
####  #   pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.log10.turnover',round(summary(pgls.strict.tropical.turnover)$coefficients[2,1],8),round(summary(pgls.strict.tropical.turnover)$coefficients[2,4],8))
####  #   
####  #   phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.log10.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####  #   phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####  #   phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####  #   phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.log10.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####  #   
####  # }
####  # 
####  
####  #   pgls.strict.tropical.lambda<-try(gls(strict.tropical.character~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.mu<-try(gls(strict.tropical.character~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.ndr<-try(gls(strict.tropical.character~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.turnover<-try(gls(strict.tropical.character~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   
####  # }
####  # 
####  # if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #   
####  #   pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.log10.lambda',round(summary(pgls.strict.tropical.lambda)$tTable[2,1],8),round(summary(pgls.strict.tropical.lambda)$tTable[2,4],8))
####  #   pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$tTable[2,1],8),round(summary(pgls.strict.tropical.mu)$tTable[2,4],8))
####  #   pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$tTable[2,1],8),round(summary(pgls.strict.tropical.ndr)$tTable[2,4],8))
####  #   pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.log10.turnover',round(summary(pgls.strict.tropical.turnover)$tTable[2,1],8),round(summary(pgls.strict.tropical.turnover)$tTable[2,4],8))
####  #   
####  #   phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.log10.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####  #   phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####  #   phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####  #   phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.log10.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####  #}
####  
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.muresults,pgls.Median.Latitude.ndrresults,pgls.Median.Latitude.turnoverresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.muresults,pgls.prop.trop.ndrresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.muresults,pgls.prop.temp.ndrresults,pgls.prop.temp.turnoverresults,pgls.strict.tropical.lambdaresults,pgls.strict.tropical.muresults,pgls.strict.tropical.ndrresults,pgls.strict.tropical.turnoverresults,phylanova.strict.tropical.lambda.results,phylanova.strict.tropical.mu.results,phylanova.strict.tropical.ndr.results,phylanova.strict.tropical.turnover.results))
####  #}else{
####  results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.muresults,pgls.Median.Latitude.ndrresults,pgls.Median.Latitude.turnoverresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.muresults,pgls.prop.trop.ndrresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.muresults,pgls.prop.temp.ndrresults,pgls.prop.temp.turnoverresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  if((min(table(phylogroups.tree.object2$data$binary.tropical))>1)&&(length(table(phylogroups.tree.object2$data$binary.tropical))>1)){
####    pgls.binary.trop.lambda<-phyloglm(binary.tropical~lambda.rpanda1, data = phylogroups.tree.object2$data, phy=phylogroups.tree2)
####    pgls.binary.trop.turnover<-phyloglm(binary.tropical~turnover.rpanda1, data = phylogroups.tree.object2$data, phy=phylogroups.tree2)
####    
####    pgls.binary.trop.lambda.results<-c('pgls.binary.trop.log10.lambda',round(summary(pgls.binary.trop.lambda)$coefficients[2,1],8),round(summary(pgls.binary.trop.lambda)$coefficients[2,4],8))
####    pgls.binary.trop.turnover.results<-c('pgls.binary.trop.log10.turnover',round(summary(pgls.binary.trop.turnover)$coefficients[2,1],8),round(summary(pgls.binary.trop.turnover)$coefficients[2,4],8))
####    results.binary.trop<-as.data.frame(rbind(pgls.binary.trop.lambda.results,pgls.binary.trop.turnover.results))
####    colnames(results.binary.trop)<-c('analysis','slope','pvalue')
####    results.df<-rbind(results.df,results.binary.trop)
####  }
####  
####  
####  
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_',strict.tropical.threshold,'.tropthreshold_pgls_table_new_GBIFsampling_log_modelave.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='log10.lambda')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='mu')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='ndr')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='log10.turnover')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #}
####  dev.off()
####}
####
run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave_weightlatdata<-function(tree,minage,maxage,mincladesize,ncores,sampling,table,strict.tropical.threshold,GBIF.sampling,loadRPANDA){
  colnames(table)[4]<-'tip'
  colnames(table)[9]<-'genus.species.sampling.fraction'
  cat('getting clades','\n')
  if(loadRPANDA==F){
    phylogroups<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'.RDS',sep=''))
    phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
    length(unique(phylogroups.species))
    cat('measuring diversification/trait evolution','\n')
    phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_phylogroup_modelave(x,table,sampling,mincladesize=mincladesize))
    
  }else if(loadRPANDA==T){
    cat('loading RPANDA phylogroups object','\n')
    phylogroups.RPANDA<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'_RPANDAmodelave.RDS',sep=''))
  }
  
  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
  df.phylogroups<-unique(df.phylogroups)
  #remove rows with no lambda.rpanda value
  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
  #count the number of lat data points per phylogroup
  phylogroup.species.lat<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
    phylogroup.species.lat[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude)])))
  }
  phylogroup.species.lat<-as.data.frame(phylogroup.species.lat,stringsAsFactors = F)
  colnames(phylogroup.species.lat)<-c('phylogroup.name','phylogroup.species.with.lat.data')
  phylogroup.species.lat$phylogroup.species.with.lat.data<-as.numeric(phylogroup.species.lat$phylogroup.species.with.lat.data)
  df.phylogroups<-merge(df.phylogroups,phylogroup.species.lat,all.x=TRUE)
  
  #get the mean lat in each phylogroup
  #mean.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,mean)
  #colnames(mean.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
  #df.phylogroups<-merge(df.phylogroups,mean.lat.phylogroups)
  
  #median is very correlated with mean
  median.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,median)
  colnames(median.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
  df.phylogroups<-merge(df.phylogroups,median.lat.phylogroups)
  
  #sd.Median.Latitude is quite big in phylogroups
  sd.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,function(x) sd(abs(x)))
  colnames(sd.lat.phylogroups)[2]<-'phylogroup.sd.Median.Latitude'
  sd.lat.phylogroups$phylogroup.sd.Median.Latitude[is.na(sd.lat.phylogroups$phylogroup.sd.Median.Latitude)]<-0
  df.phylogroups<-merge(df.phylogroups,sd.lat.phylogroups)
  
  phylogroups.counts.stricts.df<-matrix(NA,ncol=4,nrow=0)
  colnames(phylogroups.counts.stricts.df)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup','phylogroup.name')
  #get the counts of strict tropical
  for(d in 1:length(unique(df.phylogroups$phylogroup.name))){
    phylogroup.table<-df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[d],]
    phylogroup.table<-phylogroup.table[complete.cases(phylogroup.table),]
    strict.counts<-count(phylogroup.table, c("phylogroup.name","strict.tropical"))
    strict.counts<-strict.counts[,c('strict.tropical','freq')]
    #add missing 0,1,2 counts (e.g if there are not strict tropical -  1s -  add a row with strict tropical =1 and freq =0)
    if((0%in%strict.counts$strict.tropical)==FALSE){
      strict.counts<-rbind(strict.counts,c(0,0))
    }
    if((1%in%strict.counts$strict.tropical)==FALSE){
      strict.counts<-rbind(strict.counts,c(1,0))
    }
    if((2%in%strict.counts$strict.tropical)==FALSE){
      strict.counts<-rbind(strict.counts,c(2,0))
    }
    phylogroup.table.strict.counts<-as.data.frame(matrix(NA,nrow=1,ncol=3),stringsAsFactors = F)
    colnames(phylogroup.table.strict.counts)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup')
    phylogroup.table.strict.counts$n.strict.temperate.species.phylogroup<-strict.counts[strict.counts$strict.tropical==0,'freq']
    phylogroup.table.strict.counts$n.strict.tropical.species.phylogroup<-strict.counts[strict.counts$strict.tropical==1,'freq']
    phylogroup.table.strict.counts$n.strict.widespread.species.phylogroup<-strict.counts[strict.counts$strict.tropical==2,'freq']
    phylogroup.table.strict.counts$phylogroup.name<-phylogroup.table$phylogroup.name[1]
    phylogroups.counts.stricts.df<-rbind(phylogroups.counts.stricts.df,phylogroup.table.strict.counts)
  }
  count.strict.tropical.phylogroups<-aggregate(strict.tropical~phylogroup.name,data=df.phylogroups,count)
  # phylogroups.counts.stricts.df<-phylogroups.counts.stricts.df[complete.cases(phylogroups.counts.stricts.df),]
  df.phylogroups<-merge(df.phylogroups,phylogroups.counts.stricts.df,all.x=T)
  #getting vector of representatives in each phylogroup
  #then keep just one tip per phylogroup
  phylogroup.clades.species<-vector('character')
  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
  }
  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
  colnames(names.phylogroups.tree)<-'Tip'
  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
  
  df.phylogroups$turnover.rpanda1<-df.phylogroups$lambda.rpanda1+df.phylogroups$mu.rpanda1
  #another filter to check that size of clade is correct
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize,]
  #add another filter to select well sampled clades
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
  #filter by GBIF sampling
  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>GBIF.sampling,]
  
  #dataset with strict tropical and strict temperate species
  #strict tropical phylogroups =  phylogroups with > .7 of species strict.tropical
  #strict temperate phylogroups =  phylogroups with > .7 of species strict.tropical
  df.phylogroups$strict.tropical.character<-NA
  df.phylogroups[(df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-1
  df.phylogroups[(df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-0
  
  
  
  df.phylogroups$proportion.strict.tropical<-df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
  df.phylogroups$proportion.strict.temperate<-df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
  
  df.phylogroups$lambda.rpanda1<-log10(df.phylogroups$lambda.rpanda1)
  df.phylogroups$turnover.rpanda1<-log10(df.phylogroups$turnover.rpanda1)
  
  df.phylogroups.strict<-df.phylogroups[!is.na(df.phylogroups$strict.tropical.character),]
  
  
  ###############this is for lambda
  pdf(paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_new_GBIFsampling_log_modelave_weightedlat.pdf',sep=''),paper='a4')
  par(mfrow=c(2,2))
  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
  hist(df.phylogroups$phylogroup.Median.Latitude,xlab='phylogroup.Median.Latitude')
  
  if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
    wilcox.test.trop.vs.temp.lambda<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
    t.test.trop.vs.temp.lambda<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
    
    wilcox.test.trop.vs.temp.mu<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
    t.test.trop.vs.temp.mu<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
    
    wilcox.test.trop.vs.temp.ndr<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
    t.test.trop.vs.temp.ndr<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
    
    wilcox.test.trop.vs.temp.turnover<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
    t.test.trop.vs.temp.turnover<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
    
    wilcox.t.results<-as.data.frame(rbind(c('log10.lambda',wilcox.test.trop.vs.temp.lambda$p.value,t.test.trop.vs.temp.lambda$p.value),c('mu',wilcox.test.trop.vs.temp.mu$p.value,t.test.trop.vs.temp.mu$p.value),c('ndr',wilcox.test.trop.vs.temp.ndr$p.value,t.test.trop.vs.temp.ndr$p.value),c('log10.turnover',wilcox.test.trop.vs.temp.turnover$p.value,t.test.trop.vs.temp.turnover$p.value)))
    colnames(wilcox.t.results)<-c('parameter','pvalue.wilcox','pvalue.ttest')
    write.table(wilcox.t.results,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_wilcox.ttest_table_new_GBIFsampling_log_modelave_weightedlat.txt',sep=''),sep='\t',quote=F,row.names=F)
    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],main=paste('lambda strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.lambda$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.lambda.rpanda1')
    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],main=paste('mu strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.mu$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='mu.rpanda1')
    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],main=paste('ndr strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.ndr$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='ndr.rpanda1')
    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],main=paste('turnover strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.turnover$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.turnover.rpanda1')
    
  }
  
  
  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
  #phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$strict.tropical.character,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1)
  #colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$phylogroup.sd.Median.Latitude,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1,df.phylogroups$proportion.strict.tropical,df.phylogroups$proportion.strict.temperate,df.phylogroups$phylogroup.species.with.lat.data/df.phylogroups$phylogroup.size)
  colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','phylogroup.sd.Median.Latitude','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1','proportion.strict.tropical','proportion.strict.temperate','GBIF.sampling')
  
  phylogroups.tree$node.label<-NULL
  #this is for caper
  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
  #for mean family seed weight
  cat('running pgls','\n')
  row.names(phylogroups.data)<-phylogroups.data$Tip
  
  #create binary trait (trop = 1, temp =0) for phyloglm
  phylogroups.data$binary.tropical<-NA
  
  phylogroups.data[phylogroups.data$proportion.strict.temperate>strict.tropical.threshold&phylogroups.data$proportion.strict.tropical<strict.tropical.threshold,'binary.tropical']<-0
  phylogroups.data[phylogroups.data$proportion.strict.temperate<strict.tropical.threshold&phylogroups.data$proportion.strict.tropical>strict.tropical.threshold,'binary.tropical']<-1
  phylogroups.data2<-phylogroups.data[complete.cases(phylogroups.data),]
  phylogroups.tree2<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data2$Tip))
  phylogroups.tree.object2<-comparative.data(phy = phylogroups.tree2,data=phylogroups.data2,names.col='Tip',vcv=TRUE)
  
  
  
  #pgls.Median.Latitude.lambda<-try(pgls(lambda.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
  #if(class(pgls.Median.Latitude.lambda)=='try-error'){
  #  pgls.Median.Latitude.lambda<-try(gls(lambda.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
  #}
  #pgls.Median.Latitude.mu<-try(pgls(mu.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
  #if(class(pgls.Median.Latitude.mu)=='try-error'){
  #  pgls.Median.Latitude.mu<-try(gls(mu.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
  #}
  #pgls.Median.Latitude.ndr<-try(pgls(ndr.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
  #if(class(pgls.Median.Latitude.ndr)=='try-error'){
  #  pgls.Median.Latitude.ndr<-try(gls(ndr.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
  #}
  #pgls.Median.Latitude.turnover<-try(pgls(turnover.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
  #if(class(pgls.Median.Latitude.turnover)=='try-error'){
  #  pgls.Median.Latitude.turnover<-try(gls(turnover.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
  #}
  
  pgls.prop.trop.lambda<-try(gls(proportion.strict.tropical~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML",weights =~GBIF.sampling))
  if(class(pgls.prop.trop.lambda)=='try-error'){
    
   
      lam<-seq(0,1,0.1)
      fit.prop.trop.lambda<-list()
      form.prop.trop.lambda<-proportion.strict.tropical~lambda.rpanda1
      phylogroups.data<-phylogroups.data[complete.cases(phylogroups.data[,c('lambda.rpanda1')]),]
      for (i in seq_along(lam)){
        cor <- corPagel(lam[i], phy = drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data$Tip)), fixed = TRUE)
        fit.prop.trop.lambda[[i]] <- gls(form.prop.trop.lambda, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML",weights =~GBIF.sampling)
      }
      pgls.prop.trop.lambda<-fit.prop.trop.lambda[[which.min(sapply(fit.prop.trop.lambda, logLik))]]
   
  }
  pgls.prop.trop.turnover<-try(gls(proportion.strict.tropical~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML",weights =~GBIF.sampling))
  if(class(pgls.prop.trop.turnover)=='try-error'){
    
      lam<-seq(0,1,0.1)
      fit.prop.trop.turnover<-list()
      form.prop.trop.turnover<-proportion.strict.tropical~turnover.rpanda1
      phylogroups.data<-phylogroups.data[complete.cases(phylogroups.data[,c('turnover.rpanda1')]),]
      for (i in seq_along(lam)){
        cor <- corPagel(lam[i], phy = drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data$Tip)), fixed = TRUE)
        fit.prop.trop.turnover[[i]] <- gls(form.prop.trop.turnover, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML",weights =~GBIF.sampling)
      }
      pgls.prop.trop.turnover<-fit.prop.trop.turnover[[which.min(sapply(fit.prop.trop.turnover, logLik))]]
    
  }
  
  pgls.prop.temp.lambda<-try(gls(proportion.strict.temperate~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML",weights =~GBIF.sampling))
  if(class(pgls.prop.temp.lambda)=='try-error'){
    
    
    lam<-seq(0,1,0.1)
    fit.prop.temp.lambda<-list()
    form.prop.temp.lambda<-proportion.strict.temperate~lambda.rpanda1
    phylogroups.data<-phylogroups.data[complete.cases(phylogroups.data[,c('lambda.rpanda1')]),]
    for (i in seq_along(lam)){
      cor <- corPagel(lam[i], phy = drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data$Tip)), fixed = TRUE)
      fit.prop.temp.lambda[[i]] <- gls(form.prop.temp.lambda, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML",weights =~GBIF.sampling)
    }
    pgls.prop.temp.lambda<-fit.prop.temp.lambda[[which.min(sapply(fit.prop.temp.lambda, logLik))]]
    
  }
 pgls.prop.temp.turnover<-try(gls(proportion.strict.temperate~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML",weights =~GBIF.sampling))
 if(class(pgls.prop.temp.turnover)=='try-error'){
   
   lam<-seq(0,1,0.1)
   fit.prop.temp.turnover<-list()
   form.prop.temp.turnover<-proportion.strict.temperate~turnover.rpanda1
   phylogroups.data<-phylogroups.data[complete.cases(phylogroups.data[,c('turnover.rpanda1')]),]
   for (i in seq_along(lam)){
     cor <- corPagel(lam[i], phy = drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data$Tip)), fixed = TRUE)
     fit.prop.temp.turnover[[i]] <- gls(form.prop.temp.turnover, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML",weights =~GBIF.sampling)
   }
   pgls.prop.temp.turnover<-fit.prop.temp.turnover[[which.min(sapply(fit.prop.temp.turnover, logLik))]]
   
 }
 
  if(class(pgls.prop.trop.lambda)=='pgls'){
    plot(phylogroups.tree.object$data$lambda.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
    abline(pgls.prop.trop.lambda)
    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],8),round(summary(pgls.prop.trop.lambda)$coefficients[2,4],8))
  }else if (class(pgls.prop.trop.lambda)=='gls'){
    plot(phylogroups.data$lambda.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
    abline(pgls.prop.trop.lambda)
    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$tTable[2,1],8),round(summary(pgls.prop.trop.lambda)$tTable[2,4],8))
  }
  if(class(pgls.prop.trop.turnover)=='pgls'){
    plot(phylogroups.tree.object$data$turnover.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
    abline(pgls.prop.trop.turnover)
    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.log10.turnover',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],8),round(summary(pgls.prop.trop.turnover)$coefficients[2,4],8))
  }else if (class(pgls.prop.trop.turnover)=='gls'){
    plot(phylogroups.data$turnover.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
    abline(pgls.prop.trop.turnover)
    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.log10.turnover',round(summary(pgls.prop.trop.turnover)$tTable[2,1],8),round(summary(pgls.prop.trop.turnover)$tTable[2,4],8))
  }
  
  if(class(pgls.prop.temp.lambda)=='pgls'){
    plot(phylogroups.tree.object$data$lambda.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
    abline(pgls.prop.temp.lambda)
    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],8),round(summary(pgls.prop.temp.lambda)$coefficients[2,4],8))
  }else if (class(pgls.prop.temp.lambda)=='gls'){
    plot(phylogroups.data$lambda.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
    abline(pgls.prop.temp.lambda)
    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$tTable[2,1],8),round(summary(pgls.prop.temp.lambda)$tTable[2,4],8))
  }
  if(class(pgls.prop.temp.turnover)=='pgls'){
    plot(phylogroups.tree.object$data$turnover.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
    abline(pgls.prop.temp.turnover)
    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.log10.turnover',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],8),round(summary(pgls.prop.temp.turnover)$coefficients[2,4],8))
  }else if (class(pgls.prop.temp.turnover)=='gls'){
    plot(phylogroups.data$turnover.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
    abline(pgls.prop.temp.turnover)
    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.log10.turnover',round(summary(pgls.prop.temp.turnover)$tTable[2,1],8),round(summary(pgls.prop.temp.turnover)$tTable[2,4],8))
  }
  write.table(phylogroups.data$GBIF.sampling,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_',strict.tropical.threshold,'.tropthreshold_pgls_table_new_GBIFsampling_log_modelave_weightedlat_phylogroupGBIFsampling.txt',sep=''),quote=F,row.names=F,sep='\t')
  results.df<-as.data.frame(rbind(pgls.prop.trop.lambdaresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.turnoverresults))
  colnames(results.df)<-c('analysis','slope','pvalue')
  if((min(table(phylogroups.tree.object2$data$binary.tropical))>1)&&(length(table(phylogroups.tree.object2$data$binary.tropical))>1)){
    pgls.binary.trop.lambda<-phyloglm(binary.tropical~lambda.rpanda1, data = phylogroups.tree.object2$data, phy=phylogroups.tree2)
    pgls.binary.trop.turnover<-phyloglm(binary.tropical~turnover.rpanda1, data = phylogroups.tree.object2$data, phy=phylogroups.tree2)
    
    pgls.binary.trop.lambda.results<-c('pgls.binary.trop.log10.lambda',round(summary(pgls.binary.trop.lambda)$coefficients[2,1],8),round(summary(pgls.binary.trop.lambda)$coefficients[2,4],8))
    pgls.binary.trop.turnover.results<-c('pgls.binary.trop.log10.turnover',round(summary(pgls.binary.trop.turnover)$coefficients[2,1],8),round(summary(pgls.binary.trop.turnover)$coefficients[2,4],8))
    results.binary.trop<-as.data.frame(rbind(pgls.binary.trop.lambda.results,pgls.binary.trop.turnover.results))
    colnames(results.binary.trop)<-c('analysis','slope','pvalue')
    results.df<-rbind(results.df,results.binary.trop)
  }
  
  
  
  rownames(results.df)<-NULL
  write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_',strict.tropical.threshold,'.tropthreshold_pgls_table_new_GBIFsampling_log_modelave_weightedlat.txt',sep=''),quote=F,row.names=F,sep='\t')
  
  dev.off()
}

####run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave_size10<-function(tree,minage,maxage,mincladesize,ncores,sampling,table,strict.tropical.threshold,GBIF.sampling,loadRPANDA){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  if(loadRPANDA==F){
####    phylogroups<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'.RDS',sep=''))
####    phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####    length(unique(phylogroups.species))
####    cat('measuring diversification/trait evolution','\n')
####    phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_phylogroup_modelave(x,table,sampling,mincladesize=mincladesize))
####    
####  }else if(loadRPANDA==T){
####    cat('loading RPANDA phylogroups object','\n')
####    phylogroups.RPANDA<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'_RPANDAmodelave.RDS',sep=''))
####  }
####  
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of lat data points per phylogroup
####  phylogroup.species.lat<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.lat[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude)])))
####  }
####  phylogroup.species.lat<-as.data.frame(phylogroup.species.lat,stringsAsFactors = F)
####  colnames(phylogroup.species.lat)<-c('phylogroup.name','phylogroup.species.with.lat.data')
####  phylogroup.species.lat$phylogroup.species.with.lat.data<-as.numeric(phylogroup.species.lat$phylogroup.species.with.lat.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.lat,all.x=TRUE)
####  
####  #get the mean lat in each phylogroup
####  #mean.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,mean)
####  #colnames(mean.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  #df.phylogroups<-merge(df.phylogroups,mean.lat.phylogroups)
####  
####  #median is very correlated with mean
####  median.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,median)
####  colnames(median.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  df.phylogroups<-merge(df.phylogroups,median.lat.phylogroups)
####  
####  #sd.Median.Latitude is quite big in phylogroups
####  sd.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,function(x) sd(abs(x)))
####  colnames(sd.lat.phylogroups)[2]<-'phylogroup.sd.Median.Latitude'
####  sd.lat.phylogroups$phylogroup.sd.Median.Latitude[is.na(sd.lat.phylogroups$phylogroup.sd.Median.Latitude)]<-0
####  df.phylogroups<-merge(df.phylogroups,sd.lat.phylogroups)
####  
####  phylogroups.counts.stricts.df<-matrix(NA,ncol=4,nrow=0)
####  colnames(phylogroups.counts.stricts.df)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup','phylogroup.name')
####  #get the counts of strict tropical
####  for(d in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.table<-df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[d],]
####    phylogroup.table<-phylogroup.table[complete.cases(phylogroup.table),]
####    strict.counts<-count(phylogroup.table, c("phylogroup.name","strict.tropical"))
####    strict.counts<-strict.counts[,c('strict.tropical','freq')]
####    #add missing 0,1,2 counts (e.g if there are not strict tropical -  1s -  add a row with strict tropical =1 and freq =0)
####    if((0%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(0,0))
####    }
####    if((1%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(1,0))
####    }
####    if((2%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(2,0))
####    }
####    phylogroup.table.strict.counts<-as.data.frame(matrix(NA,nrow=1,ncol=3),stringsAsFactors = F)
####    colnames(phylogroup.table.strict.counts)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup')
####    phylogroup.table.strict.counts$n.strict.temperate.species.phylogroup<-strict.counts[strict.counts$strict.tropical==0,'freq']
####    phylogroup.table.strict.counts$n.strict.tropical.species.phylogroup<-strict.counts[strict.counts$strict.tropical==1,'freq']
####    phylogroup.table.strict.counts$n.strict.widespread.species.phylogroup<-strict.counts[strict.counts$strict.tropical==2,'freq']
####    phylogroup.table.strict.counts$phylogroup.name<-phylogroup.table$phylogroup.name[1]
####    phylogroups.counts.stricts.df<-rbind(phylogroups.counts.stricts.df,phylogroup.table.strict.counts)
####  }
####  count.strict.tropical.phylogroups<-aggregate(strict.tropical~phylogroup.name,data=df.phylogroups,count)
####  # phylogroups.counts.stricts.df<-phylogroups.counts.stricts.df[complete.cases(phylogroups.counts.stricts.df),]
####  df.phylogroups<-merge(df.phylogroups,phylogroups.counts.stricts.df,all.x=T)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  
####  df.phylogroups$turnover.rpanda1<-df.phylogroups$lambda.rpanda1+df.phylogroups$mu.rpanda1
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>9,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #filter by GBIF sampling
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>GBIF.sampling,]
####  
####  #dataset with strict tropical and strict temperate species
####  #strict tropical phylogroups =  phylogroups with > .7 of species strict.tropical
####  #strict temperate phylogroups =  phylogroups with > .7 of species strict.tropical
####  df.phylogroups$strict.tropical.character<-NA
####  df.phylogroups[(df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-1
####  df.phylogroups[(df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-0
####  
####  
####  
####  df.phylogroups$proportion.strict.tropical<-df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  df.phylogroups$proportion.strict.temperate<-df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  
####  df.phylogroups$lambda.rpanda1<-log10(df.phylogroups$lambda.rpanda1)
####  df.phylogroups$turnover.rpanda1<-log10(df.phylogroups$turnover.rpanda1)
####  
####  df.phylogroups.strict<-df.phylogroups[!is.na(df.phylogroups$strict.tropical.character),]
####  
####  
####  ###############this is for lambda
####  pdf(paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_new_GBIFsampling_log_modelave_size10.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Median.Latitude,xlab='phylogroup.Median.Latitude')
####  
####  if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####    wilcox.test.trop.vs.temp.lambda<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
####    t.test.trop.vs.temp.lambda<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.mu<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
####    t.test.trop.vs.temp.mu<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.ndr<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
####    t.test.trop.vs.temp.ndr<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.turnover<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    t.test.trop.vs.temp.turnover<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    
####    wilcox.t.results<-as.data.frame(rbind(c('log10.lambda',wilcox.test.trop.vs.temp.lambda$p.value,t.test.trop.vs.temp.lambda$p.value),c('mu',wilcox.test.trop.vs.temp.mu$p.value,t.test.trop.vs.temp.mu$p.value),c('ndr',wilcox.test.trop.vs.temp.ndr$p.value,t.test.trop.vs.temp.ndr$p.value),c('log10.turnover',wilcox.test.trop.vs.temp.turnover$p.value,t.test.trop.vs.temp.turnover$p.value)))
####    colnames(wilcox.t.results)<-c('parameter','pvalue.wilcox','pvalue.ttest')
####    write.table(wilcox.t.results,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_wilcox.ttest_table_new_GBIFsampling_log_modelave_size10.txt',sep=''),sep='\t',quote=F,row.names=F)
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],main=paste('lambda strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.lambda$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.lambda.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],main=paste('mu strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.mu$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='mu.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],main=paste('ndr strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.ndr$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='ndr.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],main=paste('turnover strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.turnover$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.turnover.rpanda1')
####    
####  }
####  
####  
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  #phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$strict.tropical.character,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1)
####  #colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$phylogroup.sd.Median.Latitude,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1,df.phylogroups$proportion.strict.tropical,df.phylogroups$proportion.strict.temperate)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','phylogroup.sd.Median.Latitude','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1','proportion.strict.tropical','proportion.strict.temperate')
####  
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  row.names(phylogroups.data)<-phylogroups.data$Tip
####  
####  #create binary trait (trop = 1, temp =0) for phyloglm
####  phylogroups.data$binary.tropical<-NA
####  
####  phylogroups.data[phylogroups.data$proportion.strict.temperate>strict.tropical.threshold&phylogroups.data$proportion.strict.tropical<strict.tropical.threshold,'binary.tropical']<-0
####  phylogroups.data[phylogroups.data$proportion.strict.temperate<strict.tropical.threshold&phylogroups.data$proportion.strict.tropical>strict.tropical.threshold,'binary.tropical']<-1
####  phylogroups.data2<-phylogroups.data[complete.cases(phylogroups.data),]
####  phylogroups.tree2<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data2$Tip))
####  phylogroups.tree.object2<-comparative.data(phy = phylogroups.tree2,data=phylogroups.data2,names.col='Tip',vcv=TRUE)
####  
####  
####  
####  pgls.Median.Latitude.lambda<-try(pgls(lambda.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.lambda)=='try-error'){
####    pgls.Median.Latitude.lambda<-try(gls(lambda.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.Median.Latitude.mu<-try(pgls(mu.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.mu)=='try-error'){
####    pgls.Median.Latitude.mu<-try(gls(mu.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.Median.Latitude.ndr<-try(pgls(ndr.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.ndr)=='try-error'){
####    pgls.Median.Latitude.ndr<-try(gls(ndr.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.Median.Latitude.turnover<-try(pgls(turnover.rpanda1~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.turnover)=='try-error'){
####    pgls.Median.Latitude.turnover<-try(gls(turnover.rpanda1 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  pgls.prop.trop.lambda<-try(pgls(proportion.strict.tropical~lambda.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.lambda)=='try-error'){
####    pgls.prop.trop.lambda<-try(gls(proportion.strict.tropical~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    if (class(pgls.prop.trop.lambda)=='try-error'){
####      lam<-seq(0,1,0.1)
####      fit.prop.trop.lambda<-list()
####      form.prop.trop.lambda<-proportion.strict.tropical~lambda.rpanda1
####      phylogroups.data<-phylogroups.data[complete.cases(phylogroups.data[,c('lambda.rpanda1')]),]
####      for (i in seq_along(lam)){
####        cor <- corPagel(lam[i], phy = drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data$Tip)), fixed = TRUE)
####        fit.prop.trop.lambda[[i]] <- gls(form.prop.trop.lambda, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML")
####      }
####      pgls.prop.trop.lambda<-fit.prop.trop.lambda[[which.min(sapply(fit.prop.trop.lambda, logLik))]]
####    }
####  }
####  pgls.prop.trop.mu<-try(pgls(proportion.strict.tropical~mu.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.mu)=='try-error'){
####    pgls.prop.trop.mu<-try(gls(proportion.strict.tropical~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.trop.ndr<-try(pgls(proportion.strict.tropical~ndr.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.ndr)=='try-error'){
####    pgls.prop.trop.ndr<-try(gls(proportion.strict.tropical~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.trop.turnover<-try(pgls(proportion.strict.tropical~turnover.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.turnover)=='try-error'){
####    pgls.prop.trop.turnover<-try(gls(proportion.strict.tropical~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####    if (class(pgls.prop.trop.turnover)=='try-error'){
####      lam<-seq(0,1,0.1)
####      fit.prop.trop.turnover<-list()
####      form.prop.trop.turnover<-proportion.strict.tropical~turnover.rpanda1
####      phylogroups.data<-phylogroups.data[complete.cases(phylogroups.data[,c('turnover.rpanda1')]),]
####      for (i in seq_along(lam)){
####        cor <- corPagel(lam[i], phy = drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data$Tip)), fixed = TRUE)
####        fit.prop.trop.turnover[[i]] <- gls(form.prop.trop.turnover, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML")
####      }
####      pgls.prop.trop.turnover<-fit.prop.trop.turnover[[which.min(sapply(fit.prop.trop.turnover, logLik))]]
####    }
####  }
####  
####  pgls.prop.temp.lambda<-try(pgls(proportion.strict.temperate~lambda.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.lambda)=='try-error'){
####    pgls.prop.temp.lambda<-try(gls(proportion.strict.temperate~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.mu<-try(pgls(proportion.strict.temperate~mu.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.mu)=='try-error'){
####    pgls.prop.temp.mu<-try(gls(proportion.strict.temperate~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.ndr<-try(pgls(proportion.strict.temperate~ndr.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.ndr)=='try-error'){
####    pgls.prop.temp.ndr<-try(gls(proportion.strict.temperate~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.turnover<-try(pgls(proportion.strict.temperate~turnover.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.turnover)=='try-error'){
####    pgls.prop.temp.turnover<-try(gls(proportion.strict.temperate~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  
####  #  phylogroups.tree.strict<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups.strict$Tip))
####  #  phylogroups.data.strict<-data.frame(df.phylogroups.strict$Tip,df.phylogroups.strict$phylogroup.Median.Latitude,df.phylogroups.strict$strict.tropical.character,df.phylogroups.strict$lambda.rpanda1,df.phylogroups.strict$mu.rpanda1,df.phylogroups.strict$ndr.rpanda1,df.phylogroups.strict$turnover.rpanda1)
####  #  colnames(phylogroups.data.strict)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####  #  phylogroups.tree.strict.object<-comparative.data(phy = phylogroups.tree.strict,data=phylogroups.data.strict,names.col='Tip',vcv=TRUE)
####  #  
####  #  pgls.strict.tropical.lambda<-try(pgls(strict.tropical.character~lambda.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.mu<-try(pgls(strict.tropical.character~mu.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.ndr<-try(pgls(strict.tropical.character~ndr.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.turnover<-try(pgls(strict.tropical.character~turnover.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  
####  #  
####  #  pgls.strict.tropical.lambda<-try(gls(strict.tropical.character~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.mu<-try(gls(strict.tropical.character~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.ndr<-try(gls(strict.tropical.character~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.turnover<-try(gls(strict.tropical.character~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  
####  #  strict.tropical.phylanova<-phylogroups.tree.strict.object$data$strict.tropical.character
####  #  names(strict.tropical.phylanova)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.lambda<-phylogroups.tree.strict.object$data$lambda.rpanda1
####  #  names(strict.tropical.phylanova.lambda)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.mu<-phylogroups.tree.strict.object$data$mu.rpanda1
####  #  names(strict.tropical.phylanova.mu)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.ndr<-phylogroups.tree.strict.object$data$ndr.rpanda1
####  #  names(strict.tropical.phylanova.ndr)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.turnover<-phylogroups.tree.strict.object$data$turnover.rpanda1
####  #  names(strict.tropical.phylanova.turnover)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  phylanova.strict.tropical.lambda<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.lambda,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.mu<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.mu,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.ndr<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.ndr,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.turnover<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.turnover,nsim = 100,posthoc = T)
####  #}
####  #summary(pgls.seedlambda)
####  #pgls.Median.Latitude.sigsq<-try(pgls(log10(lambda.rpanda1)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.Median.Latitude.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.log10.lambda',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.Median.Latitude.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.log10.lambda',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],8),round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.mu)
####    pgls.Median.Latitude.muresults<-c('pgls.Median.Latitude.mu',round(summary(pgls.Median.Latitude.mu)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.mu)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.mu)
####    pgls.Median.Latitude.muresults<-c('pgls.Median.Latitude.mu',round(summary(pgls.Median.Latitude.mu)$tTable[2,1],8),round(summary(pgls.Median.Latitude.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],8),round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.turnover)
####    pgls.Median.Latitude.turnoverresults<-c('pgls.Median.Latitude.log10.turnover',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.turnover)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.turnover)
####    pgls.Median.Latitude.turnoverresults<-c('pgls.Median.Latitude.log10.turnover',round(summary(pgls.Median.Latitude.turnover)$tTable[2,1],8),round(summary(pgls.Median.Latitude.turnover)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],8),round(summary(pgls.prop.trop.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$tTable[2,1],8),round(summary(pgls.prop.trop.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.mu)
####    pgls.prop.trop.muresults<-c('pgls.prop.trop.mu',round(summary(pgls.prop.trop.mu)$coefficients[2,1],8),round(summary(pgls.prop.trop.mu)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.mu)
####    pgls.prop.trop.muresults<-c('pgls.prop.trop.mu',round(summary(pgls.prop.trop.mu)$tTable[2,1],8),round(summary(pgls.prop.trop.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],8),round(summary(pgls.prop.trop.ndr)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$tTable[2,1],8),round(summary(pgls.prop.trop.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.turnover)
####    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.log10.turnover',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],8),round(summary(pgls.prop.trop.turnover)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.turnover)
####    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.log10.turnover',round(summary(pgls.prop.trop.turnover)$tTable[2,1],8),round(summary(pgls.prop.trop.turnover)$tTable[2,4],8))
####  }
####  
####  if(class(pgls.prop.temp.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],8),round(summary(pgls.prop.temp.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$tTable[2,1],8),round(summary(pgls.prop.temp.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.mu)
####    pgls.prop.temp.muresults<-c('pgls.prop.temp.mu',round(summary(pgls.prop.temp.mu)$coefficients[2,1],8),round(summary(pgls.prop.temp.mu)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.mu)
####    pgls.prop.temp.muresults<-c('pgls.prop.temp.mu',round(summary(pgls.prop.temp.mu)$tTable[2,1],8),round(summary(pgls.prop.temp.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],8),round(summary(pgls.prop.temp.ndr)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$tTable[2,1],8),round(summary(pgls.prop.temp.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.turnover)
####    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.log10.turnover',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],8),round(summary(pgls.prop.temp.turnover)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.turnover)
####    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.log10.turnover',round(summary(pgls.prop.temp.turnover)$tTable[2,1],8),round(summary(pgls.prop.temp.turnover)$tTable[2,4],8))
####  }
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #    
####  #   pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.log10.lambda',round(summary(pgls.strict.tropical.lambda)$coefficients[2,1],8),round(summary(pgls.strict.tropical.lambda)$coefficients[2,4],8))
####  #   pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$coefficients[2,1],8),round(summary(pgls.strict.tropical.mu)$coefficients[2,4],8))
####  #   pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$coefficients[2,1],8),round(summary(pgls.strict.tropical.ndr)$coefficients[2,4],8))
####  #   pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.log10.turnover',round(summary(pgls.strict.tropical.turnover)$coefficients[2,1],8),round(summary(pgls.strict.tropical.turnover)$coefficients[2,4],8))
####  #   
####  #   phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.log10.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####  #   phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####  #   phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####  #   phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.log10.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####  #   
####  # }
####  # 
####  
####  #   pgls.strict.tropical.lambda<-try(gls(strict.tropical.character~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.mu<-try(gls(strict.tropical.character~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.ndr<-try(gls(strict.tropical.character~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.turnover<-try(gls(strict.tropical.character~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   
####  # }
####  # 
####  # if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #   
####  #   pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.log10.lambda',round(summary(pgls.strict.tropical.lambda)$tTable[2,1],8),round(summary(pgls.strict.tropical.lambda)$tTable[2,4],8))
####  #   pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$tTable[2,1],8),round(summary(pgls.strict.tropical.mu)$tTable[2,4],8))
####  #   pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$tTable[2,1],8),round(summary(pgls.strict.tropical.ndr)$tTable[2,4],8))
####  #   pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.log10.turnover',round(summary(pgls.strict.tropical.turnover)$tTable[2,1],8),round(summary(pgls.strict.tropical.turnover)$tTable[2,4],8))
####  #   
####  #   phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.log10.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####  #   phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####  #   phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####  #   phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.log10.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####  #}
####  
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.muresults,pgls.Median.Latitude.ndrresults,pgls.Median.Latitude.turnoverresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.muresults,pgls.prop.trop.ndrresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.muresults,pgls.prop.temp.ndrresults,pgls.prop.temp.turnoverresults,pgls.strict.tropical.lambdaresults,pgls.strict.tropical.muresults,pgls.strict.tropical.ndrresults,pgls.strict.tropical.turnoverresults,phylanova.strict.tropical.lambda.results,phylanova.strict.tropical.mu.results,phylanova.strict.tropical.ndr.results,phylanova.strict.tropical.turnover.results))
####  #}else{
####  results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.muresults,pgls.Median.Latitude.ndrresults,pgls.Median.Latitude.turnoverresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.muresults,pgls.prop.trop.ndrresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.muresults,pgls.prop.temp.ndrresults,pgls.prop.temp.turnoverresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  if((min(table(phylogroups.tree.object2$data$binary.tropical))>1)&&(length(table(phylogroups.tree.object2$data$binary.tropical))>1)){
####    pgls.binary.trop.lambda<-phyloglm(binary.tropical~lambda.rpanda1, data = phylogroups.tree.object2$data, phy=phylogroups.tree2)
####    pgls.binary.trop.turnover<-phyloglm(binary.tropical~turnover.rpanda1, data = phylogroups.tree.object2$data, phy=phylogroups.tree2)
####    
####    pgls.binary.trop.lambda.results<-c('pgls.binary.trop.log10.lambda',round(summary(pgls.binary.trop.lambda)$coefficients[2,1],8),round(summary(pgls.binary.trop.lambda)$coefficients[2,4],8))
####    pgls.binary.trop.turnover.results<-c('pgls.binary.trop.log10.turnover',round(summary(pgls.binary.trop.turnover)$coefficients[2,1],8),round(summary(pgls.binary.trop.turnover)$coefficients[2,4],8))
####    results.binary.trop<-as.data.frame(rbind(pgls.binary.trop.lambda.results,pgls.binary.trop.turnover.results))
####    colnames(results.binary.trop)<-c('analysis','slope','pvalue')
####    results.df<-rbind(results.df,results.binary.trop)
####  }
####  
####  
####  
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_',strict.tropical.threshold,'.tropthreshold_pgls_table_new_GBIFsampling_log_modelave_size10.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='log10.lambda')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='mu')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='ndr')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='log10.turnover')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #}
####  dev.off()
####}



####run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_weighted<-function(tree,minage,maxage,mincladesize,ncores,sampling,table,strict.tropical.threshold,GBIF.sampling,loadRPANDA){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  if(loadRPANDA==F){
####    phylogroups<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'.RDS',sep=''))
####    phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####    length(unique(phylogroups.species))
####    cat('measuring diversification/trait evolution','\n')
####    phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_phylogroup(x,table,sampling,mincladesize=mincladesize))
####    
####  }else if(loadRPANDA==T){
####    cat('loading RPANDA phylogroups object','\n')
####    phylogroups.RPANDA<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'_RPANDA.RDS',sep=''))
####  }
####  
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of lat data points per phylogroup
####  phylogroup.species.lat<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.lat[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude)])))
####  }
####  phylogroup.species.lat<-as.data.frame(phylogroup.species.lat,stringsAsFactors = F)
####  colnames(phylogroup.species.lat)<-c('phylogroup.name','phylogroup.species.with.lat.data')
####  phylogroup.species.lat$phylogroup.species.with.lat.data<-as.numeric(phylogroup.species.lat$phylogroup.species.with.lat.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.lat,all.x=TRUE)
####  
####  #get the mean lat in each phylogroup
####  #mean.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,mean)
####  #colnames(mean.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  #df.phylogroups<-merge(df.phylogroups,mean.lat.phylogroups)
####  
####  #median is very correlated with mean
####  median.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,median)
####  colnames(median.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  df.phylogroups<-merge(df.phylogroups,median.lat.phylogroups)
####  
####  #sd.Median.Latitude is quite big in phylogroups
####  sd.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,function(x) sd(abs(x)))
####  colnames(sd.lat.phylogroups)[2]<-'phylogroup.sd.Median.Latitude'
####  sd.lat.phylogroups$phylogroup.sd.Median.Latitude[is.na(sd.lat.phylogroups$phylogroup.sd.Median.Latitude)]<-0
####  df.phylogroups<-merge(df.phylogroups,sd.lat.phylogroups)
####  
####  phylogroups.counts.stricts.df<-matrix(NA,ncol=4,nrow=0)
####  colnames(phylogroups.counts.stricts.df)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup','phylogroup.name')
####  #get the counts of strict tropical
####  for(d in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.table<-df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[d],]
####    phylogroup.table<-phylogroup.table[complete.cases(phylogroup.table),]
####    strict.counts<-count(phylogroup.table, c("phylogroup.name","strict.tropical"))
####    strict.counts<-strict.counts[,c('strict.tropical','freq')]
####    #add missing 0,1,2 counts (e.g if there are not strict tropical -  1s -  add a row with strict tropical =1 and freq =0)
####    if((0%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(0,0))
####    }
####    if((1%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(1,0))
####    }
####    if((2%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(2,0))
####    }
####    phylogroup.table.strict.counts<-as.data.frame(matrix(NA,nrow=1,ncol=3),stringsAsFactors = F)
####    colnames(phylogroup.table.strict.counts)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup')
####    phylogroup.table.strict.counts$n.strict.temperate.species.phylogroup<-strict.counts[strict.counts$strict.tropical==0,'freq']
####    phylogroup.table.strict.counts$n.strict.tropical.species.phylogroup<-strict.counts[strict.counts$strict.tropical==1,'freq']
####    phylogroup.table.strict.counts$n.strict.widespread.species.phylogroup<-strict.counts[strict.counts$strict.tropical==2,'freq']
####    phylogroup.table.strict.counts$phylogroup.name<-phylogroup.table$phylogroup.name[1]
####    phylogroups.counts.stricts.df<-rbind(phylogroups.counts.stricts.df,phylogroup.table.strict.counts)
####  }
####  count.strict.tropical.phylogroups<-aggregate(strict.tropical~phylogroup.name,data=df.phylogroups,count)
####  # phylogroups.counts.stricts.df<-phylogroups.counts.stricts.df[complete.cases(phylogroups.counts.stricts.df),]
####  df.phylogroups<-merge(df.phylogroups,phylogroups.counts.stricts.df,all.x=T)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  
####  df.phylogroups$turnover.rpanda1<-df.phylogroups$lambda.rpanda1+df.phylogroups$mu.rpanda1
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #filter by GBIF sampling
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>GBIF.sampling,]
####  
####  #dataset with strict tropical and strict temperate species
####  #strict tropical phylogroups =  phylogroups with > .7 of species strict.tropical
####  #strict temperate phylogroups =  phylogroups with > .7 of species strict.tropical
####  df.phylogroups$strict.tropical.character<-NA
####  df.phylogroups[(df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-1
####  df.phylogroups[(df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-0
####  
####  
####  
####  df.phylogroups$proportion.strict.tropical<-df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  df.phylogroups$proportion.strict.temperate<-df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  
####  df.phylogroups$lambda.rpanda1<-log10(df.phylogroups$lambda.rpanda1)
####  df.phylogroups$turnover.rpanda1<-log10(df.phylogroups$turnover.rpanda1)
####  
####  df.phylogroups.strict<-df.phylogroups[!is.na(df.phylogroups$strict.tropical.character),]
####  
####  
####  ###############this is for lambda
####  pdf(paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_new_GBIFsampling_log_weightsd.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Median.Latitude,xlab='phylogroup.Median.Latitude')
####  
####  if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####    wilcox.test.trop.vs.temp.lambda<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
####    t.test.trop.vs.temp.lambda<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.mu<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
####    t.test.trop.vs.temp.mu<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.ndr<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
####    t.test.trop.vs.temp.ndr<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'])
####    
####    wilcox.test.trop.vs.temp.turnover<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    t.test.trop.vs.temp.turnover<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    
####    wilcox.t.results<-as.data.frame(rbind(c('log10.lambda',wilcox.test.trop.vs.temp.lambda$p.value,t.test.trop.vs.temp.lambda$p.value),c('mu',wilcox.test.trop.vs.temp.mu$p.value,t.test.trop.vs.temp.mu$p.value),c('ndr',wilcox.test.trop.vs.temp.ndr$p.value,t.test.trop.vs.temp.ndr$p.value),c('log10.turnover',wilcox.test.trop.vs.temp.turnover$p.value,t.test.trop.vs.temp.turnover$p.value)))
####    colnames(wilcox.t.results)<-c('parameter','pvalue.wilcox','pvalue.ttest')
####    write.table(wilcox.t.results,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_wilcox.ttest_table_new_GBIFsampling_log_weightsd.txt',sep=''),sep='\t',quote=F,row.names=F)
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],main=paste('lambda strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.lambda$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.lambda.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],main=paste('mu strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.mu$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='mu.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],main=paste('ndr strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.ndr$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='ndr.rpanda1')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],main=paste('turnover strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.turnover$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.turnover.rpanda1')
####    
####  }
####  
####  
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  #phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$strict.tropical.character,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1)
####  #colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$phylogroup.sd.Median.Latitude,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1,df.phylogroups$proportion.strict.tropical,df.phylogroups$proportion.strict.temperate)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','phylogroup.sd.Median.Latitude','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1','proportion.strict.tropical','proportion.strict.temperate')
####  
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  row.names(phylogroups.data)<-phylogroups.data$Tip
####  
####  #create binary trait (trop = 1, temp =0) for phyloglm
####  phylogroups.data$binary.tropical<-NA
####  
####  phylogroups.data[phylogroups.data$proportion.strict.temperate>strict.tropical.threshold&phylogroups.data$proportion.strict.tropical<strict.tropical.threshold,'binary.tropical']<-0
####  phylogroups.data[phylogroups.data$proportion.strict.temperate<strict.tropical.threshold&phylogroups.data$proportion.strict.tropical>strict.tropical.threshold,'binary.tropical']<-1
####  phylogroups.data2<-phylogroups.data[complete.cases(phylogroups.data),]
####  phylogroups.tree2<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data2$Tip))
####  phylogroups.tree.object2<-comparative.data(phy = phylogroups.tree2,data=phylogroups.data2,names.col='Tip',vcv=TRUE)
####  
####  lam<-seq(0,1,0.1)
####  fit.Median.Latitude.lambda<-list()
####  fit.Median.Latitude.mu<-list()
####  fit.Median.Latitude.ndr<-list()
####  fit.Median.Latitude.turnover<-list()
####  form.Median.Latitude.lambda<-lambda.rpanda1 ~ abs(phylogroup.Median.Latitude)
####  form.Median.Latitude.mu<-mu.rpanda1 ~ abs(phylogroup.Median.Latitude)
####  form.Median.Latitude.ndr<-ndr.rpanda1 ~ abs(phylogroup.Median.Latitude)
####  form.Median.Latitude.turnover<-turnover.rpanda1 ~ abs(phylogroup.Median.Latitude)
####  phylogroups.data<-phylogroups.data[complete.cases(phylogroups.data[,c('lambda.rpanda1','turnover.rpanda1')]),]
####  for (i in seq_along(lam)){
####    cor <- corPagel(lam[i], phy = drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data$Tip)), fixed = TRUE)
####    fit.Median.Latitude.lambda[[i]] <- gls(form.Median.Latitude.lambda, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML",weights=~1/phylogroup.sd.Median.Latitude)
####    fit.Median.Latitude.mu[[i]] <- gls(form.Median.Latitude.mu, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML",weights=~1/phylogroup.sd.Median.Latitude)
####    fit.Median.Latitude.ndr[[i]] <- gls(form.Median.Latitude.ndr, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML",weights=~1/phylogroup.sd.Median.Latitude)
####    fit.Median.Latitude.turnover[[i]] <- gls(form.Median.Latitude.turnover, correlation = cor, data = phylogroups.data, na.action = na.exclude, method = "ML",weights=~1/phylogroup.sd.Median.Latitude)
####    
####    
####  }
####  pgls.Median.Latitude.lambda<-fit.Median.Latitude.lambda[[which.min(sapply(fit.Median.Latitude.lambda, logLik))]]
####  pgls.Median.Latitude.mu<-fit.Median.Latitude.mu[[which.min(sapply(fit.Median.Latitude.mu, logLik))]]
####  pgls.Median.Latitude.ndr<-fit.Median.Latitude.ndr[[which.min(sapply(fit.Median.Latitude.ndr, logLik))]]
####  pgls.Median.Latitude.turnover<-fit.Median.Latitude.turnover[[which.min(sapply(fit.Median.Latitude.turnover, logLik))]]
####  
####  pgls.prop.trop.lambda<-try(pgls(proportion.strict.tropical~lambda.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.lambda)=='try-error'){
####    pgls.prop.trop.lambda<-try(gls(proportion.strict.tropical~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.trop.mu<-try(pgls(proportion.strict.tropical~mu.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.mu)=='try-error'){
####    pgls.prop.trop.mu<-try(gls(proportion.strict.tropical~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.trop.ndr<-try(pgls(proportion.strict.tropical~ndr.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.ndr)=='try-error'){
####    pgls.prop.trop.ndr<-try(gls(proportion.strict.tropical~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.trop.turnover<-try(pgls(proportion.strict.tropical~turnover.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.turnover)=='try-error'){
####    pgls.prop.trop.turnover<-try(gls(proportion.strict.tropical~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  pgls.prop.temp.lambda<-try(pgls(proportion.strict.temperate~lambda.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.lambda)=='try-error'){
####    pgls.prop.temp.lambda<-try(gls(proportion.strict.temperate~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.mu<-try(pgls(proportion.strict.temperate~mu.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.mu)=='try-error'){
####    pgls.prop.temp.mu<-try(gls(proportion.strict.temperate~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.ndr<-try(pgls(proportion.strict.temperate~ndr.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.ndr)=='try-error'){
####    pgls.prop.temp.ndr<-try(gls(proportion.strict.temperate~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.turnover<-try(pgls(proportion.strict.temperate~turnover.rpanda1, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.turnover)=='try-error'){
####    pgls.prop.temp.turnover<-try(gls(proportion.strict.temperate~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  
####  #  phylogroups.tree.strict<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups.strict$Tip))
####  #  phylogroups.data.strict<-data.frame(df.phylogroups.strict$Tip,df.phylogroups.strict$phylogroup.Median.Latitude,df.phylogroups.strict$strict.tropical.character,df.phylogroups.strict$lambda.rpanda1,df.phylogroups.strict$mu.rpanda1,df.phylogroups.strict$ndr.rpanda1,df.phylogroups.strict$turnover.rpanda1)
####  #  colnames(phylogroups.data.strict)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####  #  phylogroups.tree.strict.object<-comparative.data(phy = phylogroups.tree.strict,data=phylogroups.data.strict,names.col='Tip',vcv=TRUE)
####  #  
####  #  pgls.strict.tropical.lambda<-try(pgls(strict.tropical.character~lambda.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.mu<-try(pgls(strict.tropical.character~mu.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.ndr<-try(pgls(strict.tropical.character~ndr.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.turnover<-try(pgls(strict.tropical.character~turnover.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  
####  #  
####  #  pgls.strict.tropical.lambda<-try(gls(strict.tropical.character~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.mu<-try(gls(strict.tropical.character~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.ndr<-try(gls(strict.tropical.character~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.turnover<-try(gls(strict.tropical.character~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  
####  #  strict.tropical.phylanova<-phylogroups.tree.strict.object$data$strict.tropical.character
####  #  names(strict.tropical.phylanova)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.lambda<-phylogroups.tree.strict.object$data$lambda.rpanda1
####  #  names(strict.tropical.phylanova.lambda)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.mu<-phylogroups.tree.strict.object$data$mu.rpanda1
####  #  names(strict.tropical.phylanova.mu)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.ndr<-phylogroups.tree.strict.object$data$ndr.rpanda1
####  #  names(strict.tropical.phylanova.ndr)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.turnover<-phylogroups.tree.strict.object$data$turnover.rpanda1
####  #  names(strict.tropical.phylanova.turnover)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  phylanova.strict.tropical.lambda<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.lambda,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.mu<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.mu,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.ndr<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.ndr,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.turnover<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.turnover,nsim = 100,posthoc = T)
####  #}
####  #summary(pgls.seedlambda)
####  #pgls.Median.Latitude.sigsq<-try(pgls(log10(lambda.rpanda1)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.Median.Latitude.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.log10.lambda',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.Median.Latitude.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.log10.lambda',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],8),round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.mu)
####    pgls.Median.Latitude.muresults<-c('pgls.Median.Latitude.mu',round(summary(pgls.Median.Latitude.mu)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.mu)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.mu)
####    pgls.Median.Latitude.muresults<-c('pgls.Median.Latitude.mu',round(summary(pgls.Median.Latitude.mu)$tTable[2,1],8),round(summary(pgls.Median.Latitude.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],8),round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.turnover)
####    pgls.Median.Latitude.turnoverresults<-c('pgls.Median.Latitude.log10.turnover',round(summary(pgls.Median.Latitude.turnover)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.turnover)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.turnover)
####    pgls.Median.Latitude.turnoverresults<-c('pgls.Median.Latitude.log10.turnover',round(summary(pgls.Median.Latitude.turnover)$tTable[2,1],8),round(summary(pgls.Median.Latitude.turnover)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],8),round(summary(pgls.prop.trop.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$tTable[2,1],8),round(summary(pgls.prop.trop.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.mu)
####    pgls.prop.trop.muresults<-c('pgls.prop.trop.mu',round(summary(pgls.prop.trop.mu)$coefficients[2,1],8),round(summary(pgls.prop.trop.mu)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.mu)
####    pgls.prop.trop.muresults<-c('pgls.prop.trop.mu',round(summary(pgls.prop.trop.mu)$tTable[2,1],8),round(summary(pgls.prop.trop.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],8),round(summary(pgls.prop.trop.ndr)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$tTable[2,1],8),round(summary(pgls.prop.trop.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.turnover)
####    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.log10.turnover',round(summary(pgls.prop.trop.turnover)$coefficients[2,1],8),round(summary(pgls.prop.trop.turnover)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.turnover)
####    pgls.prop.trop.turnoverresults<-c('pgls.prop.trop.log10.turnover',round(summary(pgls.prop.trop.turnover)$tTable[2,1],8),round(summary(pgls.prop.trop.turnover)$tTable[2,4],8))
####  }
####  
####  if(class(pgls.prop.temp.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],8),round(summary(pgls.prop.temp.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.lambda)=='gls'){
####    plot(phylogroups.data$lambda.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$tTable[2,1],8),round(summary(pgls.prop.temp.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.mu)=='pgls'){
####    plot(phylogroups.tree.object$data$mu.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.mu)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.mu)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.mu)
####    pgls.prop.temp.muresults<-c('pgls.prop.temp.mu',round(summary(pgls.prop.temp.mu)$coefficients[2,1],8),round(summary(pgls.prop.temp.mu)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.mu)=='gls'){
####    plot(phylogroups.data$mu.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='mu',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.mu)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.mu)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.mu)
####    pgls.prop.temp.muresults<-c('pgls.prop.temp.mu',round(summary(pgls.prop.temp.mu)$tTable[2,1],8),round(summary(pgls.prop.temp.mu)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],8),round(summary(pgls.prop.temp.ndr)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.ndr)=='gls'){
####    plot(phylogroups.data$ndr.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$tTable[2,1],8),round(summary(pgls.prop.temp.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.turnover)=='pgls'){
####    plot(phylogroups.tree.object$data$turnover.rpanda1,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.turnover)
####    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.log10.turnover',round(summary(pgls.prop.temp.turnover)$coefficients[2,1],8),round(summary(pgls.prop.temp.turnover)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.turnover)=='gls'){
####    plot(phylogroups.data$turnover.rpanda1,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.turnover',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.turnover)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.turnover)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.turnover)
####    pgls.prop.temp.turnoverresults<-c('pgls.prop.temp.log10.turnover',round(summary(pgls.prop.temp.turnover)$tTable[2,1],8),round(summary(pgls.prop.temp.turnover)$tTable[2,4],8))
####  }
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #    
####  #   pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.log10.lambda',round(summary(pgls.strict.tropical.lambda)$coefficients[2,1],8),round(summary(pgls.strict.tropical.lambda)$coefficients[2,4],8))
####  #   pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$coefficients[2,1],8),round(summary(pgls.strict.tropical.mu)$coefficients[2,4],8))
####  #   pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$coefficients[2,1],8),round(summary(pgls.strict.tropical.ndr)$coefficients[2,4],8))
####  #   pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.log10.turnover',round(summary(pgls.strict.tropical.turnover)$coefficients[2,1],8),round(summary(pgls.strict.tropical.turnover)$coefficients[2,4],8))
####  #   
####  #   phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.log10.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####  #   phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####  #   phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####  #   phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.log10.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####  #   
####  # }
####  # 
####  
####  #   pgls.strict.tropical.lambda<-try(gls(strict.tropical.character~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.mu<-try(gls(strict.tropical.character~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.ndr<-try(gls(strict.tropical.character~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.turnover<-try(gls(strict.tropical.character~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   
####  # }
####  # 
####  # if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #   
####  #   pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.log10.lambda',round(summary(pgls.strict.tropical.lambda)$tTable[2,1],8),round(summary(pgls.strict.tropical.lambda)$tTable[2,4],8))
####  #   pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$tTable[2,1],8),round(summary(pgls.strict.tropical.mu)$tTable[2,4],8))
####  #   pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$tTable[2,1],8),round(summary(pgls.strict.tropical.ndr)$tTable[2,4],8))
####  #   pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.log10.turnover',round(summary(pgls.strict.tropical.turnover)$tTable[2,1],8),round(summary(pgls.strict.tropical.turnover)$tTable[2,4],8))
####  #   
####  #   phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.log10.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####  #   phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####  #   phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####  #   phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.log10.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####  #}
####  
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.muresults,pgls.Median.Latitude.ndrresults,pgls.Median.Latitude.turnoverresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.muresults,pgls.prop.trop.ndrresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.muresults,pgls.prop.temp.ndrresults,pgls.prop.temp.turnoverresults,pgls.strict.tropical.lambdaresults,pgls.strict.tropical.muresults,pgls.strict.tropical.ndrresults,pgls.strict.tropical.turnoverresults,phylanova.strict.tropical.lambda.results,phylanova.strict.tropical.mu.results,phylanova.strict.tropical.ndr.results,phylanova.strict.tropical.turnover.results))
####  #}else{
####  results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.muresults,pgls.Median.Latitude.ndrresults,pgls.Median.Latitude.turnoverresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.muresults,pgls.prop.trop.ndrresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.muresults,pgls.prop.temp.ndrresults,pgls.prop.temp.turnoverresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  if((min(table(phylogroups.tree.object2$data$binary.tropical))>1)&&(length(table(phylogroups.tree.object2$data$binary.tropical))>1)){
####    pgls.binary.trop.lambda<-phyloglm(binary.tropical~lambda.rpanda1, data = phylogroups.tree.object2$data, phy=phylogroups.tree2)
####    pgls.binary.trop.turnover<-phyloglm(binary.tropical~turnover.rpanda1, data = phylogroups.tree.object2$data, phy=phylogroups.tree2)
####    
####    pgls.binary.trop.lambda.results<-c('pgls.binary.trop.log10.lambda',round(summary(pgls.binary.trop.lambda)$coefficients[2,1],8),round(summary(pgls.binary.trop.lambda)$coefficients[2,4],8))
####    pgls.binary.trop.turnover.results<-c('pgls.binary.trop.log10.turnover',round(summary(pgls.binary.trop.turnover)$coefficients[2,1],8),round(summary(pgls.binary.trop.turnover)$coefficients[2,4],8))
####    results.binary.trop<-as.data.frame(rbind(pgls.binary.trop.lambda.results,pgls.binary.trop.turnover.results))
####    colnames(results.binary.trop)<-c('analysis','slope','pvalue')
####    results.df<-rbind(results.df,results.binary.trop)
####  }
####  
####  
####  
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/RPANDA/RPANDA_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_',strict.tropical.threshold,'.tropthreshold_pgls_table_new_GBIFsampling_log_weightsd.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='log10.lambda')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='mu')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='ndr')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='log10.turnover')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #}
####  dev.off()
####}
####

####run_clades_MS_pgls_savedphylogroups_abslat_stricttropthreshold_log<-function(tree,minage,maxage,mincladesize,ncores,sampling,table,strict.tropical.threshold,GBIF.sampling,loadRPANDA){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  if(loadRPANDA==F){
####    phylogroups<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'.RDS',sep=''))
####    phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####    length(unique(phylogroups.species))
####    cat('measuring diversification/trait evolution','\n')
####    phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_EB_z0_phylogroup_select_lat(x,table,sampling,mincladesize=mincladesize))
####    
####  }else if(loadRPANDA==T){
####    cat('loading RPANDA phylogroups object','\n')
####    phylogroups.MS<-readRDS(file=paste('~/Dropbox/Work_in_progress/LDG_plants/phylogroups_',minage,'_',maxage,'_size',mincladesize,'_MS.RDS',sep=''))
####  }
####  
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  #count the number of lat data points per phylogroup
####  phylogroup.species.lat<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.lat[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Median.Latitude)])))
####  }
####  phylogroup.species.lat<-as.data.frame(phylogroup.species.lat,stringsAsFactors = F)
####  colnames(phylogroup.species.lat)<-c('phylogroup.name','phylogroup.species.with.lat.data')
####  phylogroup.species.lat$phylogroup.species.with.lat.data<-as.numeric(phylogroup.species.lat$phylogroup.species.with.lat.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.lat,all.x=TRUE)
####  
####  #get the mean lat in each phylogroup
####  #mean.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,mean)
####  #colnames(mean.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  #df.phylogroups<-merge(df.phylogroups,mean.lat.phylogroups)
####  
####  #median is very correlated with mean
####  median.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,median)
####  colnames(median.lat.phylogroups)[2]<-'phylogroup.Median.Latitude'
####  df.phylogroups<-merge(df.phylogroups,median.lat.phylogroups)
####  
####  #sd.Median.Latitude is quite big in phylogroups
####  sd.lat.phylogroups<-aggregate(Median.Latitude~phylogroup.name,data=df.phylogroups,function(x) sd(abs(x)))
####  colnames(sd.lat.phylogroups)[2]<-'phylogroup.sd.Median.Latitude'
####  sd.lat.phylogroups$phylogroup.sd.Median.Latitude[is.na(sd.lat.phylogroups$phylogroup.sd.Median.Latitude)]<-0
####  df.phylogroups<-merge(df.phylogroups,sd.lat.phylogroups)
####  
####  phylogroups.counts.stricts.df<-matrix(NA,ncol=4,nrow=0)
####  colnames(phylogroups.counts.stricts.df)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup','phylogroup.name')
####  #get the counts of strict tropical
####  for(d in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.table<-df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[d],]
####    phylogroup.table<-phylogroup.table[complete.cases(phylogroup.table),]
####    strict.counts<-count(phylogroup.table, c("phylogroup.name","strict.tropical"))
####    strict.counts<-strict.counts[,c('strict.tropical','freq')]
####    #add missing 0,1,2 counts (e.g if there are not strict tropical -  1s -  add a row with strict tropical =1 and freq =0)
####    if((0%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(0,0))
####    }
####    if((1%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(1,0))
####    }
####    if((2%in%strict.counts$strict.tropical)==FALSE){
####      strict.counts<-rbind(strict.counts,c(2,0))
####    }
####    phylogroup.table.strict.counts<-as.data.frame(matrix(NA,nrow=1,ncol=3),stringsAsFactors = F)
####    colnames(phylogroup.table.strict.counts)<-c('n.strict.tropical.species.phylogroup','n.strict.temperate.species.phylogroup','n.strict.widespread.species.phylogroup')
####    phylogroup.table.strict.counts$n.strict.temperate.species.phylogroup<-strict.counts[strict.counts$strict.tropical==0,'freq']
####    phylogroup.table.strict.counts$n.strict.tropical.species.phylogroup<-strict.counts[strict.counts$strict.tropical==1,'freq']
####    phylogroup.table.strict.counts$n.strict.widespread.species.phylogroup<-strict.counts[strict.counts$strict.tropical==2,'freq']
####    phylogroup.table.strict.counts$phylogroup.name<-phylogroup.table$phylogroup.name[1]
####    phylogroups.counts.stricts.df<-rbind(phylogroups.counts.stricts.df,phylogroup.table.strict.counts)
####  }
####  count.strict.tropical.phylogroups<-aggregate(strict.tropical~phylogroup.name,data=df.phylogroups,count)
####  # phylogroups.counts.stricts.df<-phylogroups.counts.stricts.df[complete.cases(phylogroups.counts.stricts.df),]
####  df.phylogroups<-merge(df.phylogroups,phylogroups.counts.stricts.df,all.x=T)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  
####  #df.phylogroups$turnover.ms<-df.phylogroups$lambda.ms09+df.phylogroups$mu.rpanda1
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #filter by GBIF sampling
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>GBIF.sampling,]
####  
####  #dataset with strict tropical and strict temperate species
####  #strict tropical phylogroups =  phylogroups with > .7 of species strict.tropical
####  #strict temperate phylogroups =  phylogroups with > .7 of species strict.tropical
####  df.phylogroups$strict.tropical.character<-NA
####  df.phylogroups[(df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-1
####  df.phylogroups[(df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data)>strict.tropical.threshold,'strict.tropical.character']<-0
####  
####  
####  
####  df.phylogroups$proportion.strict.tropical<-df.phylogroups$n.strict.tropical.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  df.phylogroups$proportion.strict.temperate<-df.phylogroups$n.strict.temperate.species.phylogroup/df.phylogroups$phylogroup.species.with.lat.data
####  
####  df.phylogroups$lambda.ms0<-log10(df.phylogroups$lambda.ms0)
####  df.phylogroups$lambda.ms05<-log10(df.phylogroups$lambda.ms05)
####  df.phylogroups$lambda.ms09<-log10(df.phylogroups$lambda.ms09)
####  #df.phylogroups$ndr.ms05<-log10(df.phylogroups$ndr.ms05)
####  #df.phylogroups$ndr.ms09<-log10(df.phylogroups$ndr.ms09)
####  
####  
####  #df.phylogroups$turnover.rpanda1<-log10(df.phylogroups$turnover.rpanda1)
####  
####  df.phylogroups.strict<-df.phylogroups[!is.na(df.phylogroups$strict.tropical.character),]
####  
####  
####  ###############this is for lambda
####  pdf(paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/MS/MS_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_new_GBIFsampling_log.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Median.Latitude,xlab='phylogroup.Median.Latitude')
####  
####  if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####    wilcox.test.trop.vs.temp.lambda<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.ms05'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.ms05'])
####    t.test.trop.vs.temp.lambda<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.ms05'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.ms05'])
####    
####    wilcox.test.trop.vs.temp.ndr<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.ms05'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.ms05'])
####    t.test.trop.vs.temp.ndr<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.ms05'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.ms05'])
####    
####    #wilcox.test.trop.vs.temp.turnover<-wilcox.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    #t.test.trop.vs.temp.turnover<-t.test(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'])
####    
####    wilcox.t.results<-as.data.frame(rbind(c('log10.lambda',wilcox.test.trop.vs.temp.lambda$p.value,t.test.trop.vs.temp.lambda$p.value),c('ndr',wilcox.test.trop.vs.temp.ndr$p.value,t.test.trop.vs.temp.ndr$p.value)))
####    colnames(wilcox.t.results)<-c('parameter','pvalue.wilcox','pvalue.ttest')
####    write.table(wilcox.t.results,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/MS/MS_EB_z0_lambda_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size',strict.tropical.threshold,'.tropthreshold_wilcox.ttest_table_new_GBIFsampling_log.txt',sep=''),sep='\t',quote=F,row.names=F)
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.ms05'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.ms05'],main=paste('lambda strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.lambda$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.lambda.ms05')
####    boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.ms05'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.ms05'],main=paste('ndr strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.ndr$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='ndr.ms05')
####    #boxplot(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],main=paste('turnover strict.temp vs strict.trop (t test pval = ',round(t.test.trop.vs.temp.turnover$p.value,3),')',sep=''),cex.main=.5,names=c(paste('temp(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,]),')',sep=''),paste('trop(',nrow(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,]),')',sep='')),notch=T,ylab='log10.turnover.rpanda1')
####    
####  }
####  
####  
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  #phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$strict.tropical.character,df.phylogroups$lambda.ms0,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$turnover.rpanda1)
####  #colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.ms0','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Median.Latitude,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$proportion.strict.tropical,df.phylogroups$proportion.strict.temperate)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Median.Latitude','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','proportion.strict.tropical','proportion.strict.temperate')
####  
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  row.names(phylogroups.data)<-phylogroups.data$Tip
####  
####  #create binary trait (trop = 1, temp =0) for phyloglm
####  phylogroups.data$binary.tropical<-NA
####  
####  phylogroups.data[phylogroups.data$proportion.strict.temperate>strict.tropical.threshold&phylogroups.data$proportion.strict.tropical<strict.tropical.threshold,'binary.tropical']<-0
####  phylogroups.data[phylogroups.data$proportion.strict.temperate<strict.tropical.threshold&phylogroups.data$proportion.strict.tropical>strict.tropical.threshold,'binary.tropical']<-1
####  phylogroups.data2<-phylogroups.data[complete.cases(phylogroups.data),]
####  phylogroups.tree2<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,phylogroups.data2$Tip))
####  phylogroups.tree.object2<-comparative.data(phy = phylogroups.tree2,data=phylogroups.data2,names.col='Tip',vcv=TRUE)
####  
####  
####  
####  pgls.Median.Latitude.lambda<-try(pgls(lambda.ms05~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.lambda)=='try-error'){
####    pgls.Median.Latitude.lambda<-try(gls(lambda.ms05 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.Median.Latitude.ndr<-try(pgls(ndr.ms05~abs(phylogroup.Median.Latitude), data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.Median.Latitude.ndr)=='try-error'){
####    pgls.Median.Latitude.ndr<-try(gls(ndr.ms05 ~ abs(phylogroup.Median.Latitude), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  pgls.prop.trop.lambda<-try(pgls(proportion.strict.tropical~lambda.ms05, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.lambda)=='try-error'){
####    pgls.prop.trop.lambda<-try(gls(proportion.strict.tropical~lambda.ms05, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.trop.ndr<-try(pgls(proportion.strict.tropical~ndr.ms05, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.trop.ndr)=='try-error'){
####    pgls.prop.trop.ndr<-try(gls(proportion.strict.tropical~ndr.ms05, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  pgls.prop.temp.lambda<-try(pgls(proportion.strict.temperate~lambda.ms05, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.lambda)=='try-error'){
####    pgls.prop.temp.lambda<-try(gls(proportion.strict.temperate~lambda.ms05, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  pgls.prop.temp.ndr<-try(pgls(proportion.strict.temperate~ndr.ms05, data = phylogroups.tree.object, lambda='ML'))
####  if(class(pgls.prop.temp.ndr)=='try-error'){
####    pgls.prop.temp.ndr<-try(gls(proportion.strict.temperate~ndr.ms05, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML"))
####  }
####  
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  
####  #  phylogroups.tree.strict<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups.strict$Tip))
####  #  phylogroups.data.strict<-data.frame(df.phylogroups.strict$Tip,df.phylogroups.strict$phylogroup.Median.Latitude,df.phylogroups.strict$strict.tropical.character,df.phylogroups.strict$lambda.rpanda1,df.phylogroups.strict$mu.rpanda1,df.phylogroups.strict$ndr.rpanda1,df.phylogroups.strict$turnover.rpanda1)
####  #  colnames(phylogroups.data.strict)<-c('Tip','phylogroup.Median.Latitude','strict.tropical.character','lambda.rpanda1','mu.rpanda1','ndr.rpanda1','turnover.rpanda1')
####  #  phylogroups.tree.strict.object<-comparative.data(phy = phylogroups.tree.strict,data=phylogroups.data.strict,names.col='Tip',vcv=TRUE)
####  #  
####  #  pgls.strict.tropical.lambda<-try(pgls(strict.tropical.character~lambda.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.mu<-try(pgls(strict.tropical.character~mu.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.ndr<-try(pgls(strict.tropical.character~ndr.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  pgls.strict.tropical.turnover<-try(pgls(strict.tropical.character~turnover.rpanda1, data = phylogroups.tree.strict.object, lambda='ML'))
####  #  
####  #  
####  #  pgls.strict.tropical.lambda<-try(gls(strict.tropical.character~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.mu<-try(gls(strict.tropical.character~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.ndr<-try(gls(strict.tropical.character~ndr.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  pgls.strict.tropical.turnover<-try(gls(strict.tropical.character~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #  
####  #  strict.tropical.phylanova<-phylogroups.tree.strict.object$data$strict.tropical.character
####  #  names(strict.tropical.phylanova)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.lambda<-phylogroups.tree.strict.object$data$lambda.rpanda1
####  #  names(strict.tropical.phylanova.lambda)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.mu<-phylogroups.tree.strict.object$data$mu.rpanda1
####  #  names(strict.tropical.phylanova.mu)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.ndr<-phylogroups.tree.strict.object$data$ndr.rpanda1
####  #  names(strict.tropical.phylanova.ndr)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  strict.tropical.phylanova.turnover<-phylogroups.tree.strict.object$data$turnover.rpanda1
####  #  names(strict.tropical.phylanova.turnover)<-rownames(phylogroups.tree.strict.object$data)
####  #  
####  #  phylanova.strict.tropical.lambda<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.lambda,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.mu<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.mu,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.ndr<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.ndr,nsim = 100,posthoc = T)
####  #  phylanova.strict.tropical.turnover<-phylANOVA(tree=phylogroups.tree.strict,x=strict.tropical.phylanova,y=strict.tropical.phylanova.turnover,nsim = 100,posthoc = T)
####  #}
####  #summary(pgls.seedlambda)
####  #pgls.Median.Latitude.sigsq<-try(pgls(log10(lambda.rpanda1)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.Median.Latitude.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.ms05~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.log10.lambda',round(summary(pgls.Median.Latitude.lambda)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.Median.Latitude.lambda)=='gls'){
####    plot(phylogroups.data$lambda.ms05~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.lambda)
####    pgls.Median.Latitude.lambdaresults<-c('pgls.Median.Latitude.log10.lambda',round(summary(pgls.Median.Latitude.lambda)$tTable[2,1],8),round(summary(pgls.Median.Latitude.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.Median.Latitude.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.ms05~abs(phylogroups.tree.object$data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$coefficients[2,1],8),round(summary(pgls.Median.Latitude.ndr)$coefficients[2,4],8))
####    
####  }else if (class(pgls.Median.Latitude.ndr)=='gls'){
####    plot(phylogroups.data$ndr.ms05~abs(phylogroups.data$phylogroup.Median.Latitude),xlab='Median.Latitude',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(0,60))
####    abline(pgls.Median.Latitude.ndr)
####    pgls.Median.Latitude.ndrresults<-c('pgls.Median.Latitude.ndr',round(summary(pgls.Median.Latitude.ndr)$tTable[2,1],8),round(summary(pgls.Median.Latitude.ndr)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.ms05,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$coefficients[2,1],8),round(summary(pgls.prop.trop.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.lambda)=='gls'){
####    plot(phylogroups.data$lambda.ms05,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.lambda)
####    pgls.prop.trop.lambdaresults<-c('pgls.prop.trop.log10.lambda',round(summary(pgls.prop.trop.lambda)$tTable[2,1],8),round(summary(pgls.prop.trop.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.trop.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.ms05,phylogroups.tree.object$data$proportion.strict.tropical,ylab='prop.trop',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$coefficients[2,1],8),round(summary(pgls.prop.trop.ndr)$coefficients[2,4],8))
####  }else if (class(pgls.prop.trop.ndr)=='gls'){
####    plot(phylogroups.data$ndr.ms05,phylogroups.data$proportion.strict.tropical,ylab='prop.trop',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.trop.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.trop.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.trop.ndr)
####    pgls.prop.trop.ndrresults<-c('pgls.prop.trop.ndr',round(summary(pgls.prop.trop.ndr)$tTable[2,1],8),round(summary(pgls.prop.trop.ndr)$tTable[2,4],8))
####  }
####  
####  if(class(pgls.prop.temp.lambda)=='pgls'){
####    plot(phylogroups.tree.object$data$lambda.ms05,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$coefficients[2,1],8),round(summary(pgls.prop.temp.lambda)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.lambda)=='gls'){
####    plot(phylogroups.data$lambda.ms05,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.lambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.lambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.lambda)
####    pgls.prop.temp.lambdaresults<-c('pgls.prop.temp.log10.lambda',round(summary(pgls.prop.temp.lambda)$tTable[2,1],8),round(summary(pgls.prop.temp.lambda)$tTable[2,4],8))
####  }
####  if(class(pgls.prop.temp.ndr)=='pgls'){
####    plot(phylogroups.tree.object$data$ndr.ms05,phylogroups.tree.object$data$proportion.strict.temperate,ylab='prop.temp',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$coefficients[2,1],8),round(summary(pgls.prop.temp.ndr)$coefficients[2,4],8))
####  }else if (class(pgls.prop.temp.ndr)=='gls'){
####    plot(phylogroups.data$ndr.ms05,phylogroups.data$proportion.strict.temperate,ylab='prop.temp',xlab='ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.prop.temp.ndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.prop.temp.ndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8)
####    abline(pgls.prop.temp.ndr)
####    pgls.prop.temp.ndrresults<-c('pgls.prop.temp.ndr',round(summary(pgls.prop.temp.ndr)$tTable[2,1],8),round(summary(pgls.prop.temp.ndr)$tTable[2,4],8))
####  }
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #    
####  #   pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.log10.lambda',round(summary(pgls.strict.tropical.lambda)$coefficients[2,1],8),round(summary(pgls.strict.tropical.lambda)$coefficients[2,4],8))
####  #   pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$coefficients[2,1],8),round(summary(pgls.strict.tropical.mu)$coefficients[2,4],8))
####  #   pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$coefficients[2,1],8),round(summary(pgls.strict.tropical.ndr)$coefficients[2,4],8))
####  #   pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.log10.turnover',round(summary(pgls.strict.tropical.turnover)$coefficients[2,1],8),round(summary(pgls.strict.tropical.turnover)$coefficients[2,4],8))
####  #   
####  #   phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.log10.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####  #   phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####  #   phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####  #   phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.log10.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####  #   
####  # }
####  # 
####  
####  #   pgls.strict.tropical.lambda<-try(gls(strict.tropical.character~lambda.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.mu<-try(gls(strict.tropical.character~mu.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.ndr<-try(gls(strict.tropical.character~ndr.ms05, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   pgls.strict.tropical.turnover<-try(gls(strict.tropical.character~turnover.rpanda1, correlation = corPagel(1,phy = phylogroups.tree.strict,fixed=FALSE),data = phylogroups.data.strict, method = "ML"))
####  #   
####  # }
####  # 
####  # if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #   
####  #   pgls.strict.tropical.lambdaresults<-c('pgls.strict.tropical.log10.lambda',round(summary(pgls.strict.tropical.lambda)$tTable[2,1],8),round(summary(pgls.strict.tropical.lambda)$tTable[2,4],8))
####  #   pgls.strict.tropical.muresults<-c('pgls.strict.tropical.mu',round(summary(pgls.strict.tropical.mu)$tTable[2,1],8),round(summary(pgls.strict.tropical.mu)$tTable[2,4],8))
####  #   pgls.strict.tropical.ndrresults<-c('pgls.strict.tropical.ndr',round(summary(pgls.strict.tropical.ndr)$tTable[2,1],8),round(summary(pgls.strict.tropical.ndr)$tTable[2,4],8))
####  #   pgls.strict.tropical.turnoverresults<-c('pgls.strict.tropical.log10.turnover',round(summary(pgls.strict.tropical.turnover)$tTable[2,1],8),round(summary(pgls.strict.tropical.turnover)$tTable[2,4],8))
####  #   
####  #   phylanova.strict.tropical.lambda.results<-c('phylanova.strict.tropical.log10.lambda',round(phylanova.strict.tropical.lambda$F,3),round(phylanova.strict.tropical.lambda$Pf,3))
####  #   phylanova.strict.tropical.mu.results<-c('phylanova.strict.tropical.mu',round(phylanova.strict.tropical.mu$F,3),round(phylanova.strict.tropical.mu$Pf,3))
####  #   phylanova.strict.tropical.ndr.results<-c('phylanova.strict.tropical.ndr',round(phylanova.strict.tropical.ndr$F,3),round(phylanova.strict.tropical.ndr$Pf,3))
####  #   phylanova.strict.tropical.turnover.results<-c('phylanova.strict.tropical.log10.turnover',round(phylanova.strict.tropical.turnover$F,3),round(phylanova.strict.tropical.turnover$Pf,3))
####  #}
####  
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.muresults,pgls.Median.Latitude.ndrresults,pgls.Median.Latitude.turnoverresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.muresults,pgls.prop.trop.ndrresults,pgls.prop.trop.turnoverresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.muresults,pgls.prop.temp.ndrresults,pgls.prop.temp.turnoverresults,pgls.strict.tropical.lambdaresults,pgls.strict.tropical.muresults,pgls.strict.tropical.ndrresults,pgls.strict.tropical.turnoverresults,phylanova.strict.tropical.lambda.results,phylanova.strict.tropical.mu.results,phylanova.strict.tropical.ndr.results,phylanova.strict.tropical.turnover.results))
####  #}else{
####  results.df<-as.data.frame(rbind(pgls.Median.Latitude.lambdaresults,pgls.Median.Latitude.ndrresults,pgls.prop.trop.lambdaresults,pgls.prop.trop.ndrresults,pgls.prop.temp.lambdaresults,pgls.prop.temp.ndrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  if((min(table(phylogroups.tree.object2$data$binary.tropical))>1)&&(length(table(phylogroups.tree.object2$data$binary.tropical))>1)){
####    pgls.binary.trop.lambda<-phyloglm(binary.tropical~lambda.ms05, data = phylogroups.tree.object2$data, phy=phylogroups.tree2)
####
####    pgls.binary.trop.lambda.results<-c('pgls.binary.trop.log10.lambda',round(summary(pgls.binary.trop.lambda)$coefficients[2,1],8),round(summary(pgls.binary.trop.lambda)$coefficients[2,4],8))
####    results.binary.trop<-as.data.frame(rbind(pgls.binary.trop.lambda.results))
####    colnames(results.binary.trop)<-c('analysis','slope','pvalue')
####    results.df<-rbind(results.df,results.binary.trop)
####  }
####  
####  
####  
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('~/Dropbox/Work_in_progress/LDG_plants/clade_analyses/MS/MS_',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_',strict.tropical.threshold,'.tropthreshold_pgls_table_new_GBIFsampling_log.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  #if((min(table(df.phylogroups.strict$strict.tropical.character))>1)&&(length(table(df.phylogroups.strict$strict.tropical.character))>1)){
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='log10.lambda')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'lambda.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='mu')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'mu.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='ndr')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'ndr.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #  
####  #  plot(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==1,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='red',pch=16,xlab='abs.median.lat',ylab='log10.turnover')
####  #  points(abs(df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'phylogroup.Median.Latitude']),df.phylogroups.strict[df.phylogroups.strict$strict.tropical.character==0,'turnover.rpanda1'],xlim=c(0,50),ylim=c(0,10),col='blue',pch=16)
####  #}
####  dev.off()
####}
###wrapper to run clade analysis with RPANDA + pgls
#this creates the phylogroups (there's another script to run if phyogroups have previously been saved)
####run_clades_RPANDA_pgls<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  #save phylogroup object to Rsave file
####  save(phylogroups,file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_BM_phylogroup_select_size(x,table,sampling,mincladesize))
####  #phylogroups.RPANDA.nomuexp<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_BM_phylogroup_select_size(x,table,sampling,mincladesize))
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############correlation with lambda
####  pdf(paste('./output/plots/clade_analyses_new/RPANDAlambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_RPANDAlambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  ###########this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/RPANDAndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$sigsq)
####  
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #####this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_RPANDAndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####}
#this will load presaved phylogroups
####run_clades_RPANDA_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  phylogroups.species<-unlist(lapply(phylogroups$phylogroups,function(x)x$tip.label))
####  length(unique(phylogroups.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.RPANDA<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_BM_phylogroup_select_size(x,table,sampling,mincladesize))
####  #phylogroups.RPANDA.nomuexp<-lapply(phylogroups$phylogroups,function(x) run_RPANDA_BM_phylogroup_select_size(x,table,sampling,mincladesize))
####  df.phylogroups <- do.call("rbind", phylogroups.RPANDA)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.rpanda value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.rpanda1),]
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############correlation with lambda
####  pdf(paste('./output/plots/clade_analyses_new/RPANDAlambda_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_RPANDAlambda_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####  ###########this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/RPANDAndr_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.rpanda1,df.phylogroups$mu.rpanda1,df.phylogroups$ndr.rpanda1,df.phylogroups$sigsq)
####  
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.rpanda','mu.rpanda','ndr.rpanda','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #####this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.rpanda)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.rpanda)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.rpanda) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.rpanda)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.rpanda) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.rpanda)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_RPANDAndr_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####}

###wrapper to run clade analysis (with Magallon & Sanderson estimator) for congeneric species only
#this creates the phylogroups (there's another script to run if phyogroups have previously been saved)
####run_clades_MSlambda_congenerics_pgls<-function(minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  #save phylogroup object to Rsave file
####  #check for phylogroups that only contain congeneric species
####  phylogroups.congenerics.positions<-unlist(lapply(phylogroups$phylogroups,function(x){genera<-length(unique(sapply(strsplit(as.character(x$tip.label),'_'),function(x) x[1]))==1)}))
####  phylogroups.congenerics<-phylogroups$phylogroups[which(phylogroups.congenerics.positions==1)]
####  phylogroups.congenerics.species<-unlist(lapply(phylogroups.congenerics,function(x)x$tip.label))
####  length(unique(phylogroups.congenerics.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_BM_phylogroup_select_size(x,table))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############this is for lambda
####  pdf(paste('./output/plots/clade_analyses_new/MSlambda_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSlambda_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  ###############this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/MSndr_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.ms0,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms0','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSndr_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####}
#####this will load presaved phylogroups
####run_clades_MSlambda_congenerics_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  #phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  #save phylogroup object to Rsave file
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  #check for phylogroups that only contain congeneric species
####  phylogroups.congenerics.positions<-unlist(lapply(phylogroups$phylogroups,function(x){genera<-length(unique(sapply(strsplit(as.character(x$tip.label),'_'),function(x) x[1]))==1)}))
####  phylogroups.congenerics<-phylogroups$phylogroups[which(phylogroups.congenerics.positions==1)]
####  phylogroups.congenerics.species<-unlist(lapply(phylogroups.congenerics,function(x)x$tip.label))
####  length(unique(phylogroups.congenerics.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups$phylogroups,function(x) run_MS_BM_phylogroup_select_size(x,table,mincladesize = mincladesize,sampling=sampling))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  ###############this is for lambda
####  pdf(paste('./output/plots/clade_analyses_new/MSlambda_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSlambda_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  ###############this is for ndr
####  pdf(paste('./output/plots/clade_analyses_new/MSndr_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$phylogroup.Log10.Seed.Weight,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/clade_analyses_new/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSndr_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
#}
#####run_clades_MSlambda_z0_congenerics_pgls_savedphylogroups<-function(tree,minage,maxage,mincladesize,ncores,sampling,table){
####  colnames(table)[4]<-'tip'
####  colnames(table)[9]<-'genus.species.sampling.fraction'
####  cat('getting clades','\n')
####  #phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=mincladesize,ncores=ncores)
####  #save phylogroup object to Rsave file
####  load(file=paste('phylogroups_',minage,'_',maxage,'.Rsave',sep=''))
####  #check for phylogroups that only contain congeneric species
####  phylogroups.congenerics.positions<-unlist(lapply(phylogroups$phylogroups,function(x){genera<-length(unique(sapply(strsplit(as.character(x$tip.label),'_'),function(x) x[1]))==1)}))
####  phylogroups.congenerics<-phylogroups$phylogroups[which(phylogroups.congenerics.positions==1)]
####  phylogroups.congenerics.species<-unlist(lapply(phylogroups.congenerics,function(x)x$tip.label))
####  length(unique(phylogroups.congenerics.species))
####  cat('measuring diversification/trait evolution','\n')
####  phylogroups.MS<-lapply(phylogroups.congenerics,function(x) run_MS_EB_z0_phylogroup_select_size(x,table,mincladesize = mincladesize,sampling=sampling))
####  df.phylogroups <- do.call("rbind", phylogroups.MS)
####  df.phylogroups<-unique(df.phylogroups)
####  #remove rows with no lambda.ms value
####  df.phylogroups<-df.phylogroups[!is.na(df.phylogroups$lambda.ms0),]
####  
####  #count the number of seed size data points per phylogroup
####  phylogroup.species.seedsize<-matrix(NA,ncol=2,nrow=length(unique(df.phylogroups$phylogroup.name)))
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.species.seedsize[i,]<-c(unique(df.phylogroups$phylogroup.name)[i],as.numeric(length(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight[!is.na(df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$Log10.Seed.Weight)])))
####  }
####  phylogroup.species.seedsize<-as.data.frame(phylogroup.species.seedsize,stringsAsFactors = F)
####  colnames(phylogroup.species.seedsize)<-c('phylogroup.name','phylogroup.species.with.seed.data')
####  phylogroup.species.seedsize$phylogroup.species.with.seed.data<-as.numeric(phylogroup.species.seedsize$phylogroup.species.with.seed.data)
####  df.phylogroups<-merge(df.phylogroups,phylogroup.species.seedsize,all.x=TRUE)
####  
####  #get the mean seed size in each phylogroup
####  mean.seedsize.phylogroups<-aggregate(Log10.Seed.Weight~phylogroup.name,data=df.phylogroups,mean)
####  colnames(mean.seedsize.phylogroups)[2]<-'phylogroup.Log10.Seed.Weight'
####  df.phylogroups<-merge(df.phylogroups,mean.seedsize.phylogroups)
####  #getting vector of representatives in each phylogroup
####  #then keep just one tip per phylogroup
####  phylogroup.clades.species<-vector('character')
####  for (i in 1:length(unique(df.phylogroups$phylogroup.name))){
####    phylogroup.clades.species<-c(phylogroup.clades.species,df.phylogroups[df.phylogroups$phylogroup.name==unique(df.phylogroups$phylogroup.name)[i],]$New_Species[1])
####  }
####  phylogroups.tree<-drop.tip(tree,setdiff(tree$tip.label,phylogroup.clades.species))
####  names.phylogroups.tree<-data.frame(phylogroups.tree$tip.label)
####  colnames(names.phylogroups.tree)<-'Tip'
####  df.phylogroups<-merge(names.phylogroups.tree,df.phylogroups,by.x='Tip',by.y='New_Species',all.x=TRUE)
####  #another filter to check that size of clade is correct
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.size>mincladesize&df.phylogroups$phylogroup.species.with.seed.data>mincladesize,]
####  #add another filter to select well sampled clades
####  df.phylogroups<-df.phylogroups[df.phylogroups$phylogroup.sampling.fraction>sampling,]
####  #AIC model selection,AICcBM-AICcEB
####  
####  df.phylogroups$diff.AICc.BM.EB<-df.phylogroups$aicc.BM-df.phylogroups$aicc.EB
####  df.phylogroups$best.model.sigsq<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'sigsq.BM']
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'best.model.sigsq']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'sigsq.EB']
####  df.phylogroups$output.best.model<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'output.best.model']<-'BM'
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'output.best.model']<-'EB'
####  output.best.model<-as.data.frame(table(df.phylogroups$output.best.model))
####  df.phylogroups$z0<-NA
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'z0']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB<=2,'z0.BM']
####  df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'z0']<-df.phylogroups[df.phylogroups$diff.AICc.BM.EB>2,'z0.EB']
####  correlation.z0.SeedWeight<-cor.test(df.phylogroups$z0,df.phylogroups$phylogroup.Log10.Seed.Weight) 
####  correlation.df<-cbind(correlation.z0.SeedWeight$estimate,correlation.z0.SeedWeight$p.value)
####  colnames(correlation.df)<-c('estimate','p.value')
####  write.table(correlation.df,file=paste('./output/clade_analyses/MS_EB_z0/MS_EB_z0_lambda_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'cor_z0_SeedWeight.txt',sep=''),quote=F,sep='\t',row.names=F)
####  write.table(output.best.model,file=paste('./output/clade_analyses/MS_EB_z0/MS_EB_z0_lambda_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'best_fit_BMEB_model.txt',sep=''),quote=F,sep='\t',row.names=F)
####  
####  
####  ###############this is for lambda
####  pdf(paste('./output/plots/MSlambda_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$z0,df.phylogroups$lambda.ms0,df.phylogroups$lambda.ms05,df.phylogroups$lambda.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$best.model.sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','lambda.ms0','lambda.ms05','lambda.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedlambda<-try(pgls(log10(lambda.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedlambda)
####  pgls.seedsigsq<-try(pgls(log10(lambda.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedlambda)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$coefficients[2,1],8),round(summary(pgls.seedlambda)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedlambda)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedlambda<-gls(log10(lambda.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedlambda)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedlambda)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedlambda)
####    pgls.seedlambdaresults<-c('pgls.seedlambda',round(summary(pgls.seedlambda)$tTable[2,1],8),round(summary(pgls.seedlambda)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$lambda.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(lambda.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$lambda.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.lambda',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratelambdaresults<-c('pgls.seedratelambda',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedlambdaresults,pgls.seedratelambdaresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/plots/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSlambda_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  ###############this is for ndr
####  pdf(paste('./output/clade_analyses/MS_EB_z0/MSndr_congenerics_analysis_clades_',minage,'to',maxage,'_',sampling,'_',mincladesize,'size.pdf',sep=''),paper='a4')
####  par(mfrow=c(2,2))
####  hist(df.phylogroups$phylogroup.crowncapture,xlab='crown.capture.prob',main=paste('clades_',minage,'to',maxage,'_',mincladesize,'_crown.prob_',median(df.phylogroups$phylogroup.crowncapture),sep=''),xaxs='i',breaks=20,cex.main=.7)
####  hist(df.phylogroups$phylogroup.Log10.Seed.Weight,xlab='phylogroup.Log10.Seed.Weight')
####  phylogroups.tree<-drop.tip(phylogroups.tree,setdiff(phylogroups.tree$tip.label,df.phylogroups$Tip))
####  phylogroups.data<-data.frame(df.phylogroups$Tip,df.phylogroups$z0,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$ndr.ms05,df.phylogroups$ndr.ms09,df.phylogroups$best.model.sigsq)
####  colnames(phylogroups.data)<-c('Tip','phylogroup.Log10.Seed.Weight','ndr.ms05','ndr.ms09','ndr.ms05','ndr.ms09','sigsq')
####  phylogroups.tree$node.label<-NULL
####  #this is for caper
####  phylogroups.tree.object<-comparative.data(phy = phylogroups.tree,data=phylogroups.data,names.col='Tip',vcv=TRUE)
####  #for mean family seed weight
####  cat('running pgls','\n')
####  pgls.seedndr<-try(pgls(log10(ndr.ms05)~phylogroup.Log10.Seed.Weight, data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedndr)
####  pgls.seedsigsq<-try(pgls(log10(ndr.ms05)~log10(sigsq), data = phylogroups.tree.object, lambda='ML'))
####  #summary(pgls.seedsigsq)
####  ####try both optimisations (caper + nlme if caper gives an error)
####  if(class(pgls.seedndr)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~phylogroups.tree.object$data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$coefficients[2,1],8),round(summary(pgls.seedndr)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedndr)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedndr<-gls(log10(ndr.ms05) ~ phylogroup.Log10.Seed.Weight, correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~phylogroups.data$phylogroup.Log10.Seed.Weight,xlab='log10.Seed.Weight',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedndr)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedndr)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-4.1,1.1),ylim=c(-2,1))
####    abline(pgls.seedndr)
####    pgls.seedndrresults<-c('pgls.seedndr',round(summary(pgls.seedndr)$tTable[2,1],8),round(summary(pgls.seedndr)$tTable[2,4],8))
####  }
####  if(class(pgls.seedsigsq)!='try-error'){
####    plot(log10(phylogroups.tree.object$data$ndr.ms05)~log10(phylogroups.tree.object$data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$coefficients[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$coefficients[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$coefficients[2,1],8),round(summary(pgls.seedsigsq)$coefficients[2,4],8))
####  }
####  if(class(pgls.seedsigsq)=='try-error'){
####    row.names(phylogroups.data)<-phylogroups.data$Tip
####    pgls.seedsigsq<-gls(log10(ndr.ms05) ~ log10(sigsq), correlation = corPagel(1,phy = phylogroups.tree,fixed=FALSE),data = phylogroups.data, method = "ML")
####    plot(log10(phylogroups.data$ndr.ms05)~log10(phylogroups.data$sigsq),xlab='log10.sigsq',ylab='log10.ndr',main=paste(sum(df.phylogroups$phylogroup.size),' sp in ',nrow(phylogroups.data), ' clades-', '(',round(summary(pgls.seedsigsq)$tTable[2,1],3),' ;pvalue=',round(summary(pgls.seedsigsq)$tTable[2,4],3),')',sep=''),cex.main=.7,pch=16,cex=.8,xlim=c(-3.5,1),ylim=c(-2,1))
####    abline(pgls.seedsigsq)
####    pgls.seedratendrresults<-c('pgls.seedratendr',round(summary(pgls.seedsigsq)$tTable[2,1],8),round(summary(pgls.seedsigsq)$tTable[2,4],8))
####  }
####  
####  dev.off()
####  results.df<-as.data.frame(rbind(pgls.seedndrresults,pgls.seedratendrresults))
####  colnames(results.df)<-c('analysis','slope','pvalue')
####  rownames(results.df)<-NULL
####  write.table(results.df,file=paste('./output/clade_analyses/MS_EB_z0/',minage,'_',maxage,'_',sampling,'_',mincladesize,'size_pgls_MSndr_congenerics_results.txt',sep=''),quote=F,row.names=F,sep='\t')
####  
####}
####
prepare_tree_and_table_dataset<-function(){
  cat('preparing dataset to run clade based analyses','\n')
  #Qian and Jin created an updated version of Zanne tree - 'Appendix S3-v2' or 'QianTree.txt'
  #reference doi: 10.1093/jpe/rtv047
  Qian<-read.tree('./raw_data/QianTree.txt')
  Qian_TPL<-read.table('./raw_data/Qian_TPL.txt',header=T,sep='\t',stringsAsFactors = F,quote='')
  Qian_TPL$Merged<-paste(Qian_TPL$New.Genus,Qian_TPL$New.Species,sep='_')
  #run through taxonlookup
  Qian_TPLunique<-Qian_TPL[!duplicated(Qian_TPL$Merged),]
  Qian_Lookup<-lookup_table(as.character(Qian_TPLunique$Merged),by_species = TRUE)
  #get table with angiosperms only
  Qian_Lookup_Angios<-subset(Qian_Lookup,Qian_Lookup$group=='Angiosperms')
  Qian_Lookup_Angios$Fullspecies<-row.names(Qian_Lookup_Angios)
  row.names(Qian_Lookup_Angios)<-NULL
  Qian_Lookup_Angios_TPL<-merge(Qian_Lookup_Angios,Qian_TPLunique,by.x='Fullspecies',by.y='Merged')
  colnames(Qian_Lookup_Angios_TPL)[1]<-'New_Species'
  Qian_Lookup_Angios_TPL$Old_Species<-paste(Qian_Lookup_Angios_TPL$Genus,Qian_Lookup_Angios_TPL$Species,sep='_')
  Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[,c(c(1:6),24)]
  #drop tips in Qian tree (not angios, duplicated, hybrids)
  remove_tips<-setdiff(Qian$tip.label,Qian_Lookup_Angios_TPL$Old_Species)
  Qian_dropped<-drop.tip(Qian,remove_tips)
  #remove "_sp" tips
  remove_sp_tips<-Qian_Lookup_Angios_TPL[grep('_sp$',Qian_Lookup_Angios_TPL$Old_Species),]$Old_Species
  Qian_dropped<-drop.tip(Qian_dropped,remove_sp_tips)
  Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[-grep('_sp$',Qian_Lookup_Angios_TPL$New_Species),]
  Qiandf<-data.frame(Qian_dropped$tip.label)
  colnames(Qiandf)<-'Tipname'
  #replace tip names in Qian Tree with TPL+taxonlookup alternative
  Qian_Lookup_Angios_TPL<-merge(Qiandf,Qian_Lookup_Angios_TPL,by.x='Tipname',by.y='Old_Species')
  Qian_Lookup_Angios_TPL$Tipname<-as.character(Qian_Lookup_Angios_TPL$Tipname)
  
  species.counts<-as.data.frame(table(Qian_Lookup_Angios_TPL$genus),stringsAsFactors = F)
  colnames(species.counts)<-c('genus','number.of.species.tree')
  Qian_Lookup_Angios_TPL<-merge(Qian_Lookup_Angios_TPL,species.counts)
  Qian_Lookup_Angios_TPL$genus.sampling.fraction<-Qian_Lookup_Angios_TPL$number.of.species.tree/Qian_Lookup_Angios_TPL$number.of.species
  Qian_Lookup_Angios_TPL<-Qian_Lookup_Angios_TPL[match(Qian_dropped$tip.label,Qian_Lookup_Angios_TPL$Tipname),]
  Qian_dropped$tip.label<-Qian_Lookup_Angios_TPL$New_Species
  
  tree<-Qian_dropped
  #get total info from TaxonLookUp
  #plant_lookup_version_current()
  PlantLookup<-plant_lookup(include_counts = TRUE)
  #get number of genera in each family
  genus.counts<-as.data.frame(table(PlantLookup$family),stringsAsFactors = F)
  colnames(genus.counts)<-c('family','number.of.genera')
  
  #run TaxonLookup on taxa from the tree
  tree.table<-lookup_table(tree$tip.label)
  genus.counts.tree<-as.data.frame(table(tree.table$family),stringsAsFactors = F)
  colnames(genus.counts.tree)<-c('family','number.of.genera.tree')
  
  #merge all info
  genus.counts<-merge(genus.counts,genus.counts.tree)
  genus.counts$family.sampling.fraction<-genus.counts$number.of.genera.tree/genus.counts$number.of.genera
  #merge with all tree info
  Qian_Lookup_Angios_TPL<-merge(Qian_Lookup_Angios_TPL,genus.counts,all.x=TRUE)
  #add info on GBIF (not discarding species WITHOUT GBIF data)
  GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
  #GBIFdata<-read.table('~/Dropbox/Work_in_progress/LDG_plants/unbiased_reduced/GBIFdata.BAMM.unbiased_reduced_1_table.txt',sep='\t',stringsAsFactors = F,header = T)
  GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
  GBIFdata$strict.tropical<-NA
  #strict tropical 1 = abs(max latitude)<23.5 & abs(min latitude)<23.5
  GBIFdata[abs(GBIFdata$Max.Latitude)<=23.5&abs(GBIFdata$Min.Latitude)<=23.5&abs(GBIFdata$Median.Latitude)<=23.5,'strict.tropical']<-1
  #strict tropical 0 (strict temperate)
  GBIFdata[abs(GBIFdata$Max.Latitude)>23.5&abs(GBIFdata$Min.Latitude)>23.5&abs(GBIFdata$Median.Latitude)>23.5,'strict.tropical']<-0
  #2: not strict tropical or strict temperate species
  GBIFdata$strict.tropical[is.na(GBIFdata$strict.tropical)]<-2
  GBIFdata<-GBIFdata[-which(duplicated(GBIFdata$binomial)),]
  GBIFdata<-unique(GBIFdata)
  
  
  #this adds GBIF data to Qian dataset BUT doesn't discard species with no seed size data
  Qian_Lookup_Angios_TPL_GBIF<-merge(Qian_Lookup_Angios_TPL,GBIFdata,by.x='New_Species',by.y='binomial',all.x=TRUE)
  #adding info on taxonomic discordance (when there are more tips in the tree than accepted species according to taxonlookup)
  #it's only 555 species(1% of the total) // 79 genera (1% of the total)
  Qian_Lookup_Angios_TPL_GBIF$taxonomic.discordance<-NA
  Qian_Lookup_Angios_TPL_GBIF[Qian_Lookup_Angios_TPL_GBIF$number.of.species<Qian_Lookup_Angios_TPL_GBIF$number.of.species.tree,]$taxonomic.discordance<-1
  Qian_Lookup_Angios_TPL_GBIF[Qian_Lookup_Angios_TPL_GBIF$number.of.species>=Qian_Lookup_Angios_TPL_GBIF$number.of.species.tree,]$taxonomic.discordance<-0
  #to simplify I make the genus.species.sampling.fraction for these species = 1. But they can be removed as well for further analyses
  #Qian_Lookup_Angios_TPL_GBIF[Qian_Lookup_Angios_TPL_GBIF$taxonomic.discordance==1,]$genus.sampling.fraction<-1
  #remove species with taxonomic discordance
  Qian_Lookup_Angios_TPL_GBIF<-Qian_Lookup_Angios_TPL_GBIF[!Qian_Lookup_Angios_TPL_GBIF$taxonomic.discordance==1,]
  Qian_dropped<-drop.tip(Qian_dropped,setdiff(Qian_dropped$tip.label,Qian_Lookup_Angios_TPL_GBIF$New_Species))
  #write.tree(Qian_dropped,file='~/Dropbox/Work_in_progress/LDG_plants/Qian_dropped.tree')
  
  #tree<-read.tree('~/Dropbox/Work_in_progress/LDG_plants/Qian_dropped.tree')
  Qian_dropped<-tree
  output.list<-list(c(Qian_dropped,Qian_Lookup_Angios_TPL_GBIF))
}

