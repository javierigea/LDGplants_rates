library(BAMMtools)
library(coda)
library(ape)
library(mvtnorm)
library(lattice)
library(plyr)
#mcmcmout:mcmcoutfile
analyse_BAMM_convergence<-function(mcmcout,burnin){
  mcmcout.div <- read.csv(mcmcout, header=T)
  #see if convergence is ok
  plot(mcmcout.div$logLik ~ mcmcout.div$generation)
  #remove initial 10% (burnin)
  burnstart.div <- floor(burnin * nrow(mcmcout.div))
  postburn.div <- mcmcout.div[burnstart.div:nrow(mcmcout.div), ]
  plot(postburn.div$logLik ~ postburn.div$generation)
  plot(postburn.div$logLik ~ postburn.div$N_shifts)
  #check ESS (should be >200)
  print(paste('ESS Lik',effectiveSize(postburn.div$logLik),sep=': '))
  print(paste('ESS Nshifts',effectiveSize(postburn.div$N_shifts),sep=': '))
  #compute posterior probabilities of models
  post_probs.div <- table(postburn.div$N_shifts) / nrow(postburn.div)
  names(post_probs.div)
  #plot the posterior distribution of Nrateshifts
  plot(names(post_probs.div),post_probs.div)
}
#this applies a burnin (proportion of generations to discard) and recalculates the generation numbers in an event data file
#outputs results as *_burnin.txt
burnin_eventdata<-function(eventdata.file,burnin){
  ngens<-system(sprintf('tail -n1 %s | cut -d "," -f1',eventdata.file),intern=TRUE)
  ngens<-as.numeric(ngens)
  ngens.new<-round_any(ngens*burnin,10000)
  csv<-read.csv(eventdata.file,header=TRUE)
  csv<-csv[csv$generation>=ngens.new,]
  csv$generation<-csv$generation-ngens.new
  write.csv(csv,file=gsub(eventdata.file,pattern='.txt',replacement='_burnin.txt'),quote=F,row.names = FALSE)
}

#this function takes a subclade.tree file and its corresponding eventfile and paste the event data into the backbone event file
#age.root.Zanne is the age of the root in Zanne tree
add_subclade_to_backbone_BAMM<-function(subclade.treefile,subclade.eventfile,backbone.treefile,backbone.eventfile,newname,age.root.Zanne){
  subclade.tree<-read.tree(subclade.treefile)
  backbone.tree<-read.tree(backbone.treefile)
  clade.species.backbone<-intersect(backbone.tree$tip.label,subclade.tree$tip.label)
  #get age of clade root
  age.clade.root<-max(branching.times(subclade.tree))
  age.clade.root<-age.root.Zanne-age.clade.root
  #load subclade analysis (eventfile)
  subclade.event<-read.csv(subclade.eventfile,stringsAsFactors = F)
  #add age of clade root to absolute time
  subclade.event$abstime<-subclade.event$abstime+age.clade.root
  #load backbone analysis (eventfile)
  backbone.event<-read.csv(backbone.eventfile,stringsAsFactors = F)
  #check for shifts in branch leading to clade 
  grep.clade.species.backbone<-grep(clade.species.backbone,c(backbone.event$leftchild,backbone.event$rightchild))
  if(length(grep.clade.species.backbone)>0){
    #drop all shifts in the branch leading to the clade
    backbone.event<-backbone.event[-grep.clade.species.backbone,]
  }
  #this method of splicing things together will place a rate shift at the exact crown node of your clipped out subtrees. There is some possibility that this will cause numerical issues if the shift occurs at the exact time of an internal node
  # It would be safest to change the time of the rate shift at the  base of the subclade by an arbitrary but tiny amount, so that the shift moves a small bit down the parent branch.
  #add a small negative number (say, -0.001) to the root event for the analysis for subclade X, which will just push this shift a tiny bit closer to the root and make sure it doesn’t sit exactly on the node.
  correction<-sample(seq(from=0.0005,to=0.0015,by=0.0001),size=1)
  cat(correction,'\n')
  subclade.event[subclade.event$abstime==min(subclade.event$abstime),]$abstime<-subclade.event[subclade.event$abstime==min(subclade.event$abstime),]$abstime-correction
  merged.event<-rbind(backbone.event,subclade.event)
  merged.event<-merged.event[order(merged.event$generation),]
  #return(merged.event)
  write.csv(merged.event,file=paste('./',newname,'.txt',sep=''),quote=F,row.names=F)
} 



#take a list of event data files + burnin (the same for all) + a path to a folder
#apply burnin to all files
#get EventData + save it
eventdatafiles_burnin_merge_clades<-function(event_data.files,burnin,name,path){
  cat('burnin to event data files','\n')
  lapply(event_data.files,function(x) burnin_eventdata(eventdata.file = x,burnin=burnin))
  burnin.event_data.files<-gsub(event_data.files,pattern='.txt',replacement='_burnin.txt')
  tree<-read.tree(paste(path,'/BAMM_Species_tree_noseed.tree',sep=''))
  age.root.Zanne<-max(branching.times(tree))
  #setwd(paste('./Zanne_clades_BAMM/',name,'/results/',sep=''))
  cat('merging event data files','\n')
  add_subclade_to_backbone_BAMM(subclade.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_1.tree',sep=''),subclade.eventfile=burnin.event_data.files[1],backbone.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_7.tree',sep=''),backbone.eventfile=burnin.event_data.files[7],newname = paste(path,name,'/results/backbone1_',name,sep=''),age.root.Zanne)
  add_subclade_to_backbone_BAMM(subclade.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_2.tree',sep=''),subclade.eventfile=burnin.event_data.files[2],backbone.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_7.tree',sep=''),backbone.eventfile=paste(path,name,'/results/backbone1_',name,'.txt',sep=''),newname = paste(path,name,'/results/backbone12_',name,sep=''),age.root.Zanne)
  add_subclade_to_backbone_BAMM(subclade.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_3.tree',sep=''),subclade.eventfile=burnin.event_data.files[3],backbone.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_7.tree',sep=''),backbone.eventfile=paste(path,name,'/results/backbone12_',name,'.txt',sep=''),newname = paste(path,name,'/results/backbone123_',name,sep=''),age.root.Zanne)
  add_subclade_to_backbone_BAMM(subclade.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_4.tree',sep=''),subclade.eventfile=burnin.event_data.files[4],backbone.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_7.tree',sep=''),backbone.eventfile=paste(path,name,'/results/backbone123_',name,'.txt',sep=''),newname = paste(path,name,'/results/backbone1234_',name,sep=''),age.root.Zanne)
  add_subclade_to_backbone_BAMM(subclade.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_5.tree',sep=''),subclade.eventfile=burnin.event_data.files[5],backbone.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_7.tree',sep=''),backbone.eventfile=paste(path,name,'/results/backbone1234_',name,'.txt',sep=''),newname = paste(path,name,'/results/backbone12345_',name,sep=''),age.root.Zanne)
  add_subclade_to_backbone_BAMM(subclade.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_6.tree',sep=''),subclade.eventfile=burnin.event_data.files[6],backbone.treefile=paste(path,name,'/GBIFdata.BAMM.',name,'_clade_7.tree',sep=''),backbone.eventfile=paste(path,name,'/results/backbone12345_',name,'.txt',sep=''),newname = paste(path,name,'/results/backbone123456_',name,sep=''),age.root.Zanne)
  
  tree<-read.tree(paste(path,'/BAMM_Species_tree_noseed.tree',sep=''))
  tables<-lapply(list.files(path=paste(path,name,sep=''),pattern='_sampling.txt'),function(x)read.table(paste(path,name,'/',x,sep=''),header=F,skip=1,sep='\t',stringsAsFactors = F))
  tables<-do.call('rbind',tables)
  tree.table<-drop.tip(tree,setdiff(tree$tip.label,unique(tables$V1)))
  cat('getting event data object','\n')
  #getEventData for combined file
  subsampled.Div_edata <- getEventData(tree.table, eventdata=paste(path,name,'/results/backbone123456_',name,'.txt',sep=''), burnin = 0, nsamples = 1000,type='diversification',verbose=T)
  saveRDS(subsampled.Div_edata,file=paste(path,name,'/results/backbone123456_',name,'_eventData.RDS',sep=''))
  
}

