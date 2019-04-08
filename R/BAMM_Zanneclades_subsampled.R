library(ape)
library(BAMMtools)
library(taxonlookup)

BAMM_Zanneclades_subsampled<-function(replicate){
  cat(paste('BAMM_Zanneclades_subsampled for replicate',replicate),'\n')
  cat('reading BAMM clade trees','\n')
  #read clade sampling fractions and trees
  tree<-read.tree('./raw_data/BAMM_Species_tree_noseed.tree')
  #load cladetrees
  clade.trees.files<-list.files(path = './raw_data/BAMM_Zanneclades/',pattern='.tree')
  clade.trees<-lapply(clade.trees.files,function(x) read.tree(paste('./raw_data/BAMM_Zanneclades/',x,sep='')))
  #these are the species in each subclade that are also present in the backbone (i.e. the "anchoring species")
  clade.trees.anchoring.species<-unlist(lapply(clade.trees[c(1:6)],function(x)intersect(x$tip.label,clade.trees[[7]]$tip.label)))
  #load table of subsamples
  table<-read.table(paste('./results/subsampled_unbiased/input_tables/GBIFdata.BAMM.subsampled_',replicate,'_table.txt',sep=''),header=T,sep='\t',stringsAsFactors = F)
  #get missing anchoring species from table
  cat('overlapping BAMM clade trees with subsampled dataset','\n')
  missing.anchoring.species<-clade.trees.anchoring.species[!clade.trees.anchoring.species%in%table$binomial]
  if(length(missing.anchoring.species)>0){
    #add the species that are missing in subsample
    GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
    GBIFdata$tropical<-0
    GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
    GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
    BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
    #subset table to BAMM species
    GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
    #TO CHECK HERE: remove duplicates (6 species with same info but two rows - one with present in garden, one with absent in garden)####
    GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
    GBIFdata.BAMM<-unique(GBIFdata.BAMM)
    GBIFdata.BAMM$latitudinal.band<-NA
    GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
    GBIFdata.BAMM.missing.anchoring.species<-GBIFdata.BAMM[GBIFdata.BAMM$binomial%in%missing.anchoring.species,]
    if(nrow(GBIFdata.BAMM.missing.anchoring.species)!=length(missing.anchoring.species)){
      cat('****all missing anchoring species are not in GBIFdata.BAMM, double check','\n')
    }else{
      table<-rbind(table,GBIFdata.BAMM.missing.anchoring.species)
    }
  }
  #add family name to table
  cat('building sampling file for BAMM','\n')
  PlantLookup<-plant_lookup(include_counts = TRUE)
  PlantLookup.2<-aggregate(number.of.species~family, data=PlantLookup, sum)
  colnames(PlantLookup.2)[2]<-'number.of.species.family'
  PlantLookup<-merge(PlantLookup, PlantLookup.2, by.x = 'family', by.y = 'family', all.x = TRUE)
  
  table<-merge(table,PlantLookup,by.x='Genus.Name',by.y='genus')
  family.count.table<-as.data.frame(table(table$family),stringsAsFactors = F)
  colnames(family.count.table)<-c('family','number.of.species.tree')
  table<-merge(table,family.count.table,by.x='family',by.y='family')
  table$family.sampling.fraction<-table$number.of.species.tree/table$number.of.species.family
  table.family<-table[,c('family','family.sampling.fraction')]
  
  #get number of species in each subclade analysis for the 30k tree
  clade.sampling.tables.files<-list.files(path = './raw_data/BAMM_Zanneclades/',pattern='_sampling.txt')
  clade.sampling.tables<-lapply(clade.sampling.tables.files,function(x)read.table(paste('./raw_data/BAMM_Zanneclades/',x,sep=''),header=F,skip=1,sep='\t',stringsAsFactors=F))
  clade.sampling.total.fractions<-unlist(lapply(clade.sampling.tables.files,function(x)as.numeric(read.table(paste('./raw_data/BAMM_Zanneclades/',x,sep=''),header=F,nrows = 1,sep='\t',stringsAsFactors=F))))
  for(i in 1:length(clade.sampling.tables)){
    clade.sampling.total.fractions[i]<-round(nrow(clade.sampling.tables[[i]])/clade.sampling.total.fractions[[i]],0)
  }
  #intersect each subclade analysis in the 30k tree with the 10k table
  clade.sampling.tables.subset<-lapply(clade.sampling.tables,function(x)x[x$V1%in%table$binomial,])
  clade.sampling.tables.subset2<-lapply(clade.sampling.tables.subset,function(x){df<-unique(merge(x,table.family,by.x='V2',by.y='family'));return(df[,c('V1','V2','family.sampling.fraction')])})
  #drop tips in cladetrees
  clade.trees.subset<-list()
  for (i in 1:length(clade.trees)){
    clade.trees.subset[[i]]<-drop.tip(clade.trees[[i]],setdiff(clade.trees[[i]]$tip.label,clade.sampling.tables.subset2[[i]]$V1))
  }
  #this intersects the initial clade trees anchoring species with the 'new' ones; should be TRUE or something is wrong
  identical.test<-identical(clade.trees.anchoring.species,unlist(lapply(clade.trees.subset[c(1:6)],function(x)intersect(x$tip.label,clade.trees.subset[[7]]$tip.label))))
  cat('test for identity of initial anchoring species vs subsampled anchoring species', identical.test,'\n')
  #create a folder for each file
  #create outputs for BAMM: sampling, tree, priors
  dir.create(path = './results/Zanne_clades_BAMM/')
  dir.create(path = paste('./results/Zanne_clades_BAMM/subsampled_',replicate,'/',sep=''))
  dir.create(path = paste('./results/Zanne_clades_BAMM/subsampled_',replicate,'/','/results/',sep=''))
  for (i in 1:length(clade.sampling.tables.subset2)){
    write.table(nrow(clade.sampling.tables.subset2[[i]])/clade.sampling.total.fractions[i],file = paste('./results/Zanne_clades_BAMM/subsampled_',replicate,'/GBIFdata.BAMM.subsampled_',replicate,'_clade_',i,'_sampling.txt',sep=''),quote=F,row.names=F,sep='\t',col.names = F)
    write.table(clade.sampling.tables.subset2[[i]],file = paste('./results/Zanne_clades_BAMM/subsampled_',replicate,'/GBIFdata.BAMM.subsampled_',replicate,'_clade_',i,'_sampling.txt',sep=''),quote=F,row.names=F,sep='\t',col.names = F,append=T)
    write.tree(clade.trees.subset[[i]],file = paste('./results/Zanne_clades_BAMM/subsampled_',replicate,'/GBIFdata.BAMM.subsampled_',replicate,'_clade_',i,'.tree',sep=''))
    setBAMMpriors(clade.trees.subset[[i]],total.taxa =clade.sampling.total.fractions[i],outfile = paste('./results/Zanne_clades_BAMM/subsampled_',replicate,'/GBIFdata.BAMM.subsampled_',replicate,'_clade_',i,'_priors.txt',sep='') )
  }
}


BAMM_Zanneclades_subsampled_extremebias<-function(replicate){
  cat(paste('BAMM_Zanneclades_subsampled for replicate',replicate),'\n')
  cat('reading BAMM clade trees','\n')
  #read clade sampling fractions and trees
  tree<-read.tree('./raw_data/BAMM_Species_tree_noseed.tree')
  #load cladetrees
  clade.trees.files<-list.files(path = './raw_data/BAMM_Zanneclades/',pattern='.tree')
  clade.trees<-lapply(clade.trees.files,function(x) read.tree(paste('./raw_data/BAMM_Zanneclades/',x,sep='')))
  #these are the species in each subclade that are also present in the backbone (i.e. the "anchoring species")
  clade.trees.anchoring.species<-unlist(lapply(clade.trees[c(1:6)],function(x)intersect(x$tip.label,clade.trees[[7]]$tip.label)))
  #load table of subsamples
  table<-read.table(paste('./results/subsampled_extremebias/input_tables/GBIFdata.BAMM.subsampled_extremebias_',replicate,'_table.txt',sep=''),header=T,sep='\t',stringsAsFactors = F)
  #get missing anchoring species from table
  cat('overlapping BAMM clade trees with subsampled dataset','\n')
  missing.anchoring.species<-clade.trees.anchoring.species[!clade.trees.anchoring.species%in%table$binomial]
  if(length(missing.anchoring.species)>0){
    #add the species that are missing in subsample
    GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
    GBIFdata$tropical<-0
    GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
    GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
    BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
    #subset table to BAMM species
    GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
    #TO CHECK HERE: remove duplicates (6 species with same info but two rows - one with present in garden, one with absent in garden)####
    GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
    GBIFdata.BAMM<-unique(GBIFdata.BAMM)
    GBIFdata.BAMM$latitudinal.band<-NA
    GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
    GBIFdata.BAMM.missing.anchoring.species<-GBIFdata.BAMM[GBIFdata.BAMM$binomial%in%missing.anchoring.species,]
    if(nrow(GBIFdata.BAMM.missing.anchoring.species)!=length(missing.anchoring.species)){
      cat('****all missing anchoring species are not in GBIFdata.BAMM, double check','\n')
    }else{
      table<-rbind(table,GBIFdata.BAMM.missing.anchoring.species)
    }
  }
  #add family name to table
  cat('building sampling file for BAMM','\n')
  PlantLookup<-plant_lookup(include_counts = TRUE)
  PlantLookup.2<-aggregate(number.of.accepted.species~family, data=PlantLookup, sum)
  colnames(PlantLookup.2)[2]<-'number.of.species.family'
  PlantLookup<-merge(PlantLookup, PlantLookup.2, by.x = 'family', by.y = 'family', all.x = TRUE)
  
  table<-merge(table,PlantLookup,by.x='Genus.Name',by.y='genus')
  family.count.table<-as.data.frame(table(table$family),stringsAsFactors = F)
  colnames(family.count.table)<-c('family','number.of.species.tree')
  table<-merge(table,family.count.table,by.x='family',by.y='family')
  table$family.sampling.fraction<-table$number.of.species.tree/table$number.of.species.family
  table.family<-table[,c('family','family.sampling.fraction')]
  
  #get number of species in each subclade analysis for the 30k tree
  clade.sampling.tables.files<-list.files(path = './raw_data/BAMM_Zanneclades/',pattern='_sampling.txt')
  clade.sampling.tables<-lapply(clade.sampling.tables.files,function(x)read.table(paste('./raw_data/BAMM_Zanneclades/',x,sep=''),header=F,skip=1,sep='\t',stringsAsFactors=F))
  clade.sampling.total.fractions<-unlist(lapply(clade.sampling.tables.files,function(x)as.numeric(read.table(paste('./raw_data/BAMM_Zanneclades/',x,sep=''),header=F,nrows = 1,sep='\t',stringsAsFactors=F))))
  for(i in 1:length(clade.sampling.tables)){
    clade.sampling.total.fractions[i]<-round(nrow(clade.sampling.tables[[i]])/clade.sampling.total.fractions[[i]],0)
  }
  #intersect each subclade analysis in the 30k tree with the 10k table
  clade.sampling.tables.subset<-lapply(clade.sampling.tables,function(x)x[x$V1%in%table$binomial,])
  clade.sampling.tables.subset2<-lapply(clade.sampling.tables.subset,function(x){df<-unique(merge(x,table.family,by.x='V2',by.y='family'));return(df[,c('V1','V2','family.sampling.fraction')])})
  #drop tips in cladetrees
  clade.trees.subset<-list()
  for (i in 1:length(clade.trees)){
    clade.trees.subset[[i]]<-drop.tip(clade.trees[[i]],setdiff(clade.trees[[i]]$tip.label,clade.sampling.tables.subset2[[i]]$V1))
  }
  #this intersects the initial clade trees anchoring species with the 'new' ones; should be TRUE or something is wrong
  identical.test<-identical(clade.trees.anchoring.species,unlist(lapply(clade.trees.subset[c(1:6)],function(x)intersect(x$tip.label,clade.trees.subset[[7]]$tip.label))))
  cat('test for identity of initial anchoring species vs subsampled anchoring species', identical.test,'\n')
  #create a folder for each file
  #create outputs for BAMM: sampling, tree, priors
  dir.create(path = './results/Zanne_clades_BAMM_extremebias/')
  dir.create(path = paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_',replicate,'/',sep=''))
  dir.create(path = paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_',replicate,'/','/results/',sep=''))
  for (i in 1:length(clade.sampling.tables.subset2)){
    write.table(nrow(clade.sampling.tables.subset2[[i]])/clade.sampling.total.fractions[i],file = paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_',replicate,'/GBIFdata.BAMM.subsampled_extremebias_',replicate,'_clade_',i,'_sampling.txt',sep=''),quote=F,row.names=F,sep='\t',col.names = F)
    write.table(clade.sampling.tables.subset2[[i]],file = paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_',replicate,'/GBIFdata.BAMM.subsampled_extremebias_',replicate,'_clade_',i,'_sampling.txt',sep=''),quote=F,row.names=F,sep='\t',col.names = F,append=T)
    write.tree(clade.trees.subset[[i]],file = paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_',replicate,'/GBIFdata.BAMM.subsampled_extremebias_',replicate,'_clade_',i,'.tree',sep=''))
    setBAMMpriors(clade.trees.subset[[i]],total.taxa =clade.sampling.total.fractions[i],outfile = paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_',replicate,'/GBIFdata.BAMM.subsampled_extremebias_',replicate,'_clade_',i,'_priors.txt',sep='') )
  }
}


BAMM_Zanneclades_subsampled_small<-function(replicate){
  cat(paste('BAMM_Zanneclades_subsampled for replicate',replicate),'\n')
  cat('reading BAMM clade trees','\n')
  #read clade sampling fractions and trees
  tree<-read.tree('./raw_data/BAMM_Species_tree_noseed.tree')
  #load cladetrees
  clade.trees.files<-list.files(path = './raw_data/BAMM_Zanneclades/',pattern='.tree')
  clade.trees<-lapply(clade.trees.files,function(x) read.tree(paste('./raw_data/BAMM_Zanneclades/',x,sep='')))
  #these are the species in each subclade that are also present in the backbone (i.e. the "anchoring species")
  clade.trees.anchoring.species<-unlist(lapply(clade.trees[c(1:6)],function(x)intersect(x$tip.label,clade.trees[[7]]$tip.label)))
  #load table of subsamples
  table<-read.table(paste('./results/subsampled_small/input_tables/GBIFdata.BAMM.subsampled_small_',replicate,'_table.txt',sep=''),header=T,sep='\t',stringsAsFactors = F)
  #get missing anchoring species from table
  cat('overlapping BAMM clade trees with subsampled dataset','\n')
  missing.anchoring.species<-clade.trees.anchoring.species[!clade.trees.anchoring.species%in%table$binomial]
  if(length(missing.anchoring.species)>0){
    #add the species that are missing in subsample
    GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
    GBIFdata$tropical<-0
    GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
    GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
    BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
    #subset table to BAMM species
    GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
    #TO CHECK HERE: remove duplicates (6 species with same info but two rows - one with present in garden, one with absent in garden)####
    GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
    GBIFdata.BAMM<-unique(GBIFdata.BAMM)
    GBIFdata.BAMM$latitudinal.band<-NA
    GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
    GBIFdata.BAMM.missing.anchoring.species<-GBIFdata.BAMM[GBIFdata.BAMM$binomial%in%missing.anchoring.species,]
    if(nrow(GBIFdata.BAMM.missing.anchoring.species)!=length(missing.anchoring.species)){
      cat('****all missing anchoring species are not in GBIFdata.BAMM, double check','\n')
    }else{
      table<-rbind(table,GBIFdata.BAMM.missing.anchoring.species)
    }
  }
  #add family name to table
  cat('building sampling file for BAMM','\n')
  PlantLookup<-plant_lookup(include_counts = TRUE)
  PlantLookup.2<-aggregate(number.of.accepted.species~family, data=PlantLookup, sum)
  colnames(PlantLookup.2)[2]<-'number.of.species.family'
  PlantLookup<-merge(PlantLookup, PlantLookup.2, by.x = 'family', by.y = 'family', all.x = TRUE)
  
  table<-merge(table,PlantLookup,by.x='Genus.Name',by.y='genus')
  family.count.table<-as.data.frame(table(table$family),stringsAsFactors = F)
  colnames(family.count.table)<-c('family','number.of.species.tree')
  table<-merge(table,family.count.table,by.x='family',by.y='family')
  table$family.sampling.fraction<-table$number.of.species.tree/table$number.of.species.family
  table.family<-table[,c('family','family.sampling.fraction')]
  
  #get number of species in each subclade analysis for the 30k tree
  clade.sampling.tables.files<-list.files(path = './raw_data/BAMM_Zanneclades/',pattern='_sampling.txt')
  clade.sampling.tables<-lapply(clade.sampling.tables.files,function(x)read.table(paste('./raw_data/BAMM_Zanneclades/',x,sep=''),header=F,skip=1,sep='\t',stringsAsFactors=F))
  clade.sampling.total.fractions<-unlist(lapply(clade.sampling.tables.files,function(x)as.numeric(read.table(paste('./raw_data/BAMM_Zanneclades/',x,sep=''),header=F,nrows = 1,sep='\t',stringsAsFactors=F))))
  for(i in 1:length(clade.sampling.tables)){
    clade.sampling.total.fractions[i]<-round(nrow(clade.sampling.tables[[i]])/clade.sampling.total.fractions[[i]],0)
  }
  #intersect each subclade analysis in the 30k tree with the 10k table
  clade.sampling.tables.subset<-lapply(clade.sampling.tables,function(x)x[x$V1%in%table$binomial,])
  clade.sampling.tables.subset2<-lapply(clade.sampling.tables.subset,function(x){df<-unique(merge(x,table.family,by.x='V2',by.y='family'));return(df[,c('V1','V2','family.sampling.fraction')])})
  #drop tips in cladetrees
  clade.trees.subset<-list()
  for (i in 1:length(clade.trees)){
    clade.trees.subset[[i]]<-drop.tip(clade.trees[[i]],setdiff(clade.trees[[i]]$tip.label,clade.sampling.tables.subset2[[i]]$V1))
  }
  #this intersects the initial clade trees anchoring species with the 'new' ones; should be TRUE or something is wrong
  identical.test<-identical(clade.trees.anchoring.species,unlist(lapply(clade.trees.subset[c(1:6)],function(x)intersect(x$tip.label,clade.trees.subset[[7]]$tip.label))))
  cat('test for identity of initial anchoring species vs subsampled anchoring species', identical.test,'\n')
  #create a folder for each file
  #create outputs for BAMM: sampling, tree, priors
  dir.create(path = './results/Zanne_clades_BAMM_small/')
  dir.create(path = paste('./results/Zanne_clades_BAMM_small/subsampled_small_',replicate,'/',sep=''))
  dir.create(path = paste('./results/Zanne_clades_BAMM_small/subsampled_small_',replicate,'/','/results/',sep=''))
  for (i in 1:length(clade.sampling.tables.subset2)){
    write.table(nrow(clade.sampling.tables.subset2[[i]])/clade.sampling.total.fractions[i],file = paste('./results/Zanne_clades_BAMM_small/subsampled_small_',replicate,'/GBIFdata.BAMM.subsampled_small_',replicate,'_clade_',i,'_sampling.txt',sep=''),quote=F,row.names=F,sep='\t',col.names = F)
    write.table(clade.sampling.tables.subset2[[i]],file = paste('./results/Zanne_clades_BAMM_small/subsampled_small_',replicate,'/GBIFdata.BAMM.subsampled_small_',replicate,'_clade_',i,'_sampling.txt',sep=''),quote=F,row.names=F,sep='\t',col.names = F,append=T)
    write.tree(clade.trees.subset[[i]],file = paste('./results/Zanne_clades_BAMM_small/subsampled_small_',replicate,'/GBIFdata.BAMM.subsampled_small_',replicate,'_clade_',i,'.tree',sep=''))
    setBAMMpriors(clade.trees.subset[[i]],total.taxa =clade.sampling.total.fractions[i],outfile = paste('./results/Zanne_clades_BAMM_small/subsampled_small_',replicate,'/GBIFdata.BAMM.subsampled_small_',replicate,'_clade_',i,'_priors.txt',sep='') )
  }
}
