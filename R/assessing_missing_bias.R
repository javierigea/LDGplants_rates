library(phytools)
library(taxonlookup)
library(BAMMtools)
library(plyr)

setwd('~/Dropbox/Work_in_progress/LDG_plants')
#this returns the BAMM sampling fraction for the species in the dataset
get_samplingfractionsBAMM<-function(){
  GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
  GBIFdata$tropical<-0
  GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
  GBIFdata$strict.tropical<-NA
  #strict tropical 1 = abs(max latitude)<23.5 & abs(min latitude)<23.5  & median.latitude <23.5
  GBIFdata[abs(GBIFdata$Max.Latitude)<=23.5&abs(GBIFdata$Min.Latitude)<=23.5&abs(GBIFdata$Median.Latitude)<=23.5,'strict.tropical']<-1
  #strict tropical 0 (strict temperate)
  GBIFdata[abs(GBIFdata$Max.Latitude)>23.5&abs(GBIFdata$Min.Latitude)>23.5&abs(GBIFdata$Median.Latitude)>23.5,'strict.tropical']<-0
  #2: not strict tropical or strict temperate species
  GBIFdata$strict.tropical[is.na(GBIFdata$strict.tropical)]<-2
  GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
  BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
  #subset table to BAMM species
  GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
  #remove duplicates in GBIFdata (6 species with same info but two rows - one with present in garden, one with absent in garden)#
  GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
  GBIFdata.BAMM<-unique(GBIFdata.BAMM)
  GBIFdata.BAMM$latitudinal.band<-NA
  GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
  ratesdf<-getTipRates(BAMM.object)
  ratesdf.netdiv<-getTipRates(BAMM.object,returnNetDiv = T)
  rates.df<-as.data.frame(cbind(names(ratesdf$lambda.avg),ratesdf$lambda.avg,ratesdf$mu.avg,ratesdf.netdiv$netdiv.avg),stringsAsFactors = F)
  row.names(rates.df)<-NULL
  colnames(rates.df)<-c('species','lambda.avg','mu.avg','netdiv.avg')
  GBIFdata.BAMM.rates<-merge(GBIFdata.BAMM,rates.df,by.x='binomial',by.y='species')
  GBIFdata.BAMM.rates$lambda.avg<-as.numeric(GBIFdata.BAMM.rates$lambda.avg)
  GBIFdata.BAMM.rates$mu.avg<-as.numeric(GBIFdata.BAMM.rates$mu.avg)
  GBIFdata.BAMM.rates$netdiv.avg<-as.numeric(GBIFdata.BAMM.rates$netdiv.avg)
  
  #select clades of the tree where sampling is >0.6
  #calculate genus level sampling fraction
  plant_lookup_table<-plant_lookup(include_counts = T)
  GBIFdata.BAMM.rates<-merge(GBIFdata.BAMM.rates,plant_lookup_table,by.x='Genus.Name',by.y='genus',all.x=T)
  GBIFdata.BAMM.genus.counts<-as.data.frame(table(GBIFdata.BAMM.rates$Genus.Name))
  colnames(GBIFdata.BAMM.genus.counts)<-c('Genus.Name','number.of.species.tree')
  GBIFdata.BAMM.rates<-merge(GBIFdata.BAMM.rates,GBIFdata.BAMM.genus.counts)
  
  GBIFdata.BAMM.rates$taxonomic.discordance<-NA
  GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$number.of.accepted.species<GBIFdata.BAMM.rates$number.of.species.tree,]$taxonomic.discordance<-1
  GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$number.of.accepted.species>=GBIFdata.BAMM.rates$number.of.species.tree,]$taxonomic.discordance<-0
  
  GBIFdata.BAMM.rates<-GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$taxonomic.discordance==0,]
  GBIFdata.BAMM.rates$genus.sampling.fraction<-GBIFdata.BAMM.rates$number.of.species.tree/GBIFdata.BAMM.rates$number.of.species
  
  #add BAMM sampling fractions
  Zanne_clades_sampling_fraction.files<-list.files('~/Documents/Work/repositories/seed_size/output/BAMM_results/Zanne_clades/input_files/',pattern='_sampling.txt')
  Zanne_clades_sampling_fraction.tables<-lapply(Zanne_clades_sampling_fraction.files,function(x)read.table(paste('~/Documents/Work/repositories/seed_size/output/BAMM_results/Zanne_clades/input_files/',x,sep=''),header=F,sep='\t',stringsAsFactors=F,skip=1))
  Zanne_clades_sampling_fraction.df<-do.call('rbind',Zanne_clades_sampling_fraction.tables)
  Zanne_clades_sampling_fraction.df<-unique(Zanne_clades_sampling_fraction.df)
  colnames(Zanne_clades_sampling_fraction.df)<-c('binomial','family','BAMM.fraction')
  GBIFdata.BAMM.rates<-merge(GBIFdata.BAMM.rates,Zanne_clades_sampling_fraction.df)
  
  median.families.df<-as.data.frame(aggregate(Median.Latitude~family,data=GBIFdata.BAMM.rates,median),stringsAsFactors=F)
  colnames(median.families.df)[2]<-'Median.Latitude.family'
  GBIFdata.BAMM.rates<-merge(GBIFdata.BAMM.rates,median.families.df,all.x=T)
  return(GBIFdata.BAMM.rates)
}
