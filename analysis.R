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
#change wd if neeeded
setwd('~/Documents/Work/repositories/LDGplants_rates/')
#setwd('~/Dropbox/Work_in_progress/repositories/LDG_plants/')
####create folder structure
dir.create('./plots/')
dir.create('./results/')
####raw data should have: GBIFdatasummary.csv (from Mounce et al Nat Plants), all_clades_50million_1000samples_latitudedata.RDS (eventdata from BAMM)
#QianTree.txt, BAMM_Species_tree_noseed.tree
####./raw_data/BAMM_Zanneclades/ should have clade_1_sampling.txt clade_1.tree (...) clade_7_sampling.txt, clade_7.tree
####./raw_data/BAMM_Zanneclades/control_files/ has all control_files.txt
####./R/run_fisse/traitDependent_functions.R should be modified to include FISSE.binary.median



####1) run STRAPP with full biased dataset (binary, continuous, latband, absolute lat band)####
####prepare dataset for STRAPP full####
dir.create('./results/full_biased/')
source('./R/strapp_functions.R')
GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
#tropical defined as abs(median latitude <23.5)
#temperate defined as abs(median latitude >23.5)
GBIFdata$tropical<-0
GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
#subset table to BAMM species
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
#remove duplicates in GBIFdata (6 species with same info but two rows - one with present in garden, one with absent in garden)#
GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
GBIFdata.BAMM<-unique(GBIFdata.BAMM)
GBIFdata.BAMM$latitudinal.band<-NA
GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
####strapp runs for full dataset####
tropical.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =GBIFdata.BAMM,mode = 'binary.tropical')
write.table(tropical.strapp,file='./results/full_biased/tropical.strapp.fulldataset.txt',sep='\t',quote=F,row.names=F)
medianlatitude.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =GBIFdata.BAMM,mode = 'continuous.medianlatitude')
write.table(medianlatitude.strapp,file='./results/full_biased/medianlatitude.strapp.fulldataset.txt',sep='\t',quote=F,row.names=F)
abs.medianlatitude.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =GBIFdata.BAMM,mode = 'continuous.abs.medianlatitude')
write.table(abs.medianlatitude.strapp,file='./results/full_biased/abs.medianlatitude.strapp.fulldataset.txt',sep='\t',quote=F,row.names=F)
latitudinalband.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =GBIFdata.BAMM,mode = 'latitudinalband')
write.table(latitudinalband.strapp,file='./results/full_biased/latitudinalband.strapp.fulldataset.txt',sep='\t',quote=F,row.names=F)
abs.latitudinalband.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =GBIFdata.BAMM,mode = 'abs.latitudinalband')
write.table(abs.latitudinalband.strapp,file='./results/full_biased/abs.latitudinalband.strapp.fulldataset.txt',sep='\t',quote=F,row.names=F)

####prepare dataset for strict tropical STRAPP full####
dir.create('./results/full_biased_strict/')
source('./R/strapp_functions.R')
GBIFdata<-read.csv('./GBIFdatasummary.csv')
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
####strapp runs with full dataset - strict.tropical as a 3-state character####
strict.tropical.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =GBIFdata.BAMM,mode = 'strict.tropical')
write.table(strict.tropical.strapp,file='./results/full_biased_strict/strict.tropical.3states.strapp.fulldataset.txt',sep='\t',quote=F,row.names=F)
####strapp runs with full dataset -  strict.tropical, only strict.tropical vs strict.temperate species####
GBIFdata.BAMM.strict.tropical<-GBIFdata.BAMM[!GBIFdata.BAMM$strict.tropical==2,]
#need to subset BAMMobject to strict trops and strict temps
BAMM.object.strict.tropical<-subtreeBAMM(ephy = BAMM.object,tips = GBIFdata.BAMM.strict.tropical$binomial)
saveRDS(BAMM.object.strict.tropical,'./results/full_biased_strict/all_clades_50million_1000samples_latitudedata_stricttroptemp.RDS')
strict.tropical.binary.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.strict.tropical,trait.table =GBIFdata.BAMM.strict.tropical,mode = 'strict.binary.tropical')
write.table(strict.tropical.binary.strapp,file='./results/full_biased_strict/strict.tropical.binary.strapp.fulldataset.txt',sep='\t',quote=F,row.names=F)

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
colnames(GBIFdata.BAMM)
GBIFdata.BAMM<-GBIFdata.BAMM[,c('binomial','Median.Latitude')]
write.table(GBIFdata.BAMM,file='./results/GBIFdata.BAMM.table.txt',sep='\t',quote=F,row.names = F)

####Figure 1 Plots####
source('./R/plot_figs.R')
#get plots for Figure 1
pdf('./plots/Fig1_map_latitudinalbands_log.pdf',paper='a4r')
plot_Fig1_map_bands()
dev.off()
plot_Fig1_boxplot_troptempbinary()
plot_Fig1_quantilerate_vs_lambdarate()
plot_Fig1_strapp_samples_from_file(correlationfile = './results/strapp_lambda_absMedianLatitude_allnew.txt',parameter = 'lambda',tablefile='./results/GBIFdata.BAMM.table.txt')
plot_Fig1_strapp_density()

####2) run STRAPP with REDUCED full biased dataset (binary, continuous, latband, absolute lat band)####
dir.create('./results/full_biased_reduced/')
source('./R/strapp_functions.R')
GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
#tropical defined as abs(median latitude <23.5)
#temperate defined as abs(median latitude >23.5)
GBIFdata$tropical<-0
GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
#subset table to BAMM species
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
#remove duplicates in GBIFdata (6 species with same info but two rows - one with present in garden, one with absent in garden)#
GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
GBIFdata.BAMM<-unique(GBIFdata.BAMM)
GBIFdata.BAMM$latitudinal.band<-NA
GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
#proportions of data points per latitudinal band
sort(round(table(GBIFdata.BAMM$latitudinal.band)/nrow(GBIFdata.BAMM),4))
####discard bands: -60,-50,60,70,80,90 (<4% of total points)####
GBIFdata.BAMM.reduced<-GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band!='-60'&GBIFdata.BAMM$latitudinal.band!='-50'&GBIFdata.BAMM$latitudinal.band!='60'&GBIFdata.BAMM$latitudinal.band!='70'&GBIFdata.BAMM$latitudinal.band!='80'&GBIFdata.BAMM$latitudinal.band!='90',]
#create new BAMM object and save
BAMM.object.lat_reduced<-subtreeBAMM(BAMM.object,tips = GBIFdata.BAMM.reduced$binomial)
saveRDS(BAMM.object.lat_reduced,'./results/full_biased_reduced/all_clades_50million_1000samples_latitudedata_reduced.RDS')
####strapp runs for full reduced dataset####
tropical.reduced.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.lat_reduced,trait.table =GBIFdata.BAMM.reduced,mode = 'binary.tropical')
write.table(tropical.reduced.strapp,file='./results/full_biased_reduced/tropical.strapp.fulldatasetreduced.txt',sep='\t',quote=F,row.names=F)
medianlatitude.reduced.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.lat_reduced,trait.table =GBIFdata.BAMM.reduced,mode = 'continuous.medianlatitude')
write.table(medianlatitude.reduced.strapp,file='./results/full_biased_reduced/medianlatitude.strapp.fulldatasetreduced.txt',sep='\t',quote=F,row.names=F)
abs.medianlatitude.reduced.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.lat_reduced,trait.table =GBIFdata.BAMM.reduced,mode = 'continuous.abs.medianlatitude')
write.table(abs.medianlatitude.reduced.strapp,file='./results/full_biased_reduced/abs.medianlatitude.strapp.fulldatasetreduced.txt',sep='\t',quote=F,row.names=F)
latitudinalband.reduced.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.lat_reduced,trait.table =GBIFdata.BAMM.reduced,mode = 'latitudinalband')
write.table(latitudinalband.reduced.strapp,file='./results/full_biased_reduced/latitudinalband.strapp.fulldatasetreduced.txt',sep='\t',quote=F,row.names=F)
abs.latitudinalband.reduced.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.lat_reduced,trait.table =GBIFdata.BAMM.reduced,mode = 'abs.latitudinalband')
write.table(abs.latitudinalband.reduced.strapp,file='./results/full_biased_reduced/abs.latitudinalband.strapp.fulldatasetreduced.txt',sep='\t',quote=F,row.names=F)

####TO ADD: PLOT FIG S1 HERE####


####TO ADD: PLOT FIG S2 HERE####

####3) run STRAPP with full biased dataset with >4 GBIF data points per species####
#prepare dataset for GBIF 4occ data
dir.create('./results/full_biased_GBIFocc/')
GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
GBIFdata.filter<-GBIFdata[GBIFdata$Georeferenced.GBIF.Records>4,]
write.csv(GBIFdata.filter,file='./results/full_biased_GBIFocc/GBIFdatasummary_filterGBIFocc.csv')
#create BAMM object for GBIF >4occ data
GBIFdata<-read.csv('./results/full_biased_GBIFocc/GBIFdatasummary_filterGBIFocc.csv')
GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
BAMM.object.GBIF4occ<-subtreeBAMM(ephy = BAMM.object,tips = GBIFdata.BAMM$binomial)
saveRDS(BAMM.object.GBIF4occ,file='./results/full_biased_GBIFocc/all_clades_50million_1000samples_latitudedata_GBIF4occ.RDS')

####run STRAPP with full biased dataset with species with >4 GBIF occ (binary, continuous, latband, absolute lat band)####
source('./R/strapp_functions.R')
GBIFdata<-read.csv('./results/full_biased_GBIFocc/GBIFdatasummary_filterGBIFocc.csv')
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
BAMM.object.GBIFocc<-readRDS(file='./results/full_biased_GBIFocc/all_clades_50million_1000samples_latitudedata_GBIF4occ.RDS')
#subset table to BAMM species
GBIFdata.BAMM.GBIFocc<-GBIFdata[GBIFdata$binomial%in%BAMM.object.GBIFocc$tip.label,]
# remove duplicates (6 species with same info but two rows - one with present in garden, one with absent in garden)#
GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
GBIFdata.BAMM<-unique(GBIFdata.BAMM)
GBIFdata.BAMM$latitudinal.band<-NA
GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
####strapp runs with GBIFOcc>4 full dataset####
tropical.GBIFocc.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =BAMM.object.GBIFocc,mode = 'binary.tropical')
write.table(tropical.strapp,file='./results/full_biased_GBIFocc/tropical.strapp.fulldataset.GBIFOcc.txt',sep='\t',quote=F,row.names=F)
medianlatitude.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =BAMM.object.GBIFocc,mode = 'continuous.medianlatitude')
write.table(medianlatitude.strapp,file='./results/full_biased_GBIFocc/medianlatitude.strapp.fulldataset.GBIFOcc.txt',sep='\t',quote=F,row.names=F)
abs.medianlatitude.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =BAMM.object.GBIFocc,mode = 'continuous.abs.medianlatitude')
write.table(abs.medianlatitude.strapp,file='./results/full_biased_GBIFocc/abs.medianlatitude.strapp.fulldataset.GBIFOcc.txt',sep='\t',quote=F,row.names=F)
latitudinalband.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =BAMM.object.GBIFocc,mode = 'latitudinalband')
write.table(latitudinalband.strapp,file='./results/full_biased_GBIFocc/latitudinalband.strapp.fulldataset.GBIFOcc.txt',sep='\t',quote=F,row.names=F)
abs.latitudinalband.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =BAMM.object.GBIFocc,mode = 'abs.latitudinalband')
write.table(abs.latitudinalband.strapp,file='./results/full_biased_GBIFocc/abs.latitudinalband.strapp.fulldataset.GBIFOcc.txt',sep='\t',quote=F,row.names=F)

####TO ADD HERE: PLOT FIG HERE, SIMILAR TO FIG1, FIGS1####

####4) geographic bias plots####
####geographic bias plot for full dataset####
source('./R/assessing_bias_dataset.R')
#this plots the proportion of trop and temp species in the GBIF dataset, and the proportion of species in each lat band
pdf('./plots/full_dataset_BAMMdataset_biases_new2.pdf',width=10,height=7)
par(mfrow=c(2,2))
par(pty='s')
plot_GBIFdata_proportions(GBIFdata.file='./raw_data/GBIFdatasummary.csv')
#this plots the proportion of trop and temp species in the GBIF BAMM dataset, and the proportion of species in each lat band
plot_GBIFdataBAMM_proportions(GBIFdata.file='./raw_data/GBIFdatasummary.csv',BAMM.object.path='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
#there's an obvious bias, with many more species in temperate regions being sequenced
dev.off()

####geographic bias plot for full reduced dataset####
source('./R/assessing_bias_dataset.R')
#this plots the proportion of trop and temp species in the GBIF dataset, and the proportion of species in each lat band
GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
GBIFdata$latitudinal.band<-NA
GBIFdata$latitudinal.band<-apply(GBIFdata,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
GBIFdata<-GBIFdata[GBIFdata$latitudinal.band!='-60'&GBIFdata$latitudinal.band!='-50'&GBIFdata$latitudinal.band!='60'&GBIFdata$latitudinal.band!='70'&GBIFdata$latitudinal.band!='80'&GBIFdata$latitudinal.band!='90',]
write.csv(GBIFdata,file='./results/full_biased_reduced/GBIFdatasummary_reduced.csv')

pdf('./plots/full_reduced_dataset_BAMMdataset_biases_new.pdf',width=10,height=7)
par(mfrow=c(2,3))
plot_GBIFdata_proportions(GBIFdata.file='./results/full_biased_reduced/GBIFdatasummary_reduced.csv')
#this plots the proportion of trop and temp species in the GBIF BAMM dataset, and the proportion of species in each lat band
plot_GBIFdataBAMM_proportions(GBIFdata.file='./results/full_biased_reduced/GBIFdatasummary_reduced.csv',BAMM.object.path='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
#there's an obvious bias, with many more species in temperate regions being sequenced
dev.off()


#####5) run STRAPP for GBIFdataBAMM subsampled datasets####
#these are 'unbiased' datasets (i.e., the proportions of tropical and temperate species are = to the proportions in original GBIFdata)
dir.create('./results/subsampled_unbiased/')
dir.create('./results/subsampled_unbiased/input_tables/')

####generate the replicate datasets####
GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
#tropical defined as abs(median latitude <23.5)
#temperate defined as abs(median latitude >23.5)
GBIFdata$tropical<-0
GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
sampling.proportions<-table(GBIFdata$tropical)/nrow(GBIFdata)
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
#subset table to BAMM species
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
#remove duplicates in GBIFdata (6 species with same info but two rows - one with present in garden, one with absent in garden)#
GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
GBIFdata.BAMM<-unique(GBIFdata.BAMM)
GBIFdata.BAMM$latitudinal.band<-NA
GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
#size.dataset is the size of the subsampled datasets
size.dataset<-10000
for (replicate in c(1:100)){
  GBIFdata.BAMM.subsampled<-GBIFdata.BAMM[c(sample(which(GBIFdata.BAMM$tropical==0),size = round(sampling.proportions[1]*size.dataset)),sample(which(GBIFdata.BAMM$tropical==1),size = round(sampling.proportions[2]*size.dataset))),]
  write.table(GBIFdata.BAMM.subsampled,file=paste('./results/subsampled_unbiased/input_tables/GBIFdata.BAMM.subsampled_',replicate,'_table.txt',sep=''),sep='\t',quote=F,row.names=F)
}
#copy the GBIFdata.BAMM.subsampled_ files in ./results/subsampled_unbiased/input_tables/ to the cluster

####run STRAPP 100 subsampled unbiased full dataset (binary, continuous, latband, absolute lat band) x replicates####
#run run_strapp_subsample_replicates_cluster.R in cluster, 100 times
#e.g.: /scripts/conscriptoR ~/ldg_plants/run_strapp_subsample_replicates_cluster.R 1
#e.g.: /scripts/conscriptoR /ldg_plants/run_strapp_subsample_replicates_cluster.R 2

#then copy the outputs in cluster to ./results/strapp_unbiased/
#and sort in folders (abs.latitudinalband, abs.medianlatitude.continuous, latitudinalband, medianlatitude.continuous, tropical.binary)
dir.create('./results/subsampled_unbiased/abs.latitudinalband/')
dir.create('./results/subsampled_unbiased/abs.medianlatitude/')
dir.create('./results/subsampled_unbiased/latitudinalband')
dir.create('./results/subsampled_unbiased/medianlatitude/')
dir.create('./results/subsampled_unbiased/tropical/')

####plot STRAPP summary results of the abs.medianlatitude, tropical.binary and abs.latitudinalband for subsampled datasets####
#absolute median latitude strapp
pdf('./plots/Fig3_strapp_results_GBIFbamm_10k_subsampledreplicates.pdf',width=10,height=7)
par(mfrow=c(2,3))
blue<-brewer.pal(11,'RdYlBu')[9]
red<-brewer.pal(11,'RdYlBu')[2]

####Figure 2 plots####
source('./R/plot_figs.R')
#get plots for Figure 2
pdf('./plots/Fig2_all_subsampled10k_quantilerate_vs_lambdarate_extremebias.pdf',paper='a4r')
par(mfrow=c(2,2))
par(pty='s')
plot_Fig2_allsubsampled_extremebias_quantilerate_vs_lambda()
dev.off()
pdf('./plots/Fig2_map_latitudinalbands_log_medians_extremebias.pdf',paper='a4r')
plot_Fig2_latitudinalband_boxplots_subsampled_extremebias()
dev.off()
pdf('./plots/Fig2_map_latitudinalbands_log_medians_extremebias_summary.pdf',paper='a4r')
plot_Fig2_latitudinalband_boxplots_subsampled_extremebias_summary()
dev.off()

pdf('./plots/Fig2_map_latitudinalbands_log_medians_unbiased_summary.pdf',paper='a4r')
plot_Fig2_latitudinalband_boxplots_subsampled_unbiased_summary()
dev.off()

pdf('./plots/Fig2_strapp_correlations.pdf')
par(pty='s')
plot_Fig2_strapp_samples_from_files(path.correlationfile='./results/subsampled_unbiased/',parameter='lambda',tablefile='./results/GBIFdata.BAMM.table.txt')
dev.off()
plot_Fig2_strapp_densities(path.strappobjects = './results/subsampled_unbiased/')
##absolute median latitude strapp
#abs.medianlat.strapp.files<-list.files('./results/subsampled_unbiased/abs.medianlatitude/',pattern='abs.medianlatitude.+.txt')
#abs.medianlat.strapp.tables<-lapply(abs.medianlat.strapp.files,function(x)read.table(paste('./results/subsampled_unbiased/abs.medianlatitude/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
#abs.medianlat.strapp.speciation<-do.call('rbind.fill',lapply(abs.medianlat.strapp.tables,function(x)x[1,]))
##plot median rho of speciation rate ~ absolute median lat
#plot(density(unlist(lapply(abs.medianlat.strapp.tables,function(x)x$estimate[1]))),xlim=c(0,0.2),xlab='Spearmans rho',ylab='',yaxt='n',main='STRAPP.abs.median.latitude.GBIF.BAMM.10k.unbiased.subsamples',yaxs='i',xaxs='i',cex.main=.8)
#plot(density(unlist(lapply(abs.medianlat.strapp.tables,function(x)x$p.value[1]))),xlim=c(0,0.1),xlab='p.value',ylab='',yaxt='n',main='STRAPP.abs.median.latitude.GBIF.BAMM.10k.unbiased.subsamples',yaxs='i',xaxs='i',cex.main=.8)
#abline(v=0.05,lty=2)
##tropical binary
#tropical.binary.strapp.files<-list.files('./results/subsampled_unbiased/tropical/',pattern='^tropical.strapp.+.txt')
#tropical.binary.strapp.tables<-lapply(tropical.binary.strapp.files,function(x)read.table(paste('./results/subsampled_unbiased/tropical/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
#tropical.binary.strapp.tables.speciation<-do.call('rbind',lapply(tropical.binary.strapp.tables,function(x)x[1,]))
##plot median of speciation rate in temp (0) and trop (1)
#boxplot(tropical.binary.strapp.tables.speciation$estimate.0,tropical.binary.strapp.tables.speciation$estimate.1,names=c('temperate','tropical'),pch=16,ylim=c(0.5,2),ylab='median.lambda',main='GBIF.BAMM.10k.unbiased.subsamples',outline=F,cex.main=.8)
#points(x=jitter(rep(1,length(tropical.binary.strapp.tables.speciation$estimate.0)),2),y=tropical.binary.strapp.tables.speciation$estimate.0,col=adjustcolor(blue,alpha.f = 0.5),pch=16)
#points(x=jitter(rep(2,length(tropical.binary.strapp.tables.speciation$estimate.1)),2),y=tropical.binary.strapp.tables.speciation$estimate.1,col=adjustcolor(red,alpha.f = 0.5),pch=16)
##and hist of pvalues
plot(density(tropical.binary.strapp.tables.speciation$p.value),xlim=c(0,0.1),xlab='p.value',ylab='',yaxt='n',main='STRAPP.tropical.GBIF.BAMM.10k.unbiased.subsamples',yaxs='i',xaxs='i',cex.main=.8)
abline(v=0.05,lty=2)
#absolute lat band strapp
abs.latband.strapp.files<-list.files('./results/subsampled_unbiased/abs.latitudinalband/',pattern='abs.latitudinalband.+.txt')
abs.latband.strapp.tables<-lapply(abs.latband.strapp.files,function(x)read.table(paste('./results/subsampled_unbiased/abs.latitudinalband/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
abs.latband.strapp.speciation<-do.call('rbind.fill',lapply(abs.latband.strapp.tables,function(x)x[1,]))
#plot median of speciation rate across absolute lat bands
boxplot(abs.latband.strapp.speciation[,c(1:(grep('p.value',colnames(abs.latband.strapp.speciation))-1))],names=c(seq(0,80,by=10)),pch=16,ylab='median.lambda',outline=F,main='GBIF.BAMM.10k.unbiased.subsamples',xlab='abs.latitudinal.band',cex.main=.8)
#and hist of pvalues
plot(density(abs.latband.strapp.speciation$p.value),xlim=c(0,0.1),xlab='p.value',ylab='',yaxt='n',main='STRAPP.abs.latitudinal.band.GBIF.BAMM.10k.unbiased.subsamples',yaxs='i',xaxs='i',cex.main=.8)
abline(v=0.05,lty=2)
dev.off()
####Figure 2 plots. alternative####
#plot Fig2 with axis as in Fig1
pdf('./plots/Fig2_all_subsampled10k_quantilerate_vs_lambdarate_full_xaxisFig1.pdf')
plot_Fig2_allsubsampled_quantilerate_vs_lambda()
dev.off()

pdf('./plots/Fig2_map_latitudinalbands_log_medians_xaxisFig1.pdf',paper='a4r')
plot_Fig2_latitudinalband_boxplots_subsampled_unbiased_xaxisFig1()
dev.off()

#plot a random replicate of the 100 unbiased STRAPP datasets with the same scale than Fig.1
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
#load a replicate
replicate<-read.table('./results/subsampled_unbiased/input_tables/GBIFdata.BAMM.subsampled_1_table.txt',header=T,sep='\t',stringsAsFactors = F)
replicate.table<-GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$binomial%in%replicate$binomial,]
pdf('./plots/FigS3_map_latitudinalbands_log_subsampled_unbiased_1.pdf',paper='a4r')
plot_Fig1_map_bands(GBIFdata.BAMM.rates=replicate.table)
dev.off()

####Figure 2 plots for extremebias####
source('./R/plot_figs.R')
#get plots for Figure 2
plot_Fig2_allsubsampled_extremebias_quantilerate_vs_lambda()
plot_Fig2_latitudinalband_boxplots_subsampled_extremebias()
pdf('./plots/Fig2_strapp_correlations_extremebias.pdf')
par(pty='s')
plot_Fig2_strapp_samples_from_files(path.correlationfile='./results/subsampled_extremebias/',parameter='lambda',tablefile='./results/GBIFdata.BAMM.table.txt')
dev.off()

pdf('./plots/Fig2_strappdensities_subsampled_extremebias.pdf')
par(pty='s')
plot_Fig2_strapp_densities(path.strappobjects = './results/subsampled_extremebias/')
dev.off()


####5B) run STRAPP for GBIFdataBAMM subsampled extremebias datasets####
#these are 'unbiased' datasets (i.e., the proportions of tropical and temperate species are = to the proportions in original GBIFdata)
dir.create('./results/subsampled_extremebias/')
dir.create('./results/subsampled_extremebias/input_tables/')
####generate the replicate datasets####
#assessing missing taxa from GBIF
#this is the total number of accepted species
plant_lookup_table<-plant_lookup(include_counts = T)
sum(plant_lookup_table$number.of.accepted.species)
#proportion of species in GBIF dataset
GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
#tropical defined as abs(median latitude <23.5)
#temperate defined as abs(median latitude >23.5)
GBIFdata$tropical<-0
GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
#these are the real numbers of tropical and temperate species on GBIF
table(GBIFdata$tropical)
sampling.proportions<-table(GBIFdata$tropical)/nrow(GBIFdata)
#the most extreme case: all the species missing from GBIF are tropical
#add missing species to tropical category
tropicality.GBIF.extreme.bias.table<-table(GBIFdata$tropical)
tropicality.GBIF.extreme.bias.table[2]<-tropicality.GBIF.extreme.bias.table[2]+(sum(plant_lookup_table$number.of.accepted.species)-nrow(GBIFdata))
#these are the proportions under the most extreme case of bias (0.29 temperate, 0.71 tropical)
sampling.proportions.extreme.bias<-tropicality.GBIF.extreme.bias.table/(sum(tropicality.GBIF.extreme.bias.table))
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
#subset table to BAMM species
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
#remove duplicates in GBIFdata (6 species with same info but two rows - one with present in garden, one with absent in garden)#
GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
GBIFdata.BAMM<-unique(GBIFdata.BAMM)
GBIFdata.BAMM$latitudinal.band<-NA
GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
#size.dataset is the size of the subsampled datasets
size.dataset<-10000
for (replicate in c(1:10)){
  GBIFdata.BAMM.subsampled.extremebias<-GBIFdata.BAMM[c(sample(which(GBIFdata.BAMM$tropical==0),size = round(sampling.proportions.extreme.bias[1]*size.dataset)),sample(which(GBIFdata.BAMM$tropical==1),size = round(sampling.proportions.extreme.bias[2]*size.dataset))),]
  write.table(GBIFdata.BAMM.subsampled.extremebias,file=paste('./results/subsampled_extremebias/input_tables/GBIFdata.BAMM.subsampled_extremebias_',replicate,'_table.txt',sep=''),sep='\t',quote=F,row.names=F)
}

#copy the GBIFdata.BAMM.subsampled_extremebias files in ./results/subsampled_extremebias/input_tables/ to the cluster

####run STRAPP 100 subsampled extremebias full dataset (binary, continuous, latband, absolute lat band) x replicates####
#run run_strapp_subsample_replicates_cluster.R in cluster, 100 times
#e.g.: /scripts/conscriptoR ~/ldg_plants/run_strapp_subsample_extremebias_replicates_cluster.R 1
#e.g.: /scripts/conscriptoR /ldg_plants/run_strapp_subsample_replicates_cluster.R 2

####5C) small datasets####
#assessing missing taxa from GBIF
#this is the total number of accepted species
plant_lookup_table<-plant_lookup(include_counts = T)
sum(plant_lookup_table$number.of.accepted.species)
#proportion of species in GBIF dataset
GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
#tropical defined as abs(median latitude <23.5)
#temperate defined as abs(median latitude >23.5)
GBIFdata$tropical<-0
GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
#these are the real numbers of tropical and temperate species on GBIF
table(GBIFdata$tropical)
sampling.proportions<-table(GBIFdata$tropical)/nrow(GBIFdata)
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
#subset table to BAMM species
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
#remove duplicates in GBIFdata (6 species with same info but two rows - one with present in garden, one with absent in garden)#
GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
GBIFdata.BAMM<-unique(GBIFdata.BAMM)
GBIFdata.BAMM$latitudinal.band<-NA
GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
#create folder to store the BAMM replicates
dir.create('./results/subsampled_small/')
dir.create('./results/subsampled_small/input_tables/')
#size.dataset is the size of the subsampled datasets
#I'll take the proportion of plant species in Zanne tree (~8%)

#get trees that have 8% of the species in Zanne and compare rates
size.dataset<-round((nrow(GBIFdata.BAMM)/sum(plant_lookup_table$number.of.accepted.species))*nrow(GBIFdata.BAMM))
for (replicate in c(1:10)){
  GBIFdata.BAMM.subsampled.small<-GBIFdata.BAMM[c(sample(which(GBIFdata.BAMM$tropical==0),size = round(sampling.proportions[1]*size.dataset)),sample(which(GBIFdata.BAMM$tropical==1),size = round(sampling.proportions[2]*size.dataset))),]
  write.table(GBIFdata.BAMM.subsampled.small,file=paste('./results/subsampled_small/input_tables/GBIFdata.BAMM.subsampled_small_',replicate,'_table.txt',sep=''),sep='\t',quote=F,row.names=F)
}




####6) compare BAMM rate estimates of 28k tree vs 10k subsampled trees####
####split BAMM_subsampled_replicates into clades for Zanne tree analyses in BAMM####
#generate BAMM input files for 10 replicates (control files are in './raw_data/BAMM_Zanneclades/control_files/)
source('./R/BAMM_Zanneclades_subsampled.R')
lapply(c(1:10),function(x) BAMM_Zanneclades_subsampled(replicate=x))

#copy priors from priors file into control file, run BAMM clades on hydrogen (see './R/csmit_BAMM_run_example.txt')
#each clade may take around 3-4 to run on 4 processors
#most will require a second 50 million run (cont1.txt control files included)

#copy mcmc_out and event_data files into the corresponding /results/ folder (e.g './results/Zanne_clades_BAMM/subsampled_1/results/')

####analyse the BAMM runs, check convergence and merge all clades into a single analysis for each replicate####
source('./R/BAMMoutput_subsampled_replicates_analyse.R')
dir.create('./results/Zanne_clades_BAMM/event_data/')

#check mcmc convergence in subsampled 1
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_1/results/mcmc_out_GBIFdata.BAMM.subsampled_1_clade_1_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_1/results/mcmc_out_GBIFdata.BAMM.subsampled_1_clade_2_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_1/results/mcmc_out_GBIFdata.BAMM.subsampled_1_clade_3_50.txt',burnin=0.5)
#analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_1/results/mcmc_out_GBIFdata.BAMM.subsampled_1_clade_4_50.txt',burnin=0.8)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_1/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_1_clade_4_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_1/results/mcmc_out_GBIFdata.BAMM.subsampled_1_clade_5_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_1/results/mcmc_out_GBIFdata.BAMM.subsampled_1_clade_6_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_1/results/mcmc_out_GBIFdata.BAMM.subsampled_1_clade_7_50.txt',burnin=0.5)
#merging all analyses into one for subsampled 1 and saving eventdata to event_data folder
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM/subsampled_1/results/event_data_GBIFdata.BAMM.subsampled_1_clade_1_50.txt','./results/Zanne_clades_BAMM/subsampled_1/results/event_data_GBIFdata.BAMM.subsampled_1_clade_2_50.txt','./results/Zanne_clades_BAMM/subsampled_1/results/event_data_GBIFdata.BAMM.subsampled_1_clade_3_50.txt','./results/Zanne_clades_BAMM/subsampled_1/results/run2_event_data_GBIFdata.BAMM.subsampled_1_clade_4_50.txt','./results/Zanne_clades_BAMM/subsampled_1/results/event_data_GBIFdata.BAMM.subsampled_1_clade_5_50.txt','./results/Zanne_clades_BAMM/subsampled_1/results/event_data_GBIFdata.BAMM.subsampled_1_clade_6_50.txt','./results/Zanne_clades_BAMM/subsampled_1/results/event_data_GBIFdata.BAMM.subsampled_1_clade_7_50.txt'),burnin = 0.5,name='subsampled_1')

#check mcmc convergence in subsampled 2
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_2/results/mcmc_out_GBIFdata.BAMM.subsampled_2_clade_1_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_2/results/mcmc_out_GBIFdata.BAMM.subsampled_2_clade_2_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_2/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_2_clade_3_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_2/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_2_clade_4_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_2/results/mcmc_out_GBIFdata.BAMM.subsampled_2_clade_5_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_2/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_2_clade_6_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_2/results/mcmc_out_GBIFdata.BAMM.subsampled_2_clade_7_50.txt',burnin=0.4)
#merging all analyses into one for subsampled 2 and saving eventdata to event_data folder
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM/subsampled_2/results/event_data_GBIFdata.BAMM.subsampled_2_clade_1_50.txt','./results/Zanne_clades_BAMM/subsampled_2/results/event_data_GBIFdata.BAMM.subsampled_2_clade_2_50.txt','./results/Zanne_clades_BAMM/subsampled_2/results/run2_event_data_GBIFdata.BAMM.subsampled_2_clade_3_50.txt','./results/Zanne_clades_BAMM/subsampled_2/results/run2_event_data_GBIFdata.BAMM.subsampled_2_clade_4_50.txt','./results/Zanne_clades_BAMM/subsampled_2/results/event_data_GBIFdata.BAMM.subsampled_2_clade_5_50.txt','./results/Zanne_clades_BAMM/subsampled_2/results/run2_event_data_GBIFdata.BAMM.subsampled_2_clade_6_50.txt','./results/Zanne_clades_BAMM/subsampled_2/results/event_data_GBIFdata.BAMM.subsampled_2_clade_7_50.txt'),burnin = 0.4,name='subsampled_2')


#check mcmc convergence in subsampled 3
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_3/results/mcmc_out_GBIFdata.BAMM.subsampled_3_clade_1_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_3/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_3_clade_1_50.txt',burnin=0)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_3/results/mcmc_out_GBIFdata.BAMM.subsampled_3_clade_2_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_3/results/mcmc_out_GBIFdata.BAMM.subsampled_3_clade_3_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_3/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_3_clade_2_50.txt',burnin=0)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_3/results/mcmc_out_GBIFdata.BAMM.subsampled_3_clade_4_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_3/results/mcmc_out_GBIFdata.BAMM.subsampled_3_clade_5_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_3/results/mcmc_out_GBIFdata.BAMM.subsampled_3_clade_6_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_3/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_3_clade_6_50.txt',burnin=0)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_3/results/mcmc_out_GBIFdata.BAMM.subsampled_3_clade_7_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_3/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_3_clade_7_50.txt',burnin=0)
#merging all analyses into one for subsampled 3 and saving eventdata to event_data folder
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM//subsampled_3/results/run2_event_data_GBIFdata.BAMM.subsampled_3_clade_1_50.txt','./results/Zanne_clades_BAMM//subsampled_3/results/run2_event_data_GBIFdata.BAMM.subsampled_3_clade_2_50.txt','./results/Zanne_clades_BAMM//subsampled_3/results/run2_event_data_GBIFdata.BAMM.subsampled_3_clade_2_50.txt','./results/Zanne_clades_BAMM//subsampled_3/results/event_data_GBIFdata.BAMM.subsampled_3_clade_4_50.txt','./results/Zanne_clades_BAMM//subsampled_3/results/event_data_GBIFdata.BAMM.subsampled_3_clade_5_50.txt','./results/Zanne_clades_BAMM//subsampled_3/results/run2_event_data_GBIFdata.BAMM.subsampled_3_clade_6_50.txt','./results/Zanne_clades_BAMM//subsampled_3/results/run2_event_data_GBIFdata.BAMM.subsampled_3_clade_7_50.txt'),burnin = 0.5,name='subsampled_3')

#check mcmc convergence in subsampled 4
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/mcmc_out_GBIFdata.BAMM.subsampled_4_clade_1_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_4_clade_1_50.txt',burnin=0)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/mcmc_out_GBIFdata.BAMM.subsampled_4_clade_2_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_4_clade_2_50.txt',burnin=0)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/mcmc_out_GBIFdata.BAMM.subsampled_4_clade_3_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_4_clade_3_50.txt',burnin=0)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/mcmc_out_GBIFdata.BAMM.subsampled_4_clade_4_50.txt',burnin=0.8)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_4_clade_4_50.txt',burnin=0)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/mcmc_out_GBIFdata.BAMM.subsampled_4_clade_5_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/mcmc_out_GBIFdata.BAMM.subsampled_4_clade_6_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_4_clade_6_50.txt',burnin=0)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_4/results/mcmc_out_GBIFdata.BAMM.subsampled_4_clade_7_50.txt',burnin=0.2)
#merging all analyses into one for subsampled 4 and saving eventdata to event_data folder
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM//subsampled_4/results/run2_event_data_GBIFdata.BAMM.subsampled_4_clade_1_50.txt','./results/Zanne_clades_BAMM//subsampled_4/results/run2_event_data_GBIFdata.BAMM.subsampled_4_clade_2_50.txt','./results/Zanne_clades_BAMM//subsampled_4/results/run2_event_data_GBIFdata.BAMM.subsampled_4_clade_3_50.txt','./results/Zanne_clades_BAMM//subsampled_4/results/run2_event_data_GBIFdata.BAMM.subsampled_4_clade_4_50.txt','./results/Zanne_clades_BAMM//subsampled_4/results/event_data_GBIFdata.BAMM.subsampled_4_clade_5_50.txt','./results/Zanne_clades_BAMM//subsampled_4/results/run2_event_data_GBIFdata.BAMM.subsampled_4_clade_6_50.txt','./results/Zanne_clades_BAMM//subsampled_4/results/event_data_GBIFdata.BAMM.subsampled_4_clade_7_50.txt'),burnin = 0.4,name='subsampled_4')

#check mcmc convergence in subsampled 5
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_5/results/mcmc_out_GBIFdata.BAMM.subsampled_5_clade_1_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_5/results/mcmc_out_GBIFdata.BAMM.subsampled_5_clade_2_50.txt',burnin=0.2)
#analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_5/results/mcmc_out_GBIFdata.BAMM.subsampled_5_clade_3_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_5/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_5_clade_3_50.txt',burnin=0)
#analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_5/results/mcmc_out_GBIFdata.BAMM.subsampled_5_clade_4_50.txt',burnin=0.8)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_5/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_5_clade_4_50.txt',burnin=0)
#analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_5/results/mcmc_out_GBIFdata.BAMM.subsampled_5_clade_5_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_5/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_5_clade_5_50.txt',burnin=0)
#analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_5/results/mcmc_out_GBIFdata.BAMM.subsampled_5_clade_6_50.txt',burnin=0.3)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_5/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_5_clade_6_50.txt',burnin=0)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_5/results/mcmc_out_GBIFdata.BAMM.subsampled_5_clade_7_50.txt',burnin=0.5)
#merging all analyses into one for subsampled 5 and saving eventdata to event_data folder
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM//subsampled_5/results/event_data_GBIFdata.BAMM.subsampled_5_clade_1_50.txt','./results/Zanne_clades_BAMM//subsampled_5/results/event_data_GBIFdata.BAMM.subsampled_5_clade_2_50.txt','./results/Zanne_clades_BAMM//subsampled_5/results/run2_event_data_GBIFdata.BAMM.subsampled_5_clade_3_50.txt','./results/Zanne_clades_BAMM//subsampled_5/results/run2_event_data_GBIFdata.BAMM.subsampled_5_clade_4_50.txt','./results/Zanne_clades_BAMM//subsampled_5/results/run2_event_data_GBIFdata.BAMM.subsampled_5_clade_5_50.txt','./results/Zanne_clades_BAMM//subsampled_5/results/run2_event_data_GBIFdata.BAMM.subsampled_5_clade_6_50.txt','./results/Zanne_clades_BAMM//subsampled_5/results/event_data_GBIFdata.BAMM.subsampled_5_clade_7_50.txt'),burnin=0.5,name='subsampled_5')

#check mcmc convergence in subsampled 6
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_6/results/mcmc_out_GBIFdata.BAMM.subsampled_6_clade_1_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_6/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_6_clade_1_50.txt',burnin=0)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_6/results/mcmc_out_GBIFdata.BAMM.subsampled_6_clade_2_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_6/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_6_clade_2_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_6/results/mcmc_out_GBIFdata.BAMM.subsampled_6_clade_3_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_6/results/mcmc_out_GBIFdata.BAMM.subsampled_6_clade_4_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_6/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_6_clade_4_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_6/results/mcmc_out_GBIFdata.BAMM.subsampled_6_clade_5_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_6/results/mcmc_out_GBIFdata.BAMM.subsampled_6_clade_6_50.txt',burnin=0.3)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_6/results/mcmc_out_GBIFdata.BAMM.subsampled_6_clade_7_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_6/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_6_clade_7_50.txt',burnin=0.2)
#merging all analyses into one for subsampled 6 and saving eventdata to event_data folder
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM//subsampled_6/results/run2_event_data_GBIFdata.BAMM.subsampled_6_clade_1_50.txt','./results/Zanne_clades_BAMM//subsampled_6/results/run2_event_data_GBIFdata.BAMM.subsampled_6_clade_2_50.txt','./results/Zanne_clades_BAMM//subsampled_6/results/event_data_GBIFdata.BAMM.subsampled_6_clade_3_50.txt','./results/Zanne_clades_BAMM//subsampled_6/results/run2_event_data_GBIFdata.BAMM.subsampled_6_clade_4_50.txt','./results/Zanne_clades_BAMM//subsampled_6/results/event_data_GBIFdata.BAMM.subsampled_6_clade_5_50.txt','./results/Zanne_clades_BAMM//subsampled_6/results/event_data_GBIFdata.BAMM.subsampled_6_clade_6_50.txt','./results/Zanne_clades_BAMM//subsampled_6/results/run2_event_data_GBIFdata.BAMM.subsampled_6_clade_7_50.txt'),burnin = 0.4,name='subsampled_6')

#check mcmc convergence in subsampled 7
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_7/results/mcmc_out_GBIFdata.BAMM.subsampled_7_clade_1_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_7/results/mcmc_out_GBIFdata.BAMM.subsampled_7_clade_2_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_7/results/mcmc_out_GBIFdata.BAMM.subsampled_7_clade_3_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_7/results/mcmc_out_GBIFdata.BAMM.subsampled_7_clade_4_50.txt',burnin=0.6)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_7/results/mcmc_out_GBIFdata.BAMM.subsampled_7_clade_5_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_7/results/mcmc_out_GBIFdata.BAMM.subsampled_7_clade_6_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_7/results/mcmc_out_GBIFdata.BAMM.subsampled_7_clade_7_50.txt',burnin=0.2)
#merging all analyses into one for subsampled 7 and saving eventdata to event_data folder
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM//subsampled_7/results/event_data_GBIFdata.BAMM.subsampled_7_clade_1_50.txt','./results/Zanne_clades_BAMM//subsampled_7/results/event_data_GBIFdata.BAMM.subsampled_7_clade_2_50.txt','./results/Zanne_clades_BAMM//subsampled_7/results/event_data_GBIFdata.BAMM.subsampled_7_clade_3_50.txt','./results/Zanne_clades_BAMM//subsampled_7/results/event_data_GBIFdata.BAMM.subsampled_7_clade_4_50.txt','./results/Zanne_clades_BAMM//subsampled_7/results/event_data_GBIFdata.BAMM.subsampled_7_clade_5_50.txt','./results/Zanne_clades_BAMM//subsampled_7/results/event_data_GBIFdata.BAMM.subsampled_7_clade_6_50.txt','./results/Zanne_clades_BAMM//subsampled_7/results/event_data_GBIFdata.BAMM.subsampled_7_clade_7_50.txt'),burnin=0.6,name='subsampled_7')

#check mcmc convergence in subsampled 8
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM//subsampled_8/results/mcmc_out_GBIFdata.BAMM.subsampled_8_clade_1_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM//subsampled_8/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_8_clade_1_50.txt',burnin=0)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM//subsampled_8/results/mcmc_out_GBIFdata.BAMM.subsampled_8_clade_2_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM//subsampled_8/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_8_clade_2_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM//subsampled_8/results/mcmc_out_GBIFdata.BAMM.subsampled_8_clade_3_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM//subsampled_8/results/mcmc_out_GBIFdata.BAMM.subsampled_8_clade_4_50.txt',burnin=0.6)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM//subsampled_8/results/mcmc_out_GBIFdata.BAMM.subsampled_8_clade_5_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM//subsampled_8/results/mcmc_out_GBIFdata.BAMM.subsampled_8_clade_6_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM//subsampled_8/results/mcmc_out_GBIFdata.BAMM.subsampled_8_clade_7_50.txt',burnin=0.3)
#merging all analyses into one for subsampled 8 and saving eventdata to event_data folder
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM//subsampled_8/results/run2_event_data_GBIFdata.BAMM.subsampled_8_clade_1_50.txt','./results/Zanne_clades_BAMM//subsampled_8/results/run2_event_data_GBIFdata.BAMM.subsampled_8_clade_2_50.txt','./results/Zanne_clades_BAMM//subsampled_8/results/event_data_GBIFdata.BAMM.subsampled_8_clade_3_50.txt','./results/Zanne_clades_BAMM//subsampled_8/results/event_data_GBIFdata.BAMM.subsampled_8_clade_4_50.txt','./results/Zanne_clades_BAMM//subsampled_8/results/event_data_GBIFdata.BAMM.subsampled_8_clade_5_50.txt','./results/Zanne_clades_BAMM//subsampled_8/results/event_data_GBIFdata.BAMM.subsampled_8_clade_6_50.txt','./results/Zanne_clades_BAMM//subsampled_8/results/event_data_GBIFdata.BAMM.subsampled_8_clade_7_50.txt'),burnin=0.6,name='subsampled_8')

#check mcmc convergence in subsampled 9
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_9/results/mcmc_out_GBIFdata.BAMM.subsampled_9_clade_1_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_9/results/mcmc_out_GBIFdata.BAMM.subsampled_9_clade_2_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_9/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_9_clade_2_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_9/results/mcmc_out_GBIFdata.BAMM.subsampled_9_clade_3_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_9/results/mcmc_out_GBIFdata.BAMM.subsampled_9_clade_4_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_9/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_9_clade_4_50.txt',burnin=0.6)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_9/results/mcmc_out_GBIFdata.BAMM.subsampled_9_clade_5_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_9/results/mcmc_out_GBIFdata.BAMM.subsampled_9_clade_6_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_9/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_9_clade_6_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_9/results/mcmc_out_GBIFdata.BAMM.subsampled_9_clade_7_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_9/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_9_clade_7_50.txt',burnin=0.2)
#merging all analyses into one for subsampled 9 and saving eventdata to event_data folder
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM//subsampled_9/results/event_data_GBIFdata.BAMM.subsampled_9_clade_1_50.txt','./results/Zanne_clades_BAMM//subsampled_9/results/run2_event_data_GBIFdata.BAMM.subsampled_9_clade_2_50.txt','./results/Zanne_clades_BAMM//subsampled_9/results/event_data_GBIFdata.BAMM.subsampled_9_clade_3_50.txt','./results/Zanne_clades_BAMM//subsampled_9/results/run2_event_data_GBIFdata.BAMM.subsampled_9_clade_4_50.txt','./results/Zanne_clades_BAMM//subsampled_9/results/event_data_GBIFdata.BAMM.subsampled_9_clade_5_50.txt','./results/Zanne_clades_BAMM//subsampled_9/results/run2_event_data_GBIFdata.BAMM.subsampled_9_clade_6_50.txt','./results/Zanne_clades_BAMM//subsampled_9/results/run2_event_data_GBIFdata.BAMM.subsampled_9_clade_7_50.txt'),burnin=0.4,name='subsampled_9')

#check mcmc convergence in subsampled 10
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_10/results/mcmc_out_GBIFdata.BAMM.subsampled_10_clade_1_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_10/results/mcmc_out_GBIFdata.BAMM.subsampled_10_clade_2_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_10/results/mcmc_out_GBIFdata.BAMM.subsampled_10_clade_3_50.txt',burnin=0.4)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_10/results/mcmc_out_GBIFdata.BAMM.subsampled_10_clade_4_50.txt',burnin=0.6)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_10/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_10_clade_4_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_10/results/mcmc_out_GBIFdata.BAMM.subsampled_10_clade_5_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_10/results/mcmc_out_GBIFdata.BAMM.subsampled_10_clade_6_50.txt',burnin=0.5)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_10/results/run2_mcmc_out_GBIFdata.BAMM.subsampled_10_clade_6_50.txt',burnin=0.2)
analyse_BAMM_convergence(mcmcout = './results/Zanne_clades_BAMM/subsampled_10/results/mcmc_out_GBIFdata.BAMM.subsampled_10_clade_7_50.txt',burnin=0.2)
#merging all analyses into one for subsampled 10 and saving eventdata to event_data folder
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM//subsampled_10/results/event_data_GBIFdata.BAMM.subsampled_10_clade_1_50.txt','./results/Zanne_clades_BAMM//subsampled_10/results/event_data_GBIFdata.BAMM.subsampled_10_clade_2_50.txt','./results/Zanne_clades_BAMM//subsampled_10/results/event_data_GBIFdata.BAMM.subsampled_10_clade_3_50.txt','./results/Zanne_clades_BAMM//subsampled_10/results/run2_event_data_GBIFdata.BAMM.subsampled_10_clade_4_50.txt','./results/Zanne_clades_BAMM//subsampled_10/results/event_data_GBIFdata.BAMM.subsampled_10_clade_5_50.txt','./results/Zanne_clades_BAMM//subsampled_10/results/run2_event_data_GBIFdata.BAMM.subsampled_10_clade_6_50.txt','./results/Zanne_clades_BAMM//subsampled_10/results/event_data_GBIFdata.BAMM.subsampled_10_clade_7_50.txt'),burnin=0.5,name='subsampled_10')

####get data frame with rates from full 28k tip dataset + 10 replicates of unbiased 10k datasets####
source('./R/compare_BAMMestimates_full_vs_subsampled.R')
rates.df.merge<-data_frame_compare_rates_full_vs_subsampled(folder.subsampled.event_data='./results/Zanne_clades_BAMM/event_data/')

full28k.vs.10k.results<-replicates_phylolm_rho_trop_vs_temp(rates.df.merge=rates.df.merge)
#this has the results of phylolm
full28k.vs.10k.lms<-full28k.vs.10k.results[[1]]
#check the regression coefficients for tropicality here
lapply(full28k.vs.10k.lms,function(x)summary(x)$coefficients)
lapply(full28k.vs.10k.lms,function(x)confint(x))
#calculate grand means and sds
grand.mean<-function(M, N) {weighted.mean(M, N)}
grand.sd<-function(S, M, N) {sqrt(weighted.mean(S^2+M^2,N)-weighted.mean(M, N)^2)}
sizes<-unlist(lapply(full28k.vs.10k.lms,function(x)x$n))
sds<-unlist(lapply(full28k.vs.10k.lms,function(x)summary(x)$coefficients[4,2]))
means<-unlist(lapply(full28k.vs.10k.lms,function(x)summary(x)$coefficients[4,1]))
mean(means)
grand.mean(M=means,N=sizes)
grand.sd(S=sds,M=means,N=sizes)
#ths has the results of rhos
full28k.vs.10k.rhos<-full28k.vs.10k.results[[2]]
median(unlist(lapply(full28k.vs.10k.rhos,function(x)x$estimate)))
median(unlist(lapply(full28k.vs.10k.rhos,function(x)x$p.value)))
plot_replicates_phylolm_trop_vs_temp(rates.df.merge=rates.df.merge)
####plot rate comparisons in all 10 subsampled datasets####
pdf('./plots/FigS4_lm_BAMMrates_28ktree_vs_10subsamples_10ktrees.pdf',paper='a4r')
par(mfrow=c(2,5))
plot_replicates_phylolm_trop_vs_temp(rates.df.merge=rates.df.merge)
dev.off()

####
####6A) compare BAMM rate estimates of full dataset and extremebias dataset####
source('./R/BAMM_Zanneclades_subsampled.R')
lapply(c(1:10),function(x) BAMM_Zanneclades_subsampled(replicate=x))

#copy priors from priors file into control file, run BAMM clades on hydrogen (see './R/csmit_BAMM_run_example.txt')
#each clade may take around 3-4 to run on 4 processors
#most will require a second 50 million run (cont1.txt control files included)

#copy mcmc_out and event_data files into the corresponding /results/ folder (e.g './results/Zanne_clades_BAMM/subsampled_1/results/')

####analyse the BAMM runs, check convergence and merge all clades into a single analysis for each replicate####
source('./R/BAMMoutput_subsampled_replicates_analyse.R')
dir.create('./results/Zanne_clades_BAMM/event_data/')

subsampled1.extremebias.files<-list.files('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_1/results/',pattern='run1_mcmc_out_')
lapply(subsampled1.extremebias.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_1/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_1/results/run1_run1_event_data_GBIFdata.BAMM.subsampled_extremebias_1_clade_1_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_1/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_1_clade_2_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_1/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_1_clade_3_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_1/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_1_clade_4_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_1/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_1_clade_5_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_1/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_1_clade_6_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_1/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_1_clade_7_50.txt'),burnin = 0.5,name='subsampled_extremebias_1',path = './results/Zanne_clades_BAMM_extremebias/')

subsampled2.extremebias.files<-list.files('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_2/results/',pattern='run1_mcmc_out_')
lapply(subsampled2.extremebias.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_2/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_2/results/run1_run1_event_data_GBIFdata.BAMM.subsampled_extremebias_2_clade_1_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_2/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_2_clade_2_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_2/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_2_clade_3_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_2/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_2_clade_4_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_2/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_2_clade_5_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_2/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_2_clade_6_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_2/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_2_clade_7_50.txt'),burnin = 0.5,name='subsampled_extremebias_2',path = './results/Zanne_clades_BAMM_extremebias/')

subsampled3.extremebias.files<-list.files('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_3/results/',pattern='run1_mcmc_out_')
lapply(subsampled3.extremebias.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_3/results/',x,sep=''),burnin=0.7))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_3/results/run1_run1_event_data_GBIFdata.BAMM.subsampled_extremebias_3_clade_1_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_3/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_3_clade_2_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_3/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_3_clade_3_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_3/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_3_clade_4_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_3/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_3_clade_5_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_3/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_3_clade_6_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_3/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_3_clade_7_50.txt'),burnin = 0.7,name='subsampled_extremebias_3',path = './results/Zanne_clades_BAMM_extremebias/')

subsampled4.extremebias.files<-list.files('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_4/results/',pattern='run1_mcmc_out_')
lapply(subsampled4.extremebias.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_4/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_4/results/run1_run1_event_data_GBIFdata.BAMM.subsampled_extremebias_4_clade_1_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_4/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_4_clade_2_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_4/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_4_clade_3_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_4/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_4_clade_4_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_4/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_4_clade_5_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_4/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_4_clade_6_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_4/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_4_clade_7_50.txt'),burnin = 0.5,name='subsampled_extremebias_4',path = './results/Zanne_clades_BAMM_extremebias/')

subsampled5.extremebias.files<-list.files('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_5/results/',pattern='run1_mcmc_out_')
lapply(subsampled5.extremebias.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_5/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_5/results/run1_run1_event_data_GBIFdata.BAMM.subsampled_extremebias_5_clade_1_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_5/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_5_clade_2_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_5/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_5_clade_3_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_5/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_5_clade_4_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_5/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_5_clade_5_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_5/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_5_clade_6_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_5/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_5_clade_7_50.txt'),burnin = 0.5,name='subsampled_extremebias_5',path = './results/Zanne_clades_BAMM_extremebias/')

subsampled6.extremebias.files<-list.files('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_6/results/',pattern='run1_mcmc_out_')
lapply(subsampled6.extremebias.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_6/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_6/results/run1_run1_event_data_GBIFdata.BAMM.subsampled_extremebias_6_clade_1_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_6/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_6_clade_2_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_6/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_6_clade_3_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_6/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_6_clade_4_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_6/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_6_clade_5_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_6/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_6_clade_6_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_6/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_6_clade_7_50.txt'),burnin = 0.5,name='subsampled_extremebias_6',path = './results/Zanne_clades_BAMM_extremebias/')

subsampled7.extremebias.files<-list.files('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_7/results/',pattern='run1_mcmc_out_')
lapply(subsampled7.extremebias.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_7/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_7/results/run1_run1_event_data_GBIFdata.BAMM.subsampled_extremebias_7_clade_1_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_7/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_7_clade_2_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_7/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_7_clade_3_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_7/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_7_clade_4_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_7/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_7_clade_5_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_7/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_7_clade_6_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_7/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_7_clade_7_50.txt'),burnin = 0.5,name='subsampled_extremebias_7',path = './results/Zanne_clades_BAMM_extremebias/')

subsampled8.extremebias.files<-list.files('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_8/results/',pattern='run1_mcmc_out_')
lapply(subsampled8.extremebias.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_8/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_8/results/run1_run1_event_data_GBIFdata.BAMM.subsampled_extremebias_8_clade_1_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_8/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_8_clade_2_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_8/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_8_clade_3_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_8/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_8_clade_4_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_8/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_8_clade_5_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_8/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_8_clade_6_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_8/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_8_clade_7_50.txt'),burnin = 0.5,name='subsampled_extremebias_8',path = './results/Zanne_clades_BAMM_extremebias/')

subsampled9.extremebias.files<-list.files('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_9/results/',pattern='run1_mcmc_out_')
lapply(subsampled9.extremebias.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_9/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_9/results/run1_run1_event_data_GBIFdata.BAMM.subsampled_extremebias_9_clade_1_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_9/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_9_clade_2_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_9/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_9_clade_3_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_9/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_9_clade_4_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_9/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_9_clade_5_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_9/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_9_clade_6_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_9/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_9_clade_7_50.txt'),burnin = 0.5,name='subsampled_extremebias_9',path = './results/Zanne_clades_BAMM_extremebias/')

subsampled10.extremebias.files<-list.files('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_10/results/',pattern='run1_mcmc_out_')
lapply(subsampled10.extremebias.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_10/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_10/results/run1_run1_event_data_GBIFdata.BAMM.subsampled_extremebias_10_clade_1_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_10/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_10_clade_2_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_10/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_10_clade_3_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_10/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_10_clade_4_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_10/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_10_clade_5_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_10/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_10_clade_6_50.txt','./results/Zanne_clades_BAMM_extremebias/subsampled_extremebias_10/results/run1_event_data_GBIFdata.BAMM.subsampled_extremebias_10_clade_7_50.txt'),burnin = 0.5,name='subsampled_extremebias_10',path = './results/Zanne_clades_BAMM_extremebias/')

source('./R/compare_BAMMestimates_full_vs_subsampled.R')
rates.df.merge.extremebias<-data_frame_compare_rates_full_vs_subsampled(folder.subsampled.event_data='./results/Zanne_clades_BAMM_extremebias/event_data/')

full28k.vs.10kextremebias.results<-replicates_phylolm_rho_trop_vs_temp_extremebias(rates.df.merge=rates.df.merge.extremebias)
#this has the results of phylolm
full28k.vs.10k.extremebias.lms<-full28k.vs.10kextremebias.results[[1]]
#check the regression coefficients for tropicality here
lapply(full28k.vs.10k.extremebias.lms,function(x)summary(x)$coefficients)
lapply(full28k.vs.10k.extremebias.lms,function(x)confint(x))
#calculate grand means and sds
grand.mean<-function(M, N) {weighted.mean(M, N)}
grand.sd<-function(S, M, N) {sqrt(weighted.mean(S^2+M^2,N)-weighted.mean(M, N)^2)}
sizes<-unlist(lapply(full28k.vs.10k.extremebias.lms,function(x)x$n))
sds<-unlist(lapply(full28k.vs.10k.extremebias.lms,function(x)summary(x)$coefficients[4,2]))
means<-unlist(lapply(full28k.vs.10k.extremebias.lms,function(x)summary(x)$coefficients[4,1]))
mean(means)
grand.mean(M=means,N=sizes)
grand.sd(S=sds,M=means,N=sizes)
#ths has the results of rhos
full28k.vs.10k.extremebias.rhos<-full28k.vs.10kextremebias.results[[2]]
median(unlist(lapply(full28k.vs.10k.extremebias.rhos,function(x)x$estimate)))
median(unlist(lapply(full28k.vs.10k.extremebias.rhos,function(x)x$p.value)))
plot_replicates_phylolm_trop_vs_temp_extremebias(rates.df.merge=rates.df.merge.extremebias)
####plot rate comparisons in all 10 subsampled extremebiased_datasets####
pdf('./plots/FigS4_lm_BAMMrates_28ktree_vs_10subsamples_10ktrees_extremebiased.pdf',paper='a4r')
par(mfrow=c(2,5))
plot_replicates_phylolm_trop_vs_temp_extremebias(rates.df.merge=rates.df.merge.extremebias)
dev.off()

####6B) compare BAMM rate estimates with "small" datasets####
#subsampled_small
subsampled1.small.files<-list.files('./results/Zanne_clades_BAMM_small/subsampled_small_1/results/',pattern='_out_')
lapply(subsampled1.small.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_small/subsampled_small_1/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_small/subsampled_small_1/results/event_data_GBIFdata.BAMM.subsampled_small_1_clade_1_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_1/results/event_data_GBIFdata.BAMM.subsampled_small_1_clade_2_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_1/results/event_data_GBIFdata.BAMM.subsampled_small_1_clade_3_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_1/results/event_data_GBIFdata.BAMM.subsampled_small_1_clade_4_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_1/results/event_data_GBIFdata.BAMM.subsampled_small_1_clade_5_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_1/results/event_data_GBIFdata.BAMM.subsampled_small_1_clade_6_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_1/results/event_data_GBIFdata.BAMM.subsampled_small_1_clade_7_50.txt'),burnin = 0.5,name='subsampled_small_1',path = './results/Zanne_clades_BAMM_small/')

subsampled2.small.files<-list.files('./results/Zanne_clades_BAMM_small/subsampled_small_2/results/',pattern='_out_')
lapply(subsampled2.small.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_small/subsampled_small_2/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_small/subsampled_small_2/results/event_data_GBIFdata.BAMM.subsampled_small_2_clade_1_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_2/results/event_data_GBIFdata.BAMM.subsampled_small_2_clade_2_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_2/results/event_data_GBIFdata.BAMM.subsampled_small_2_clade_3_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_2/results/event_data_GBIFdata.BAMM.subsampled_small_2_clade_4_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_2/results/event_data_GBIFdata.BAMM.subsampled_small_2_clade_5_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_2/results/event_data_GBIFdata.BAMM.subsampled_small_2_clade_6_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_2/results/event_data_GBIFdata.BAMM.subsampled_small_2_clade_7_50.txt'),burnin = 0.5,name='subsampled_small_2',path = './results/Zanne_clades_BAMM_small/')

subsampled3.small.files<-list.files('./results/Zanne_clades_BAMM_small/subsampled_small_3/results/',pattern='_out_')
lapply(subsampled3.small.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_small/subsampled_small_3/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_small/subsampled_small_3/results/event_data_GBIFdata.BAMM.subsampled_small_3_clade_1_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_3/results/event_data_GBIFdata.BAMM.subsampled_small_3_clade_2_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_3/results/event_data_GBIFdata.BAMM.subsampled_small_3_clade_3_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_3/results/event_data_GBIFdata.BAMM.subsampled_small_3_clade_4_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_3/results/event_data_GBIFdata.BAMM.subsampled_small_3_clade_5_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_3/results/event_data_GBIFdata.BAMM.subsampled_small_3_clade_6_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_3/results/event_data_GBIFdata.BAMM.subsampled_small_3_clade_7_50.txt'),burnin = 0.5,name='subsampled_small_3',path = './results/Zanne_clades_BAMM_small/')

subsampled4.small.files<-list.files('./results/Zanne_clades_BAMM_small/subsampled_small_4/results/',pattern='_out_')
lapply(subsampled4.small.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_small/subsampled_small_4/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_small/subsampled_small_4/results/event_data_GBIFdata.BAMM.subsampled_small_4_clade_1_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_4/results/event_data_GBIFdata.BAMM.subsampled_small_4_clade_2_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_4/results/event_data_GBIFdata.BAMM.subsampled_small_4_clade_3_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_4/results/event_data_GBIFdata.BAMM.subsampled_small_4_clade_4_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_4/results/event_data_GBIFdata.BAMM.subsampled_small_4_clade_5_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_4/results/event_data_GBIFdata.BAMM.subsampled_small_4_clade_6_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_4/results/event_data_GBIFdata.BAMM.subsampled_small_4_clade_7_50.txt'),burnin = 0.5,name='subsampled_small_4',path = './results/Zanne_clades_BAMM_small/')

subsampled5.small.files<-list.files('./results/Zanne_clades_BAMM_small/subsampled_small_5/results/',pattern='_out_')
lapply(subsampled5.small.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_small/subsampled_small_5/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_small/subsampled_small_5/results/event_data_GBIFdata.BAMM.subsampled_small_5_clade_1_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_5/results/event_data_GBIFdata.BAMM.subsampled_small_5_clade_2_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_5/results/event_data_GBIFdata.BAMM.subsampled_small_5_clade_3_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_5/results/event_data_GBIFdata.BAMM.subsampled_small_5_clade_4_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_5/results/event_data_GBIFdata.BAMM.subsampled_small_5_clade_5_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_5/results/event_data_GBIFdata.BAMM.subsampled_small_5_clade_6_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_5/results/event_data_GBIFdata.BAMM.subsampled_small_5_clade_7_50.txt'),burnin = 0.5,name='subsampled_small_5',path = './results/Zanne_clades_BAMM_small/')

subsampled6.small.files<-list.files('./results/Zanne_clades_BAMM_small/subsampled_small_6/results/',pattern='_out_')
lapply(subsampled6.small.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_small/subsampled_small_6/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_small/subsampled_small_6/results/event_data_GBIFdata.BAMM.subsampled_small_6_clade_1_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_6/results/event_data_GBIFdata.BAMM.subsampled_small_6_clade_2_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_6/results/event_data_GBIFdata.BAMM.subsampled_small_6_clade_3_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_6/results/event_data_GBIFdata.BAMM.subsampled_small_6_clade_4_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_6/results/event_data_GBIFdata.BAMM.subsampled_small_6_clade_5_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_6/results/event_data_GBIFdata.BAMM.subsampled_small_6_clade_6_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_6/results/event_data_GBIFdata.BAMM.subsampled_small_6_clade_7_50.txt'),burnin = 0.5,name='subsampled_small_6',path = './results/Zanne_clades_BAMM_small/')

subsampled7.small.files<-list.files('./results/Zanne_clades_BAMM_small/subsampled_small_7/results/',pattern='_out_')
lapply(subsampled7.small.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_small/subsampled_small_7/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_small/subsampled_small_7/results/event_data_GBIFdata.BAMM.subsampled_small_7_clade_1_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_7/results/event_data_GBIFdata.BAMM.subsampled_small_7_clade_2_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_7/results/event_data_GBIFdata.BAMM.subsampled_small_7_clade_3_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_7/results/event_data_GBIFdata.BAMM.subsampled_small_7_clade_4_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_7/results/event_data_GBIFdata.BAMM.subsampled_small_7_clade_5_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_7/results/event_data_GBIFdata.BAMM.subsampled_small_7_clade_6_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_7/results/event_data_GBIFdata.BAMM.subsampled_small_7_clade_7_50.txt'),burnin = 0.5,name='subsampled_small_7',path = './results/Zanne_clades_BAMM_small/')

subsampled8.small.files<-list.files('./results/Zanne_clades_BAMM_small/subsampled_small_8/results/',pattern='_out_')
lapply(subsampled8.small.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_small/subsampled_small_8/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_small/subsampled_small_8/results/event_data_GBIFdata.BAMM.subsampled_small_8_clade_1_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_8/results/event_data_GBIFdata.BAMM.subsampled_small_8_clade_2_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_8/results/event_data_GBIFdata.BAMM.subsampled_small_8_clade_3_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_8/results/event_data_GBIFdata.BAMM.subsampled_small_8_clade_4_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_8/results/event_data_GBIFdata.BAMM.subsampled_small_8_clade_5_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_8/results/event_data_GBIFdata.BAMM.subsampled_small_8_clade_6_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_8/results/event_data_GBIFdata.BAMM.subsampled_small_8_clade_7_50.txt'),burnin = 0.5,name='subsampled_small_8',path = './results/Zanne_clades_BAMM_small/')

subsampled9.small.files<-list.files('./results/Zanne_clades_BAMM_small/subsampled_small_9/results/',pattern='_out_')
lapply(subsampled9.small.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_small/subsampled_small_9/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_small/subsampled_small_9/results/event_data_GBIFdata.BAMM.subsampled_small_9_clade_1_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_9/results/event_data_GBIFdata.BAMM.subsampled_small_9_clade_2_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_9/results/event_data_GBIFdata.BAMM.subsampled_small_9_clade_3_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_9/results/event_data_GBIFdata.BAMM.subsampled_small_9_clade_4_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_9/results/event_data_GBIFdata.BAMM.subsampled_small_9_clade_5_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_9/results/event_data_GBIFdata.BAMM.subsampled_small_9_clade_6_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_9/results/event_data_GBIFdata.BAMM.subsampled_small_9_clade_7_50.txt'),burnin = 0.5,name='subsampled_small_9',path = './results/Zanne_clades_BAMM_small/')

subsampled10.small.files<-list.files('./results/Zanne_clades_BAMM_small/subsampled_small_10/results/',pattern='_out_')
lapply(subsampled10.small.files,function(x)analyse_BAMM_convergence(paste('./results/Zanne_clades_BAMM_small/subsampled_small_10/results/',x,sep=''),burnin=0.5))
eventdatafiles_burnin_merge_clades(event_data.files = c('./results/Zanne_clades_BAMM_small/subsampled_small_10/results/event_data_GBIFdata.BAMM.subsampled_small_10_clade_1_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_10/results/event_data_GBIFdata.BAMM.subsampled_small_10_clade_2_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_10/results/event_data_GBIFdata.BAMM.subsampled_small_10_clade_3_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_10/results/event_data_GBIFdata.BAMM.subsampled_small_10_clade_4_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_10/results/event_data_GBIFdata.BAMM.subsampled_small_10_clade_5_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_10/results/event_data_GBIFdata.BAMM.subsampled_small_10_clade_6_50.txt','./results/Zanne_clades_BAMM_small/subsampled_small_10/results/event_data_GBIFdata.BAMM.subsampled_small_10_clade_7_50.txt'),burnin = 0.5,name='subsampled_small_10',path = './results/Zanne_clades_BAMM_small/')

source('./R/compare_BAMMestimates_full_vs_subsampled.R')
rates.df.merge.small<-data_frame_compare_rates_full_vs_subsampled(folder.subsampled.event_data='./results/Zanne_clades_BAMM_small/event_data/')

full28k.vs.10ksmall.results<-replicates_phylolm_rho_trop_vs_temp_small(rates.df.merge=rates.df.merge.small)
#this has the results of phylolm
full28k.vs.10k.small.lms<-full28k.vs.10ksmall.results[[1]]
#check the regression coefficients for tropicality here
lapply(full28k.vs.10k.small.lms,function(x)summary(x)$coefficients)
lapply(full28k.vs.10k.small.lms,function(x)confint(x))
#calculate grand means and sds
grand.mean<-function(M, N) {weighted.mean(M, N)}
grand.sd<-function(S, M, N) {sqrt(weighted.mean(S^2+M^2,N)-weighted.mean(M, N)^2)}
sizes<-unlist(lapply(full28k.vs.10k.small.lms,function(x)x$n))
sds<-unlist(lapply(full28k.vs.10k.small.lms,function(x)summary(x)$coefficients[4,2]))
means<-unlist(lapply(full28k.vs.10k.small.lms,function(x)summary(x)$coefficients[4,1]))
mean(means)
grand.mean(M=means,N=sizes)
grand.sd(S=sds,M=means,N=sizes)
#ths has the results of rhos
full28k.vs.10k.small.rhos<-full28k.vs.10ksmall.results[[2]]
median(unlist(lapply(full28k.vs.10k.small.rhos,function(x)x$estimate)))
median(unlist(lapply(full28k.vs.10k.small.rhos,function(x)x$p.value)))
plot_replicates_phylolm_trop_vs_temp_small(rates.df.merge=rates.df.merge.small)
####plot rate comparisons in all 10 subsampled smalled_datasets####
pdf('./plots/FigS4_lm_BAMMrates_28ktree_vs_10subsamples_10ktrees_smalldataset.pdf',paper='a4r')
par(mfrow=c(2,5))
plot_replicates_phylolm_trop_vs_temp_small(rates.df.merge=rates.df.merge.small)
dev.off()

####7) FiSSE analyses####
####FiSSE for full biased dataset (binary tropical)####
source('./R/run_fisse/traitDependent_functions.R')
GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
GBIFdata$tropical<-0
GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
#subset table to BAMM species
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
tree<-read.tree("./raw_data/Qian_GBIFdata.tree")
traits<-GBIFdata.BAMM$tropical
names(traits)<-GBIFdata.BAMM$binomial

#FISSE.binary commented out (uses means instead of medians so not great for our distribution of rates)
#system.time(fisse.binary.tropical.full<- FISSE.binary(tree, traits,reps=100))
#lambda0 = 1.757793; lambda1 = 1.030903; p-value = 0.099
#saveRDS(fisse.binary.tropical.full,file='./results/fisse.binary.tropical.full.100results.RDS')

#this takes 1 hour
#FiSSE modified as in Faurby and Antonelli, uses medians instead of means
dir.create('./results/full_biased/fisse/')
system.time(fisse.binary.tropical.full.median<- FISSE.binary.median(tree, traits,reps=100))
saveRDS(fisse.binary.tropical.full.median,file='./results/full_biased/fisse/fisse.binary.tropical.full.median.100results.RDS')
#RESULT: lambda0 = 0.2884; lambda1 = 0.1821; p-value = 0.01

####FiSSE for subsampled dataset####
####run FiSSE subsampled datasets on cluster####
dir.create('./results/subsampled_unbiased/fisse/')
#copy './results/subsampled_unbiased/input_tables/
#mkdir './results/subsampled_unbiased/fisse/ on cluster
#run './R/run_fisse_unbiased_replicates_cluster.R on cluster
#create 100 conscriptoR commands to rrun_fisse_unbiased_replicates_cluster on hydrogen
for (i in seq(from=1,to=97,by=4)){
  write(x=paste('/scripts/conscriptoR ~/ldg_plants/run_fisse_unbiased_replicates_cluster.R',i,i+3),append=T,file='./R/conscriptor_fisse_unbiased_replicates_cluster_100_by4.txt')
}
#copy ./results/subsampled_unbiased/fisse/fisse.binary.subsampled_*_median.RDS to local ./results/subsampled_unbiased/fisse/

####check results of FiSSE for subsampled dataset####
#check results of fisse local subsampled datasets median
fisse.binary.subsampled.files<-list.files('./results/subsampled_unbiased/fisse/',pattern='fisse.binary.subsampled_.+_median.RDS')
fisse.binary.subsampled.list<-lapply(fisse.binary.subsampled.files,function(x)readRDS(file=paste('./results/subsampled_unbiased/fisse/',x,sep='')))
#no significant differences
hist(unlist(lapply(fisse.binary.subsampled.list,function(x)x$pval)))
median(unlist(lapply(fisse.binary.subsampled.list,function(x)x$pval)))
#median pvalue=0.168
#lambda0 and lambda1 look different, but the difference is not significantly different from the neutral
plot(density(unlist(lapply(fisse.binary.subsampled.list,function(x)x$lambda0))),xlim=c(0.1,0.3),main='fisse.unbiased',col='blue',ylab='',xlab='lambda',yaxt='n',ylim=c(0,200))
lines(density(unlist(lapply(fisse.binary.subsampled.list,function(x)x$lambda1))),col='red',ylab='',xlab='',yaxt='n')

plot(density(unlist(lapply(fisse.binary.subsampled.list,function(x)x$lambda1-x$lambda0))),main='fisse.unbiased',ylab='',xlab='delta.lambda(lambda1-lambda0)',xlim=c(-0.1,0.1),ylim=c(0,200),yaxt='n')
lines(density(unlist(lapply(fisse.binary.subsampled.list,function(x)x$null_mean_diff))),main='fisse.unbiased',ylab='',xlab='median.delta.lambda(lambda1-lambda0)',yaxt='n',lty=2)

####8)DR rates comparison: 28k (biased) tree vs 10k (unbiased) trees####
dir.create('./results/DRrates_compare_full_subsampled/')
source('./R/run_fisse/traitDependent_functions.R')

####get DRrates for 28k tree####
tree.full<-read.tree('./raw_data/BAMM_Species_tree_noseed.tree')
DR.full.tree <- DR_statistic(tree.full)
DR.full.tree.df<-as.data.frame(names(DR.full.tree),DR.full.tree)
rownames(DR.full.tree.df)<-NULL
write.table(DR.full.tree.df,file='./results/DRrates_compare_full_subsampled/DR.30ktree_table.txt',sep='\t',quote=F,row.names=F)

####get rates for subsampled replicates (for 10)####
#this takes a few hours, maybe worth running on the cluster
table.list<-list.files('./results/subsampled_unbiased/input_tables/',pattern='GBIFdata.BAMM.subsampled_')
tables<-lapply(c(1:10),function(x) read.table(paste('./results/subsampled_unbiased/input_tables/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
tree.full<-read.tree('./raw_data/BAMM_Species_tree_noseed.tree')
tree.subsampled.list<-lapply(tables,function(x) drop.tip(tree.full,setdiff(tree.full$tip.label,x$binomial)))
DR.subsampled.tree.list <- lapply(tree.subsampled.list,function(x) DR_statistic(x))
DR.subsampled.tree.list.df<-lapply(DR.subsampled.tree.list,function(x) as.data.frame(names(x),x))
DR.subsampled.tree.list.df<-lapply(DR.subsampled.tree.list.df,function(x) {rownames(x)<-NULL;return(x)})
for(i in c(1:10)){
  write.table(DR.subsampled.tree.list.df[[i]],file=paste('./results/DRrates_compare_full_subsampled/DR.10ktree_',i,'_table.txt',sep=''),sep='\t',quote=F,row.names=F)  
}

####plot DR comparisons Fig S6, Fig 6####
source('./R/compare_DRrates_full_subsampled.R')

####9) clade analyses####
source('./R/clade_based_analyses.R')
dir.create('./results/clade_analyses/')
####generate tree + table for input####
inputlist<-prepare_tree_and_table_dataset()
#inputlist[[1]]<-tree (Qian_dropped)
#inputlist[[2]]<-table (Qian_Lookup_Angios_TPL_GBIF)
tree<-inputlist[[1]]
Qian_Lookup_Angios_TPL_GBIF<-inputlist[[2]]

####get and save phylogroups####
#for 0-18 by 2
for (i in seq(from=0,to=18,by=2)){
  cat(i,'\n')
  minage<-i
  maxage<-i+2
  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=5,ncores=1)
  saveRDS(phylogroups,file=paste('./results/clade_analyses/phylogroups_',minage,'_',maxage,'_3size.RDS',sep=''))
}

#for 0-20 by 5
for (i in seq(from=0,to=20,by=5)){
  cat(i,'\n')
  minage<-i
  maxage<-i+5
  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=5,ncores=1)
  saveRDS(phylogroups,file=paste('./results/clade_analyses/phylogroups_',minage,'_',maxage,'_5size.RDS',sep=''))
}

#for 0-20 by 4
for (i in seq(from=0,to=20,by=4)){
  cat(i,'\n')
  minage<-i
  maxage<-i+4
  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=5,ncores=1)
  saveRDS(phylogroups,file=paste('./results/clade_analyses//phylogroups_',minage,'_',maxage,'_5size.RDS',sep=''))
}

####run RPANDA + save results####
#for 0-20 by 5
for (i in seq(from=0,to=20,by=5)){
  cat(i,'\n')
  minage<-i
  maxage<-i+5
  run_clades_RPANDAmodelave_pgls_savedphylogroups_and_save(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.5,GBIF.sampling = 0.5)
}

#for 0-20 by 4
for (i in seq(from=0,to=20,by=4)){
  cat(i,'\n')
  minage<-i
  maxage<-i+4
  run_clades_RPANDAmodelave_pgls_savedphylogroups_and_save(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.5,GBIF.sampling = 0.5)
}

#for 0-18 by 2
for (i in seq(from=0,to=18,by=2)){
  cat(i,'\n')
  minage<-i
  maxage<-i+2
  run_clades_RPANDAmodelave_pgls_savedphylogroups_and_save(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.5,GBIF.sampling = 0.5)
}

####plot results####
#for 0-20 by 4
for (i in seq(from=0,to=20,by=4)){
  cat(i,'\n')
  minage<-i
  maxage<-i+4
  run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave_weightlatdata(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.5,GBIF.sampling = 0.5,loadRPANDA = T)
  run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave_weightlatdata(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.7,GBIF.sampling = 0.5,loadRPANDA = T)
  run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave_weightlatdata(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.9,GBIF.sampling = 0.5,loadRPANDA = T)
}
#for 0-20 by 5
for (i in seq(from=0,to=20,by=5)){
  cat(i,'\n')
  minage<-i
  maxage<-i+5
  run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave_weightlatdata(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.5,GBIF.sampling = 0.5,loadRPANDA = T)
  run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave_weightlatdata(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.7,GBIF.sampling = 0.5,loadRPANDA = T)
  run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave_weightlatdata(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.9,GBIF.sampling = 0.5,loadRPANDA = T)
}
#for 0-20 by 2
for (i in seq(from=0,to=18,by=2)){
  cat(i,'\n')
  minage<-i
  maxage<-i+2
  run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave_weightlatdata(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.5,GBIF.sampling = 0.5,loadRPANDA = T)
  run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave_weightlatdata(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.7,GBIF.sampling = 0.5,loadRPANDA = T)
  run_clades_RPANDA_pgls_savedphylogroups_abslat_stricttropthreshold_log_modelave_weightlatdata(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=Qian_Lookup_Angios_TPL_GBIF,strict.tropical.threshold = 0.9,GBIF.sampling = 0.5,loadRPANDA = T)
}

####10) assess the relationship of missing taxa and latitude####
source('./R/assessing_missing_bias.R')
GBIFdata.BAMM.rates<-get_samplingfractionsBAMM()
#can the sampling fraction be predicted by the median latitude?
pdf('./plots/FigS8_medianlatitude_vs_BAMMsamplingfraction.pdf',paper='a4r')
par(mfrow=(c(2,2)))
plot(abs(GBIFdata.BAMM.rates$Median.Latitude),GBIFdata.BAMM.rates$BAMM.fraction*100,pch=16,main='species',ylab='BAMM sampling fraction (%)',xlab='species median latitude',cex=.7,col=rgb(0,0,0,0.05),yaxt='n')
axis(2,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),las=2)
lm.lat.fraction.species<-lm(BAMM.fraction*100~abs(Median.Latitude),data=GBIFdata.BAMM.rates)
abline(lm.lat.fraction.species,col='red',lwd=2,lty=2)
dev.off()
summary(lm.lat.fraction.species)
confint(lm.lat.fraction.species)
#slope=0.05903 ***
#for every 10 degree increase in latitude = 0.5903% increase in sampling fraction

#is lambda predicted by sampling fraction?
GBIFdata.BAMM.rates$BAMM.fraction.percentage<-GBIFdata.BAMM.rates$BAMM.fraction*100
pdf('./plots/FigS9_lambda_vs_BAMMsamplingfraction.pdf',paper='a4r')
par(mfrow=(c(2,2)))
plot(GBIFdata.BAMM.rates$BAMM.fraction.percentage,log(GBIFdata.BAMM.rates$lambda.avg),pch=16,main='species',ylab='lambda (lineages/myr)',xlab='BAMM sampling fraction (%)',cex=.7,col=rgb(0,0,0,0.05),yaxt='n')
axis(2,at=log(c(0.05,0.25,1,5)),labels=c(0.05,0.25,1,5),las=2)
#looks like higher lambdas at lower sampling fractions
#linear regression of lambda and median.latitude+fraction
#lm.lambda.lat.fraction.species<-lm(log(lambda.avg)~abs(Median.Latitude)+BAMM.fraction.percentage,data=GBIFdata.BAMM.rates)
lm.lambda.lat.fraction.species<-lm(log(lambda.avg)~BAMM.fraction.percentage,data=GBIFdata.BAMM.rates)
lm.lambda.lat.fraction.species.poly2<-lm(log(lambda.avg)~poly(abs(BAMM.fraction.percentage),2,raw=T),data=GBIFdata.BAMM.rates)
abline(lm.lambda.lat.fraction.species,col='red',lwd=2,lty=2)
#abline(log(lm.lat.fraction.species$coefficients[1]),log(lm.lat.fraction.species$coefficients[2]))
summary(lm.lambda.lat.fraction.species)
dev.off()
#slope for abs.Median.Lattude  = 0.0159
#slope for BAMM.fraction.percentage = -0.057
#for every 10% increase in sampling = -0.56 (decrease) in speciation rate







