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
dir.create('./output/')
dir.create('./output/trees/')

####raw data should have: GBIFdatasummary.csv (from Mounce et al Nat Plants), all_clades_50million_1000samples_latitudedata.RDS (eventdata from BAMM)
#QianTree.txt, BAMM_Species_tree_noseed.tree
####./raw_data/BAMM_Zanneclades/ should have clade_1_sampling.txt clade_1.tree (...) clade_7_sampling.txt, clade_7.tree
####./raw_data/BAMM_Zanneclades/control_files/ has all control_files.txt
####./R/run_fisse/traitDependent_functions.R should be modified to include FISSE.binary.median

################
####1) prepare BAMM files for GBOTB tree####
source('./R/GBOBT_BAMM.R')
#this will create GBOTB trees in raw_data
####TO DO:add the Div_edata here and clean below####

####subtreeBAMM on BAMM object
treeGBOTB<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
GBIFdata$tropical<-0
GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
GBIFdata<-GBIFdata[-which(duplicated(GBIFdata$binomial)),]
GBIFdata<-unique(GBIFdata)
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%treeGBOTB$tip.label,]

subtreeGBIF<-subtreeBAMM(subsampled.Div_edata, tips = GBIFdata.BAMM$binomial, node = NULL)
saveRDS(subtreeGBIF,file=paste0(path,name,'/results/all_clades_50million_1000samples_latitudedata_GBOTB.RDS'))
#####
#####generate the table of results, will be stored in ./output/tables/####
source('./R/get_GBIFBAMMrates_table.R')

####GBOTB ANALYSES (MAIN)####

####2) run STRAPP with full biased dataset (binary, continuous, latband, absolute lat band)####
####prepare dataset for STRAPP full####
folder.strapp.GBOTB.full<-'./results/full_biased_GBOTB/'
dir.create(folder.strapp.GBOTB.full)
source('./R/strapp_functions.R')
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata_GBOTB.RDS')
GBIFdata.BAMM<-read.table('./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',sep='\t',header=T,stringsAsFactors = F)
####strapp runs for full dataset####
tropical.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =GBIFdata.BAMM,mode = 'binary.tropical')
write.table(tropical.strapp,file=paste0(folder.strapp.GBOTB.full,'tropical.strapp.fulldataset.GBOTB.txt'),sep='\t',quote=F,row.names=F)
#medianlatitude.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =GBIFdata.BAMM,mode = 'continuous.medianlatitude')
#write.table(medianlatitude.strapp,paste0(folder.strapp.GBOTB.full,'medianlatitude.strapp.fulldataset.GBOTB.txt'),sep='\t',quote=F,row.names=F)
abs.medianlatitude.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =GBIFdata.BAMM,mode = 'continuous.abs.medianlatitude')
write.table(abs.medianlatitude.strapp,file=paste0(folder.strapp.GBOTB.full,'abs.medianlatitude.strapp.fulldataset.GBOTB.txt'),sep='\t',quote=F,row.names=F)
#latitudinalband.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =GBIFdata.BAMM,mode = 'latitudinalband')
#write.table(latitudinalband.strapp,file=paste0(folder.strapp.GBOTB.full,'latitudinalband.strapp.fulldataset.GBOTB.txt'),sep='\t',quote=F,row.names=F)
abs.latitudinalband.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object,trait.table =GBIFdata.BAMM,mode = 'abs.latitudinalband')
write.table(abs.latitudinalband.strapp,file=paste0(folder.strapp.GBOTB.full,'abs.latitudinalband.strapp.fulldataset.GBOTB.txt'),sep='\t',quote=F,row.names=F)

####prepare dataset for strict tropical STRAPP full
folder.strapp.GBOTB.strict<-'./results/full_biased_strict_GBOTB/'

dir.create(folder.strapp.GBOTB.strict)
source('./R/strapp_functions.R')
GBIFdata.BAMM<-read.table('./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',sep='\t',header=T,stringsAsFactors = F)
GBIFdata.BAMM$strict.tropical<-NA
#strict tropical 1 = abs(max latitude)<23.5 & abs(min latitude)<23.5  & median.latitude <23.5
GBIFdata.BAMM[abs(GBIFdata.BAMM$Max.Latitude)<=23.5&abs(GBIFdata.BAMM$Min.Latitude)<=23.5&abs(GBIFdata.BAMM$Median.Latitude)<=23.5,'strict.tropical']<-1
#strict tropical 0 (strict temperate)
GBIFdata.BAMM[abs(GBIFdata.BAMM$Max.Latitude)>23.5&abs(GBIFdata.BAMM$Min.Latitude)>23.5&abs(GBIFdata.BAMM$Median.Latitude)>23.5,'strict.tropical']<-0
#2: not strict tropical or strict temperate species
GBIFdata.BAMM$strict.tropical[is.na(GBIFdata.BAMM$strict.tropical)]<-2
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata_GBOTB.RDS')
####strapp runs with full dataset -  strict.tropical, only strict.tropical vs strict.temperate species####
GBIFdata.BAMM.strict.tropical<-GBIFdata.BAMM[!GBIFdata.BAMM$strict.tropical==2,]
#need to subset BAMMobject to strict trops and strict temps
BAMM.object.strict.tropical<-subtreeBAMM(ephy = BAMM.object,tips = GBIFdata.BAMM.strict.tropical$binomial)
saveRDS(BAMM.object.strict.tropical,'./results/full_biased_strict_GBOTB/all_clades_50million_1000samples_latitudedata_GBOTB_stricttroptemp.RDS')
strict.tropical.binary.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.strict.tropical,trait.table =GBIFdata.BAMM.strict.tropical,mode = 'strict.binary.tropical')
write.table(strict.tropical.binary.strapp,file='./results/full_biased_strict_GBOTB/strict.tropical.binary.strapp.GBOTB.fulldataset.txt',sep='\t',quote=F,row.names=F)

####3) Figure 1 Plots####
source('./R/plot_figs.R')
GBIFdata.BAMM.rates<-read.table('./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',header=T,sep='\t',stringsAsFactors = F)
pdf('./plots/Fig1_map_latitudinalbands_log_GBOTB.pdf',paper='a4r')
plot_Fig1_map_altbands(GBIFdata.BAMM.rates)
dev.off()
pdf('./plots/Fig1_boxplots_quantile_log_GBOTB.pdf',paper='a4r')
plot_Fig1_boxplot_troptempbinary_GBOTB(GBIFdata.BAMM.rates)
plot_Fig1_quantilerate_vs_lambdarate(GBIFdata.BAMM.rates)
dev.off()

pdf('./plots/Fig1_strapp_correlations_full_GBOTB.pdf')
par(pty='s')
plot_Fig1_strapp_samples_from_file(correlationfile = './results/full_biased_GBOTB/strapp_lambda_absMedianLatitude_GBOTB_all.txt',parameter = 'lambda',tablefile='./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt')
dev.off()

pdf('./plots/Fig1_strappdensity_full_GBOTB.pdf',paper='a4') 
par(pty='s')
plot_Fig1_strapp_density()
dev.off()

####4) TABLE S1####
###run STRAPP with REDUCED full biased dataset (binary, continuous, latband, absolute lat band)
folder.strapp.GBOTB.reduced<-'./results/full_biased_reduced_GBOTB/'

dir.create(folder.strapp.GBOTB.reduced)
source('./R/strapp_functions.R')

GBIFdata.BAMM<-read.table('./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',sep='\t',header=T,stringsAsFactors = F)
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata_GBOTB.RDS')
#proportions of data points per latitudinal band
sort(round(table(GBIFdata.BAMM$latitudinal.band)/nrow(GBIFdata.BAMM),4))
####discard bands: -60,-50,70,80,90 (<1% of total points)####
GBIFdata.BAMM.reduced<-GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band!='-60'&GBIFdata.BAMM$latitudinal.band!='-50'&GBIFdata.BAMM$latitudinal.band!='70'&GBIFdata.BAMM$latitudinal.band!='80'&GBIFdata.BAMM$latitudinal.band!='90',]
#create new BAMM object and save
BAMM.object.lat_reduced<-subtreeBAMM(BAMM.object,tips = GBIFdata.BAMM.reduced$binomial)
saveRDS(BAMM.object.lat_reduced,file=paste0(folder.strapp.GBOTB.reduced,'all_clades_50million_1000samples_latitudedata_GBOTB_reduced.RDS'))
####strapp runs for full reduced dataset####
tropical.reduced.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.lat_reduced,trait.table =GBIFdata.BAMM.reduced,mode = 'binary.tropical')
write.table(tropical.reduced.strapp,file=paste0(folder.strapp.GBOTB.reduced,'tropical.strapp.GBOTB.fulldatasetreduced.txt'),sep='\t',quote=F,row.names=F)
#medianlatitude.reduced.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.lat_reduced,trait.table =GBIFdata.BAMM.reduced,mode = 'continuous.medianlatitude')
#write.table(medianlatitude.reduced.strapp,file=paste0(folder.strapp.GBOTB.reduced,'medianlatitude.strapp.fulldatasetreduced.txt'),sep='\t',quote=F,row.names=F)
abs.medianlatitude.reduced.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.lat_reduced,trait.table =GBIFdata.BAMM.reduced,mode = 'continuous.abs.medianlatitude')
write.table(abs.medianlatitude.reduced.strapp,file=paste0(folder.strapp.GBOTB.reduced,'abs.medianlatitude.strapp.fulldatasetreduced.txt'),sep='\t',quote=F,row.names=F)
#latitudinalband.reduced.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.lat_reduced,trait.table =GBIFdata.BAMM.reduced,mode = 'latitudinalband')
#write.table(latitudinalband.reduced.strapp,file=paste0(folder.strapp.GBOTB.reduced,'latitudinalband.strapp.fulldatasetreduced.txt'),sep='\t',quote=F,row.names=F)
abs.latitudinalband.reduced.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.lat_reduced,trait.table =GBIFdata.BAMM.reduced,mode = 'abs.latitudinalband')
write.table(abs.latitudinalband.reduced.strapp,file=paste0(folder.strapp.GBOTB.reduced,'abs.latitudinalband.strapp.fulldatasetreduced.txt'),sep='\t',quote=F,row.names=F)

#####5) Figure S1####
source('./R/assessing_bias_dataset.R')
#this plots the proportion of trop and temp species in the GBIF dataset, and the proportion of species in each lat band
pdf('./plots/FigS1_full_dataset_BAMMdataset_biases_new_GBOTB.pdf',width=10,height=7)
par(mfrow=c(2,2))
par(pty='s')
plot_GBIFdata_proportions(GBIFdata.file='./raw_data/GBIFdatasummary.csv')
#this plots the proportion of trop and temp species in the GBIF BAMM dataset, and the proportion of species in each lat band
plot_GBIFdataBAMM_proportions(GBIFdata.file='./raw_data/GBIFdatasummary.csv',BAMM.object.path='./raw_data/all_clades_50million_1000samples_latitudedata_GBOTB.RDS')
#there's an obvious bias, with many more species in temperate regions being sequenced
dev.off()

####6) GBIFdataBAMM subsampled unbiased datasets####
####generate subsampled_unbiased datasets 30k
dir.create('./results/subsampled_unbiased_GBOTB/')
dir.create('./results/subsampled_unbiased_GBOTB/input_tables/')
GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
GBIFdata$tropical<-0
GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
sampling.proportions.GBIF<-table(GBIFdata$tropical)/nrow(GBIFdata)

#load tree
GBOTB.bif.angios_BAMM.tree<-read.tree('GBOTB.bif.angios_BAMM.tree')
#subset table to BAMM species
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%GBOTB.bif.angios_BAMM.tree$tip.label,]
GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]

#size.dataset is the size of the subsampled datasets
size.dataset<-30000
for (replicate in c(1:100)){
  GBIFdata.BAMM.subsampled<-GBIFdata.BAMM[c(sample(which(GBIFdata.BAMM$tropical==0),size = round(sampling.proportions.GBIF[1]*size.dataset)),sample(which(GBIFdata.BAMM$tropical==1),size = round(sampling.proportions.GBIF[2]*size.dataset))),]
  #output full table
  write.table(GBIFdata.BAMM.subsampled,file=paste('./results/subsampled_unbiased_GBOTB/input_tables/GBIFdata.BAMM.GBOTB.subsampled_',replicate,'_table.txt',sep=''),sep='\t',quote=F,row.names=F)
}
#copy the GBIFdata.BAMM.GBOTB.subsampled_ files in ./results/subsampled_unbiased_GBOTB/input_tables/ to the cluster
####run STRAPP 100 subsampled unbiased full dataset (binary, continuous, latband, absolute lat band) x replicates####
#run run_strapp_subsample_replicates_cluster.R in cluster, 100 times
#e.g.: /scripts/conscriptoR ~/ldg_plants/run_strapp_subsample_replicates_cluster.R 1
#e.g.: /scripts/conscriptoR /ldg_plants/run_strapp_subsample_replicates_cluster.R 2

#then copy the outputs in cluster to ./results/strapp_unbiased/
#and sort in folders (abs.latitudinalband, abs.medianlatitude.continuous, latitudinalband, medianlatitude.continuous, tropical.binary)
dir.create('./results/subsampled_unbiased_GBOTB/abs.latitudinalband/')
dir.create('./results/subsampled_unbiased_GBOTB/abs.medianlatitude/')
dir.create('./results/subsampled_unbiased_GBOTB/tropical/')

####plot STRAPP summary results of the abs.medianlatitude, tropical.binary and abs.latitudinalband for subsampled datasets####
#absolute median latitude strapp
pdf('./plots/Fig3_strapp_results_GBIFbamm_10k_subsampledreplicates.pdf',width=10,height=7)
par(mfrow=c(2,3))
blue<-brewer.pal(11,'RdYlBu')[9]
red<-brewer.pal(11,'RdYlBu')[2]



####7)Figure 3 plots####
source('./R/plot_figs.R')
#for GBOTB#
pdf('./plots/Fig3_map_all_subsampled30k_latitudinalbands_log_medians_unbiased_summary_GBOTB.pdf',paper='a4r')
plot_Fig2_latitudinalband_boxplots_subsampled_unbiased_summary_GBOTB()
dev.off()
pdf('./plots/Fig3_all_subsampled30k_quantilerate_vs_lambdarate_full_GBOTB.pdf')
plot_Fig2_allsubsampled_quantilerate_vs_lambda_GBOTB()
dev.off()
pdf('./plots/Fig3_subsampled_unbiased_GBOTB_strapp_densities.pdf')
par(pty='s')
plot_Fig2_strapp_densities(path.strappobjects = './results/subsampled_unbiased_GBOTB/')
dev.off()
pdf('./plots/Fig3_subsampled_unbiased_GBOTB_strapp_correlations.pdf')
par(pty='s')
plot_Fig2_strapp_samples_from_files(path.correlationfile='./results/subsampled_unbiased_GBOTB/',parameter='lambda',tablefile='./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt')
dev.off()

#####8) compare BAMM rate estimates of GBOTB tree vs unbiased trees####
source('./R/BAMM_Zanneclades_subsampled.R')
lapply(c(1:10),function(x) BAMM_GBOTBclades_subsampled_unbiased(replicate=x))
###run BAMM 10 replicates and store in ./results/subsampled_unbiased_1/results/ ./results/subsampled_unbiased_2/results/ etc
###analyse the BAMM runs, check convergence and merge all clades into a single analysis for each replicate
###this will print PDFs to check convergence of the runs, see burnin settings etc in BAMM_check_GBOTB_clades_unbiased.R
source('./R/BAMMoutput_subsampled_replicates_analyse.R')
source('./R/BAMM_check_GBOTB_clades_unbiased.R')

#compare rate differences full vs unbiased
source('./R/compare_BAMMestimates_full_vs_subsampled.R')
rates.df.merge.unbiased<-data_frame_compare_rates_full_vs_subsampled_GBOTB(folder.subsampled.event_data='./results/GBOTB_clades_BAMM_unbiased/event_data/')
write.table(rates.df.merge.unbiased,file='./results/GBOTB_clades_BAMM_unbiased/rates_df_merge_unbiased_table.txt',quote=F,row.names=F,sep='\t')

rates.df.merge.unbiased<-read.table(file='./results/subsampled_unbiased_GBOTB//rates_df_merge_unbiased_table.txt',header=T,sep='\t')
full60k.vs.30kunbiased.results<-replicates_phylolm_rho_trop_vs_temp_unbiased_GBOTB(rates.df.merge=rates.df.merge.unbiased)
#this has the results of phylolm
full60k.vs.30k.unbiased.lms<-full60k.vs.30kunbiased.results[[1]]
#check the regression coefficients for tropicality here
lapply(full60k.vs.30k.unbiased.lms,function(x)summary(x)$coefficients)
lapply(full60k.vs.30k.unbiased.lms,function(x)confint(x))
#calculate grand means and sds
grand.mean<-function(M, N) {weighted.mean(M, N)}
grand.sd<-function(S, M, N) {sqrt(weighted.mean(S^2+M^2,N)-weighted.mean(M, N)^2)}
sizes<-unlist(lapply(full60k.vs.30k.unbiased.lms,function(x)x$n))
sds<-unlist(lapply(full60k.vs.30k.unbiased.lms,function(x)summary(x)$coefficients[4,2]))
means<-unlist(lapply(full60k.vs.30k.unbiased.lms,function(x)summary(x)$coefficients[4,1]))
mean(means)
grand.mean(M=means,N=sizes)
grand.sd(S=sds,M=means,N=sizes)
#ths has the results of rhos
full60k.vs.30k.unbiased.rhos<-full60k.vs.30kunbiased.results[[2]]
median(unlist(lapply(full60k.vs.30k.unbiased.rhos,function(x)x$estimate)))
median(unlist(lapply(full60k.vs.30k.unbiased.rhos,function(x)x$p.value)))
plot_replicates_phylolm_trop_vs_temp_unbiased_GBOTB(rates.df.merge=rates.df.merge.unbiased)

####9) Figure 2 plots: plot rate comparisons in all 10 subsampled unbiaseded_datasets####
pdf('./plots/FigF2_lm_BAMMrates_60ktree_vs_10subsamples_10ktrees_unbiased_log_GBOTB.pdf',paper='a4r')
par(mfrow=c(2,5))
plot_replicates_phylolm_trop_vs_temp_unbiased_log_GBOTB(rates.df.merge=rates.df.merge.unbiased)
dev.off()

####10) Supplementary Note 2: DR analyses####
source('./R/DRmetric.R')
#get DR values for GBOTB_GBIF tree and write to tablle
tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
GBIFdata.BAMM<-read.table('./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',sep='\t',header=T,stringsAsFactors = F)
tree<-drop.tip(tree,setdiff(tree$tip.label,GBIFdata.BAMM$binomial))
DR<-DR_statistic(tree)
DR.df<-as.data.frame(DR)
DR.df$species<-rownames(DR.df)
rownames(DR.df)<-NULL
GBIFdata.BAMM.DR<-merge(GBIFdata.BAMM,DR.df,by.x='binomial',by.y='species')
write.table(GBIFdata.BAMM.DR,file='./output/GBIFdata_GBOTB_BAMM_DR_rates_table.txt',quote=F,sep='\t',row.names = F)

####11) Fig S9: DR figure####
GBIFdata.BAMM.DR.rates<-read.table('./GBIFdata_GBOTB_BAMM_DR_rates_table.txt',header=T,sep='\t',stringsAsFactors = F)
source('./R/plot_figs.R')
pdf('./plots/FigS9_DR_map_latitudinalbands_log_GBOTB.pdf',paper='a4r')
plot_Fig1_map_altbands_DR(GBIFdata.BAMM.DR.rates)
dev.off()
pdf('./plots/Fig1_DR_boxplots_quantile_log_GBOTB.pdf',paper='a4r')
plot_Fig1_DRboxplot_troptempbinary_GBOTB(GBIFdata.BAMM.DR.rates)
plot_Fig1_quantilerate_vs_DRrate(GBIFdata.BAMM.DR.rates)
dev.off()

####12) Fig S4####
###assess the relationship of missing taxa and latitude for GBOTB
source('./R/assessing_missing_bias.R')
GBIFdata.BAMM.rates<-read.table('./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',header=T,sep='\t',stringsAsFactors = F)
#can the sampling fraction be predicted by the median latitude?
pdf('./plots/FigS4_medianlatitude_vs_BAMMsamplingfraction_GBOTB.pdf',paper='a4r')
par(mfrow=(c(2,2)))
plot(abs(GBIFdata.BAMM.rates$Median.Latitude),GBIFdata.BAMM.rates$BAMM.fraction*100,pch=16,main='species',ylab='BAMM sampling fraction (%)',xlab='species median latitude',cex=.7,col=rgb(0,0,0,0.05),yaxt='n')
axis(2,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),las=2)
lm.lat.fraction.species<-lm(BAMM.fraction*100~abs(Median.Latitude),data=GBIFdata.BAMM.rates)
abline(lm.lat.fraction.species,col='red',lwd=2,lty=2)
dev.off()
summary(lm.lat.fraction.species)
confint(lm.lat.fraction.species)
#slope=0.200979 ***
#for every 10 degree increase in latitude = 2.0098% increase in sampling fraction

####13) Fig S5####
###assess the relationship of missing taxa and lambda estimates for GBOTB
#is lambda predicted by sampling fraction?
source('./R/assessing_missing_bias.R')
GBIFdata.BAMM.rates<-read.table('./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',header=T,sep='\t',stringsAsFactors = F)
GBIFdata.BAMM.rates$BAMM.fraction.percentage<-GBIFdata.BAMM.rates$BAMM.fraction*100
pdf('./plots/FigS5_lambda_vs_BAMMsamplingfraction_GBOTB.pdf',paper='a4r')
par(mfrow=(c(2,2)))
plot(GBIFdata.BAMM.rates$BAMM.fraction.percentage,log(GBIFdata.BAMM.rates$lambda.avg),pch=16,main='species',ylab='lambda (lineages/myr)',xlab='BAMM sampling fraction (%)',cex=.7,col=rgb(0,0,0,0.02),yaxt='n')
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
#slope for BAMM.fraction.percentage = -0.0207
#for every 10% increase in sampling = -0.21 (decrease) in speciation rate


####14) GBIFdataBAMM subsampled extremebias datasets####
####generate the replicate datasets
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
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata_GBOTB.RDS')
#subset table to BAMM species
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
#remove duplicates in GBIFdata (6 species with same info but two rows - one with present in garden, one with absent in garden)#
GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
GBIFdata.BAMM<-unique(GBIFdata.BAMM)
GBIFdata.BAMM$latitudinal.band<-NA
GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))

#create folder to store the BAMM replicates
dir.create('./results/subsampled_extremebias_GBOTB/')
dir.create('./results/subsampled_extremebias_GBOTB/input_tables/')

#size.dataset is the size of the subsampled datasets
size.dataset<-30000
for (replicate in c(1:10)){
  GBIFdata.BAMM.subsampled.extremebias<-GBIFdata.BAMM[c(sample(which(GBIFdata.BAMM$tropical==0),size = round(sampling.proportions.extreme.bias[1]*size.dataset)),sample(which(GBIFdata.BAMM$tropical==1),size = round(sampling.proportions.extreme.bias[2]*size.dataset))),]
  write.table(GBIFdata.BAMM.subsampled.extremebias,file=paste('./results/subsampled_extremebias_GBOTB/input_tables/GBIFdata.BAMM.GBOTB.subsampled_extremebias_',replicate,'_table.txt',sep=''),sep='\t',quote=F,row.names=F)
}
#copy the GBIFdata.BAMM.GBOTB.subsampled_ files in ./results/subsampled_extremebias_GBOTB/input_tables/ to the cluster
####run STRAPP 100 subsampled extremebias full dataset (binary, continuous, latband, absolute lat band) x replicates####
#run run_strapp_subsample_replicates_cluster.R in cluster, 100 times
#e.g.: /scripts/conscriptoR ~/ldg_plants/run_strapp_subsample_replicates_cluster.R 1
#e.g.: /scripts/conscriptoR /ldg_plants/run_strapp_subsample_replicates_cluster.R 2

#then copy the outputs in cluster to ./results/strapp_extremebias/
#and sort in folders (abs.latitudinalband, abs.medianlatitude.continuous, latitudinalband, medianlatitude.continuous, tropical.binary)
dir.create('./results/subsampled_extremebias_GBOTB/abs.latitudinalband/')
dir.create('./results/subsampled_extremebias_GBOTB/abs.medianlatitude/')
dir.create('./results/subsampled_extremebias_GBOTB/tropical/')

####plot STRAPP summary results of the abs.medianlatitude, tropical.binary and abs.latitudinalband for subsampled datasets####
#absolute median latitude strapp
pdf('./plots/Fig3_strapp_results_GBIFbamm_10k_subsampledreplicates.pdf',width=10,height=7)
par(mfrow=c(2,3))
blue<-brewer.pal(11,'RdYlBu')[9]
red<-brewer.pal(11,'RdYlBu')[2]



####15)Figure S3 plots####
source('./R/plot_figs.R')
#for GBOTB#
pdf('./plots/FigS3_map_all_subsampled30k_latitudinalbands_log_medians_extremebias_summary_GBOTB.pdf',paper='a4r')
plot_Fig2_latitudinalband_boxplots_subsampled_extremebias_summary_GBOTB()
dev.off()
pdf('./plots/FigS3_all_subsampled30k_quantilerate_vs_lambdarate_full_GBOTB.pdf')
plot_Fig2_allsubsampled_quantilerate_vs_lambda_GBOTB()
dev.off()
pdf('./plots/FigS3_subsampled_extremebias_GBOTB_strapp_densities.pdf')
par(pty='s')
plot_Fig2_strapp_densities(path.strappobjects = './results/subsampled_extremebias_GBOTB/')
dev.off()
pdf('./plots/FigS3_subsampled_extremebias_GBOTB_strapp_correlations.pdf')
par(pty='s')
plot_Fig2_strapp_samples_from_files(path.correlationfile='./results/subsampled_extremebias_GBOTB/',parameter='lambda',tablefile='./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt')
dev.off()

#####16) compare BAMM rate estimates of GBOTB tree vs extremebias trees####
source('./R/BAMM_Zanneclades_subsampled.R')
lapply(c(1:10),function(x) BAMM_GBOTBclades_subsampled_extremebias(replicate=x))
###run BAMM 10 replicates and store in ./results/subsampled_extremebias_1/results/ ./results/subsampled_extremebias_2/results/ etc
###analyse the BAMM runs, check convergence and merge all clades into a single analysis for each replicate
###this will print PDFs to check convergence of the runs, see burnin settings etc in BAMM_check_GBOTB_clades_extremebias.R
source('./R/BAMMoutput_subsampled_replicates_analyse.R')
source('./R/BAMM_check_GBOTB_clades_extremebias.R')

#compare rate differences full vs extremebias
source('./R/compare_BAMMestimates_full_vs_subsampled.R')
rates.df.merge.extremebias<-data_frame_compare_rates_full_vs_subsampled_GBOTB(folder.subsampled.event_data='./results/GBOTB_clades_BAMM_extremebias/event_data/')
write.table(rates.df.merge.extremebias,file='./results/GBOTB_clades_BAMM_extremebias/rates_df_merge_extremebias_table.txt',quote=F,row.names=F,sep='\t')

rates.df.merge.extremebias<-read.table(file='./results/subsampled_extremebias_GBOTB//rates_df_merge_extremebias_table.txt',header=T,sep='\t')
full60k.vs.30kextremebias.results<-replicates_phylolm_rho_trop_vs_temp_extremebias_GBOTB(rates.df.merge=rates.df.merge.extremebias)
#this has the results of phylolm
full60k.vs.30k.extremebias.lms<-full60k.vs.30kextremebias.results[[1]]
#check the regression coefficients for tropicality here
lapply(full60k.vs.30k.extremebias.lms,function(x)summary(x)$coefficients)
lapply(full60k.vs.30k.extremebias.lms,function(x)confint(x))
#calculate grand means and sds
grand.mean<-function(M, N) {weighted.mean(M, N)}
grand.sd<-function(S, M, N) {sqrt(weighted.mean(S^2+M^2,N)-weighted.mean(M, N)^2)}
sizes<-unlist(lapply(full60k.vs.30k.extremebias.lms,function(x)x$n))
sds<-unlist(lapply(full60k.vs.30k.extremebias.lms,function(x)summary(x)$coefficients[4,2]))
means<-unlist(lapply(full60k.vs.30k.extremebias.lms,function(x)summary(x)$coefficients[4,1]))
mean(means)
grand.mean(M=means,N=sizes)
grand.sd(S=sds,M=means,N=sizes)
#ths has the results of rhos
full60k.vs.30k.extremebias.rhos<-full60k.vs.30kextremebias.results[[2]]
median(unlist(lapply(full60k.vs.30k.extremebias.rhos,function(x)x$estimate)))
median(unlist(lapply(full60k.vs.30k.extremebias.rhos,function(x)x$p.value)))
plot_replicates_phylolm_trop_vs_temp_extremebias_GBOTB(rates.df.merge=rates.df.merge.extremebias)

####18) Figure S2 plots: plot rate comparisons in all 10 subsampled extremebiased_datasets####
pdf('./plots/FigS2_lm_BAMMrates_60ktree_vs_10subsamples_10ktrees_extremebias_log_GBOTB.pdf',paper='a4r')
par(mfrow=c(2,5))
plot_replicates_phylolm_trop_vs_temp_extremebias_log_GBOTB(rates.df.merge=rates.df.merge.extremebias)
dev.off()

####19) prepare "small" datasets####
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
BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata_GBOTB.RDS')
#subset table to BAMM species
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
#remove duplicates in GBIFdata (6 species with same info but two rows - one with present in garden, one with absent in garden)#
GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
GBIFdata.BAMM<-unique(GBIFdata.BAMM)
GBIFdata.BAMM$latitudinal.band<-NA
GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
#create folder to store the BAMM replicates
dir.create('./results/subsampled_small_GBOTB/')
dir.create('./results/subsampled_small_GBOTB/input_tables/')
#size.dataset is the size of the subsampled datasets
#I'll take the proportion of plant species in GBOTB tree
#get trees that have 8% of the species in GBOTB and compare rates
size.dataset<-round((nrow(GBIFdata.BAMM)/sum(plant_lookup_table$number.of.accepted.species))*nrow(GBIFdata.BAMM))
for (replicate in c(1:10)){
  GBIFdata.BAMM.subsampled.small<-GBIFdata.BAMM[c(sample(which(GBIFdata.BAMM$tropical==0),size = round(sampling.proportions[1]*size.dataset)),sample(which(GBIFdata.BAMM$tropical==1),size = round(sampling.proportions[2]*size.dataset))),]
  write.table(GBIFdata.BAMM.subsampled.small,file=paste('./results/subsampled_small_GBOTB/input_tables/GBIFdata.BAMM.GBOTB.subsampled_small_',replicate,'_table.txt',sep=''),sep='\t',quote=F,row.names=F)
}
####split BAMM_subsampled_replicates into clades for Zanne tree analyses in BAMM
#generate BAMM input files for 10 replicates (control files are in './raw_data/BAMM_Zanneclades/control_files/)
source('./R/BAMM_Zanneclades_subsampled.R')
lapply(c(1:10),function(x) BAMM_GBOTBclades_subsampled_small(replicate=x))

source('./R/BAMMoutput_subsampled_replicates_analyse.R')
source('./R/BAMM_check_GBOTB_clades_small.R')

#compare rate differences full vs small
source('./R/compare_BAMMestimates_full_vs_subsampled.R')
rates.df.merge.small<-data_frame_compare_rates_full_vs_subsampled_GBOTB(folder.subsampled.event_data='./results/GBOTB_clades_BAMM_small/event_data/')
write.table(rates.df.merge.small,file='./results/GBOTB_clades_BAMM_small/rates_df_merge_small_table.txt',quote=F,row.names=F,sep='\t')

rates.df.merge.small<-read.table(file='./results/subsampled_small_GBOTB//rates_df_merge_small_table.txt',header=T,sep='\t')
full60k.vs.small.results<-replicates_phylolm_rho_trop_vs_temp_unbiased_GBOTB(rates.df.merge=rates.df.merge.unbiased)
#this has the results of phylolm
full60k.vs.small.lms<-full60k.vs.small.results[[1]]
#check the regression coefficients for tropicality here
lapply(full60k.vs.small.lms,function(x)summary(x)$coefficients)
lapply(full60k.vs.small.lms,function(x)confint(x))
#calculate grand means and sds
grand.mean<-function(M, N) {weighted.mean(M, N)}
grand.sd<-function(S, M, N) {sqrt(weighted.mean(S^2+M^2,N)-weighted.mean(M, N)^2)}
sizes<-unlist(lapply(full60k.vs.small.lms,function(x)x$n))
sds<-unlist(lapply(full60k.vs.small.lms,function(x)summary(x)$coefficients[4,2]))
means<-unlist(lapply(full60k.vs.small.lms,function(x)summary(x)$coefficients[4,1]))
mean(means)
grand.mean(M=means,N=sizes)
grand.sd(S=sds,M=means,N=sizes)
#ths has the results of rhos
full60k.vs.small.rhos<-full60k.vs.small.results[[2]]
median(unlist(lapply(full60k.vs.small.rhos,function(x)x$estimate)))
median(unlist(lapply(full60k.vs.small.rhos,function(x)x$p.value)))
plot_replicates_phylolm_trop_vs_temp_small_GBOTB(rates.df.merge=rates.df.merge.unbiased)

####20) Figure S6 plots: plot rate comparisons in all 10 subsampled small_datasets####
pdf('./plots/FigS6_lm_BAMMrates_60ktree_vs_10subsamples_10ktrees_small_log_GBOTB.pdf',paper='a4r')
par(mfrow=c(2,5))
plot_replicates_phylolm_trop_vs_temp_small_GBOTB(rates.df.merge=rates.df.merge.small)
dev.off()

####21) analyses with Zanne tree####
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
write.table(GBIFdata.BAMM,file='./results/GBIFdata.BAMM.Zanne_table.txt',sep='\t',quote=F,row.names = F)

####22)Figure S7 Plots####
source('./R/plot_figs.R')
#get plots for Figure S7 = Figure 1 with Zanne
pdf('./plots/FigS7_mainfigureZanne_map_latitudinalbands_log.pdf',paper='a4r')
plot_Fig1_map_bands()
dev.off()
plot_Fig1_boxplot_troptempbinary()
plot_Fig1_quantilerate_vs_lambdarate()
plot_Fig1_strapp_samples_from_file(correlationfile = './results/strapp_lambda_absMedianLatitude_allnew.txt',parameter = 'lambda',tablefile='./results/GBIFdata.BAMM.table.txt')
plot_Fig1_strapp_density()
dev.off()

####23) clade-based analyses GBOTB (RPANDA)####
source('./R/clade_based_analyses.R')
dir.create('./results/clade_analyses_GBOTB/')
####generate tree + table for input####
inputlist<-prepare_tree_and_table_dataset_GBOTB()
#inputlist[[1]]<-tree (Qian_dropped)
#inputlist[[2]]<-table (Qian_Lookup_Angios_TPL_GBIF)
tree<-inputlist[[1]]
GBOTB_Lookup_Angios_GBIF<-inputlist[[2]]
#get and save phylogroups
#for 0-20 by 4
for (i in seq(from=0,to=20,by=4)){
  cat(i,'\n')
  minage<-i
  maxage<-i+4
  phylogroups<-get.phylogroups(tree,minage=minage,maxage=maxage,mincladesize=5,ncores=1)
  saveRDS(phylogroups,file=paste('./results/clade_analyses_GBOTB//phylogroups_',minage,'_',maxage,'_5size.RDS',sep=''))
}
#run RPANDA and save results
#for 0-20 by 4
for (i in seq(from=0,to=20,by=4)){
  cat(i,'\n')
  minage<-i
  maxage<-i+4
  run_clades_RPANDAmodelave_pgls_savedphylogroups_and_save_GBOTB(tree=Qian_dropped,minage=minage,maxage=maxage,mincladesize=5,sampling=0.3,ncores=2,table=GBOTB_Lookup_Angios_GBIF,strict.tropical.threshold = 0.5,GBIF.sampling = 0.5)
}
####23)Figure S8 Plots####
source('./R/plot_correlations_throughtimeslices.R')

####24) GeoHiSSE####
###run fGeoHisse_model_run_cluster_ldgplants_mod.R on cluster
###warning: some models take 1-2 months to finish!

