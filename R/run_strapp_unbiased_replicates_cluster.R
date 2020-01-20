.libPaths('~/R/x86_64-pc-linux-gnu-library/3.5')
setwd('~/ldg_plants/')

library(BAMMtools)
library(plyr)

assign_equalwidth_band<-function(number){
  vector.limits<-c(rev(seq(from=-4.7,to=-62,by=-9.4)),seq(from=4.7,to=90,by=9.4))
  #this checks if the median.latitude is identical to any of the bounds
  if(number%in%vector.limits){
    return(vector.limits[which(vector.limits==number)])
  }
  sorted.vector.limits<-sort(c(number,vector.limits))
  position<-which(sorted.vector.limits==number)
  if(length(position)>1){
    cat(number,'\n')
  }
  if(position==1){
    return(vector.limits[1])
  }else if (position==length(sorted.vector.limits)){
    return(vector.limits[length(vector.limits)-1])
  }else{
    return(vector.limits[position-1]) 
  }
}

run_strapp_BAMM.object<-function(BAMM.object,trait.table,mode){
  if(mode=='binary.tropical'){
    tropical<-trait.table$tropical
    names(tropical)<-trait.table$binomial
    cat('running strapp lambda','\n')
    lambda.tropical<-traitDependentBAMM(BAMM.object, tropical, 1000, rate = "speciation",return.full = FALSE, method = "mann-whitney", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp mu','\n')
    mu.tropical<-traitDependentBAMM(BAMM.object, tropical, 1000, rate = "extinction",return.full = FALSE, method = "mann-whitney", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp netdiv','\n')
    netdiv.tropical<-traitDependentBAMM(BAMM.object, tropical, 1000, rate = "net diversification",return.full = FALSE, method = "mann-whitney", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    strapp.tropical.df<-as.data.frame(rbind(c(unlist(lambda.tropical$estimate),lambda.tropical$p.value,lambda.tropical$method,lambda.tropical$two.tailed,lambda.tropical$rate),c(unlist(mu.tropical$estimate),mu.tropical$p.value,mu.tropical$method,mu.tropical$two.tailed,mu.tropical$rate),c(unlist(netdiv.tropical$estimate),netdiv.tropical$p.value,netdiv.tropical$method,netdiv.tropical$two.tailed,netdiv.tropical$rate)))
    colnames(strapp.tropical.df)<-c(paste0('estimate.',names(unlist(lambda.tropical$estimate))),'p.value','method','two.tailed','rate')
    return(strapp.tropical.df)
  }else if(mode=='strict.binary.tropical'){
    tropical<-trait.table$strict.tropical
    names(tropical)<-trait.table$binomial
    cat('running strapp lambda','\n')
    lambda.tropical<-traitDependentBAMM(BAMM.object, tropical, 1000, rate = "speciation",return.full = FALSE, method = "mann-whitney", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp mu','\n')
    mu.tropical<-traitDependentBAMM(BAMM.object, tropical, 1000, rate = "extinction",return.full = FALSE, method = "mann-whitney", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp netdiv','\n')
    netdiv.tropical<-traitDependentBAMM(BAMM.object, tropical, 1000, rate = "net diversification",return.full = FALSE, method = "mann-whitney", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    strapp.tropical.df<-as.data.frame(rbind(c(unlist(lambda.tropical$estimate),lambda.tropical$p.value,lambda.tropical$method,lambda.tropical$two.tailed,lambda.tropical$rate),c(unlist(mu.tropical$estimate),mu.tropical$p.value,mu.tropical$method,mu.tropical$two.tailed,mu.tropical$rate),c(unlist(netdiv.tropical$estimate),netdiv.tropical$p.value,netdiv.tropical$method,netdiv.tropical$two.tailed,netdiv.tropical$rate)))
    colnames(strapp.tropical.df)<-c(paste0('estimate.',names(unlist(lambda.tropical$estimate))),'p.value','method','two.tailed','rate')
    return(strapp.tropical.df)
  }else if(mode=='strict.tropical'){
    tropical<-trait.table$strict.tropical
    names(tropical)<-trait.table$binomial
    cat('running strapp lambda','\n')
    lambda.tropical<-traitDependentBAMM(BAMM.object, tropical, 1000, rate = "speciation",return.full = FALSE, method = "k", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp mu','\n')
    mu.tropical<-traitDependentBAMM(BAMM.object, tropical, 1000, rate = "extinction",return.full = FALSE, method = "k", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp netdiv','\n')
    netdiv.tropical<-traitDependentBAMM(BAMM.object, tropical, 1000, rate = "net diversification",return.full = FALSE, method = "k", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    strapp.tropical.df<-as.data.frame(rbind(c(unlist(lambda.tropical$estimate),lambda.tropical$p.value,lambda.tropical$method,lambda.tropical$two.tailed,lambda.tropical$rate),c(unlist(mu.tropical$estimate),mu.tropical$p.value,mu.tropical$method,mu.tropical$two.tailed,mu.tropical$rate),c(unlist(netdiv.tropical$estimate),netdiv.tropical$p.value,netdiv.tropical$method,netdiv.tropical$two.tailed,netdiv.tropical$rate)))
    colnames(strapp.tropical.df)<-c(paste0('estimate.',names(unlist(lambda.tropical$estimate))),'p.value','method','two.tailed','rate')
    return(strapp.tropical.df)
  }else if (mode=='continuous.medianlatitude'){
    median.lat<-trait.table$Median.Latitude
    names(median.lat)<-trait.table$binomial
    cat('running strapp lambda','\n')
    lambda.median.lat<-traitDependentBAMM(BAMM.object, median.lat, 1000, rate = "speciation",return.full = FALSE, method = "spearman", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp mu','\n')
    mu.median.lat<-traitDependentBAMM(BAMM.object, median.lat, 1000, rate = "extinction",return.full = FALSE, method = "spearman", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp netdiv','\n')
    netdiv.median.lat<-traitDependentBAMM(BAMM.object, median.lat, 1000, rate = "net diversification",return.full = FALSE, method = "spearman", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    strapp.median.lat.df<-as.data.frame(rbind(c(lambda.median.lat$estimate,lambda.median.lat$p.value,lambda.median.lat$method,lambda.median.lat$two.tailed,lambda.median.lat$rate),c(mu.median.lat$estimate,mu.median.lat$p.value,mu.median.lat$method,mu.median.lat$two.tailed,mu.median.lat$rate),c(netdiv.median.lat$estimate,netdiv.median.lat$p.value,netdiv.median.lat$method,netdiv.median.lat$two.tailed,netdiv.median.lat$rate)))
    colnames(strapp.median.lat.df)<-c('estimate','p.value','method','two.tailed','rate')
    return(strapp.median.lat.df)
  }else if (mode=='continuous.abs.medianlatitude'){
    median.lat<-abs(trait.table$Median.Latitude)
    names(median.lat)<-trait.table$binomial
    cat('running strapp lambda','\n')
    lambda.median.lat<-traitDependentBAMM(BAMM.object, median.lat, 1000, rate = "speciation",return.full = FALSE, method = "spearman", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp mu','\n')
    mu.median.lat<-traitDependentBAMM(BAMM.object, median.lat, 1000, rate = "extinction",return.full = FALSE, method = "spearman", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp netdiv','\n')
    netdiv.median.lat<-traitDependentBAMM(BAMM.object, median.lat, 1000, rate = "net diversification",return.full = FALSE, method = "spearman", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    strapp.median.lat.df<-as.data.frame(rbind(c(lambda.median.lat$estimate,lambda.median.lat$p.value,lambda.median.lat$method,lambda.median.lat$two.tailed,lambda.median.lat$rate),c(mu.median.lat$estimate,mu.median.lat$p.value,mu.median.lat$method,mu.median.lat$two.tailed,mu.median.lat$rate),c(netdiv.median.lat$estimate,netdiv.median.lat$p.value,netdiv.median.lat$method,netdiv.median.lat$two.tailed,netdiv.median.lat$rate)))
    colnames(strapp.median.lat.df)<-c('estimate','p.value','method','two.tailed','rate')
    return(strapp.median.lat.df)
  }else if (mode=='latitudinalband'){
  	 lat.band<-trait.table$latitudinal.band
    names(lat.band)<-trait.table$binomial
    cat('running strapp lambda','\n')
    lambda.lat.band<-traitDependentBAMM(BAMM.object, lat.band, 1000, rate = "speciation",return.full = FALSE, method = "kruskal", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp mu','\n')
    mu.lat.band<-traitDependentBAMM(BAMM.object, lat.band, 1000, rate = "extinction",return.full = FALSE, method = "kruskal", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp netdiv','\n')
    netdiv.lat.band<-traitDependentBAMM(BAMM.object, lat.band, 1000, rate = "net diversification",return.full = FALSE, method = "kruskal", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    strapp.latband.df<-as.data.frame(rbind(c(unlist(lambda.lat.band$estimate),lambda.lat.band$p.value,lambda.lat.band$method,lambda.lat.band$two.tailed,lambda.lat.band$rate),c(unlist(mu.lat.band$estimate),mu.lat.band$p.value,mu.lat.band$method,mu.lat.band$two.tailed,mu.lat.band$rate),c(unlist(netdiv.lat.band$estimate),netdiv.lat.band$p.value,netdiv.lat.band$method,netdiv.lat.band$two.tailed,netdiv.lat.band$rate)))
    strapp.latband.df<-strapp.latband.df[c(match(as.character(seq(from=min(lat.band),to=max(lat.band),by=10)),colnames(strapp.latband.df)),grep('V',colnames(strapp.latband.df)))]
    colnames(strapp.latband.df)[grep('V',colnames(strapp.latband.df))]<-c('p.value','method','two.tailed','rate')
    return(strapp.latband.df)
  }else if (mode=='abs.latitudinalband'){
    lat.band.abs<-abs(trait.table$latitudinal.band)
    names(lat.band.abs)<-trait.table$binomial
    cat('running strapp lambda','\n')
    lambda.lat.band.abs<-traitDependentBAMM(BAMM.object, lat.band.abs, 1000, rate = "speciation",return.full = FALSE, method = "kruskal", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp mu','\n')
    mu.lat.band.abs<-traitDependentBAMM(BAMM.object, lat.band.abs, 1000, rate = "extinction",return.full = FALSE, method = "kruskal", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    cat('running strapp netdiv','\n')
    netdiv.lat.band.abs<-traitDependentBAMM(BAMM.object, lat.band.abs, 1000, rate = "net diversification",return.full = FALSE, method = "kruskal", logrates = TRUE,two.tailed = TRUE, traitorder = NA, nthreads = 1)
    strapp.latband.abs.df<-as.data.frame(rbind(c(unlist(lambda.lat.band.abs$estimate),lambda.lat.band.abs$p.value,lambda.lat.band.abs$method,lambda.lat.band.abs$two.tailed,lambda.lat.band.abs$rate),c(unlist(mu.lat.band.abs$estimate),mu.lat.band.abs$p.value,mu.lat.band.abs$method,mu.lat.band.abs$two.tailed,mu.lat.band.abs$rate),c(unlist(netdiv.lat.band.abs$estimate),netdiv.lat.band.abs$p.value,netdiv.lat.band.abs$method,netdiv.lat.band.abs$two.tailed,netdiv.lat.band.abs$rate)))
    strapp.latband.abs.df<-strapp.latband.abs.df[c(match(as.character(seq(from=min(lat.band.abs),to=max(lat.band.abs),by=10)),colnames(strapp.latband.abs.df)),grep('V',colnames(strapp.latband.abs.df)))]
    colnames(strapp.latband.abs.df)[grep('V',colnames(strapp.latband.abs.df))]<-c('p.value','method','two.tailed','rate')
    return(strapp.latband.abs.df)
  }
  
}
args<-commandArgs(trailingOnly = TRUE)
print(args)
replicate<-as.numeric(args[1])
#this checks whether the GBIFdata.BAMM.subsampled table and the subsetted BAMM object exist and creates them if they don't
if((file.exists(paste('./results/subsampled_unbiased_GBOTB/input_tables/GBIFdata.BAMM.GBOTB.subsampled_',replicate,'_table.txt',sep=''))==F)||(file.exists(paste('./results/subsampled_unbiased_GBOTB/input_tables/all_clades_50million_1000samples_latitudedata_GBOTB_unbiased_',replicate,'.RDS',sep=''))==F)){
	cat('files not found, generating','\n')
  	GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
	#tropical defined as abs(median latitude <23.5)
	#temperate defined as abs(median latitude >23.5)
	GBIFdata$tropical<-0
	GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
	GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
	sampling.proportions<-table(GBIFdata[GBIFdata$strict.tropical!=2,]$strict.tropical)/nrow(GBIFdata[GBIFdata$strict.tropical!=2,])
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
	GBIFdata.BAMM.subsampled<-GBIFdata.BAMM[c(sample(which(GBIFdata.BAMM$tropical==0),size = round(sampling.proportions[1]*size.dataset)),sample(which(GBIFdata.BAMM$tropical==1),size = round(sampling.proportions[2]*size.dataset))),]
	write.table(GBIFdata.BAMM.subsampled,file=paste('./results/subsampled_unbiased/input_tables/GBIFdata.BAMM.subsampled_',replicate,'_table.txt',sep=''),sep='\t',quote=F,row.names=F)
	BAMM.object.subsampled<-subtreeBAMM(BAMM.object,tips=GBIFdata.BAMM.subsampled$binomial)
	saveRDS(BAMM.object.subsampled,file=paste('./results/subsampled_unbiased/all_clades_50million_1000samples_latitudedata_unbiased_',replicate,'.RDS',sep=''))

}else{
	BAMM.object.subsampled<-readRDS(file=paste('./results/subsampled_unbiased_GBOTB/input_tables/all_clades_50million_1000samples_latitudedata_GBOTB_unbiased_',replicate,'.RDS',sep=''))
	GBIFdata.BAMM.subsampled<-read.table(file=paste('./results/subsampled_unbiased_GBOTB/input_tables/GBIFdata.BAMM.GBOTB.subsampled_',replicate,'_table.txt',sep=''),sep='\t',header=T,stringsAsFactors=F)
	GBIFdata.BAMM.subsampled$latitudinal.band.new<-apply(GBIFdata.BAMM.subsampled,1,function(x) assign_equalwidth_band(number=as.numeric(x['Median.Latitude'])))
	#change the names of the bands to multiples of ten
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='-61.1','latitudinal.band.new']<-(-60)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='-51.7','latitudinal.band.new']<-(-50)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='-42.3','latitudinal.band.new']<-(-40)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='-32.9','latitudinal.band.new']<-(-30)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='-23.5','latitudinal.band.new']<-(-20)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='-14.1','latitudinal.band.new']<-(-10)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='-4.7','latitudinal.band.new']<-(0)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='4.7','latitudinal.band.new']<-(10)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='14.1','latitudinal.band.new']<-(20)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='23.5','latitudinal.band.new']<-(30)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='32.9','latitudinal.band.new']<-(40)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='42.3','latitudinal.band.new']<-(50)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='51.7','latitudinal.band.new']<-(60)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='61.1','latitudinal.band.new']<-(70)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='70.5','latitudinal.band.new']<-(80)
	GBIFdata.BAMM.subsampled[GBIFdata.BAMM.subsampled$latitudinal.band.new=='79.9','latitudinal.band.new']<-(90)
	#use 9.4 degree width bands as the default latitudinal bands
	GBIFdata.BAMM.subsampled$latitudinal.band<-GBIFdata.BAMM.subsampled$latitudinal.band.new
}
tropical.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.subsampled,trait.table =GBIFdata.BAMM.subsampled,mode = 'binary.tropical')
write.table(tropical.strapp,file=paste('./results/subsampled_unbiased_GBOTB/tropical.strapp.unbiased_',replicate,'.txt',sep=''),sep='\t',quote=F,row.names=F)

#medianlatitude.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.subsampled,trait.table =GBIFdata.BAMM.subsampled,mode = 'continuous.medianlatitude')
#write.table(medianlatitude.strapp,file=paste('./results/subsampled_unbiased/medianlatitude.strapp.unbiased_',replicate,'.txt',sep=''),sep='\t',quote=F,row.names=F)

abs.medianlatitude.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.subsampled,trait.table =GBIFdata.BAMM.subsampled,mode = 'continuous.abs.medianlatitude')
write.table(abs.medianlatitude.strapp,file=paste('./results/subsampled_unbiased_GBOTB/abs.medianlatitude.strapp.unbiased_',replicate,'.txt',sep=''),sep='\t',quote=F,row.names=F)

#latitudinalband.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.subsampled,trait.table =GBIFdata.BAMM.subsampled,mode = 'latitudinalband')
#write.table(latitudinalband.strapp,file=paste('./results/strapp_unbiased/latitudinalband.strapp.unbiased_',replicate,'.txt',sep=''),sep='\t',quote=F,row.names=F)

abs.latitudinalband.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.subsampled,trait.table =GBIFdata.BAMM.subsampled,mode = 'abs.latitudinalband')
write.table(abs.latitudinalband.strapp,file=paste('./results/subsampled_unbiased_GBOTB/abs.latitudinalband.strapp.unbiased_',replicate,'.txt',sep=''),sep='\t',quote=F,row.names=F)



