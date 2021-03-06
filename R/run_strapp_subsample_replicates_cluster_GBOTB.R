setwd('/home/ji247/ldg_plants/')

library(BAMMtools)
library(plyr)

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
    colnames(strapp.tropical.df)<-c('estimate.0','estimate.1','p.value','method','two.tailed','rate')
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
if(file.exists(paste('./results/subsampled_unbiased_GBOTB/input_tables/GBIFdata.BAMM.GBOTB.subsampled_',replicate,'_table.txt',sep=''))==F){
	return('file not found')
}else if((file.exists(paste0('./results/subsampled_unbiased_GBOTB/input_tables/GBIFdata.BAMM.GBOTB.subsampled_',replicate,'_table.txt'))==T)&&(file.exists(paste0('./results/subsampled_unbiased_GBOTB/input_tables/all_clades_50million_1000samples_latitudedata_GBOTB_unbiased_',replicate,'.RDS'))==F)){
	GBIFdata.BAMM.subsampled<-read.table(file=paste0('./results/subsampled_unbiased_GBOTB/input_tables/GBIFdata.BAMM.GBOTB.subsampled_',replicate,'_table.txt'),header=T,stringsAsFactors=F)
	BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata_GBOTB.RDS')
	cat('subsampling BAMM object','\n')
	BAMM.object.subsampled<-subtreeBAMM(BAMM.object,tips=GBIFdata.BAMM.subsampled$binomial)
	saveRDS(BAMM.object.subsampled,file=paste0('./results/subsampled_unbiased_GBOTB/input_tables/all_clades_50million_1000samples_latitudedata_GBOTB_unbiased_',replicate,'.RDS'))
}else{
	BAMM.object.subsampled<-readRDS(file=paste0('./results/subsampled_unbiased_GBOTB/input_tables/all_clades_50million_1000samples_latitudedata_GBOTB_unbiased_',replicate,'.RDS'))
	GBIFdata.BAMM.subsampled<-read.table(file=paste0('./results/subsampled_unbiased_GBOTB/input_tables/GBIFdata.BAMM.GBOTB.subsampled_',replicate,'_table.txt'),header=T,stringsAsFactors=F)
}
GBIFdata.BAMM.subsampled$latitudinal.band<-NA
GBIFdata.BAMM.subsampled$latitudinal.band<-apply(GBIFdata.BAMM.subsampled,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))

tropical.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.subsampled,trait.table =GBIFdata.BAMM.subsampled,mode = 'binary.tropical')
write.table(tropical.strapp,file=paste('./results/subsampled_unbiased_GBOTB/tropical.strapp.subsampled_',replicate,'.txt',sep=''),sep='\t',quote=F,row.names=F)

medianlatitude.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.subsampled,trait.table =GBIFdata.BAMM.subsampled,mode = 'continuous.medianlatitude')
write.table(medianlatitude.strapp,file=paste('./results/subsampled_unbiased_GBOTB/medianlatitude.strapp.subsampled_',replicate,'.txt',sep=''),sep='\t',quote=F,row.names=F)

abs.medianlatitude.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.subsampled,trait.table =GBIFdata.BAMM.subsampled,mode = 'continuous.abs.medianlatitude')
write.table(abs.medianlatitude.strapp,file=paste('./results/subsampled_unbiased_GBOTB/abs.medianlatitude.strapp.subsampled_',replicate,'.txt',sep=''),sep='\t',quote=F,row.names=F)

latitudinalband.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.subsampled,trait.table =GBIFdata.BAMM.subsampled,mode = 'latitudinalband')
write.table(latitudinalband.strapp,file=paste('./results/subsampled_unbiased_GBOTB/latitudinalband.strapp.subsampled_',replicate,'.txt',sep=''),sep='\t',quote=F,row.names=F)

abs.latitudinalband.strapp<-run_strapp_BAMM.object(BAMM.object = BAMM.object.subsampled,trait.table =GBIFdata.BAMM.subsampled,mode = 'abs.latitudinalband')
write.table(abs.latitudinalband.strapp,file=paste('./results/subsampled_unbiased_GBOTB/abs.latitudinalband.strapp.subsampled_',replicate,'.txt',sep=''),sep='\t',quote=F,row.names=F)

