library(plyr)
library(BAMMtools)
library(RColorBrewer)
library(nlme)
library(caper)
library(phylolm)

blue<-brewer.pal(11,'RdYlBu')[9]
red<-brewer.pal(11,'RdYlBu')[2]

data_frame_compare_rates_full_vs_subsampled<-function(folder.subsampled.event_data){
  #get all eventData RDS of subsampled replicates
  cat('getting all eventdatas in input folder','\n')
  eventdata.files<-list.files(folder.subsampled.event_data,pattern='.RDS')
  eventdata.list<-lapply(eventdata.files,function(x)readRDS(paste(folder.subsampled.event_data,x,sep='')))
  names(eventdata.list)<-gsub(gsub(eventdata.files,pattern='.*_subsampled',replacement='subsampled'),pattern='_eventData.RDS',replacement='')
  #get TipRates of all subsampled dataset
  cat('getting rates of input eventdatas','\n')
  ratesdf.subsampled.list<-lapply(eventdata.list,function(x)getTipRates(x))
  ratesdf.subsampled.list<-lapply(ratesdf.subsampled.list,function(x)as.data.frame(cbind(names(x$lambda.avg),x$lambda.avg),stringsAsFactors = F))
  ratesdf.subsampled.list<-lapply(ratesdf.subsampled.list,function(x){rownames(x)<-NULL;return(x)})
  names(ratesdf.subsampled.list)<-names(eventdata.list)
  ratesdf.subsampled.list<-lapply(c(1:length(ratesdf.subsampled.list)),function(x){colnames(ratesdf.subsampled.list[[x]])<-c('species',paste(names(ratesdf.subsampled.list)[x],'_lambda.avg',sep=''));return(ratesdf.subsampled.list[[x]])})
  #load the 28,057 BAMM object
  cat('getting rates of 28k tip dataset','\n')
  BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
  #and get TipRates
  ratesdf<-getTipRates(BAMM.object)
  rates.df<-as.data.frame(cbind(names(ratesdf$lambda.avg),ratesdf$lambda.avg),stringsAsFactors = F)
  row.names(rates.df)<-NULL
  colnames(rates.df)<-c('species','lambda.avg.28k')
  #combine all tip estimates across subsamples with the full dataset estimates
  cat('merging all rate estimates into a single dataframe','\n')
  rates.df.merge<-rates.df
  for (i in 1:length(ratesdf.subsampled.list)){
    rates.df.merge<-merge(rates.df.merge,ratesdf.subsampled.list[[i]],by.x='species',by.y='species',all.x=T)
  }
  #add tropical info + other GBIFstuff
  cat('adding GBIF info etc','\n')
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
  GBIFdata<-GBIFdata[-which(duplicated(GBIFdata$binomial)),]
  GBIFdata<-unique(GBIFdata)
  GBIFdata$latitudinal.band<-NA
  GBIFdata$latitudinal.band<-apply(GBIFdata,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
  #merge with rates.df.merge
  rates.df.merge<-merge(rates.df.merge,GBIFdata,by.x='species',by.y='binomial',all.x=T)
  return(rates.df.merge)
}

replicates_phylolm_rho_trop_vs_temp_extremebias_GBOTB<-function(rates.df.merge){
  #run linear models for trop and temp for each subsample
  lm.both<-list()
  tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
  cor.spearman.log<-list()
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_extremebias_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    cor.spearman.log[[i]]<-cor.test(log(as.numeric(rates.df.merge.subset$lambda.avg.28k)),log(as.numeric(rates.df.merge.subset$lambda.avg.10k)),method='s')
    
  }
  return(list(lm.both,cor.spearman.log))
}
replicates_phylolm_rho_trop_vs_temp_unbiased_GBOTB<-function(rates.df.merge){
  #run linear models for trop and temp for each subsample
  lm.both<-list()
  tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
  cor.spearman.log<-list()
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_unbiased_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    cor.spearman.log[[i]]<-cor.test(log(as.numeric(rates.df.merge.subset$lambda.avg.28k)),log(as.numeric(rates.df.merge.subset$lambda.avg.10k)),method='s')
    
  }
  return(list(lm.both,cor.spearman.log))
}
replicates_phylolm_rho_trop_vs_temp_small_GBOTB<-function(rates.df.merge){
  #run linear models for trop and temp for each subsample
  lm.both<-list()
  tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
  cor.spearman.log<-list()
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_small_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    cor.spearman.log[[i]]<-cor.test(log(as.numeric(rates.df.merge.subset$lambda.avg.28k)),log(as.numeric(rates.df.merge.subset$lambda.avg.10k)),method='s')
    
  }
  return(list(lm.both,cor.spearman.log))
}
replicates_phylolm_rho_trop_vs_temp_extremebias_GBOTB<-function(rates.df.merge){
  #run linear models for trop and temp for each subsample
  lm.both<-list()
  tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
  cor.spearman.log<-list()
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_extremebias_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    cor.spearman.log[[i]]<-cor.test(log(as.numeric(rates.df.merge.subset$lambda.avg.28k)),log(as.numeric(rates.df.merge.subset$lambda.avg.10k)),method='s')
    
  }
  return(list(lm.both,cor.spearman.log))
}
replicates_phylolm_rho_trop_vs_temp_small_GBOTB<-function(rates.df.merge){
  #run linear models for trop and temp for each subsample
  lm.both<-list()
  tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
  cor.spearman.log<-list()
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_small_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    cor.spearman.log[[i]]<-cor.test(log(as.numeric(rates.df.merge.subset$lambda.avg.28k)),log(as.numeric(rates.df.merge.subset$lambda.avg.10k)),method='s')
    
  }
  return(list(lm.both,cor.spearman.log))
}
data_frame_compare_rates_full_vs_subsampled_GBOTB<-function(folder.subsampled.event_data){
  #get all eventData RDS of subsampled replicates
  cat('getting all eventdatas in input folder','\n')
  eventdata.files<-list.files(folder.subsampled.event_data,pattern='.RDS')
  eventdata.list<-lapply(eventdata.files,function(x)readRDS(paste(folder.subsampled.event_data,x,sep='')))
  names(eventdata.list)<-gsub(gsub(eventdata.files,pattern='.*_subsampled',replacement='subsampled'),pattern='_eventData.RDS',replacement='')
  #get TipRates of all subsampled dataset
  cat('getting rates of input eventdatas','\n')
  ratesdf.subsampled.list<-lapply(eventdata.list,function(x)getTipRates(x))
  ratesdf.subsampled.list<-lapply(ratesdf.subsampled.list,function(x)as.data.frame(cbind(names(x$lambda.avg),x$lambda.avg),stringsAsFactors = F))
  ratesdf.subsampled.list<-lapply(ratesdf.subsampled.list,function(x){rownames(x)<-NULL;return(x)})
  names(ratesdf.subsampled.list)<-names(eventdata.list)
  ratesdf.subsampled.list<-lapply(c(1:length(ratesdf.subsampled.list)),function(x){colnames(ratesdf.subsampled.list[[x]])<-c('species',paste(names(ratesdf.subsampled.list)[x],'_lambda.avg',sep=''));return(ratesdf.subsampled.list[[x]])})
  #load the 28,057 BAMM object
  cat('getting rates of 60k tip dataset','\n')
  BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata_GBOTB.RDS')
  #and get TipRates
  ratesdf<-getTipRates(BAMM.object)
  rates.df<-as.data.frame(cbind(names(ratesdf$lambda.avg),ratesdf$lambda.avg),stringsAsFactors = F)
  row.names(rates.df)<-NULL
  colnames(rates.df)<-c('species','lambda.avg.28k')
  #combine all tip estimates across subsamples with the full dataset estimates
  cat('merging all rate estimates into a single dataframe','\n')
  rates.df.merge<-rates.df
  for (i in 1:length(ratesdf.subsampled.list)){
    rates.df.merge<-merge(rates.df.merge,ratesdf.subsampled.list[[i]],by.x='species',by.y='species',all.x=T)
  }
  #add tropical info + other GBIFstuff
  cat('adding GBIF info etc','\n')
  GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
  GBIFdata$tropical<-0
  GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
  GBIFdata$strict.tropical<-NA
  #strict tropical 1 = abs(max latitude)<23.5 & abs(min latitude)<23.5  & median.latitude <23.5
  GBIFdata[abs(GBIFdata$Max.Latitude)<=23.5&abs(GBIFdata$Min.Latitude)<=23.5&abs(GBIFdata$Median.Latitude)<=23.5,'strict.tropical']<-1
  #strict tropical 0 (strict temperate)
  GBIFdata[abs(GBIFdata$Max.Latitude)>23.5&abs(GBIFdata$Min.Latitude)>23.5&abs(GBIFdata$Median.Latitude)>23.5,'strict.tropical']<-0
  #2: not strict tropical or strict temperate species
  GBIFdata$strict.tropical[is.na(GBIFdata$strict.tropical)]<-2
  GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
  GBIFdata<-GBIFdata[-which(duplicated(GBIFdata$binomial)),]
  GBIFdata<-unique(GBIFdata)
  GBIFdata$latitudinal.band<-NA
  GBIFdata$latitudinal.band<-apply(GBIFdata,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
  #merge with rates.df.merge
  rates.df.merge<-merge(rates.df.merge,GBIFdata,by.x='species',by.y='binomial',all.x=T)
  return(rates.df.merge)
}

plot_replicates_lm_trop_vs_temp<-function(rates.df.merge){
  rates.df.merge$colour<-NA
  rates.df.merge[rates.df.merge$tropical==0,'colour']<-blue
  rates.df.merge[rates.df.merge$tropical==1,'colour']<-red
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(pty="s") 
  lm.trop<-list()
  lm.temp<-list()
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  for (i in c(1:10)){
    lm.trop[[i]]<-lm(as.numeric(rates.df.merge[rates.df.merge$tropical==1,'lambda.avg.28k'])~as.numeric(rates.df.merge[rates.df.merge$tropical==1,i+2]))
    lm.temp[[i]]<-lm(as.numeric(rates.df.merge[rates.df.merge$tropical==0,'lambda.avg.28k'])~as.numeric(rates.df.merge[rates.df.merge$tropical==0,i+2]))
    
    plot(as.numeric(rates.df.merge[,2])~as.numeric(rates.df.merge[,i+2]),pch=16,xlab=paste('rates.',gsub(colnames(rates.df.merge)[i+2],pattern='_lambda.avg',replacement=''),sep=''),ylab='rates.28k',col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1),axes=F)
    #points(as.numeric(rates.df.merge[rates.df.merge$tropical==1,2])~as.numeric(rates.df.merge[rates.df.merge$tropical==1,i+2]),pch=16,col=adjustcolor('red',alpha.f = 0.1))
    abline(lm.trop[[i]],col=red)
    abline(lm.temp[[i]],col=blue)
    
    mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
    if (i %in% c(6,7,8,9,10)){
      axis(1, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
      
    }else if (i %in% c(1,6)){
      axis(2, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
    }
    box(col = "grey60")
    
  }
  mtext("rates.10k.subsampled.tree", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
  mtext("rates.28k.tree", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
}

########THIS PLOTS FIG4
plot_replicates_lm_trop_vs_temp_mean_sd<-function(rates.df.merge){
  #run linear models for trop and temp for each subsample
  lm.trop<-list()
  lm.temp<-list()
  for (i in c(1:10)){
    lm.trop[[i]]<-lm(as.numeric(rates.df.merge[rates.df.merge$tropical==1,'lambda.avg.28k'])~as.numeric(rates.df.merge[rates.df.merge$tropical==1,i+2]))
    lm.temp[[i]]<-lm(as.numeric(rates.df.merge[rates.df.merge$tropical==0,'lambda.avg.28k'])~as.numeric(rates.df.merge[rates.df.merge$tropical==0,i+2]))
  }
  #this takes the dimensions of the plot
  u<-par('usr')
  plot(c(1,1),xlim=c(0,15),ylim=c(0,15),type='n',xlab='rates.subsampled.10k.trees',ylab='rates.28k.tree',main='')
  #plot tropical mean slope + intercept
  abline(a=mean(unlist(lapply(lm.trop,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.trop,function(x)x$coefficients[2]))),col=red,lwd=3)
  #abline(a=mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.trop,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.trop,function(x)x$coefficients[2]))),col='red',lty=2,add/)
  #abline(a=mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.trop,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.trop,function(x)x$coefficients[2]))),col='red',lty=2)
  
  #get +1sd -1sd for slope and intercept, build a curve to get points to fill a polygon
  lower.sd.trop <- curve(mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))) + x*(mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.trop,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha.f = 0.5))
  upper.sd.trop <- curve(mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[1])))  + x*(mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.trop,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha.f = 0.5))
  polygon(c(lower.sd.trop$x,rev(upper.sd.trop$x)), c(lower.sd.trop$y, rev(upper.sd.trop$y)),col=adjustcolor(red,alpha.f = 0.5), border=NA)
  
  #plot tropical mean slope + intercept
  abline(a=mean(unlist(lapply(lm.temp,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.temp,function(x)x$coefficients[2]))),col=blue,lwd=3)
  #abline(a=mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[2]))),col='blue',lty=2)
  #abline(a=mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[2]))),col='blue',lty=2)
  #get +1sd -1sd for slope and intercept, build a curve to get points to fill a polygon
  lower.sd.temp <- curve(mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))) + x*(mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha.f = 0.5))
  upper.sd.temp <- curve(mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[1])))  + x*(mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha.f = 0.5))
  polygon(c(lower.sd.temp$x,rev(upper.sd.temp$x)), c(lower.sd.temp$y, rev(upper.sd.temp$y)),col=adjustcolor(blue,alpha.f = 0.5), border=NA)
  legend('topleft',legend=c('tropical','temperate'),col=c(red,blue),lty=1,cex=.7,box.lty=0)
  
}


#this function takes a rates.df.merge dataframe ('species',') and returns a set of phylolms
#as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical
replicates_phylolm_rho_trop_vs_temp<-function(rates.df.merge){
  #run linear models for trop and temp for each subsample
  lm.both<-list()
  tree<-read.tree('./raw_data/BAMM_Species_tree_noseed.tree')
  cor.spearman.log<-list()
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    cor.spearman.log[[i]]<-cor.test(log(as.numeric(rates.df.merge.subset$lambda.avg.28k)),log(as.numeric(rates.df.merge.subset$lambda.avg.10k)),method='s')
    
  }
  return(list(lm.both,cor.spearman.log))
}

replicates_phylolm_rho_trop_vs_temp_extremebias<-function(rates.df.merge){
  #run linear models for trop and temp for each subsample
  lm.both<-list()
  tree<-read.tree('./raw_data/BAMM_Species_tree_noseed.tree')
  cor.spearman.log<-list()
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_extremebias_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    cor.spearman.log[[i]]<-cor.test(log(as.numeric(rates.df.merge.subset$lambda.avg.28k)),log(as.numeric(rates.df.merge.subset$lambda.avg.10k)),method='s')
    
  }
  return(list(lm.both,cor.spearman.log))
}

plot_replicates_phylolm_trop_vs_temp_extremebias_GBOTB<-function(rates.df.merge){
  par(mfrow=c(2,5))
  lm.both<-list()
  tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
  rates.df.merge$colour<-NA
  rates.df.merge[rates.df.merge$tropical==0,'colour']<-blue
  rates.df.merge[rates.df.merge$tropical==1,'colour']<-red
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(pty="s") 
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  par(mfrow=c(2,5))
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_extremebias_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    plot(as.numeric(rates.df.merge[,grep('lambda.avg.28k',colnames(rates.df.merge))])~as.numeric(rates.df.merge[,grep(paste('subsampled_extremebias_',i,'_lambda.avg',sep=''),colnames(rates.df.merge))]),pch=16,xlab=paste('rates.',gsub(colnames(rates.df.merge)[i+2],pattern='_lambda.avg',replacement=''),sep=''),ylab='rates.28k',col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1),axes=F,main=i)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2],col=blue)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2]+lm.both[[1]]$coefficients[3],col=red)
    u<-par('usr')
    lower.confint.lm.temp<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    upper.confint.lm.temp<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    lower.confint.lm.trop<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1]+x*confint(lm.both[[i]])[3,1]+x*confint(lm.both[[i]])[4,1], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    upper.confint.lm.trop<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2]+x*confint(lm.both[[i]])[3,2]+x*confint(lm.both[[i]])[4,2], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    polygon(c(lower.confint.lm.trop$x,rev(upper.confint.lm.trop$x)), c(lower.confint.lm.trop$y, rev(upper.confint.lm.trop$y)),col=adjustcolor(red,alpha=0.3),border=adjustcolor(red,alpha=0.3))
    polygon(c(lower.confint.lm.temp$x,rev(upper.confint.lm.temp$x)), c(lower.confint.lm.temp$y, rev(upper.confint.lm.temp$y)),col=adjustcolor(blue,alpha=0.3),border=adjustcolor(blue,alpha=0.3))
    mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
    if (i %in% c(6,7,8,9,10)){
      axis(1, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
      
    }else if (i %in% c(1,6)){
      axis(2, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
    }
    box(col = "grey60")
    mtext("λ.extremebiased(lineages/myr)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
    mtext("λ.full (lineages/myr)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
  }
}

plot_replicates_phylolm_trop_vs_temp_unbiased_GBOTB<-function(rates.df.merge){
  par(mfrow=c(2,5))
  lm.both<-list()
  tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
  rates.df.merge$colour<-NA
  rates.df.merge[rates.df.merge$tropical==0,'colour']<-blue
  rates.df.merge[rates.df.merge$tropical==1,'colour']<-red
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(pty="s") 
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  par(mfrow=c(2,5))
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_unbiased_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    plot(as.numeric(rates.df.merge[,grep('lambda.avg.28k',colnames(rates.df.merge))])~as.numeric(rates.df.merge[,grep(paste('subsampled_unbiased_',i,'_lambda.avg',sep=''),colnames(rates.df.merge))]),pch=16,xlab=paste('rates.',gsub(colnames(rates.df.merge)[i+2],pattern='_lambda.avg',replacement=''),sep=''),ylab='rates.28k',col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1),axes=F,main=i)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2],col=blue)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2]+lm.both[[1]]$coefficients[3],col=red)
    u<-par('usr')
    lower.confint.lm.temp<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    upper.confint.lm.temp<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    lower.confint.lm.trop<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1]+x*confint(lm.both[[i]])[3,1]+x*confint(lm.both[[i]])[4,1], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    upper.confint.lm.trop<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2]+x*confint(lm.both[[i]])[3,2]+x*confint(lm.both[[i]])[4,2], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    polygon(c(lower.confint.lm.trop$x,rev(upper.confint.lm.trop$x)), c(lower.confint.lm.trop$y, rev(upper.confint.lm.trop$y)),col=adjustcolor(red,alpha=0.3),border=adjustcolor(red,alpha=0.3))
    polygon(c(lower.confint.lm.temp$x,rev(upper.confint.lm.temp$x)), c(lower.confint.lm.temp$y, rev(upper.confint.lm.temp$y)),col=adjustcolor(blue,alpha=0.3),border=adjustcolor(blue,alpha=0.3))
    mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
    if (i %in% c(6,7,8,9,10)){
      axis(1, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
      
    }else if (i %in% c(1,6)){
      axis(2, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
    }
    box(col = "grey60")
    mtext("λ.unbiased(lineages/myr)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
    mtext("λ.full (lineages/myr)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
  }
}

plot_replicates_phylolm_trop_vs_temp_small_GBOTB<-function(rates.df.merge){
  par(mfrow=c(2,5))
  lm.both<-list()
  tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
  rates.df.merge$colour<-NA
  rates.df.merge[rates.df.merge$tropical==0,'colour']<-blue
  rates.df.merge[rates.df.merge$tropical==1,'colour']<-red
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(pty="s") 
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  par(mfrow=c(2,5))
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_small_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    plot(as.numeric(rates.df.merge[,grep('lambda.avg.28k',colnames(rates.df.merge))])~as.numeric(rates.df.merge[,grep(paste('subsampled_small_',i,'_lambda.avg',sep=''),colnames(rates.df.merge))]),pch=16,xlab=paste('rates.',gsub(colnames(rates.df.merge)[i+2],pattern='_lambda.avg',replacement=''),sep=''),ylab='rates.28k',col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1),axes=F,main=i)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2],col=blue)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2]+lm.both[[1]]$coefficients[3],col=red)
    u<-par('usr')
    lower.confint.lm.temp<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    upper.confint.lm.temp<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    lower.confint.lm.trop<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1]+x*confint(lm.both[[i]])[3,1]+x*confint(lm.both[[i]])[4,1], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    upper.confint.lm.trop<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2]+x*confint(lm.both[[i]])[3,2]+x*confint(lm.both[[i]])[4,2], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    polygon(c(lower.confint.lm.trop$x,rev(upper.confint.lm.trop$x)), c(lower.confint.lm.trop$y, rev(upper.confint.lm.trop$y)),col=adjustcolor(red,alpha=0.3),border=adjustcolor(red,alpha=0.3))
    polygon(c(lower.confint.lm.temp$x,rev(upper.confint.lm.temp$x)), c(lower.confint.lm.temp$y, rev(upper.confint.lm.temp$y)),col=adjustcolor(blue,alpha=0.3),border=adjustcolor(blue,alpha=0.3))
    mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
    if (i %in% c(6,7,8,9,10)){
      axis(1, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
      
    }else if (i %in% c(1,6)){
      axis(2, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
    }
    box(col = "grey60")
    mtext("λ.small(lineages/myr)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
    mtext("λ.full (lineages/myr)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
  }
}

plot_replicates_phylolm_trop_vs_temp_extremebias_log_GBOTB<-function(rates.df.merge){
  par(mfrow=c(2,5))
  lm.both<-list()
  tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
  rates.df.merge[,grep('lambda',colnames(rates.df.merge))]<-log(rates.df.merge[,grep('lambda',colnames(rates.df.merge))])
  rates.df.merge$colour<-NA
  rates.df.merge[rates.df.merge$tropical==0,'colour']<-blue
  rates.df.merge[rates.df.merge$tropical==1,'colour']<-red
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(pty="s") 
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  par(mfrow=c(2,5))
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_extremebias_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    plot(as.numeric(rates.df.merge[,grep('lambda.avg.28k',colnames(rates.df.merge))])~as.numeric(rates.df.merge[,grep(paste('subsampled_extremebias_',i,'_lambda.avg',sep=''),colnames(rates.df.merge))]),pch=16,xlab=paste('rates.',gsub(colnames(rates.df.merge)[i+2],pattern='_lambda.avg',replacement=''),sep=''),ylab='rates.28k',col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1),axes=F,main=i)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2],col=blue)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2]+lm.both[[1]]$coefficients[3],col=red)
    u<-par('usr')
    lower.confint.lm.temp<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    upper.confint.lm.temp<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    lower.confint.lm.trop<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1]+x*confint(lm.both[[i]])[3,1]+x*confint(lm.both[[i]])[4,1], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    upper.confint.lm.trop<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2]+x*confint(lm.both[[i]])[3,2]+x*confint(lm.both[[i]])[4,2], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    polygon(c(lower.confint.lm.trop$x,rev(upper.confint.lm.trop$x)), c(lower.confint.lm.trop$y, rev(upper.confint.lm.trop$y)),col=adjustcolor(red,alpha=0.3),border=adjustcolor(red,alpha=0.3))
    polygon(c(lower.confint.lm.temp$x,rev(upper.confint.lm.temp$x)), c(lower.confint.lm.temp$y, rev(upper.confint.lm.temp$y)),col=adjustcolor(blue,alpha=0.3),border=adjustcolor(blue,alpha=0.3))
    mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
    if (i %in% c(6,7,8,9,10)){
      axis(1, col = "grey40", col.axis = "grey20", at = log(c(0.05,1,5)),labels=c(0.05,1,5))
      
    }else if (i %in% c(1,6)){
      axis(2, col = "grey40", col.axis = "grey20", at = log(c(0.05,1,5)),labels=c(0.05,1,5))
    }
    box(col = "grey60")
    mtext("λ.extremebiased(lineages/myr)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
    mtext("λ.full (lineages/myr)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
  }
}
plot_replicates_phylolm_trop_vs_temp_unbiased_log_GBOTB<-function(rates.df.merge){
  par(mfrow=c(2,5))
  lm.both<-list()
  tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
  rates.df.merge[,grep('lambda',colnames(rates.df.merge))]<-log(rates.df.merge[,grep('lambda',colnames(rates.df.merge))])
  rates.df.merge$colour<-NA
  rates.df.merge[rates.df.merge$tropical==0,'colour']<-blue
  rates.df.merge[rates.df.merge$tropical==1,'colour']<-red
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(pty="s") 
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  par(mfrow=c(2,5))
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_unbiased_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    plot(as.numeric(rates.df.merge[,grep('lambda.avg.28k',colnames(rates.df.merge))])~as.numeric(rates.df.merge[,grep(paste('subsampled_unbiased_',i,'_lambda.avg',sep=''),colnames(rates.df.merge))]),pch=16,xlab=paste('rates.',gsub(colnames(rates.df.merge)[i+2],pattern='_lambda.avg',replacement=''),sep=''),ylab='rates.28k',col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1),axes=F,main=i)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2],col=blue)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2]+lm.both[[1]]$coefficients[3],col=red)
    u<-par('usr')
    lower.confint.lm.temp<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    upper.confint.lm.temp<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    lower.confint.lm.trop<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1]+x*confint(lm.both[[i]])[3,1]+x*confint(lm.both[[i]])[4,1], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    upper.confint.lm.trop<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2]+x*confint(lm.both[[i]])[3,2]+x*confint(lm.both[[i]])[4,2], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    polygon(c(lower.confint.lm.trop$x,rev(upper.confint.lm.trop$x)), c(lower.confint.lm.trop$y, rev(upper.confint.lm.trop$y)),col=adjustcolor(red,alpha=0.3),border=adjustcolor(red,alpha=0.3))
    polygon(c(lower.confint.lm.temp$x,rev(upper.confint.lm.temp$x)), c(lower.confint.lm.temp$y, rev(upper.confint.lm.temp$y)),col=adjustcolor(blue,alpha=0.3),border=adjustcolor(blue,alpha=0.3))
    mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
    if (i %in% c(6,7,8,9,10)){
      axis(1, col = "grey40", col.axis = "grey20", at = log(c(0.05,1,5)),labels=c(0.05,1,5))
      
    }else if (i %in% c(1,6)){
      axis(2, col = "grey40", col.axis = "grey20", at = log(c(0.05,1,5)),labels=c(0.05,1,5))
    }
    box(col = "grey60")
    mtext("λ.unbiased(lineages/myr)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
    mtext("λ.full (lineages/myr)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
  }
}
plot_replicates_phylolm_trop_vs_temp_small_log_GBOTB<-function(rates.df.merge){
  par(mfrow=c(2,5))
  lm.both<-list()
  tree<-read.tree('./raw_data/GBOTB_extended_bif_angios.tree')
  rates.df.merge[,grep('lambda',colnames(rates.df.merge))]<-log(rates.df.merge[,grep('lambda',colnames(rates.df.merge))])
  rates.df.merge$colour<-NA
  rates.df.merge[rates.df.merge$tropical==0,'colour']<-blue
  rates.df.merge[rates.df.merge$tropical==1,'colour']<-red
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(pty="s") 
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  par(mfrow=c(2,5))
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_small_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    plot(as.numeric(rates.df.merge[,grep('lambda.avg.28k',colnames(rates.df.merge))])~as.numeric(rates.df.merge[,grep(paste('subsampled_small_',i,'_lambda.avg',sep=''),colnames(rates.df.merge))]),pch=16,xlab=paste('rates.',gsub(colnames(rates.df.merge)[i+2],pattern='_lambda.avg',replacement=''),sep=''),ylab='rates.28k',col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1),axes=F,main=i)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2],col=blue)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2]+lm.both[[i]]$coefficients[3],col=red)
    u<-par('usr')
    lower.confint.lm.temp<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    upper.confint.lm.temp<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    lower.confint.lm.trop<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1]+x*confint(lm.both[[i]])[3,1]+x*confint(lm.both[[i]])[4,1], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    upper.confint.lm.trop<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2]+x*confint(lm.both[[i]])[3,2]+x*confint(lm.both[[i]])[4,2], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    polygon(c(lower.confint.lm.trop$x,rev(upper.confint.lm.trop$x)), c(lower.confint.lm.trop$y, rev(upper.confint.lm.trop$y)),col=adjustcolor(red,alpha=0.3),border=adjustcolor(red,alpha=0.3))
    polygon(c(lower.confint.lm.temp$x,rev(upper.confint.lm.temp$x)), c(lower.confint.lm.temp$y, rev(upper.confint.lm.temp$y)),col=adjustcolor(blue,alpha=0.3),border=adjustcolor(blue,alpha=0.3))
    mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
    if (i %in% c(6,7,8,9,10)){
      axis(1, col = "grey40", col.axis = "grey20", at = log(c(0.05,1,5)),labels=c(0.05,1,5))
      
    }else if (i %in% c(1,6)){
      axis(2, col = "grey40", col.axis = "grey20", at = log(c(0.05,1,5)),labels=c(0.05,1,5))
    }
    box(col = "grey60")
    mtext("λ.small(lineages/myr)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
    mtext("λ.full (lineages/myr)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
  }
}

replicates_phylolm_rho_trop_vs_temp_small<-function(rates.df.merge){
  #run linear models for trop and temp for each subsample
  lm.both<-list()
  tree<-read.tree('./raw_data/BAMM_Species_tree_noseed.tree')
  cor.spearman.log<-list()
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_small_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    cor.spearman.log[[i]]<-cor.test(log(as.numeric(rates.df.merge.subset$lambda.avg.28k)),log(as.numeric(rates.df.merge.subset$lambda.avg.10k)),method='s')
    
  }
  return(list(lm.both,cor.spearman.log))
}
plot_replicates_phylolm_trop_vs_temp_extremebias<-function(rates.df.merge){
  par(mfrow=c(2,5))
  lm.both<-list()
  tree<-read.tree('./raw_data/BAMM_Species_tree_noseed.tree')
  rates.df.merge$colour<-NA
  rates.df.merge[rates.df.merge$tropical==0,'colour']<-blue
  rates.df.merge[rates.df.merge$tropical==1,'colour']<-red
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(pty="s") 
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  par(mfrow=c(2,5))
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_extremebias_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    plot(as.numeric(rates.df.merge[,grep('lambda.avg.28k',colnames(rates.df.merge))])~as.numeric(rates.df.merge[,grep(paste('subsampled_extremebias_',i,'_lambda.avg',sep=''),colnames(rates.df.merge))]),pch=16,xlab=paste('rates.',gsub(colnames(rates.df.merge)[i+2],pattern='_lambda.avg',replacement=''),sep=''),ylab='rates.28k',col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1),axes=F,main=i)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2],col=blue)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2]+lm.both[[1]]$coefficients[3],col=red)
    u<-par('usr')
    lower.confint.lm.temp<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    upper.confint.lm.temp<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    lower.confint.lm.trop<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1]+x*confint(lm.both[[i]])[3,1]+x*confint(lm.both[[i]])[4,1], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    upper.confint.lm.trop<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2]+x*confint(lm.both[[i]])[3,2]+x*confint(lm.both[[i]])[4,2], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    polygon(c(lower.confint.lm.trop$x,rev(upper.confint.lm.trop$x)), c(lower.confint.lm.trop$y, rev(upper.confint.lm.trop$y)),col=adjustcolor(red,alpha=0.3),border=adjustcolor(red,alpha=0.3))
    polygon(c(lower.confint.lm.temp$x,rev(upper.confint.lm.temp$x)), c(lower.confint.lm.temp$y, rev(upper.confint.lm.temp$y)),col=adjustcolor(blue,alpha=0.3),border=adjustcolor(blue,alpha=0.3))
    mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
    if (i %in% c(6,7,8,9,10)){
      axis(1, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
      
    }else if (i %in% c(1,6)){
      axis(2, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
    }
    box(col = "grey60")
    mtext("λ.extremebiased(lineages/myr)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
    mtext("λ.full (lineages/myr)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
  }
}
plot_replicates_phylolm_trop_vs_temp_small<-function(rates.df.merge){
  par(mfrow=c(2,5))
  lm.both<-list()
  tree<-read.tree('./raw_data/BAMM_Species_tree_noseed.tree')
  rates.df.merge$colour<-NA
  rates.df.merge[rates.df.merge$tropical==0,'colour']<-blue
  rates.df.merge[rates.df.merge$tropical==1,'colour']<-red
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(pty="s") 
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  par(mfrow=c(2,5))
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_small_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    plot(as.numeric(rates.df.merge[,grep('lambda.avg.28k',colnames(rates.df.merge))])~as.numeric(rates.df.merge[,grep(paste('subsampled_small_',i,'_lambda.avg',sep=''),colnames(rates.df.merge))]),pch=16,xlab=paste('rates.',gsub(colnames(rates.df.merge)[i+2],pattern='_lambda.avg',replacement=''),sep=''),ylab='rates.28k',col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1),axes=F,main=i)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2],col=blue)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2]+lm.both[[1]]$coefficients[3],col=red)
    u<-par('usr')
    lower.confint.lm.temp<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    upper.confint.lm.temp<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    lower.confint.lm.trop<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1]+x*confint(lm.both[[i]])[3,1]+x*confint(lm.both[[i]])[4,1], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    upper.confint.lm.trop<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2]+x*confint(lm.both[[i]])[3,2]+x*confint(lm.both[[i]])[4,2], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    polygon(c(lower.confint.lm.trop$x,rev(upper.confint.lm.trop$x)), c(lower.confint.lm.trop$y, rev(upper.confint.lm.trop$y)),col=adjustcolor(red,alpha=0.3),border=adjustcolor(red,alpha=0.3))
    polygon(c(lower.confint.lm.temp$x,rev(upper.confint.lm.temp$x)), c(lower.confint.lm.temp$y, rev(upper.confint.lm.temp$y)),col=adjustcolor(blue,alpha=0.3),border=adjustcolor(blue,alpha=0.3))
    mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
    if (i %in% c(6,7,8,9,10)){
      axis(1, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
      
    }else if (i %in% c(1,6)){
      axis(2, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
    }
    box(col = "grey60")
    mtext("λ.small(lineages/myr)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
    mtext("λ.full (lineages/myr)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
  }
}

plot_replicates_phylolm_trop_vs_temp<-function(rates.df.merge){
  par(mfrow=c(2,5))
  lm.both<-list()
  tree<-read.tree('./raw_data/BAMM_Species_tree_noseed.tree')
  rates.df.merge$colour<-NA
  rates.df.merge[rates.df.merge$tropical==0,'colour']<-blue
  rates.df.merge[rates.df.merge$tropical==1,'colour']<-red
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
  par(pty="s") 
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  par(mfrow=c(2,5))
  #this is 10 because there are 10 replicates
  for (i in c(1:10)){
    selected.columns<-c('species','lambda.avg.28k',paste('subsampled_',i,'_lambda.avg',sep=''),'tropical')
    cat(i,'\n')
    #select columns with 28k rates and subsampled rate of replicate i
    rates.df.merge.subset<-rates.df.merge[,selected.columns]
    rates.df.merge.subset<-rates.df.merge.subset[complete.cases(rates.df.merge.subset),]
    colnames(rates.df.merge.subset)[3]<-'lambda.avg.10k'
    rownames(rates.df.merge.subset)<-rates.df.merge.subset$species
    #create a tree with species in subsampled dataset only
    replicate.tree<-drop.tip(tree,setdiff(tree$tip.label,rates.df.merge.subset$species))
    #phylolm of full rates ~ subsampled rates * tropical
    lm.both[[i]]<-phylolm(as.numeric(lambda.avg.28k)~as.numeric(lambda.avg.10k)*tropical,data=rates.df.merge.subset,phy = replicate.tree,model = "lambda")
    plot(as.numeric(rates.df.merge[,grep('lambda.avg.28k',colnames(rates.df.merge))])~as.numeric(rates.df.merge[,grep(paste('subsampled_',i,'_lambda.avg',sep=''),colnames(rates.df.merge))]),pch=16,xlab=paste('rates.',gsub(colnames(rates.df.merge)[i+2],pattern='_lambda.avg',replacement=''),sep=''),ylab='rates.28k',col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1),axes=F,main=i)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2],col=blue)
    abline(lm.both[[i]]$coefficients[1],lm.both[[i]]$coefficients[2]+lm.both[[1]]$coefficients[3],col=red)
    u<-par('usr')
    lower.confint.lm.temp<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    upper.confint.lm.temp<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2], from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha=0.3))
    lower.confint.lm.trop<-curve(confint(lm.both[[i]])[1,1] + x*confint(lm.both[[i]])[2,1]+x*confint(lm.both[[i]])[3,1]+x*confint(lm.both[[i]])[4,1], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    upper.confint.lm.trop<-curve(confint(lm.both[[i]])[1,2] + x*confint(lm.both[[i]])[2,2]+x*confint(lm.both[[i]])[3,2]+x*confint(lm.both[[i]])[4,2], from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha=0.3))
    polygon(c(lower.confint.lm.trop$x,rev(upper.confint.lm.trop$x)), c(lower.confint.lm.trop$y, rev(upper.confint.lm.trop$y)),col=adjustcolor(red,alpha=0.3),border=adjustcolor(red,alpha=0.3))
    polygon(c(lower.confint.lm.temp$x,rev(upper.confint.lm.temp$x)), c(lower.confint.lm.temp$y, rev(upper.confint.lm.temp$y)),col=adjustcolor(blue,alpha=0.3),border=adjustcolor(blue,alpha=0.3))
    mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
    if (i %in% c(6,7,8,9,10)){
      axis(1, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
      
    }else if (i %in% c(1,6)){
      axis(2, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
    }
    box(col = "grey60")
    mtext("λ.unbiased (lineages/myr)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
    mtext("λ.full (lineages/myr)", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
  }
}

###check proportional differences in rates 30k vs rates 10k (1- rate sub/rate full)
##mean.prop.diff.full.10k.temp<-list()
##mean.prop.diff.full.10k.trop<-list()
##library(DescTools)
##cor.test.full.10k.temp<-list()
##cor.test.full.10k.trop<-list()
##for (i in c(3:12)){
##  mean.prop.diff.full.10k.temp[[i-2]]<-abs(1-(as.numeric(rates.df.merge[rates.df.merge$tropical==0,i])/as.numeric(rates.df.merge[rates.df.merge$tropical==0,'lambda.avg.28k'])))
##  mean.prop.diff.full.10k.trop[[i-2]]<-abs(1-(as.numeric(rates.df.merge[rates.df.merge$tropical==1,i])/as.numeric(rates.df.merge[rates.df.merge$tropical==1,'lambda.avg.28k'])))
##  #cor.test.full.10k.temp[[i-2]]<-cor.test(as.numeric(rates.df.merge[rates.df.merge$tropical==0,i]),as.numeric(rates.df.merge[rates.df.merge$tropical==0,'lambda.avg.28k']),method='s')
##  #cor.test.full.10k.trop[[i-2]]<-cor.test(as.numeric(rates.df.merge[rates.df.merge$tropical==1,i]),as.numeric(rates.df.merge[rates.df.merge$tropical==1,'lambda.avg.28k']),method='s')
##  cor.test.full.10k.temp[[i-2]]<-SpearmanRho(as.numeric(rates.df.merge[rates.df.merge$tropical==0,i]),as.numeric(rates.df.merge[rates.df.merge$tropical==0,'lambda.avg.28k']),use='complete.obs',conf.level=0.95)
##  cor.test.full.10k.trop[[i-2]]<-SpearmanRho(as.numeric(rates.df.merge[rates.df.merge$tropical==1,i]),as.numeric(rates.df.merge[rates.df.merge$tropical==1,'lambda.avg.28k']),use='complete.obs',conf.level=0.95)
##}
###plot(unlist(lapply(cor.test.full.10k.temp,function(x)x$estimate)),unlist(lapply(cor.test.full.10k.trop,function(x)x$estimate)),xlim=c(0.5,1),ylim=c(0.5,1))
###abline(0,1)
##t.test.full.10k.subsampled.temp.trop<-list()
##for (i in c(1:length(mean.prop.diff.full.10k.trop))){
##  t.test.full.10k.subsampled.temp.trop[[i]]<-t.test(mean.prop.diff.full.10k.temp[[i]],mean.prop.diff.full.10k.trop[[i]],na.rm=T)
##}
##
##par(mfrow=c(2,5))
##for (i in c(1:length(mean.prop.diff.full.10k.temp))){
##  boxplot(mean.prop.diff.full.10k.temp[[i]],mean.prop.diff.full.10k.trop[[i]],outline=F,notch=T)
##}
##rates.df.merge$colour<-NA
##rates.df.merge[rates.df.merge$tropical==0,'colour']<-blue
##rates.df.merge[rates.df.merge$tropical==1,'colour']<-red
###this plots Fig S4
##pdf('./BAMMrates_28ktree_vs_10subsamples_10ktrees.pdf',paper='a4r')
##par(mfrow=c(2,5))
##par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
##par(pty="s") 
##lm.trop<-list()
##lm.temp<-list()
##par(tcl = -0.25)
##par(mgp = c(2, 0.6, 0))
##for (i in c(1:length(mean.prop.diff.full.10k.temp))){
##  lm.trop[[i]]<-lm(as.numeric(rates.df.merge[rates.df.merge$tropical==1,'lambda.avg.28k'])~as.numeric(rates.df.merge[rates.df.merge$tropical==1,i+2]))
##  lm.temp[[i]]<-lm(as.numeric(rates.df.merge[rates.df.merge$tropical==0,'lambda.avg.28k'])~as.numeric(rates.df.merge[rates.df.merge$tropical==0,i+2]))
##  
##  plot(as.numeric(rates.df.merge[,2])~as.numeric(rates.df.merge[,i+2]),pch=16,col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1),xaxt='n',yaxt='n',xlab='',ylab='')
##  #points(as.numeric(rates.df.merge[rates.df.merge$tropical==1,2])~as.numeric(rates.df.merge[rates.df.merge$tropical==1,i+2]),pch=16,col=adjustcolor('red',alpha.f = 0.1))
##  abline(lm.trop[[i]],col=red)
##  abline(lm.temp[[i]],col=blue)
##  
##  mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
##  if (i %in% c(6,7,8,9,10)){
##    axis(1, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
##    
##  }else if (i %in% c(1,6)){
##    axis(2, col = "grey40", col.axis = "grey20", at = seq(0,15,5))
##  }
##  box(col = "grey60")
##  
##}
##mtext("rates.10k.subsampled.tree", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
##mtext("rates.28k.tree", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
##dev.off()
##for (i in c(1:length(mean.prop.diff.full.10k.temp))){
##  lm.trop[[i]]<-lm(as.numeric(rates.df.merge[rates.df.merge$tropical==1,'lambda.avg.28k'])~as.numeric(rates.df.merge[rates.df.merge$tropical==1,i+2]))
##  lm.temp[[i]]<-lm(as.numeric(rates.df.merge[rates.df.merge$tropical==0,'lambda.avg.28k'])~as.numeric(rates.df.merge[rates.df.merge$tropical==0,i+2]))
##  
##  plot(as.numeric(rates.df.merge[,2])~as.numeric(rates.df.merge[,i+2]),pch=16,xlab=paste('rates.',gsub(colnames(rates.df.merge)[i+2],pattern='_lambda.avg',replacement=''),sep=''),ylab='rates.28k',col=adjustcolor(rates.df.merge$colour,alpha.f = 0.1))
##  #points(as.numeric(rates.df.merge[rates.df.merge$tropical==1,2])~as.numeric(rates.df.merge[rates.df.merge$tropical==1,i+2]),pch=16,col=adjustcolor('red',alpha.f = 0.1))
##  abline(lm.trop[[i]],col=red)
##  abline(lm.temp[[i]],col=blue)
##  
##}
##dev.off()
##
##
########THIS PLOTS FIG4
##u<-par('usr')
##pdf('./Fig4.pdf')
##
##plot(c(1,1),xlim=c(0,15),ylim=c(0,15),type='n',xlab='rates.subsampled.10k.trees',ylab='rates.28k.tree',main='')
###plot tropical mean slope and intercept
##abline(a=mean(unlist(lapply(lm.trop,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.trop,function(x)x$coefficients[2]))),col=red,lwd=3)
###abline(a=mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.trop,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.trop,function(x)x$coefficients[2]))),col='red',lty=2,add/)
###abline(a=mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.trop,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.trop,function(x)x$coefficients[2]))),col='red',lty=2)
##
###get +1sd -1sd for slope and intercept, build a curve to get points to fill a polygon
##lower.sd.trop <- curve(mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))) + x*(mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.trop,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha.f = 0.5))
##upper.sd.trop <- curve(mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[1])))  + x*(mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.trop,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha.f = 0.5))
##polygon(c(lower.sd.trop$x,rev(upper.sd.trop$x)), c(lower.sd.trop$y, rev(upper.sd.trop$y)),col=adjustcolor(red,alpha.f = 0.5), border=NA)
##
###plot tropical mean slope and intercept
##abline(a=mean(unlist(lapply(lm.temp,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.temp,function(x)x$coefficients[2]))),col=blue,lwd=3)
###abline(a=mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[2]))),col='blue',lty=2)
###abline(a=mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[2]))),col='blue',lty=2)
###get +1sd -1sd for slope and intercept, build a curve to get points to fill a polygon
##lower.sd.temp <- curve(mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))) + x*(mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha.f = 0.5))
##upper.sd.temp <- curve(mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[1])))  + x*(mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha.f = 0.5))
##polygon(c(lower.sd.temp$x,rev(upper.sd.temp$x)), c(lower.sd.temp$y, rev(upper.sd.temp$y)),col=adjustcolor(blue,alpha.f = 0.5), border=NA)
##legend('topleft',legend=c('tropical','temperate'),col=c(red,blue),lty=1,cex=.7,box.lty=0)
##dev.off()
##
##
##
##
##summary(mod1)
##plot(x=c(1:length(cor.test.full.10k.trop)),unlist(lapply(cor.test.full.10k.trop,function(x)x['rho'])),col='red',ylim=c(0,1),pch=16,xlab='replicate',ylab='Spearman rho',main='correlation of rates.28kbiased.vs.10kunbiased',xaxt='n')
##axis(1,at=c(1:10),labels=c(1:10))
##segments(x0=c(1:length(cor.test.full.10k.trop)),y0=unlist(lapply(cor.test.full.10k.trop,function(x)x['lwr.ci'])),y1=unlist(lapply(cor.test.full.10k.trop,function(x)x['upr.ci'])),col='red')
##points(x=c(1:length(cor.test.full.10k.temp)),unlist(lapply(cor.test.full.10k.temp,function(x)x['rho'])),col='blue',ylim=c(0.5,1),pch=16)
##segments(x0=c(1:length(cor.test.full.10k.temp)),y0=unlist(lapply(cor.test.full.10k.temp,function(x)x['lwr.ci'])),y1=unlist(lapply(cor.test.full.10k.temp,function(x)x['upr.ci'])),col='blue')
##abline(1,1)
##
##
##plot(x=c(1:length(t.test.full.10k.subsampled.temp.trop)),unlist(lapply(t.test.full.10k.subsampled.temp.trop,function(x)x$estimate[1])),col='blue',pch=16,xlab='replicate',ylab='mean.diff',main='mean diff of rates.28kbiased.vs.10kunbiased',xaxt='n',ylim=c(0,5))
##points(x=c(1:length(t.test.full.10k.subsampled.temp.trop)),unlist(lapply(t.test.full.10k.subsampled.temp.trop,function(x)x$estimate[2])),col='red',pch=16,xlab='replicate',ylab='mean.diff',main='mean diff of rates.28kbiased.vs.10kunbiased',xaxt='n')
##
##
##axis(1,at=c(1:10),labels=c(1:10))
##segments(x0=c(1:length(cor.test.full.10k.trop)),y0=unlist(lapply(cor.test.full.10k.trop,function(x)x['lwr.ci'])),y1=unlist(lapply(cor.test.full.10k.trop,function(x)x['upr.ci'])),col='red')
##points(x=c(1:length(cor.test.full.10k.temp)),unlist(lapply(cor.test.full.10k.temp,function(x)x['rho'])),col='blue',ylim=c(0.5,1),pch=16)
##segments(x0=c(1:length(cor.test.full.10k.temp)),y0=unlist(lapply(cor.test.full.10k.temp,function(x)x['lwr.ci'])),y1=unlist(lapply(cor.test.full.10k.temp,function(x)x['upr.ci'])),col='blue')
##abline(1,1)
##
##
##plot(density(mean.prop.diff.full.10k.temp[[1]],na.rm=T),col='blue')
##lines(density(mean.prop.diff.full.10k.trop[[1]],na.rm=T),col='red')
##
##
#######OLD, DELETE
##Div_edata_1<-readRDS('./Zanne_clades_BAMM/subsampled_1/results/backbone123456_subsampled_1_eventData.RDS')
##Div_edata_2<-readRDS('./Zanne_clades_BAMM/subsampled_2/results/backbone123456_subsampled_2_eventData.RDS')
##Div_edata_5<-readRDS('./Zanne_clades_BAMM/subsampled_5/results/backbone123456_subsampled_5_eventData.RDS')
##
##ratesdf.subsample1<-getTipRates(Div_edata_1)
##ratesdf.subsample2<-getTipRates(Div_edata_2)
##ratesdf.subsample5<-getTipRates(Div_edata_5)
##ratesdf.subsample.df1<-as.data.frame(cbind(names(ratesdf.subsample1$lambda.avg),ratesdf.subsample1$lambda.avg),stringsAsFactors = F)
##ratesdf.subsample.df2<-as.data.frame(cbind(names(ratesdf.subsample2$lambda.avg),ratesdf.subsample2$lambda.avg),stringsAsFactors = F)
##ratesdf.subsample.df5<-as.data.frame(cbind(names(ratesdf.subsample5$lambda.avg),ratesdf.subsample5$lambda.avg),stringsAsFactors = F)
##row.names(ratesdf.subsample.df1)<-NULL
##row.names(ratesdf.subsample.df2)<-NULL
##row.names(ratesdf.subsample.df5)<-NULL
##colnames(ratesdf.subsample.df1)<-c('species','lambda.avg.subsample1')
##colnames(ratesdf.subsample.df2)<-c('species','lambda.avg.subsample2')
##colnames(ratesdf.subsample.df5)<-c('species','lambda.avg.subsample5')
##
##
##
##source('./R/plots.R')
##GBIFdata<-read.csv('~/Dropbox/Work_in_progress/LDG_plants/GBIFdatasummary.csv')
##GBIFdata$tropical<-0
##GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
##GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
##BAMM.object<-readRDS(file='/Users/javier/Documents/Work/repositories/LDG_plants/all_clades_50million_1000samples_latitudedata.RDS')
###subset table to BAMM species
##GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
###TO CHECK HERE: remove duplicates (6 species with same info but two rows - one with present in garden, one with absent in garden)####
##GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
##ratesdf<-getTipRates(BAMM.object)
##ratesdf.netdiv<-getTipRates(BAMM.object,returnNetDiv = T)
##rates.df<-as.data.frame(cbind(names(ratesdf$lambda.avg),ratesdf$lambda.avg,ratesdf$mu.avg,ratesdf.netdiv$netdiv.avg),stringsAsFactors = F)
##row.names(rates.df)<-NULL
##colnames(rates.df)<-c('species','lambda.avg','mu.avg','netdiv.avg')
##GBIFdata.BAMM.rates<-merge(GBIFdata.BAMM,rates.df,by.x='binomial',by.y='species')
##GBIFdata.BAMM.rates$lambda.avg<-as.numeric(GBIFdata.BAMM.rates$lambda.avg)
##GBIFdata.BAMM.rates$mu.avg<-as.numeric(GBIFdata.BAMM.rates$mu.avg)
##GBIFdata.BAMM.rates$netdiv.avg<-as.numeric(GBIFdata.BAMM.rates$netdiv.avg)
##GBIFdata.BAMM.rates.subsample1<-merge(GBIFdata.BAMM.rates,ratesdf.subsample.df1,by.x='binomial',by.y='species')
##GBIFdata.BAMM.rates.subsample2<-merge(GBIFdata.BAMM.rates,ratesdf.subsample.df2,by.x='binomial',by.y='species')
##GBIFdata.BAMM.rates.subsample5<-merge(GBIFdata.BAMM.rates,ratesdf.subsample.df5,by.x='binomial',by.y='species')
##
##GBIFdata.BAMM.rates.subsample1$lambda.avg.subsample1<-as.numeric(GBIFdata.BAMM.rates.subsample1$lambda.avg.subsample1)
##GBIFdata.BAMM.rates.subsample2$lambda.avg.subsample2<-as.numeric(GBIFdata.BAMM.rates.subsample2$lambda.avg.subsample2)
##GBIFdata.BAMM.rates.subsample5$lambda.avg.subsample5<-as.numeric(GBIFdata.BAMM.rates.subsample5$lambda.avg.subsample5)
##
##rho.subsample1.temp<-cor.test(GBIFdata.BAMM.rates.subsample1[GBIFdata.BAMM.rates.subsample1$tropical==0,'lambda.avg'],GBIFdata.BAMM.rates.subsample1[GBIFdata.BAMM.rates.subsample1$tropical==0,'lambda.avg.subsample1'])
##rho.subsample1.trop<-cor.test(GBIFdata.BAMM.rates.subsample1[GBIFdata.BAMM.rates.subsample1$tropical==1,'lambda.avg'],GBIFdata.BAMM.rates.subsample1[GBIFdata.BAMM.rates.subsample1$tropical==1,'lambda.avg.subsample1'])
##
##rho.subsample2.temp<-cor.test(GBIFdata.BAMM.rates.subsample2[GBIFdata.BAMM.rates.subsample2$tropical==0,'lambda.avg'],GBIFdata.BAMM.rates.subsample2[GBIFdata.BAMM.rates.subsample2$tropical==0,'lambda.avg.subsample2'])
##rho.subsample2.trop<-cor.test(GBIFdata.BAMM.rates.subsample2[GBIFdata.BAMM.rates.subsample2$tropical==1,'lambda.avg'],GBIFdata.BAMM.rates.subsample2[GBIFdata.BAMM.rates.subsample2$tropical==1,'lambda.avg.subsample2'])
##
##rho.subsample5.temp<-cor.test(GBIFdata.BAMM.rates.subsample5[GBIFdata.BAMM.rates.subsample5$tropical==0,'lambda.avg'],GBIFdata.BAMM.rates.subsample5[GBIFdata.BAMM.rates.subsample5$tropical==0,'lambda.avg.subsample5'])
##rho.subsample5.trop<-cor.test(GBIFdata.BAMM.rates.subsample5[GBIFdata.BAMM.rates.subsample5$tropical==1,'lambda.avg'],GBIFdata.BAMM.rates.subsample5[GBIFdata.BAMM.rates.subsample5$tropical==1,'lambda.avg.subsample5'])
##
##
##GBIFdata.BAMM.rates.subsample1$meanproperror<-1-(GBIFdata.BAMM.rates.subsample1$lambda.avg.subsample1/GBIFdata.BAMM.rates.subsample1$lambda.avg)
##t.test(GBIFdata.BAMM.rates.subsample1[GBIFdata.BAMM.rates.subsample1$tropical==0,'meanproperror'],GBIFdata.BAMM.rates.subsample1[GBIFdata.BAMM.rates.subsample1$tropical==1,'meanproperror'])
##
##GBIFdata.BAMM.rates.subsample2$meanproperror<-1-(GBIFdata.BAMM.rates.subsample2$lambda.avg.subsample2/GBIFdata.BAMM.rates.subsample2$lambda.avg)
##t.test(GBIFdata.BAMM.rates.subsample2[GBIFdata.BAMM.rates.subsample2$tropical==0,'meanproperror'],GBIFdata.BAMM.rates.subsample2[GBIFdata.BAMM.rates.subsample2$tropical==1,'meanproperror'])
##
##GBIFdata.BAMM.rates.subsample5$meanproperror<-1-(GBIFdata.BAMM.rates.subsample5$lambda.avg.subsample5/GBIFdata.BAMM.rates.subsample5$lambda.avg)
##t.test(GBIFdata.BAMM.rates.subsample5[GBIFdata.BAMM.rates.subsample5$tropical==0,'meanproperror'],GBIFdata.BAMM.rates.subsample5[GBIFdata.BAMM.rates.subsample5$tropical==1,'meanproperror'])
##
##lm.rates.subsample.1<-lm(lambda.avg.subsample1~lambda.avg*as.factor(tropical),data=GBIFdata.BAMM.rates.subsample1)
##summary(lm.rates.subsample.1)
##lm.rates.subsample.2<-lm(lambda.avg.subsample2~lambda.avg*as.factor(tropical),data=GBIFdata.BAMM.rates.subsample2)
##summary(lm.rates.subsample.2)
##lm.rates.subsample.5<-lm(lambda.avg.subsample5~lambda.avg*as.factor(tropical),data=GBIFdata.BAMM.rates.subsample5)
##summary(lm.rates.subsample.5)
