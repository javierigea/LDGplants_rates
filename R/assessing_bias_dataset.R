

plot_GBIFdata_proportions<-function(GBIFdata.file){
  blue<-brewer.pal(11,'RdYlBu')[9]
  red<-brewer.pal(11,'RdYlBu')[2]
  #get GBIF data and build tropical variable (=median latitude)
  GBIFdata<-read.csv(GBIFdata.file)
  GBIFdata$tropical<-0
  GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
  GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
  table(GBIFdata$tropical)/nrow(GBIFdata)
  pie(table(GBIFdata$tropical),labels = paste(c('temperate','tropical'),table(GBIFdata$tropical),sep=' '), col=c(adjustcolor(blue,0.3),adjustcolor(red,0.3)),main="GBIF data")
  #group into latitudinal bands
  GBIFdata$latitudinal.band<-NA
  GBIFdata$latitudinal.band<-apply(GBIFdata,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
  plot(names(table(GBIFdata$latitudinal.band)),table(GBIFdata$latitudinal.band)/nrow(GBIFdata),xlab='median latitude',ylab='proportion',ylim=c(0,0.2),main='proportions.GBIFdata',xaxt='n',yaxt='n')
  lines(names(table(GBIFdata$latitudinal.band)),table(GBIFdata$latitudinal.band)/nrow(GBIFdata))
  axis(1,at=seq(from=-60,to=90,by=10),labels=seq(from=-60,to=90,by=10))
  axis(2,at=seq(from=0,to=0.2,by=0.05),labels=seq(from=0,to=0.2,by=0.05),las=2)
  abline(v=-23.5,lty=2)
  abline(v=23.5,lty=2)
  #absolute latitudinal_band
  GBIFdata$abs.latitudinal.band<-abs(GBIFdata$latitudinal.band)
  plot(names(table(GBIFdata$abs.latitudinal.band)),table(GBIFdata$abs.latitudinal.band)/nrow(GBIFdata),xlab='median absolute latitude',ylab='proportion',ylim=c(0,0.3),main='proportions.GBIFdata',xaxt='n',yaxt='n')
  lines(names(table(GBIFdata$abs.latitudinal.band)),table(GBIFdata$abs.latitudinal.band)/nrow(GBIFdata))
  axis(1,at=seq(from=-60,to=90,by=10),labels=seq(from=-60,to=90,by=10))
  axis(2,at=seq(from=0,to=0.3,by=0.05),labels=seq(from=0,to=0.3,by=0.05),las=2)
  abline(v=23.5,lty=2)
  #plot number of species in each latitudinal band
  plot(as.numeric(unname(table(GBIFdata$latitudinal.band))),as.numeric(names(table(GBIFdata$latitudinal.band))),xlab='number.of.species',ylab='median.latitude',ylim=c(-70,90),main='GBIFdata.nspecies',xaxt='n',yaxt='n')
  lines(as.numeric(unname(table(GBIFdata$latitudinal.band))),as.numeric(names(table(GBIFdata$latitudinal.band))))
  axis(1,at=seq(from=0,to=30000,by=5000),labels=seq(from=0,to=30000,by=5000),las=1)
  axis(2,at=seq(from=-70,to=90,by=10),labels=seq(from=-70,to=90,by=10),las=2)
  abline(h=-23.5,lty=2)
  abline(h=23.5,lty=2)
  
}


plot_GBIFdataBAMM_proportions<-function(GBIFdata.file,BAMM.object.path){
  blue<-brewer.pal(11,'RdYlBu')[9]
  red<-brewer.pal(11,'RdYlBu')[2]
  GBIFdata<-read.csv(GBIFdata.file)
  GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
  Div_edata.spp.lat<-readRDS(BAMM.object.path)
  #get GBIF data and build tropical variable (=median latitude)
  GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%Div_edata.spp.lat$tip.label,]
  GBIFdata.BAMM$tropical<-0
  GBIFdata.BAMM[(GBIFdata.BAMM$Median.Latitude<23.5)&(GBIFdata.BAMM$Median.Latitude>(-23.5)),'tropical']<-1
  GBIFdata.BAMM$binomial<-paste(GBIFdata.BAMM$Genus.Name,GBIFdata.BAMM$Species.Name,sep='_')
  GBIFdata.BAMM$latitudinal.band<-NA
  GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
  GBIFdata.BAMM$abs.latitudinal.band<-abs(GBIFdata.BAMM$latitudinal.band)
  
  table(GBIFdata.BAMM$tropical)/nrow(GBIFdata.BAMM)
  pie(table(GBIFdata.BAMM$tropical),labels = paste(c('temperate','tropical'),table(GBIFdata.BAMM$tropical),sep=' '), col=c(adjustcolor(blue,0.3),adjustcolor(red,0.3)),main="GBIFdata.BAMM")
  #group into latitudinal bands
  plot(names(table(GBIFdata.BAMM$latitudinal.band)),table(GBIFdata.BAMM$latitudinal.band)/nrow(GBIFdata.BAMM),xlab='median latitude',ylab='proportion',ylim=c(0,0.2),main='proportions.GBIFdata.BAMM',xaxt='n',yaxt='n')
  lines(names(table(GBIFdata.BAMM$latitudinal.band)),table(GBIFdata.BAMM$latitudinal.band)/nrow(GBIFdata.BAMM))
  axis(1,at=seq(from=-60,to=90,by=10),labels=seq(from=-60,to=90,by=10))
  axis(2,at=seq(from=0,to=0.2,by=0.05),labels=seq(from=0,to=0.2,by=0.05),las=2)
  abline(v=-23.5,lty=2)
  abline(v=23.5,lty=2)
  
  #absolute latitudinal_band
  plot(names(table(GBIFdata.BAMM$abs.latitudinal.band)),table(GBIFdata.BAMM$abs.latitudinal.band)/nrow(GBIFdata.BAMM),xlab='median absolute latitude',ylab='proportion',ylim=c(0,0.3),main='proportions.GBIFdata.BAMM',xaxt='n',yaxt='n')
  lines(names(table(GBIFdata.BAMM$abs.latitudinal.band)),table(GBIFdata.BAMM$abs.latitudinal.band)/nrow(GBIFdata.BAMM))
  axis(1,at=seq(from=-60,to=90,by=10),labels=seq(from=-60,to=90,by=10))
  axis(2,at=seq(from=0,to=0.3,by=0.05),labels=seq(from=0,to=0.3,by=0.05),las=2)
  abline(v=23.5,lty=2)
  
  #plot number of species in each latitudinal band
  plot(as.numeric(unname(table(GBIFdata.BAMM$latitudinal.band))),as.numeric(names(table(GBIFdata.BAMM$latitudinal.band))),xlab='number.of.species',ylab='median.latitude',ylim=c(-70,90),main='GBIFdata.BAMM.nspecies',xaxt='n',yaxt='n')
  lines(as.numeric(unname(table(GBIFdata.BAMM$latitudinal.band))),as.numeric(names(table(GBIFdata.BAMM$latitudinal.band))))
  axis(1,at=seq(from=0,to=5000,by=1000),labels=seq(from=0,to=5000,by=1000),las=1)
  axis(2,at=seq(from=-70,to=90,by=10),labels=seq(from=-70,to=90,by=10),las=2)
  abline(h=-23.5,lty=2)
  abline(h=23.5,lty=2)
  
}


plot_GBIFdata_proportions_subsample<-function(GBIFdata.file,threshold){
  blue<-brewer.pal(11,'RdYlBu')[9]
  red<-brewer.pal(11,'RdYlBu')[2]
  #get GBIF data and build tropical variable (=median latitude)
  GBIFdata<-read.csv(GBIFdata.file)
  GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
  #group into latitudinal bands
  GBIFdata$latitudinal.band<-NA
  GBIFdata$latitudinal.band<-apply(GBIFdata,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
  GBIFdata<-GBIFdata[abs(GBIFdata$Median.Latitude)<=threshold,]
  GBIFdata$tropical<-0
  GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
  
  table(GBIFdata$tropical)/nrow(GBIFdata)
  pie(table(GBIFdata$tropical),labels = paste(c('temperate','tropical'),table(GBIFdata$tropical),sep=' '), col=c(adjustcolor(blue,0.3),adjustcolor('red',0.3)),main="GBIF data")
  
  plot(names(table(GBIFdata$latitudinal.band)),table(GBIFdata$latitudinal.band)/nrow(GBIFdata),xlab='median latitude',ylab='proportion',ylim=c(0,0.2),main='proportions.GBIFdata',xaxt='n',yaxt='n')
  lines(names(table(GBIFdata$latitudinal.band)),table(GBIFdata$latitudinal.band)/nrow(GBIFdata))
  axis(1,at=seq(from=min(GBIFdata$latitudinal.band),to=max(GBIFdata$latitudinal.band),by=10),labels=seq(from=min(GBIFdata$latitudinal.band),to=max(GBIFdata$latitudinal.band),by=10))
  axis(2,at=seq(from=0,to=0.2,by=0.05),labels=seq(from=0,to=0.2,by=0.05),las=2)
  abline(v=-23.5,lty=2)
  abline(v=23.5,lty=2)
  #absolute latitudinal_band
  GBIFdata$abs.latitudinal.band<-abs(GBIFdata$latitudinal.band)
  plot(names(table(GBIFdata$abs.latitudinal.band)),table(GBIFdata$abs.latitudinal.band)/nrow(GBIFdata),xlab='median absolute latitude',ylab='proportion',ylim=c(0,0.3),main='proportions.GBIFdata',xaxt='n',yaxt='n')
  lines(names(table(GBIFdata$abs.latitudinal.band)),table(GBIFdata$abs.latitudinal.band)/nrow(GBIFdata))
  axis(1,at=seq(from=-60,to=90,by=10),labels=seq(from=-60,to=90,by=10))
  axis(2,at=seq(from=0,to=0.3,by=0.05),labels=seq(from=0,to=0.3,by=0.05),las=2)
  abline(v=23.5,lty=2)
}


plot_GBIFdataBAMM_proportions_subsample<-function(GBIFdata.file,BAMM.object.path,threshold){
  blue<-brewer.pal(11,'RdYlBu')[9]
  red<-brewer.pal(11,'RdYlBu')[2]
  
  GBIFdata<-read.csv(GBIFdata.file)
  GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
  Div_edata.spp.lat<-readRDS(BAMM.object.path)
  #get GBIF data and build tropical variable (=median latitude)
  GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%Div_edata.spp.lat$tip.label,]
  GBIFdata.BAMM$tropical<-0
  GBIFdata.BAMM[(GBIFdata.BAMM$Median.Latitude<23.5)&(GBIFdata.BAMM$Median.Latitude>(-23.5)),'tropical']<-1
  GBIFdata.BAMM$binomial<-paste(GBIFdata.BAMM$Genus.Name,GBIFdata.BAMM$Species.Name,sep='_')
  GBIFdata.BAMM$latitudinal.band<-NA
  GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
  GBIFdata.BAMM<-GBIFdata.BAMM[abs(GBIFdata.BAMM$Median.Latitude)<=threshold,]
  GBIFdata.BAMM$abs.latitudinal.band<-abs(GBIFdata.BAMM$latitudinal.band)
  
  table(GBIFdata.BAMM$tropical)/nrow(GBIFdata.BAMM)
  pie(table(GBIFdata.BAMM$tropical),labels = paste(c('temperate','tropical'),table(GBIFdata.BAMM$tropical),sep=' '), col=c(adjustcolor(blue,0.3),adjustcolor(red,0.3)),main="GBIFdata BAMM")
  #group into latitudinal bands
  plot(names(table(GBIFdata.BAMM$latitudinal.band)),table(GBIFdata.BAMM$latitudinal.band)/nrow(GBIFdata.BAMM),xlab='median latitude',ylab='proportion',ylim=c(0,0.2),main='proportions.GBIFdata.BAMM',xaxt='n',yaxt='n')
  lines(names(table(GBIFdata.BAMM$latitudinal.band)),table(GBIFdata.BAMM$latitudinal.band)/nrow(GBIFdata.BAMM))
  axis(1,at=seq(from=min(GBIFdata.BAMM$latitudinal.band),to=max(GBIFdata.BAMM$latitudinal.band),by=10),labels=seq(from=min(GBIFdata.BAMM$latitudinal.band),to=max(GBIFdata.BAMM$latitudinal.band),by=10))
  axis(2,at=seq(from=0,to=0.2,by=0.05),labels=seq(from=0,to=0.2,by=0.05),las=2)
  abline(v=-23.5,lty=2)
  abline(v=23.5,lty=2)
  
  #absolute latitudinal_band
  plot(names(table(GBIFdata.BAMM$abs.latitudinal.band)),table(GBIFdata.BAMM$abs.latitudinal.band)/nrow(GBIFdata.BAMM),xlab='median absolute latitude',ylab='proportion',ylim=c(0,0.3),main='proportions.GBIFdata.BAMM',xaxt='n',yaxt='n')
  lines(names(table(GBIFdata.BAMM$abs.latitudinal.band)),table(GBIFdata.BAMM$abs.latitudinal.band)/nrow(GBIFdata.BAMM))
  axis(1,at=seq(from=min(GBIFdata.BAMM$abs.latitudinal.band),to=max(GBIFdata.BAMM$abs.latitudinal.band),by=10),labels=seq(from=min(GBIFdata.BAMM$abs.latitudinal.band),to=max(GBIFdata.BAMM$abs.latitudinal.band),by=10))
  axis(2,at=seq(from=0,to=0.3,by=0.05),labels=seq(from=0,to=0.3,by=0.05),las=2)
  abline(v=23.5,lty=2)
}
