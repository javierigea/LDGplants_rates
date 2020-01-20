
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
GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
GBIFdata$tropical<-0
GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
GBIFdata<-GBIFdata[-which(duplicated(GBIFdata$binomial)),]
GBIFdata<-unique(GBIFdata)
BAMM.object<-readRDS('./raw_data/all_clades_50million_1000samples_latitudedata_GBOTB.RDS')
GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
GBIFdata.BAMM$latitudinal.band<-NA
GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))

GBIFdata.BAMM$latitudinal.band.new<-apply(GBIFdata.BAMM,1,function(x) assign_equalwidth_band(number=as.numeric(x['Median.Latitude'])))
#change the names of the bands to multiples of ten
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='-61.1','latitudinal.band.new']<-(-60)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='-51.7','latitudinal.band.new']<-(-50)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='-42.3','latitudinal.band.new']<-(-40)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='-32.9','latitudinal.band.new']<-(-30)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='-23.5','latitudinal.band.new']<-(-20)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='-14.1','latitudinal.band.new']<-(-10)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='-4.7','latitudinal.band.new']<-(0)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='4.7','latitudinal.band.new']<-(10)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='14.1','latitudinal.band.new']<-(20)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='23.5','latitudinal.band.new']<-(30)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='32.9','latitudinal.band.new']<-(40)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='42.3','latitudinal.band.new']<-(50)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='51.7','latitudinal.band.new']<-(60)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='61.1','latitudinal.band.new']<-(70)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='70.5','latitudinal.band.new']<-(80)
GBIFdata.BAMM[GBIFdata.BAMM$latitudinal.band.new=='79.9','latitudinal.band.new']<-(90)
#use 9.4 degree width bands as the default latitudinal bands
GBIFdata.BAMM$latitudinal.band<-GBIFdata.BAMM$latitudinal.band.new
ratesdf<-getTipRates(BAMM.object)
ratesdf.netdiv<-getTipRates(BAMM.object,returnNetDiv = T)
rates.df<-as.data.frame(cbind(names(ratesdf$lambda.avg),ratesdf$lambda.avg,ratesdf$mu.avg,ratesdf.netdiv$netdiv.avg),stringsAsFactors = F)
row.names(rates.df)<-NULL
colnames(rates.df)<-c('species','lambda.avg','mu.avg','netdiv.avg')
GBIFdata.BAMM.rates<-merge(GBIFdata.BAMM,rates.df,by.x='binomial',by.y='species')
GBIFdata.BAMM.rates$lambda.avg<-as.numeric(GBIFdata.BAMM.rates$lambda.avg)
GBIFdata.BAMM.rates$mu.avg<-as.numeric(GBIFdata.BAMM.rates$mu.avg)
GBIFdata.BAMM.rates$netdiv.avg<-as.numeric(GBIFdata.BAMM.rates$netdiv.avg)
write.table(GBIFdata.BAMM.rates,file='./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',sep='\t',quote=F,row.names=F)


