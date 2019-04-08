.libPaths('~/R/x86_64-pc-linux-gnu-library/3.5')
setwd('~/ldg_plants/')

library(BAMMtools)
library(plyr)

args<-commandArgs(trailingOnly = TRUE)
print(args)
replicate<-as.numeric(args[1])



Div_edata<-readRDS(paste('./results/extremebias/all_clades_50million_1000samples_latitudedata_extremebias_',replicate,'.RDS',sep=''))
trait_data<-read.csv('./raw_data/GBIFdatasummary.csv')
trait_data$binomial<-paste(trait_data$Genus.Name,trait_data$Species.Name,sep='_')
trait_data<-trait_data[-which(duplicated(trait_data$binomial)),]
trait_data<-unique(trait_data)

trait_data<-trait_data[trait_data$binomial%in%Div_edata$tip.label,]
trait_data<-trait_data[,c('binomial','Median.Latitude')]
trait_data$Median.Latitude<-abs(trait_data$Median.Latitude)
colnames(trait_data)<-c('V1','V2')
parameter<-'speciation'
spearman.regressions.df<-matrix(nrow=1000,ncol=2)
#this is set to 1000 -  for 1000 strapp permutations, change if different
#for all strapp permutations
#for (i in 1:10){
for (i in 1:1000){
  cat(i, '\n')	 
  #run through posterior
  gen<-i
  #get parameters from eventobject
  strapp.eventdata <- cbind(Div_edata$edge, Div_edata$eventVectors[[gen]])
  strapp.eventdata <- strapp.eventdata[strapp.eventdata[,2] <= length(Div_edata$tip.label), ]
  strapp.eventdata <- data.frame(Div_edata$tip.label, strapp.eventdata[,3], Div_edata$tipLambda[[gen]],Div_edata$tipMu[[gen]],(Div_edata$tipLambda[[gen]]-Div_edata$tipMu[[gen]]))
  colnames(strapp.eventdata)<-c('tip.label','event','tip.lambda','tip.mu','tip.div')
  
  #add trait data to parameters
  strapp.eventdata<-merge(strapp.eventdata,trait_data,by.x='tip.label','V1')
  strapp.eventdata$V2<-strapp.eventdata$V2
  #this vector determines the number of correlations to test (we need to get the value where the correlation equals 0)
  #see http://stats.stackexchange.com/questions/64938/if-linear-regression-is-related-to-pearsons-correlation-are-there-any-regressi/110112#110112
  x<-seq(-10,10,by=0.005)
  beta<-vector('numeric')
  if(parameter=='speciation'){
    for (a in 1:length(x)){
      #cat(a,'\n')
      beta[a]<-cor(log10(strapp.eventdata$tip.lambda)-(x[a]*strapp.eventdata$V2),strapp.eventdata$V2,method='s')  
    }
  }
  if(parameter=='extinction'){
    for (a in 1:length(x)){
      beta[a]<-cor(log10(strapp.eventdata$tip.mu)-(x[a]*strapp.eventdata$V2),strapp.eventdata$V2,method='spearman')  
    }
  }
  if(parameter=='net diversification'){
    #this removes any negative values from diversification (they're only 0.3% of the rows approx)
    strapp.eventdata<-strapp.eventdata[strapp.eventdata$tip.div>0,]
    for (a in 1:length(x)){
      beta[a]<-cor(log10(strapp.eventdata$tip.div)-(x[a]*strapp.eventdata$V2),strapp.eventdata$V2,method='spearman')  
    }
  }
  
  df<-data.frame(x,beta)
  df<-df[complete.cases(df),]
  
  #get the value of beta that's closest to 0
  value<-which.min(abs(df$beta - 0))
  slope<-df[df$beta==df$beta[value],]$x
  if (parameter=='speciation'){
    model<-lm(log10(strapp.eventdata$tip.lambda)-(slope*strapp.eventdata$V2)~0)  
  }
  if (parameter=='extinction'){
    model<-lm(log10(strapp.eventdata$tip.mu)-(slope*strapp.eventdata$V2)~0)  
  }
  if (parameter=='net diversification'){
    model<-lm(log10(strapp.eventdata$tip.div)-(slope*strapp.eventdata$V2)~0)  
  }
  intercept<-median(model$residuals)
  spearman.regressions.df[i,1]<-intercept
  spearman.regressions.df[i,2]<-slope 
}    
spearman.regressions.df<-data.frame(spearman.regressions.df)
if (parameter=='speciation'){
  output.parameter<-'lambda'
}
if (parameter=='extinction'){
  output.parameter<-'mu'
}
if (parameter=='net diversification'){
  output.parameter<-'div'
}
correlationfile<-paste('./results/extremebias/strapp_',output.parameter,'_absMedianLatitude_',replicate,'.txt',sep='')
write.table(spearman.regressions.df,correlationfile,row.names=F,sep='\t',quote=F)
