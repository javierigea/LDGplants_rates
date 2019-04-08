library(plyr)
library(BAMMtools)
library(RColorBrewer)

blue<-brewer.pal(11,'RdYlBu')[9]
red<-brewer.pal(11,'RdYlBu')[2]

#get DR estimates of full 30k tree
DR30k.df<-read.table('./results/DRrates_compare_full_subsampled/DR.30ktree_table.txt',header=T,sep='\t',stringsAsFactors = F)
colnames(DR30k.df)<-c('species','DR.30ktree')
#get DR estimates of 10k trees
DR10k.df.list.files<-list.files(path = './results/DRrates_compare_full_subsampled/',pattern = '*10k*')
DR10k.df.list<-lapply(DR10k.df.list.files,function(x)read.table(paste('./results/DRrates_compare_full_subsampled/',x,sep=''),sep='\t',header=T,stringsAsFactors = F))
for(i in 1:length(DR10k.df.list)){
  cat(i,'\n')
  colnames(DR10k.df.list[[i]])<-c('species',gsub(DR10k.df.list.files,pattern='_table.txt',replacement='')[i])
}
#merge DR 30k estimates with DR 10k estimates
DR30k.df.merge<-DR30k.df
for (i in 1:length(DR10k.df.list)){
  DR30k.df.merge<-merge(DR30k.df.merge,DR10k.df.list[[i]],by.x='species',by.y='species',all.x=T)
}

#add tropical info + other GBIFstuff
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
GBIFdata<-GBIFdata[-which(duplicated(GBIFdata$binomial)),]
GBIFdata<-unique(GBIFdata)
GBIFdata$latitudinal.band<-NA
GBIFdata$latitudinal.band<-apply(GBIFdata,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))

DR30k.df.merge<-merge(DR30k.df.merge,GBIFdata,by.x='species',by.y='binomial',all.x=T)
DR30k.df.merge$colour<-NA
DR30k.df.merge[DR30k.df.merge$tropical==0,'colour']<-blue
DR30k.df.merge[DR30k.df.merge$tropical==1,'colour']<-red

#this plots Fig S4
pdf('./plots/DRrates_28ktree_vs_10subsamples_10ktrees.pdf',paper='a4r')
par(mfrow=c(2,5))
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(pty="s") 
lm.trop<-list()
lm.temp<-list()
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
for (i in c(1:length(DR10k.df.list))){
  lm.trop[[i]]<-lm(as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==1,'DR.30ktree'])~as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==1,i+2]))
  lm.temp[[i]]<-lm(as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==0,'DR.30ktree'])~as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==0,i+2]))
  
  plot(as.numeric(DR30k.df.merge[,2])~as.numeric(DR30k.df.merge[,i+2]),pch=16,col=adjustcolor(DR30k.df.merge$colour,alpha.f = 0.3),xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(0,100),ylim=c(0,100))
  #points(as.numeric(rates.df.merge[rates.df.merge$tropical==1,2])~as.numeric(rates.df.merge[rates.df.merge$tropical==1,i+2]),pch=16,col=adjustcolor('red',alpha.f = 0.1))
  abline(lm.trop[[i]],col=red)
  abline(lm.temp[[i]],col=blue)
  
  mtext(paste(letters[i],')',sep=''), side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")
  if (i %in% c(6,7,8,9,10)){
    axis(1, col = "grey40", col.axis = "grey20", at = seq(0,100,25))
    
  }else if (i %in% c(1,6)){
    axis(2, col = "grey40", col.axis = "grey20", at = seq(0,100,25))
  }
  box(col = "grey60")
  
}
mtext("DR.rates.10k.subsampled.tree", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
mtext("DR.rates.28k.tree", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
dev.off()


####this plots Fig5
pdf('./plots/Fig6.pdf')

plot(c(1,1),type='n',xlab='DR.rates.subsampled.10k.trees',ylab='DR.rates.28k.tree',main='',xlim=c(0,100),ylim=c(0,100))
u<-par('usr')
#plot tropical mean slope and intercept
abline(a=mean(unlist(lapply(lm.trop,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.trop,function(x)x$coefficients[2]))),col=red,lwd=3)
#abline(a=mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.trop,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.trop,function(x)x$coefficients[2]))),col='red',lty=2,add/)
#abline(a=mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.trop,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.trop,function(x)x$coefficients[2]))),col='red',lty=2)

#get +1sd -1sd for slope and intercept, build a curve to get points to fill a polygon
lower.sd.trop <- curve(mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))) + x*(mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.trop,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha.f = 0.5))
upper.sd.trop <- curve(mean(unlist(lapply(lm.trop,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[1])))  + x*(mean(unlist(lapply(lm.trop,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.trop,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(red,alpha.f = 0.5))
polygon(c(lower.sd.trop$x,rev(upper.sd.trop$x)), c(lower.sd.trop$y, rev(upper.sd.trop$y)),col=adjustcolor(red,alpha.f = 0.5), border=NA)

#plot tropical mean slope and intercept
abline(a=mean(unlist(lapply(lm.temp,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.temp,function(x)x$coefficients[2]))),col=blue,lwd=3)
#abline(a=mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[2]))),col='blue',lty=2)
#abline(a=mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))),mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[2]))),col='blue',lty=2)
#get +1sd -1sd for slope and intercept, build a curve to get points to fill a polygon
lower.sd.temp <- curve(mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[1]))) + x*(mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))-sd(unlist(lapply(lm.temp,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha.f = 0.5))
upper.sd.temp <- curve(mean(unlist(lapply(lm.temp,function(x)x$coefficients[1])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[1])))  + x*(mean(unlist(lapply(lm.temp,function(x)x$coefficients[2])))+sd(unlist(lapply(lm.temp,function(x)x$coefficients[2])))), from=u[1], to=u[2],add=T,col=adjustcolor(blue,alpha.f = 0.5))
polygon(c(lower.sd.temp$x,rev(upper.sd.temp$x)), c(lower.sd.temp$y, rev(upper.sd.temp$y)),col=adjustcolor(blue,alpha.f = 0.5), border=NA)
legend('topleft',legend=c('tropical','temperate'),col=c(red,blue),lty=1,cex=.7,box.lty=0)
dev.off()











###check proportional differences in rates 30k vs rates 10k (1- rate sub/rate full)
##mean.prop.diff.full.10k.temp<-list()
##mean.prop.diff.full.10k.trop<-list()
##library(DescTools)
##cor.test.full.10k.temp<-list()
##cor.test.full.10k.trop<-list()
##for (i in c(3:12)){
##  mean.prop.diff.full.10k.temp[[i-2]]<-abs(1-(as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==0,i])/as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==0,'DR.30ktree'])))
##  mean.prop.diff.full.10k.trop[[i-2]]<-abs(1-(as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==1,i])/as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==1,'DR.30ktree'])))
##  cor.test.full.10k.temp[[i-2]]<-SpearmanRho(as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==0,i]),as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==0,'DR.30ktree']),use='complete.obs',conf.level=0.95)
##  cor.test.full.10k.trop[[i-2]]<-SpearmanRho(as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==1,i]),as.numeric(DR30k.df.merge[DR30k.df.merge$tropical==1,'DR.30ktree']),use='complete.obs',conf.level=0.95)
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
##DR30k.df.merge$colour<-NA
##DR30k.df.merge[DR30k.df.merge$tropical==0,'colour']<-blue
##DR30k.df.merge[DR30k.df.merge$tropical==1,'colour']<-red

##DR30k.df.merge$mean.10k<-rowMeans(DR30k.df.merge[,c(2:12)],na.rm = T)
##DR30k.df.merge$median.10k<-apply(DR30k.df.merge[,c(2:12)],1,function(x)median(x,na.rm=T))
##cor.test(DR30k.df.merge$DR.30ktree,DR30k.df.merge$median.10k)
##cor.test(DR30k.df.merge$median.10k,DR30k.df.merge$mean.10k)
##DR30k.df.merge<-merge(DR30k.df.merge,GBIFdata.BAMM,by.x='species',by.y='binomial',all.x=T)
##plot(density(DR30k.df.merge[DR30k.df.merge$tropical==0,'DR.30ktree']),xlim=c(0,5),col='blue',ylim=c(0,10))
##lines(density(DR30k.df.merge[DR30k.df.merge$tropical==1,'DR.30ktree']),col='red')
##
##cor.test(DR30k.df.merge[DR30k.df.merge$tropical==0,]$DR.30ktree,DR30k.df.merge[DR30k.df.merge$tropical==0,]$median.10k)
##cor.test(DR30k.df.merge[DR30k.df.merge$tropical==1,]$DR.30ktree,DR30k.df.merge[DR30k.df.merge$tropical==1,]$median.10k)
##
##cor.test(DR30k.df.merge[DR30k.df.merge$tropical==0,]$DR.30ktree,DR30k.df.merge[DR30k.df.merge$tropical==0,]$mean.10k)
##cor.test(DR30k.df.merge[DR30k.df.merge$tropical==1,]$DR.30ktree,DR30k.df.merge[DR30k.df.merge$tropical==1,]$mean.10k)
##
##cor.test.trop<-list()
##cor.test.temp<-list()
##for (i in c(3:12)){
##  cor.test.temp[[i-2]]<-cor.test(DR30k.df.merge[DR30k.df.merge$tropical==0,'DR.30ktree'],DR30k.df.merge[DR30k.df.merge$tropical==0,i],method='s')
##  cor.test.trop[[i-2]]<-cor.test(DR30k.df.merge[DR30k.df.merge$tropical==1,'DR.30ktree'],DR30k.df.merge[DR30k.df.merge$tropical==1,i],method='s')
##}
##plot(x=c(1:10),unlist(lapply(cor.test.trop,function(x)x$estimate)),col='red',ylim=c(0.7,1),pch=16,xlab='replicate',ylab='Pearsons r')
##segments(x0=c(1:10),y0=unlist(lapply(cor.test.trop,function(x)x$conf.int[1])),y1=unlist(lapply(cor.test.trop,function(x)x$conf.int[2])),col='red')
##points(x=c(1:10),unlist(lapply(cor.test.temp,function(x)x$estimate)),col='blue',ylim=c(0.7,1),pch=16)
##segments(x0=c(1:10),y0=unlist(lapply(cor.test.temp,function(x)x$conf.int[1])),y1=unlist(lapply(cor.test.temp,function(x)x$conf.int[2])),col='blue')
##abline(1,1)
##
##plot(density(DR30k.df.merge[(DR30k.df.merge$DR.30ktree<5)&(DR30k.df.merge$tropical==0),]$DR.30ktree),col='black',ylim=c(0,5))
##for(i in c(3:12)){
##  lines(density(DR30k.df.merge[which(DR30k.df.merge[,i]<5&DR30k.df.merge$tropical==0),i]),col='blue')
##}
##
##lines(density(DR30k.df.merge[(DR30k.df.merge$DR.30ktree<5)&(DR30k.df.merge$tropical==1),]$DR.30ktree),col='red',ylim=c(0,5))
##for(i in c(3:12)){
##  lines(density(DR30k.df.merge[which(DR30k.df.merge[,i]<5&DR30k.df.merge$tropical==1),i]),col='red')
##}
##
##cor.test.trop.rho<-list()
##cor.test.temp.rho<-list()
##for (i in c(3:12)){
##  cor.test.temp.rho[[i-2]]<-SpearmanRho(DR30k.df.merge[DR30k.df.merge$tropical==0,'DR.30ktree'],DR30k.df.merge[DR30k.df.merge$tropical==0,i],use='complete',conf.level = 0.95)
##  cor.test.trop.rho[[i-2]]<-SpearmanRho(DR30k.df.merge[DR30k.df.merge$tropical==1,'DR.30ktree'],DR30k.df.merge[DR30k.df.merge$tropical==1,i],use='complete',conf.level = 0.95)
##}
##library(DescTools)
##spearman.ci(var1, var2, nrep = 1000, conf.level = 0.95)
##
##mean.prop.diff.full.10k.temp<-list()
##mean.prop.diff.full.10k.trop<-list()
##
##for (i in c(3:12)){
##  mean.prop.diff.full.10k.temp[[i-2]]<-1-(DR30k.df.merge[DR30k.df.merge$tropical==0,'DR.30ktree']/DR30k.df.merge[DR30k.df.merge$tropical==0,i])
##  mean.prop.diff.full.10k.trop[[i-2]]<-1-(DR30k.df.merge[DR30k.df.merge$tropical==1,'DR.30ktree']/DR30k.df.merge[DR30k.df.merge$tropical==1,i])
##}
##
##