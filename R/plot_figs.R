library(RColorBrewer)
library(plyr)
library(BAMMtools)
library(rgdal)
library(rgeos)

blue<-brewer.pal(11,'RdYlBu')[9]
red<-brewer.pal(11,'RdYlBu')[2]

plot_Fig1_map_bands<-function(GBIFdata.BAMM.rates){
  blue<-brewer.pal(11,'RdYlBu')[9]
  red<-brewer.pal(11,'RdYlBu')[2]
  
  #GBIFdata.BAMM.rates is a table with species latitudinal.band,lambda.avg
  #####CHANGE PATH OF REALMS LAYER####
  realms.layer<-readOGR('./raw_data/realms/','newRealms')
  realms.layer<-gUnaryUnion(realms.layer)
  #pdf('./plots/Fig1_map_latitudinalbands_log.pdf',paper='a4r')
  plot(realms.layer,lwd=0.5)
  #for (i in seq(from=85,to=-55,by=-10)){
  #  abline(h=i,col='grey70')
  #}
  abline(v=0)
  abline(h=23.5,col='black')
  abline(h=-23.5,col='black')
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-55,-55,-23.5,-23.5,-55),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-23.5,-23.5,23.5,23.5,-23.5),col=adjustcolor(red,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(85,85,23.5,23.5,85),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  u<-par('usr')
  #plot(c(0,0),type='n',xlim=c(u[1]/30,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #plot(c(0,0),type='n',xlim=c(u[1]/50,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  plot(c(0,0),type='n',xlim=c(-3,3),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  axis(1,c(log(c(0.25,1,2,5,10))),labels=c(0.25,1,2,5,10))
  axis(4,seq(from=-50/10,to=80/10,by=10/10),labels=seq(from=-50,to=80,by=10),las=2,cex.axis=.5)
  #axis(1,c(log(0),log(3),log(6),log(9),log(12),log(15)),labels=c(0,3,6,9,12,15))
  abline(v=0)
  #par(new=T)
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='10','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='0','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=0/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-10','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-20','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-30','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-40','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-50','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=-50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='20','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='30','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='40','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='50','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='60','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=60/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='70','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=70/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='80','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=80/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  #dev.off()
  
}

plot_Fig1_map_altbands<-function(GBIFdata.BAMM.rates){
  blue<-brewer.pal(11,'RdYlBu')[9]
  red<-brewer.pal(11,'RdYlBu')[2]
  
  #GBIFdata.BAMM.rates is a table with species latitudinal.band,lambda.avg
  #####CHANGE PATH OF REALMS LAYER####
  realms.layer<-readOGR('./raw_data/realms/','newRealms')
  realms.layer<-gUnaryUnion(realms.layer)
  #pdf('./plots/Fig1_map_latitudinalbands_log.pdf',paper='a4r')
  plot(realms.layer,lwd=0.5)
  #for (i in seq(from=85,to=-55,by=-10)){
  #  abline(h=i,col='grey70')
  #}
  abline(v=0)
  abline(h=23.5,col='black')
  abline(h=-23.5,col='black')
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-55,-55,-23.5,-23.5,-55),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-23.5,-23.5,23.5,23.5,-23.5),col=adjustcolor(red,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(85,85,23.5,23.5,85),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  u<-par('usr')
  #plot(c(0,0),type='n',xlim=c(u[1]/30,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #plot(c(0,0),type='n',xlim=c(u[1]/50,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  plot(c(0,0),type='n',xlim=c(-2,1.5),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  axis(1,c(log(c(0.25,1,2,5))),labels=c(0.25,1,2,5))
  axis(4,seq(from=-50/10,to=80/10,by=10/10),labels=seq(from=-50,to=80,by=10),las=2,cex.axis=.5)
  #axis(1,c(log(0),log(3),log(6),log(9),log(12),log(15)),labels=c(0,3,6,9,12,15))
  abline(v=0)
  #par(new=T)
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='10','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='0','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=0/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-10','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-20','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-30','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-40','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-50','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=-50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='20','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='30','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='40','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='50','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='60','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=60/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='70','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=70/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='80','lambda.avg']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=80/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  #dev.off()
  
}

plot_Fig1_map_altbands_DR<-function(GBIFdata.BAMM.rates){
  blue<-brewer.pal(11,'RdYlBu')[9]
  red<-brewer.pal(11,'RdYlBu')[2]
  
  #GBIFdata.BAMM.rates is a table with species latitudinal.band,lambda.avg
  #####CHANGE PATH OF REALMS LAYER####
  realms.layer<-readOGR('./raw_data/realms/','newRealms')
  realms.layer<-gUnaryUnion(realms.layer)
  #pdf('./plots/Fig1_map_latitudinalbands_log.pdf',paper='a4r')
  plot(realms.layer,lwd=0.5)
  #for (i in seq(from=85,to=-55,by=-10)){
  #  abline(h=i,col='grey70')
  #}
  abline(v=0)
  abline(h=23.5,col='black')
  abline(h=-23.5,col='black')
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-55,-55,-23.5,-23.5,-55),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-23.5,-23.5,23.5,23.5,-23.5),col=adjustcolor(red,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(85,85,23.5,23.5,85),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  u<-par('usr')
  #plot(c(0,0),type='n',xlim=c(u[1]/30,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #plot(c(0,0),type='n',xlim=c(u[1]/50,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  plot(c(0,0),type='n',xlim=c(-3,1),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  axis(1,c(log(c(0.1,0.25,0.5,1,2))),labels=c(0.1,0.25,0.5,1,2))
  axis(4,seq(from=-50/10,to=80/10,by=10/10),labels=seq(from=-50,to=80,by=10),las=2,cex.axis=.5)
  #axis(1,c(log(0),log(3),log(6),log(9),log(12),log(15)),labels=c(0,3,6,9,12,15))
  abline(v=0)
  #par(new=T)
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='10','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='0','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=0/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-10','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-20','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-30','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-40','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='-50','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=-50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='20','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='30','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='40','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='50','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='60','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=60/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='70','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=70/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$latitudinal.band=='80','DR']),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=80/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  #dev.off()
  
}
plot_Fig1_boxplot_troptempbinary<-function(GBIFdata.BAMM.rates){

  blue<-brewer.pal(11,'RdYlBu')[9]
  red<-brewer.pal(11,'RdYlBu')[2]
  #boxplot(lambda.avg~tropical,data=GBIFdata.BAMM.rates,ylim=c(0,5),outline=F,whisklty=0,staplelty=0,notch=T,names=c('temperate','tropical'),col=c(blue,red),ylab='BAMM.tip.lambda',main='BAMM.allspecies.28057.lambda')
  boxplot(log(lambda.avg)~tropical,data=GBIFdata.BAMM.rates,ylim=c(log(0.15),log(10)),outline=F,whisklty=0,staplelty=0,notch=T,names=c('temperate','tropical'),col=c(blue,red),ylab='BAMM.tip.lambda',main='BAMM.allspecies.28057.lambda',yaxt='n')
  axis(2,at=c(log(c(0.25,1,2,5))),labels=c(0.25,1,2,5))

}

plot_Fig1_boxplot_troptempbinary_GBOTB<-function(GBIFdata.BAMM.rates){
  
  blue<-brewer.pal(11,'RdYlBu')[9]
  red<-brewer.pal(11,'RdYlBu')[2]
  #boxplot(lambda.avg~tropical,data=GBIFdata.BAMM.rates,ylim=c(0,5),outline=F,whisklty=0,staplelty=0,notch=T,names=c('temperate','tropical'),col=c(blue,red),ylab='BAMM.tip.lambda',main='BAMM.allspecies.28057.lambda')
  boxplot(log(lambda.avg)~tropical,data=GBIFdata.BAMM.rates,ylim=c(-2,1.5),outline=F,whisklty=0,staplelty=0,notch=T,names=c('temperate','tropical'),col=c(blue,red),ylab='BAMM.tip.lambda',main='BAMM.allspecies.28057.lambda',yaxt='n')
  axis(2,at=c(log(c(0.25,1,2,5))),labels=c(0.25,1,2,5))
  
}

plot_Fig1_DRboxplot_troptempbinary_GBOTB<-function(GBIFdata.BAMM.rates){
  
  blue<-brewer.pal(11,'RdYlBu')[9]
  red<-brewer.pal(11,'RdYlBu')[2]
  #boxplot(lambda.avg~tropical,data=GBIFdata.BAMM.rates,ylim=c(0,5),outline=F,whisklty=0,staplelty=0,notch=T,names=c('temperate','tropical'),col=c(blue,red),ylab='BAMM.tip.lambda',main='BAMM.allspecies.28057.lambda')
  boxplot(log(DR)~tropical,data=GBIFdata.BAMM.rates,ylim=c(-3,1),outline=F,whisklty=0,staplelty=0,notch=T,names=c('temperate','tropical'),col=c(blue,red),ylab='DR',main='DR.all.species',yaxt='n')
  axis(2,at=c(log(c(0.1,0.25,0.5,1,2,5))),labels=c(0.1,0.25,0.5,1,2,5),las=2)
  
}


#plot Fig1 quantile rate distrib vs loglambdarate
plot_Fig1_quantilerate_vs_lambdarate<-function(GBIFdata.BAMM.rates){
  #GBIFdata<-read.csv('./raw_data/GBIFdatasummary.csv')
  #GBIFdata$tropical<-0
  #GBIFdata[(GBIFdata$Median.Latitude<23.5)&(GBIFdata$Median.Latitude>(-23.5)),'tropical']<-1
  #GBIFdata$strict.tropical<-NA
  ##strict tropical 1 = abs(max latitude)<23.5 & abs(min latitude)<23.5  & median.latitude <23.5
  #GBIFdata[abs(GBIFdata$Max.Latitude)<=23.5&abs(GBIFdata$Min.Latitude)<=23.5&abs(GBIFdata$Median.Latitude)<=23.5,'strict.tropical']<-1
  ##strict tropical 0 (strict temperate)
  #GBIFdata[abs(GBIFdata$Max.Latitude)>23.5&abs(GBIFdata$Min.Latitude)>23.5&abs(GBIFdata$Median.Latitude)>23.5,'strict.tropical']<-0
  ##2: not strict tropical or strict temperate species
  #GBIFdata$strict.tropical[is.na(GBIFdata$strict.tropical)]<-2
  #GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
  #BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
  ##subset table to BAMM species
  #GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%BAMM.object$tip.label,]
  ##remove duplicates in GBIFdata (6 species with same info but two rows - one with present in garden, one with absent in garden)#
  #GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
  #GBIFdata.BAMM<-unique(GBIFdata.BAMM)
  #GBIFdata.BAMM$latitudinal.band<-NA
  #GBIFdata.BAMM$latitudinal.band<-apply(GBIFdata.BAMM,1,function(x) round_any(as.numeric(x['Median.Latitude']),10))
  #ratesdf<-getTipRates(BAMM.object)
  #ratesdf.netdiv<-getTipRates(BAMM.object,returnNetDiv = T)
  #rates.df<-as.data.frame(cbind(names(ratesdf$lambda.avg),ratesdf$lambda.avg,ratesdf$mu.avg,ratesdf.netdiv$netdiv.avg),stringsAsFactors = F)
  #row.names(rates.df)<-NULL
  #colnames(rates.df)<-c('species','lambda.avg','mu.avg','netdiv.avg')
  #GBIFdata.BAMM.rates<-merge(GBIFdata.BAMM,rates.df,by.x='binomial',by.y='species')
  #GBIFdata.BAMM.rates$lambda.avg<-as.numeric(GBIFdata.BAMM.rates$lambda.avg)
  #GBIFdata.BAMM.rates$mu.avg<-as.numeric(GBIFdata.BAMM.rates$mu.avg)
  #GBIFdata.BAMM.rates$netdiv.avg<-as.numeric(GBIFdata.BAMM.rates$netdiv.avg)
  
  
  GBIFdata.BAMM.rates$rank.lambda<-0
  GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==0,]$rank.lambda<-rank(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==0,]$lambda.avg)/length(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==0,]$lambda.avg)
  GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==1,]$rank.lambda<-rank(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==1,]$lambda.avg)/length(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==1,]$lambda.avg)
  #pdf('./plots/Fig1_quantilerate_vs_lambdarate_full.pdf')
  #par(mfrow=c(2,2))
  plot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==0,]$lambda.avg),GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==0,]$rank.lambda,pch=16,col=blue,main='lambda rates',xlab='log.lambda.rate',ylab='quantile of lambda distribution',cex=.5,xaxt='n',yaxt='n')
  points(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==1,]$lambda.avg),GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==1,]$rank.lambda,pch=16,col=red,cex=.5)
  axis(1,log(c(0.05,0.25,1,5)),labels=c(0.05,0.25,1,5))
  axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1),las=2)
  legend ('topleft',c('temperate','tropical'),col = c(blue,red), lty=1,lwd=2, bty="n", cex=.6)
  #dev.off()
  
}


plot_Fig1_quantilerate_vs_DRrate<-function(GBIFdata.BAMM.rates){
  
  GBIFdata.BAMM.rates$rank.DR<-0
  GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==0,]$rank.DR<-rank(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==0,]$DR)/length(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==0,]$DR)
  GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==1,]$rank.DR<-rank(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==1,]$DR)/length(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==1,]$DR)
  plot(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==0,]$DR),GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==0,]$rank.DR,pch=16,col=blue,main='DR',xlab='log.DR',ylab='quantile of DR distribution',cex=.5,xaxt='n',yaxt='n')
  points(log(GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==1,]$DR),GBIFdata.BAMM.rates[GBIFdata.BAMM.rates$tropical==1,]$rank.DR,pch=16,col=red,cex=.5)
  axis(1,log(c(0.1,0.5,2,5)),labels=c(0.1,0.5,2,5))
  axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1),las=2)
  legend ('topleft',c('temperate','tropical'),col = c(blue,red), lty=1,lwd=2, bty="n", cex=.6)
  #dev.off()
  
}


plot_Fig1_strapp_samples_from_file<-function(correlationfile,parameter,tablefile){
  spearman.regressions.df<-read.table(correlationfile,header=T,sep='\t')
  table<-read.table(tablefile,header=T,stringsAsFactors = F,sep='\t')
  table<-table[,c('binomial','Median.Latitude')]
  colnames(table)<-c('species','Median.Latitude')
  #pdf(pdffile,paper='a4')
  #min.y<-min(c(table$mean.lambda.rates,table$mean.mu.rates))
  #max.y<-max(c(table$mean.lambda.rates,table$mean.mu.rates))
  #min.x<-min(c(table$log.Seed.Weight))
  #max.x<-max(c(table$log.Seed.Weight))
  par(pty='s')
  min.x<-min(abs(table$Median.Latitude))
  max.x<-max(abs(table$Median.Latitude))
  #min.y<-min(unlist(BAMM.object$tipLambda))
  #max.y<-max(unlist(BAMM.object$tipLambda))
  min.y<-0.1
  max.y<-5
  if(parameter=='lambda'){
    plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='abs.median.latitude',ylab='speciation rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
    ticks.x<-seq(from=0,to=90,by=10)
    #ticks.x<-log10(seq(min.x, 80, by = 10))
    axis(1, at = ticks.x, labels=seq(from=0,to=90,by=10), las=1)
    #ticks.y<-c(log10(c(0.05,0.25,1,5)))
    #axis(2, at = ticks.y, labels=10^(ticks.y), las=1)
    axis(2, at = log10(c(0.25,1,2,5)), labels=c(0.25,1,2,5), las=1)
    #col.lines<-rgb(150,192,231,maxColorValue = 255)
    col.lines<-'grey80'
  }
  if(parameter=='mu'){
    plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='seed mass (g)',ylab='extinction rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
    ticks.x<-seq(-6, 4, by = 2)
    axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
    ticks.y<-seq(-3, 1, by = 1)
    axis(2, at = ticks.y, labels=10^(ticks.y), las=1)      
    #col.lines<-rgb(150,192,231,maxColorValue = 255)
    col.lines<-'grey80'
  }
  if(parameter=='div'){
    plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='seed mass (g)',ylab='diversification rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
    ticks.x<-seq(-6, 4, by = 2)
    axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
    ticks.y<-seq(-3, 1, by = 1)
    axis(2, at = ticks.y, labels=10^(ticks.y), las=1)      
    #col.lines<-rgb(150,192,231,maxColorValue = 255)
    col.lines<-'grey80'
  }
  
  apply(spearman.regressions.df,1,function(x) abline(x[1],x[2],col=adjustcolor('grey60',alpha=0.3)))
  median.slope<-which.min(abs(median(spearman.regressions.df$X2)-spearman.regressions.df$X2))
  abline(spearman.regressions.df$X1[median.slope],spearman.regressions.df$X2[median.slope])  
}


plot_Fig1_strapp_density<-function(){
  BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata_GBOTB.RDS')
  GBIFdata.BAMM<-read.table('./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',sep='\t',header=T,stringsAsFactors = F)
  trait.table =GBIFdata.BAMM
  names(median.lat)<-trait.table$binomial
  median.latitude<-abs(GBIFdata.BAMM.rates$Median.Latitude)
  names(median.latitude)<-GBIFdata.BAMM$binomial
  cat('running STRAPP','\n')
  strapp.lambda<-traitDependentBAMM(BAMM.object, median.latitude, 1000, rate = "speciation",return.full = T, method = "spearman", logrates = TRUE,two.tailed = TRUE)

  print('Plotting strapp histograms')
  pvalue<-strapp.lambda$p.value
  estimate<-strapp.lambda$estimate
  strapp.diff<-abs(strapp.lambda$obs.corr)-abs(strapp.lambda$null)
  title<-paste('rho = ',estimate,'(',pvalue,')',sep='')
  d<-density(strapp.diff)
  plot(d,xlim=c(-0.2,0.2),xaxs='i',yaxs='i',main=title,yaxt='n',ylab='',bty='n',xlab='rho[observed-null]')
  abline(v=0,lty=2)
  
}

plot_Fig2_allsubsampled_quantilerate_vs_lambda<-function(){
  #read all subsampled input tables
  BAMM.subsampled.files.list<-list.files('./results/strapp_unbiased/input_tables/',pattern='_table.txt')
  BAMM.subsampled.list<-lapply(BAMM.subsampled.files.list,function(x)read.table(paste('./results/strapp_unbiased/input_tables/',x,sep=''),stringsAsFactors = F,sep='\t',header=T))
  
  
  BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
  ratesdf<-getTipRates(BAMM.object)
  rates.df<-as.data.frame(cbind(names(ratesdf$lambda.avg),ratesdf$lambda.avg),stringsAsFactors = F)
  row.names(rates.df)<-NULL
  colnames(rates.df)<-c('species','lambda.avg')
  
  BAMM.subsampled.list.rates<-lapply(BAMM.subsampled.list,function(x)merge(x,rates.df,by.x='binomial',by.y='species',all.x=T))
  BAMM.subsampled.list.rates<-lapply(BAMM.subsampled.list.rates,function(x){x$lambda.avg<-as.numeric(x$lambda.avg);x$rank.lambda<-0;x[x$tropical==0,]$rank.lambda<-rank(x[x$tropical==0,]$lambda.avg)/length(x[x$tropical==0,]$lambda.avg);x[x$tropical==1,]$rank.lambda<-rank(x[x$tropical==1,]$lambda.avg)/length(x[x$tropical==1,]$lambda.avg);x$rank.lambda<-as.numeric(x$rank.lambda);return(x)})
  
  par(mfrow=c(2,2))
  par(pty='s')
  
  plot(log(as.numeric(BAMM.subsampled.list.rates[[1]][BAMM.subsampled.list.rates[[1]]$tropical==0,]$lambda.avg)),as.numeric(BAMM.subsampled.list.rates[[1]][BAMM.subsampled.list.rates[[1]]$tropical==0,]$rank.lambda),pch=16,col=blue,main='lambda rates',xlab='log.lambda.rate',ylab='quantile of lambda distribution',cex=.5,type='n',xaxt='n',yaxt='n',xlim=log(c(0.1,6)))
  lapply(BAMM.subsampled.list.rates,function(x)lines(log10(x[x$tropical==0,]$lambda.avg[order(log10(x[x$tropical==0,]$lambda.avg))]),x[x$tropical==0,]$rank.lambda[order(x[x$tropical==0,]$lambda.avg)],col=blue))
  lapply(BAMM.subsampled.list.rates,function(x)lines(log10(x[x$tropical==1,]$lambda.avg[order(log10(x[x$tropical==1,]$lambda.avg))]),x[x$tropical==1,]$rank.lambda[order(x[x$tropical==1,]$lambda.avg)],col=red))
  axis(1,log(c(0.05,0.25,1,5)),labels=c(0.05,0.25,1,5))
  axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1),las=2)
  legend ('topleft',c('temperate','tropical'),col = c(blue,red), lty=1,lwd=2, bty="n", cex=.6)
  
}

plot_Fig2_allsubsampled_quantilerate_vs_lambda_GBOTB<-function(){
  #read all subsampled input tables
  BAMM.subsampled.files.list<-list.files('./results/subsampled_unbiased_GBOTB/input_tables/',pattern='_table.txt')
  BAMM.subsampled.list<-lapply(BAMM.subsampled.files.list,function(x)read.table(paste('./results/subsampled_unbiased_GBOTB/input_tables/',x,sep=''),stringsAsFactors = F,sep='\t',header=T))
  
  GBIFdata.BAMM.rates<-read.table(file='./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',sep='\t',header=T,stringsAsFactors = F)
  GBIFdata.BAMM.rates<-GBIFdata.BAMM.rates[,c('binomial','lambda.avg')]  
  
  BAMM.subsampled.list.rates<-lapply(BAMM.subsampled.list,function(x)merge(x,GBIFdata.BAMM.rates,by.x='binomial',by.y='binomial',all.x=T))
  BAMM.subsampled.list.rates<-lapply(BAMM.subsampled.list.rates,function(x){x$lambda.avg<-as.numeric(x$lambda.avg);x$rank.lambda<-0;x[x$tropical==0,]$rank.lambda<-rank(x[x$tropical==0,]$lambda.avg)/length(x[x$tropical==0,]$lambda.avg);x[x$tropical==1,]$rank.lambda<-rank(x[x$tropical==1,]$lambda.avg)/length(x[x$tropical==1,]$lambda.avg);x$rank.lambda<-as.numeric(x$rank.lambda);return(x)})
  
  par(mfrow=c(2,2))
  par(pty='s')
  
  plot(log(as.numeric(BAMM.subsampled.list.rates[[1]][BAMM.subsampled.list.rates[[1]]$tropical==0,]$lambda.avg)),as.numeric(BAMM.subsampled.list.rates[[1]][BAMM.subsampled.list.rates[[1]]$tropical==0,]$rank.lambda),pch=16,col=blue,main='lambda rates',xlab='log.lambda.rate',ylab='quantile of lambda distribution',cex=.5,type='n',xaxt='n',yaxt='n',xlim=log(c(0.1,6)))
  lapply(BAMM.subsampled.list.rates,function(x)lines(log(x[x$tropical==0,]$lambda.avg[order(log(x[x$tropical==0,]$lambda.avg))]),x[x$tropical==0,]$rank.lambda[order(x[x$tropical==0,]$lambda.avg)],col=blue))
  lapply(BAMM.subsampled.list.rates,function(x)lines(log(x[x$tropical==1,]$lambda.avg[order(log(x[x$tropical==1,]$lambda.avg))]),x[x$tropical==1,]$rank.lambda[order(x[x$tropical==1,]$lambda.avg)],col=red))
  axis(1,log(c(0.05,0.25,1,5)),labels=c(0.05,0.25,1,5))
  axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1),las=2)
  legend ('topleft',c('temperate','tropical'),col = c(blue,red), lty=1,lwd=2, bty="n", cex=.6)
  
}
plot_FigS5_allsubsampled_extremebias_quantilerate_vs_lambda_GBOTB<-function(){
  #read all subsampled input tables
  BAMM.subsampled.files.list<-list.files('./results/subsampled_extremebias_GBOTB//input_tables/',pattern='_table.txt')
  BAMM.subsampled.list<-lapply(BAMM.subsampled.files.list,function(x)read.table(paste('./results/subsampled_extremebias_GBOTB/input_tables/',x,sep=''),stringsAsFactors = F,sep='\t',header=T))
  
  GBIFdata.BAMM.rates<-read.table(file='./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',sep='\t',header=T,stringsAsFactors = F)
  GBIFdata.BAMM.rates<-GBIFdata.BAMM.rates[,c('binomial','lambda.avg')]  
  
  BAMM.subsampled.list.rates<-lapply(BAMM.subsampled.list,function(x)merge(x,GBIFdata.BAMM.rates,by.x='binomial',by.y='binomial',all.x=T))
  BAMM.subsampled.list.rates<-lapply(BAMM.subsampled.list.rates,function(x){x$lambda.avg<-as.numeric(x$lambda.avg);x$rank.lambda<-0;x[x$tropical==0,]$rank.lambda<-rank(x[x$tropical==0,]$lambda.avg)/length(x[x$tropical==0,]$lambda.avg);x[x$tropical==1,]$rank.lambda<-rank(x[x$tropical==1,]$lambda.avg)/length(x[x$tropical==1,]$lambda.avg);x$rank.lambda<-as.numeric(x$rank.lambda);return(x)})
  
  par(mfrow=c(2,2))
  par(pty='s')
  
  plot(log(as.numeric(BAMM.subsampled.list.rates[[1]][BAMM.subsampled.list.rates[[1]]$tropical==0,]$lambda.avg)),as.numeric(BAMM.subsampled.list.rates[[1]][BAMM.subsampled.list.rates[[1]]$tropical==0,]$rank.lambda),pch=16,col=blue,main='lambda rates',xlab='log.lambda.rate',ylab='quantile of lambda distribution',cex=.5,type='n',xaxt='n',yaxt='n',xlim=log(c(0.1,6)))
  lapply(BAMM.subsampled.list.rates,function(x)lines(log(x[x$tropical==0,]$lambda.avg[order(log(x[x$tropical==0,]$lambda.avg))]),x[x$tropical==0,]$rank.lambda[order(x[x$tropical==0,]$lambda.avg)],col=blue))
  lapply(BAMM.subsampled.list.rates,function(x)lines(log(x[x$tropical==1,]$lambda.avg[order(log(x[x$tropical==1,]$lambda.avg))]),x[x$tropical==1,]$rank.lambda[order(x[x$tropical==1,]$lambda.avg)],col=red))
  axis(1,log(c(0.05,0.25,1,5)),labels=c(0.05,0.25,1,5))
  axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1),las=2)
  legend ('topleft',c('temperate','tropical'),col = c(blue,red), lty=1,lwd=2, bty="n", cex=.6)
  
}
plot_Fig2_latitudinalband_boxplots_subsampled_unbiased<-function(){
  latband.strapp.files<-list.files('./results/subsampled_unbiased/latitudinalband/',pattern='latitudinalband.+.txt')
  latband.strapp.tables<-lapply(latband.strapp.files,function(x)read.table(paste('./results/subsampled_unbiased/latitudinalband/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
  latband.strapp.speciation<-do.call('rbind.fill',lapply(latband.strapp.tables,function(x)x[1,]))
  #plot median of speciation rate across absolute lat bands
  latband.strapp.speciation<-latband.strapp.speciation[,-c(1,ncol(latband.strapp.speciation))]
  #####CHANGE PATH OF REALMS LAYER####
  realms.layer<-readOGR('/Users/javier/Documents/Work/HOTSPOTS/repositories/raw_data/GIS/realms/','newRealms')
  realms.layer<-gUnaryUnion(realms.layer)
  pdf('./plots/Fig2_map_latitudinalbands_log_medians.pdf',paper='a4r')
  plot(realms.layer,lwd=0.5,border='grey80')
  #for (i in seq(from=85,to=-55,by=-10)){
  #  abline(h=i,col='grey70')
  #}
  abline(v=0)
  abline(h=23.5,lty=2,col='grey80')
  abline(h=-23.5,lty=2,col='grey80')
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-55,-55,-23.5,-23.5,-55),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-23.5,-23.5,23.5,23.5,-23.5),col=adjustcolor(red,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(85,85,23.5,23.5,85),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  u<-par('usr')
  #plot(c(0,0),type='n',xlim=c(u[1]/30,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #plot(c(0,0),type='n',xlim=c(u[1]/50,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  plot(c(0,0),type='n',xlim=c(-log(6),log(6)),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  axis(1,c(log(c(0.5,1,2,5))),labels=c(0.5,1,2,5))
  axis(4,seq(from=-50/10,to=80/10,by=10/10),labels=seq(from=-50,to=80,by=10),las=2,cex.axis=.5)
  #axis(1,c(log(0),log(3),log(6),log(9),log(12),log(15)),labels=c(0,3,6,9,12,15))
  abline(v=0)
  #par(new=T)
  boxplot(log(latband.strapp.speciation[,7]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,6]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=0/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,5]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,4]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,3]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,2]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,1]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=-50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,8]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,9]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,10]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,11]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,12]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=60/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,13]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=70/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,14]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=80/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  dev.off()
}

plot_Fig2_latitudinalband_boxplots_subsampled_unbiased_xaxisFig1<-function(){
  latband.strapp.files<-list.files('./results/subsampled_unbiased/latitudinalband/',pattern='latitudinalband.+.txt')
  latband.strapp.tables<-lapply(latband.strapp.files,function(x)read.table(paste('./results/subsampled_unbiased/latitudinalband/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
  latband.strapp.speciation<-do.call('rbind.fill',lapply(latband.strapp.tables,function(x)x[1,]))
  #plot median of speciation rate across absolute lat bands
  latband.strapp.speciation<-latband.strapp.speciation[,-c(1,ncol(latband.strapp.speciation))]
  #####CHANGE PATH OF REALMS LAYER####
  realms.layer<-readOGR('/Users/javier/Documents/Work/HOTSPOTS/repositories/raw_data/GIS/realms/','newRealms')
  realms.layer<-gUnaryUnion(realms.layer)
  plot(realms.layer,lwd=0.5,border='grey80')
  #for (i in seq(from=85,to=-55,by=-10)){
  #  abline(h=i,col='grey70')
  #}
  abline(v=0)
  abline(h=23.5,lty=2,col='grey80')
  abline(h=-23.5,lty=2,col='grey80')
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-55,-55,-23.5,-23.5,-55),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-23.5,-23.5,23.5,23.5,-23.5),col=adjustcolor(red,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(85,85,23.5,23.5,85),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  u<-par('usr')
  #plot(c(0,0),type='n',xlim=c(u[1]/30,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #plot(c(0,0),type='n',xlim=c(u[1]/50,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  plot(c(0,0),type='n',xlim=c(-3,3),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  axis(1,c(log(c(0.25,1,2,5,10))),labels=c(0.25,1,2,5,10))
  axis(4,seq(from=-50/10,to=80/10,by=10/10),labels=seq(from=-50,to=80,by=10),las=2,cex.axis=.5)
  #axis(1,c(log(0),log(3),log(6),log(9),log(12),log(15)),labels=c(0,3,6,9,12,15))
  abline(v=0)
  #par(new=T)
  boxplot(log(latband.strapp.speciation[,7]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,6]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=0/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,5]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,4]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,3]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,2]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,1]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=-50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,8]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,9]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,10]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,11]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,12]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=60/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,13]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=70/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,14]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=80/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
}

plot_Fig2_map_latitudinalband_boxplots_subsampled_unbiased<-function(){
  latband.strapp.files<-list.files('./results/subsampled_unbiased/latitudinalband/',pattern='latitudinalband.+.txt')
  latband.strapp.tables<-lapply(latband.strapp.files,function(x)read.table(paste('./results/subsampled_unbiased/latitudinalband/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
  latband.strapp.speciation<-do.call('rbind.fill',lapply(latband.strapp.tables,function(x)x[1,]))
  #plot median of speciation rate across absolute lat bands
  latband.strapp.speciation<-latband.strapp.speciation[,-c(1,ncol(latband.strapp.speciation))]
  bands<-seq(-50,80,by=10)
  #calculate IQRs for ylim
  min.ylim<-sort(apply(latband.strapp.speciation[,c(1:14)],2,function(x) quantile(x,c(0.25,0.75),na.rm = T)))[1]
  max.ylim<-sort(apply(latband.strapp.speciation[,c(1:14)],2,function(x) quantile(x,c(0.25,0.75),na.rm = T)),decreasing = T)[1]
  pdf('./plots/Fig2_map_latitudinalband_boxplots_subsampled_unbiased.pdf')
  par(pty='s')
  plot(c(1,1),xaxt='n',yaxt='n',xlab='latitude (degrees)',ylab='median.lambda',main='unbiased subsamples',type='n',xlim=c(0,15),ylim=log(c(0.5,6)))
  axis(1,at=c(1:14),labels=seq(-50,80,by=10),cex.axis=.8)
  axis(2,log(c(0.5,1,2,5)),labels=c(0.5,1,2,5),las=2)
  for(i in 1:length(bands)){
    if(abs(bands[i])<30){
      col<-red
    }else{
      col<-blue
    }
    if(bands[i]==80){
      boxplot(log(latband.strapp.speciation[,i]),at=i,outline=F,whisklty=0,staplelty=0,frame=F,col=adjustcolor(col,alpha.f = 0.25),medcol=col,border=adjustcolor(col,alpha.f=0.5),notch=F,add=T,yaxt='n')
    }else{
      boxplot(log(latband.strapp.speciation[,i]),at=i,outline=F,whisklty=0,staplelty=0,frame=F,col=adjustcolor(col,alpha.f = 0.25),medcol=col,border=adjustcolor(col,alpha.f=0.5),notch=T,add=T,yaxt='n')
    }
    
  }
  dev.off()
  
}

plot_Fig2_strapp_densities<-function(path.strappobjects){
  strapp.abs.latband.object.files<-list.files(path=path.strappobjects,pattern='abs.medianlatitude.strapp.*.RDS')
  strapp.abs.latband.object<-lapply(strapp.abs.latband.object.files,function(x)readRDS(paste(path.strappobjects,x,sep='')))
  strapp.abs.latband.object.diff<-lapply(strapp.abs.latband.object,function(x)abs(x$obs.corr)-abs(x$null))
  
  strapp.abs.latband.object.diff.densities<-lapply(strapp.abs.latband.object.diff,function(x)density(x,from=-0.2,to=0.2,n=1000))
  strapp.abs.latband.object.diff.densities.y<-lapply(strapp.abs.latband.object.diff.densities,function(x)x$y)
  strapp.abs.latband.object.diff.densities.all<-lapply(c(1:1000),function(x) sapply(strapp.abs.latband.object.diff.densities.y,'[',x))
  strapp.abs.latband.object.diff.densities.all.median<-lapply(strapp.abs.latband.object.diff.densities.all,function(x)median(x))
  strapp.abs.latband.object.diff.densities.all.CIup<-lapply(strapp.abs.latband.object.diff.densities.all,function(x)quantile(x,c(0.975)))
  strapp.abs.latband.object.diff.densities.all.CIdw<-lapply(strapp.abs.latband.object.diff.densities.all,function(x)quantile(x,c(0.025)))
  
  plot(c(0.1,0.1),type='n',xlim=c(-0.2,0.2),ylim=c(0,25),xaxs='i',yaxs='i',xlab='rho',ylab='',yaxt='n',bty='n')
  polygon(x=c(strapp.abs.latband.object.diff.densities[[1]]$x,rev(strapp.abs.latband.object.diff.densities[[1]]$x)),y=c(strapp.abs.latband.object.diff.densities.all.CIup,rev(strapp.abs.latband.object.diff.densities.all.CIdw)),col='grey80',border=NA)
  abline(v=0,lty=2)
  lines(x=strapp.abs.latband.object.diff.densities[[1]]$x,y=strapp.abs.latband.object.diff.densities.all.median,lwd=1)
  
}

plot_Fig2_strapp_samples_from_files<-function(path.correlationfile,parameter,tablefile){
  correlationfiles<-list.files(path.correlationfile,pattern='strapp_lambda_absMedianLatitude+.*')
  correlation.tables<-lapply(correlationfiles,function(x)read.table(paste(path.correlationfile,x,sep=''),header=T,sep='\t',stringsAsFactors = F))
  correlation.tables.df<-do.call('rbind',correlation.tables)
  
  table<-read.table(tablefile,header=T,stringsAsFactors = F,sep='\t')
  table<-table[,c('binomial','Median.Latitude')]
  colnames(table)<-c('species','Median.Latitude')
  #pdf(pdffile,paper='a4')
  #min.y<-min(c(table$mean.lambda.rates,table$mean.mu.rates))
  #max.y<-max(c(table$mean.lambda.rates,table$mean.mu.rates))
  #min.x<-min(c(table$log.Seed.Weight))
  #max.x<-max(c(table$log.Seed.Weight))
  par(pty='s')
  min.x<-min(abs(table$Median.Latitude))
  max.x<-max(abs(table$Median.Latitude))
  #min.y<-min(unlist(BAMM.object$tipLambda))
  #max.y<-max(unlist(BAMM.object$tipLambda))
  min.y<-0.1
  max.y<-5
  if(parameter=='lambda'){
    plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='abs.median.latitude',ylab='speciation rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
    ticks.x<-seq(from=0,to=90,by=10)
    #ticks.x<-log10(seq(min.x, 80, by = 10))
    axis(1, at = ticks.x, labels=seq(from=0,to=90,by=10), las=1)
    #ticks.y<-c(log10(c(0.05,0.25,1,5)))
    #axis(2, at = ticks.y, labels=10^(ticks.y), las=1)
    axis(2, at = log10(c(0.25,1,2,5)), labels=c(0.25,1,2,5), las=1)
    #col.lines<-rgb(150,192,231,maxColorValue = 255)
    col.lines<-'grey80'
  }
  if(parameter=='mu'){
    plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='seed mass (g)',ylab='extinction rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
    ticks.x<-seq(-6, 4, by = 2)
    axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
    ticks.y<-seq(-3, 1, by = 1)
    axis(2, at = ticks.y, labels=10^(ticks.y), las=1)      
    #col.lines<-rgb(150,192,231,maxColorValue = 255)
    col.lines<-'grey80'
  }
  if(parameter=='div'){
    plot(c(min.x,max.x),c(log10(min.y),log10(max.y)),xlab='seed mass (g)',ylab='diversification rate (lineages myr-1)',xaxt='n',yaxt='n',type='n')
    ticks.x<-seq(-6, 4, by = 2)
    axis(1, at = ticks.x, labels=10^(ticks.x), las=1)
    ticks.y<-seq(-3, 1, by = 1)
    axis(2, at = ticks.y, labels=10^(ticks.y), las=1)      
    #col.lines<-rgb(150,192,231,maxColorValue = 255)
    col.lines<-'grey80'
  }
  u<-par('usr')
  
  
  
  
  
  lower.confint<-curve(quantile(correlation.tables.df[,1],c(0.025)) + x*quantile(correlation.tables.df[,2],c(0.025)), from=u[1], to=u[2],add=T,col=adjustcolor('grey60',alpha=0.3))
  upper.confint<-curve(quantile(correlation.tables.df[,1],c(0.975)) + x*quantile(correlation.tables.df[,2],c(0.975)), from=u[1], to=u[2],add=T,col=adjustcolor('grey60',alpha=0.3))
  polygon(c(lower.confint$x,rev(upper.confint$x)), c(lower.confint$y, rev(upper.confint$y)),col=adjustcolor('grey60',alpha=0.3),border='grey60')
  abline(quantile(correlation.tables.df[,1],c(0.5)),quantile(correlation.tables.df[,2],c(0.5)))
  
  #abline(quantile(correlation.tables.df[,1],c(0.025)),quantile(correlation.tables.df[,2],c(0.025)),lty=2)
  #abline(quantile(correlation.tables.df[,1],c(0.975)),quantile(correlation.tables.df[,2],c(0.975)),lty=2)
  
}


plot_Fig2_latitudinalband_boxplots_subsampled_unbiased_axisFig1<-function(){
  latband.strapp.files<-list.files('./results/subsampled_unbiased/latitudinalband/',pattern='latitudinalband.+.txt')
  latband.strapp.tables<-lapply(latband.strapp.files,function(x)read.table(paste('./results/subsampled_unbiased/latitudinalband/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
  latband.strapp.speciation<-do.call('rbind.fill',lapply(latband.strapp.tables,function(x)x[1,]))
  #plot median of speciation rate across absolute lat bands
  latband.strapp.speciation<-latband.strapp.speciation[,-c(1,ncol(latband.strapp.speciation))]
  #####CHANGE PATH OF REALMS LAYER####
  realms.layer<-readOGR('/Users/javier/Documents/Work/HOTSPOTS/repositories/raw_data/GIS/realms/','newRealms')
  realms.layer<-gUnaryUnion(realms.layer)
  pdf('./plots/Fig2_map_latitudinalbands_log_medians.pdf',paper='a4r')
  plot(realms.layer,lwd=0.5,border='grey80')
  #for (i in seq(from=85,to=-55,by=-10)){
  #  abline(h=i,col='grey70')
  #}
  abline(v=0)
  abline(h=23.5,lty=2,col='grey80')
  abline(h=-23.5,lty=2,col='grey80')
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-55,-55,-23.5,-23.5,-55),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-23.5,-23.5,23.5,23.5,-23.5),col=adjustcolor(red,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(85,85,23.5,23.5,85),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  u<-par('usr')
  #plot(c(0,0),type='n',xlim=c(u[1]/30,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #plot(c(0,0),type='n',xlim=c(u[1]/50,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  plot(c(0,0),type='n',xlim=c(-log(6),log(6)),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  axis(1,c(log(c(0.5,1,2,5))),labels=c(0.5,1,2,5))
  axis(4,seq(from=-50/10,to=80/10,by=10/10),labels=seq(from=-50,to=80,by=10),las=2,cex.axis=.5)
  #axis(1,c(log(0),log(3),log(6),log(9),log(12),log(15)),labels=c(0,3,6,9,12,15))
  abline(v=0)
  #par(new=T)
  boxplot(log(latband.strapp.speciation[,7]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,6]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=0/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,5]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,4]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,3]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,2]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,1]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=-50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,8]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,9]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,10]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,11]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,12]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=60/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,13]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=70/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,14]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=80/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  dev.off()
}



plot_Fig2_allsubsampled_extremebias_quantilerate_vs_lambda<-function(){
  #read all subsampled input tables
  BAMM.subsampled.files.list<-list.files('./results/subsampled_extremebias/input_tables/',pattern='_table.txt')
  BAMM.subsampled.list<-lapply(BAMM.subsampled.files.list,function(x)read.table(paste('./results/subsampled_extremebias/input_tables/',x,sep=''),stringsAsFactors = F,sep='\t',header=T))
  
  
  BAMM.object<-readRDS(file='./raw_data/all_clades_50million_1000samples_latitudedata.RDS')
  ratesdf<-getTipRates(BAMM.object)
  rates.df<-as.data.frame(cbind(names(ratesdf$lambda.avg),ratesdf$lambda.avg),stringsAsFactors = F)
  row.names(rates.df)<-NULL
  colnames(rates.df)<-c('species','lambda.avg')
  
  BAMM.subsampled.list.rates<-lapply(BAMM.subsampled.list,function(x)merge(x,rates.df,by.x='binomial',by.y='species',all.x=T))
  BAMM.subsampled.list.rates<-lapply(BAMM.subsampled.list.rates,function(x){x$lambda.avg<-as.numeric(x$lambda.avg);x$rank.lambda<-0;x[x$tropical==0,]$rank.lambda<-rank(x[x$tropical==0,]$lambda.avg)/length(x[x$tropical==0,]$lambda.avg);x[x$tropical==1,]$rank.lambda<-rank(x[x$tropical==1,]$lambda.avg)/length(x[x$tropical==1,]$lambda.avg);x$rank.lambda<-as.numeric(x$rank.lambda);return(x)})
  
  par(mfrow=c(2,2))
  par(pty='s')
  
  plot(log(as.numeric(BAMM.subsampled.list.rates[[1]][BAMM.subsampled.list.rates[[1]]$tropical==0,]$lambda.avg)),as.numeric(BAMM.subsampled.list.rates[[1]][BAMM.subsampled.list.rates[[1]]$tropical==0,]$rank.lambda),pch=16,col=blue,main='lambda rates',xlab='log.lambda.rate',ylab='quantile of lambda distribution',cex=.5,type='n',xaxt='n',yaxt='n',xlim=log(c(0.1,6)))
  lapply(BAMM.subsampled.list.rates,function(x)lines(log10(x[x$tropical==0,]$lambda.avg[order(log10(x[x$tropical==0,]$lambda.avg))]),x[x$tropical==0,]$rank.lambda[order(x[x$tropical==0,]$lambda.avg)],col=blue))
  lapply(BAMM.subsampled.list.rates,function(x)lines(log10(x[x$tropical==1,]$lambda.avg[order(log10(x[x$tropical==1,]$lambda.avg))]),x[x$tropical==1,]$rank.lambda[order(x[x$tropical==1,]$lambda.avg)],col=red))
  axis(1,log(c(0.05,0.25,1,5)),labels=c(0.05,0.25,1,5))
  axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1),las=2)
  legend ('topleft',c('temperate','tropical'),col = c(blue,red), lty=1,lwd=2, bty="n", cex=.6)
  
}

plot_Fig2_latitudinalband_boxplots_subsampled_extremebias<-function(){
  latband.strapp.files<-list.files('./results/subsampled_extremebias/latitudinalband/',pattern='latitudinalband.+.txt')
  latband.strapp.tables<-lapply(latband.strapp.files,function(x)read.table(paste('./results/subsampled_extremebias/latitudinalband/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
  latband.strapp.speciation<-do.call('rbind.fill',lapply(latband.strapp.tables,function(x)x[1,]))
  #plot median of speciation rate across absolute lat bands
  latband.strapp.speciation<-latband.strapp.speciation[,-c(1,ncol(latband.strapp.speciation))]
  #####CHANGE PATH OF REALMS LAYER####
  realms.layer<-readOGR('/Users/javier/Documents/Work/HOTSPOTS/repositories/raw_data/GIS/realms/','newRealms')
  realms.layer<-gUnaryUnion(realms.layer)
  plot(realms.layer,lwd=0.5,border='grey80')
  #for (i in seq(from=85,to=-55,by=-10)){
  #  abline(h=i,col='grey70')
  #}
  abline(v=0)
  abline(h=23.5,lty=2,col='grey80')
  abline(h=-23.5,lty=2,col='grey80')
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-55,-55,-23.5,-23.5,-55),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-23.5,-23.5,23.5,23.5,-23.5),col=adjustcolor(red,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(85,85,23.5,23.5,85),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  u<-par('usr')
  #plot(c(0,0),type='n',xlim=c(u[1]/30,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #plot(c(0,0),type='n',xlim=c(u[1]/50,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  plot(c(0,0),type='n',xlim=c(-log(6),log(6)),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  axis(1,c(log(c(0.5,1,2,5))),labels=c(0.5,1,2,5))
  axis(4,seq(from=-50/10,to=80/10,by=10/10),labels=seq(from=-50,to=80,by=10),las=2,cex.axis=.5)
  #axis(1,c(log(0),log(3),log(6),log(9),log(12),log(15)),labels=c(0,3,6,9,12,15))
  abline(v=0)
  #par(new=T)
  boxplot(log(latband.strapp.speciation[,7]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,6]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=0/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,5]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-10/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,4]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,3]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,2]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=-40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,1]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=-50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,8]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=20/10,add=T,col=adjustcolor(red,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,9]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=30/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,10]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=40/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,11]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=50/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,12]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=60/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,13]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=T,at=70/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
  boxplot(log(latband.strapp.speciation[,14]),horizontal=T,boxwex=1.2,outline=F,whisklty=0,staplelty=0,notch=F,at=80/10,add=T,col=adjustcolor(blue,alpha=.7),frame=F,xaxt='n')
}

plot_Fig2_latitudinalband_boxplots_subsampled_extremebias_summary<-function(){
  #get rates from BAMM
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
  GBIFdata.BAMM.rates<-GBIFdata.BAMM.rates[,c('binomial','lambda.avg')]
  latband.strapp.files<-list.files(path='./results/subsampled_extremebias/input_tables/',pattern='GBIFdata.BAMM.subsampled_extremebias_.+.txt')
  latband.strapp.tables<-lapply(latband.strapp.files,function(x)read.table(paste('./results/subsampled_extremebias/input_tables/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
  latband.strapp.tables<-lapply(latband.strapp.tables,function(x)merge(x,GBIFdata.BAMM.rates))
  #get quartile 25,50,75 for each table across all lat bands
  latband.strapp.tables.quartiles<-lapply(latband.strapp.tables,function(x)aggregate(lambda.avg~latitudinal.band,data=x,function(x)quantile(x,c(0.25,0.5,0.75))))
  #get averages of quartiles 25,50,75 across all tables for each latitudinal band
  latband.average.quartiles<-lapply(seq(-60,90,10),function(x)lapply(latband.strapp.tables.quartiles,function(y){if(nrow(y[y$latitudinal.band==x,])==0){NA}else{y[y$latitudinal.band==x,]}}))
  latband.average.quartiles.all.tables<-lapply(seq(-60,90,10),function(x)lapply(latband.strapp.tables.quartiles,function(y){if(nrow(y[y$latitudinal.band==x,])==0){NA}else{y[y$latitudinal.band==x,]}}))
  latband.average.quartiles.all.tables.clean<-lapply(latband.average.quartiles.all.tables,function(x) x[!unlist(lapply(x,function(y){all(is.na(y))}))])
  latband.average.quartiles.all.tables.clean.df<-lapply(latband.average.quartiles.all.tables.clean,function(x)do.call('rbind',x))
  latband.average.quartiles.all.tables.clean.IQR.median<-lapply(latband.average.quartiles.all.tables.clean.df,function(x)unlist(lapply(c(1,2,3),function(y)mean(x[,2][,y]))))
  #####CHANGE PATH OF REALMS LAYER####
  realms.layer<-readOGR('/Users/javier/Documents/Work/HOTSPOTS/repositories/raw_data/GIS/realms/','newRealms')
  realms.layer<-gUnaryUnion(realms.layer)
  plot(realms.layer,lwd=0.5,border='grey80')
  #for (i in seq(from=85,to=-55,by=-10)){
  #  abline(h=i,col='grey70')
  #}
  abline(v=0)
  abline(h=23.5,lty=2,col='grey80')
  abline(h=-23.5,lty=2,col='grey80')
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-55,-55,-23.5,-23.5,-55),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-23.5,-23.5,23.5,23.5,-23.5),col=adjustcolor(red,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(85,85,23.5,23.5,85),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  u<-par('usr')
  #plot(c(0,0),type='n',xlim=c(u[1]/30,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #plot(c(0,0),type='n',xlim=c(u[1]/50,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  #plot(c(0,0),type='n',xlim=c(-log(6),log(6)),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  plot(c(0,0),type='n',xlim=c(-3,3),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  axis(1,c(log(c(0.25,1,2,5,10))),labels=c(0.25,1,2,5,10))
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  #axis(1,c(log(c(0.5,1,2,5))),labels=c(0.5,1,2,5))
  axis(4,seq(from=-50/10,to=80/10,by=10/10),labels=seq(from=-50,to=80,by=10),las=2,cex.axis=.5)
  #axis(1,c(log(0),log(3),log(6),log(9),log(12),log(15)),labels=c(0,3,6,9,12,15))
  abline(v=0)
  #par(new=T)
  
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][3]),ytop = -50/10,ybottom=(-50/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][2])),y=c(-50/10,(-50/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][3]),ytop = -40/10,ybottom=(-40/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][2])),y=c(-40/10,(-40/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][3]),ytop = -30/10,ybottom=(-30/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][2])),y=c(-30/10,(-30/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][3]),ytop = -20/10,ybottom=(-20/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][2])),y=c(-20/10,(-20/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][3]),ytop = -10/10,ybottom=(-10/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][2])),y=c(-10/10,(-10/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][3]),ytop = 0/10,ybottom=(0/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][2])),y=c(0/10,(0/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][3]),ytop = 10/10,ybottom=(10/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][2])),y=c(10/10,(10/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][3]),ytop = 20/10,ybottom=(20/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][2])),y=c(20/10,(20/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][3]),ytop = 30/10,ybottom=(30/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][2])),y=c(30/10,(30/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][3]),ytop = 40/10,ybottom=(40/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][2])),y=c(40/10,(40/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][3]),ytop = 50/10,ybottom=(50/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][2])),y=c(50/10,(50/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][3]),ytop = 60/10,ybottom=(60/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][2])),y=c(60/10,(60/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][3]),ytop = 70/10,ybottom=(70/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][2])),y=c(70/10,(70/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][3]),ytop = 80/10,ybottom=(80/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][2])),y=c(80/10,(80/10)+(1.2/3)),col='black',lwd=1)
  
  
  
  
}

plot_Fig2_latitudinalband_boxplots_subsampled_unbiased_summary<-function(){
  #get rates from BAMM
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
  GBIFdata.BAMM.rates<-GBIFdata.BAMM.rates[,c('binomial','lambda.avg')]
  latband.strapp.files<-list.files(path='./results/subsampled_unbiased/input_tables/',pattern='GBIFdata.BAMM.subsampled_.+.txt')
  latband.strapp.tables<-lapply(latband.strapp.files,function(x)read.table(paste('./results/subsampled_unbiased/input_tables/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
  latband.strapp.tables<-lapply(latband.strapp.tables,function(x)merge(x,GBIFdata.BAMM.rates))
  #get quartile 25,50,75 for each table across all lat bands
  latband.strapp.tables.quartiles<-lapply(latband.strapp.tables,function(x)aggregate(lambda.avg~latitudinal.band,data=x,function(x)quantile(x,c(0.25,0.5,0.75))))
  #get averages of quartiles 25,50,75 across all tables for each latitudinal band
  latband.average.quartiles<-lapply(seq(-60,90,10),function(x)lapply(latband.strapp.tables.quartiles,function(y){if(nrow(y[y$latitudinal.band==x,])==0){NA}else{y[y$latitudinal.band==x,]}}))
  latband.average.quartiles.all.tables<-lapply(seq(-60,90,10),function(x)lapply(latband.strapp.tables.quartiles,function(y){if(nrow(y[y$latitudinal.band==x,])==0){NA}else{y[y$latitudinal.band==x,]}}))
  latband.average.quartiles.all.tables.clean<-lapply(latband.average.quartiles.all.tables,function(x) x[!unlist(lapply(x,function(y){all(is.na(y))}))])
  latband.average.quartiles.all.tables.clean.df<-lapply(latband.average.quartiles.all.tables.clean,function(x)do.call('rbind',x))
  latband.average.quartiles.all.tables.clean.IQR.median<-lapply(latband.average.quartiles.all.tables.clean.df,function(x)unlist(lapply(c(1,2,3),function(y)mean(x[,2][,y]))))
  #####CHANGE PATH OF REALMS LAYER####
  realms.layer<-readOGR('/Users/javier/Documents/Work/HOTSPOTS/repositories/raw_data/GIS/realms/','newRealms')
  realms.layer<-gUnaryUnion(realms.layer)
  plot(realms.layer,lwd=0.5,border='grey80')
  #for (i in seq(from=85,to=-55,by=-10)){
  #  abline(h=i,col='grey70')
  #}
  abline(v=0)
  abline(h=23.5,lty=2,col='grey80')
  abline(h=-23.5,lty=2,col='grey80')
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-55,-55,-23.5,-23.5,-55),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-23.5,-23.5,23.5,23.5,-23.5),col=adjustcolor(red,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(85,85,23.5,23.5,85),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  u<-par('usr')
  #plot(c(0,0),type='n',xlim=c(u[1]/30,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #plot(c(0,0),type='n',xlim=c(u[1]/50,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  #plot(c(0,0),type='n',xlim=c(-log(6),log(6)),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  plot(c(0,0),type='n',xlim=c(-3,3),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  axis(1,c(log(c(0.25,1,2,5,10))),labels=c(0.25,1,2,5,10))
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  #axis(1,c(log(c(0.5,1,2,5))),labels=c(0.5,1,2,5))
  axis(4,seq(from=-50/10,to=80/10,by=10/10),labels=seq(from=-50,to=80,by=10),las=2,cex.axis=.5)
  #axis(1,c(log(0),log(3),log(6),log(9),log(12),log(15)),labels=c(0,3,6,9,12,15))
  abline(v=0)
  #par(new=T)
  
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][3]),ytop = -50/10,ybottom=(-50/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][2])),y=c(-50/10,(-50/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][3]),ytop = -40/10,ybottom=(-40/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][2])),y=c(-40/10,(-40/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][3]),ytop = -30/10,ybottom=(-30/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][2])),y=c(-30/10,(-30/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][3]),ytop = -20/10,ybottom=(-20/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][2])),y=c(-20/10,(-20/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][3]),ytop = -10/10,ybottom=(-10/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][2])),y=c(-10/10,(-10/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][3]),ytop = 0/10,ybottom=(0/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][2])),y=c(0/10,(0/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][3]),ytop = 10/10,ybottom=(10/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][2])),y=c(10/10,(10/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][3]),ytop = 20/10,ybottom=(20/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][2])),y=c(20/10,(20/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][3]),ytop = 30/10,ybottom=(30/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][2])),y=c(30/10,(30/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][3]),ytop = 40/10,ybottom=(40/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][2])),y=c(40/10,(40/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][3]),ytop = 50/10,ybottom=(50/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][2])),y=c(50/10,(50/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][3]),ytop = 60/10,ybottom=(60/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][2])),y=c(60/10,(60/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][3]),ytop = 70/10,ybottom=(70/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][2])),y=c(70/10,(70/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][3]),ytop = 80/10,ybottom=(80/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][2])),y=c(80/10,(80/10)+(1.2/3)),col='black',lwd=1)
  
  
  
  
}


plot_Fig2_latitudinalband_boxplots_subsampled_unbiased_summary_GBOTB<-function(){
  GBIFdata.BAMM.rates<-read.table('./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',header=T,sep='\t',stringsAsFactors = F)
  GBIFdata.BAMM.rates<-GBIFdata.BAMM.rates[,c('binomial','lambda.avg')]
  latband.strapp.files<-list.files(path='./results/subsampled_unbiased_GBOTB/input_tables/',pattern='GBIFdata.BAMM.GBOTB.subsampled_.+.txt')
  latband.strapp.tables<-lapply(latband.strapp.files,function(x)read.table(paste('./results/subsampled_unbiased_GBOTB/input_tables/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
  latband.strapp.tables<-lapply(latband.strapp.tables,function(x)merge(x,GBIFdata.BAMM.rates))
  #get quartile 25,50,75 for each table across all lat bands
  latband.strapp.tables.quartiles<-lapply(latband.strapp.tables,function(x)aggregate(lambda.avg~latitudinal.band,data=x,function(x)quantile(x,c(0.25,0.5,0.75))))
  #get averages of quartiles 25,50,75 across all tables for each latitudinal band
  latband.average.quartiles<-lapply(seq(-60,90,10),function(x)lapply(latband.strapp.tables.quartiles,function(y){if(nrow(y[y$latitudinal.band==x,])==0){NA}else{y[y$latitudinal.band==x,]}}))
  latband.average.quartiles.all.tables<-lapply(seq(-60,90,10),function(x)lapply(latband.strapp.tables.quartiles,function(y){if(nrow(y[y$latitudinal.band==x,])==0){NA}else{y[y$latitudinal.band==x,]}}))
  latband.average.quartiles.all.tables.clean<-lapply(latband.average.quartiles.all.tables,function(x) x[!unlist(lapply(x,function(y){all(is.na(y))}))])
  latband.average.quartiles.all.tables.clean.df<-lapply(latband.average.quartiles.all.tables.clean,function(x)do.call('rbind',x))
  latband.average.quartiles.all.tables.clean.IQR.median<-lapply(latband.average.quartiles.all.tables.clean.df,function(x)unlist(lapply(c(1,2,3),function(y)mean(x[,2][,y]))))
  #####CHANGE PATH OF REALMS LAYER####
  realms.layer<-readOGR('./raw_data/realms/','newRealms')
  realms.layer<-gUnaryUnion(realms.layer)
  plot(realms.layer,lwd=0.5,border='grey80')
  #for (i in seq(from=85,to=-55,by=-10)){
  #  abline(h=i,col='grey70')
  #}
  abline(v=0)
  abline(h=23.5,lty=2,col='grey80')
  abline(h=-23.5,lty=2,col='grey80')
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-55,-55,-23.5,-23.5,-55),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-23.5,-23.5,23.5,23.5,-23.5),col=adjustcolor(red,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(85,85,23.5,23.5,85),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  u<-par('usr')
  #plot(c(0,0),type='n',xlim=c(u[1]/30,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #plot(c(0,0),type='n',xlim=c(u[1]/50,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  #plot(c(0,0),type='n',xlim=c(-log(6),log(6)),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  plot(c(0,0),type='n',xlim=c(-2,1.5),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  axis(1,c(log(c(0.25,1,2,5))),labels=c(0.25,1,2,5))
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  #axis(1,c(log(c(0.5,1,2,5))),labels=c(0.5,1,2,5))
  axis(4,seq(from=-50/10,to=80/10,by=10/10),labels=seq(from=-50,to=80,by=10),las=2,cex.axis=.5)
  #axis(1,c(log(0),log(3),log(6),log(9),log(12),log(15)),labels=c(0,3,6,9,12,15))
  abline(v=0)
  #par(new=T)
  
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][3]),ytop = -50/10,ybottom=(-50/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][2])),y=c(-50/10,(-50/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][3]),ytop = -40/10,ybottom=(-40/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][2])),y=c(-40/10,(-40/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][3]),ytop = -30/10,ybottom=(-30/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][2])),y=c(-30/10,(-30/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][3]),ytop = -20/10,ybottom=(-20/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][2])),y=c(-20/10,(-20/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][3]),ytop = -10/10,ybottom=(-10/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][2])),y=c(-10/10,(-10/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][3]),ytop = 0/10,ybottom=(0/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][2])),y=c(0/10,(0/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][3]),ytop = 10/10,ybottom=(10/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][2])),y=c(10/10,(10/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][3]),ytop = 20/10,ybottom=(20/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][2])),y=c(20/10,(20/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][3]),ytop = 30/10,ybottom=(30/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][2])),y=c(30/10,(30/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][3]),ytop = 40/10,ybottom=(40/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][2])),y=c(40/10,(40/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][3]),ytop = 50/10,ybottom=(50/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][2])),y=c(50/10,(50/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][3]),ytop = 60/10,ybottom=(60/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][2])),y=c(60/10,(60/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][3]),ytop = 70/10,ybottom=(70/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][2])),y=c(70/10,(70/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][3]),ytop = 80/10,ybottom=(80/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][2])),y=c(80/10,(80/10)+(1.2/3)),col='black',lwd=1)
  
  
  
  
}

plot_FigS5_latitudinalband_boxplots_extremebias_summary_GBOTB<-function(){  
  GBIFdata.BAMM.rates<-read.table('./output/tables/GBIFdata_GBOTB_BAMM_rates_table.txt',header=T,sep='\t',stringsAsFactors = F)
  GBIFdata.BAMM.rates<-GBIFdata.BAMM.rates[,c('binomial','lambda.avg')]
  latband.strapp.files<-list.files(path='./results/subsampled_extremebias_GBOTB//input_tables/',pattern='GBIFdata.BAMM.GBOTB.subsampled_extremebias.+.txt')
  latband.strapp.tables<-lapply(latband.strapp.files,function(x)read.table(paste('./results/subsampled_extremebias_GBOTB//input_tables/',x,sep=''),header=T,sep='\t',stringsAsFactors = F))
  latband.strapp.tables<-lapply(latband.strapp.tables,function(x)merge(x,GBIFdata.BAMM.rates))
  #get quartile 25,50,75 for each table across all lat bands
  latband.strapp.tables.quartiles<-lapply(latband.strapp.tables,function(x)aggregate(lambda.avg~latitudinal.band,data=x,function(x)quantile(x,c(0.25,0.5,0.75))))
  #get averages of quartiles 25,50,75 across all tables for each latitudinal band
  latband.average.quartiles<-lapply(seq(-60,90,10),function(x)lapply(latband.strapp.tables.quartiles,function(y){if(nrow(y[y$latitudinal.band==x,])==0){NA}else{y[y$latitudinal.band==x,]}}))
  latband.average.quartiles.all.tables<-lapply(seq(-60,90,10),function(x)lapply(latband.strapp.tables.quartiles,function(y){if(nrow(y[y$latitudinal.band==x,])==0){NA}else{y[y$latitudinal.band==x,]}}))
  latband.average.quartiles.all.tables.clean<-lapply(latband.average.quartiles.all.tables,function(x) x[!unlist(lapply(x,function(y){all(is.na(y))}))])
  latband.average.quartiles.all.tables.clean.df<-lapply(latband.average.quartiles.all.tables.clean,function(x)do.call('rbind',x))
  latband.average.quartiles.all.tables.clean.IQR.median<-lapply(latband.average.quartiles.all.tables.clean.df,function(x)unlist(lapply(c(1,2,3),function(y)mean(x[,2][,y]))))
  #####CHANGE PATH OF REALMS LAYER####
  realms.layer<-readOGR('./raw_data/realms/','newRealms')
  realms.layer<-gUnaryUnion(realms.layer)
  plot(realms.layer,lwd=0.5,border='grey80')
  #for (i in seq(from=85,to=-55,by=-10)){
  #  abline(h=i,col='grey70')
  #}
  abline(v=0)
  abline(h=23.5,lty=2,col='grey80')
  abline(h=-23.5,lty=2,col='grey80')
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-55,-55,-23.5,-23.5,-55),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(-23.5,-23.5,23.5,23.5,-23.5),col=adjustcolor(red,alpha.f = 0.3),border=NA)
  #polygon(x=c(min(u.box[c(1,2)]),max(u.box[c(1,2)]),max(u.box[c(1,2)]),min(u.box[c(1,2)]),min(u.box[c(1,2)])),y=c(85,85,23.5,23.5,85),col=adjustcolor(blue,alpha.f = 0.3),border=NA)
  u<-par('usr')
  #plot(c(0,0),type='n',xlim=c(u[1]/30,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #plot(c(0,0),type='n',xlim=c(u[1]/50,u[2]/30),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  #plot(c(0,0),type='n',xlim=c(-log(6),log(6)),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  plot(c(0,0),type='n',xlim=c(-2,1.5),ylim=c(u[3]/10,u[4]/10),xlab='',ylab='',xaxt='n',yaxt='n')
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  axis(1,c(log(c(0.25,1,2,5))),labels=c(0.25,1,2,5))
  
  #axis(1,seq(from=-3,to=3,by=1),labels=seq(from=-3,to=3,by=1))
  #axis(1,c(log(c(0.5,1,2,5))),labels=c(0.5,1,2,5))
  axis(4,seq(from=-50/10,to=80/10,by=10/10),labels=seq(from=-50,to=80,by=10),las=2,cex.axis=.5)
  #axis(1,c(log(0),log(3),log(6),log(9),log(12),log(15)),labels=c(0,3,6,9,12,15))
  abline(v=0)
  #par(new=T)
  
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][3]),ytop = -50/10,ybottom=(-50/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[2]][2])),y=c(-50/10,(-50/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][3]),ytop = -40/10,ybottom=(-40/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[3]][2])),y=c(-40/10,(-40/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][3]),ytop = -30/10,ybottom=(-30/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[4]][2])),y=c(-30/10,(-30/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][3]),ytop = -20/10,ybottom=(-20/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[5]][2])),y=c(-20/10,(-20/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][3]),ytop = -10/10,ybottom=(-10/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[6]][2])),y=c(-10/10,(-10/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][3]),ytop = 0/10,ybottom=(0/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[7]][2])),y=c(0/10,(0/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][3]),ytop = 10/10,ybottom=(10/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[8]][2])),y=c(10/10,(10/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][3]),ytop = 20/10,ybottom=(20/10)+(1.2/3),col=adjustcolor(red,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[9]][2])),y=c(20/10,(20/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][3]),ytop = 30/10,ybottom=(30/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[10]][2])),y=c(30/10,(30/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][3]),ytop = 40/10,ybottom=(40/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[11]][2])),y=c(40/10,(40/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][3]),ytop = 50/10,ybottom=(50/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[12]][2])),y=c(50/10,(50/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][3]),ytop = 60/10,ybottom=(60/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[13]][2])),y=c(60/10,(60/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][3]),ytop = 70/10,ybottom=(70/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[14]][2])),y=c(70/10,(70/10)+(1.2/3)),col='black',lwd=1)
  
  rect(xleft=log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][1]),xright=log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][3]),ytop = 80/10,ybottom=(80/10)+(1.2/3),col=adjustcolor(blue,alpha=.7),xaxt='n')
  lines(x=c(log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][2]),log(latband.average.quartiles.all.tables.clean.IQR.median[[15]][2])),y=c(80/10,(80/10)+(1.2/3)),col='black',lwd=1)
  
  
  
  
}
