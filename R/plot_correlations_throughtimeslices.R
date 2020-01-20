plot_correlations_table_4Ma<-function(pgls.table,cladesizefile,name,parameter){
  clade.sizes<-read.table(cladesizefile,header=T,sep='\t',stringsAsFactors = F)
  pgls.table$Slice<-rownames(pgls.table)
  rownames(pgls.table)<-NULL
  pgls.table$Nclades<-clade.sizes$Nclades
  
  slice.table<-pgls.table
  #for seed size
  slice.table$point.size<-0
  if (nrow(slice.table[slice.table$Nclades<150,])>0){
    slice.table[(slice.table$Nclades<150),]$point.size<-1
  }
  if (nrow(slice.table[slice.table$Nclades>=150&slice.table$Nclades<300,])>0){
    slice.table[(slice.table$Nclades>=150)&(slice.table$Nclades<300),]$point.size<-2
  }
  if (nrow(slice.table[slice.table$Nclades>=300,])>0){
    slice.table[(slice.table$Nclades>=300),]$point.size<-3
  }
  
  slice.table$point.colour<-'light grey'
  if (nrow(slice.table[slice.table$pvalue<0.05,])>0){
    slice.table[slice.table$pvalue<0.05,]$point.colour<-'black'
  }
  if (nrow(slice.table[slice.table$pvalue>0.05&slice.table$pvalue<0.1,])>0){
    slice.table[slice.table$pvalue>0.05&slice.table$pvalue<0.1,]$point.colour<-'dark grey'
  }
  if (nrow(slice.table[slice.table$pvalue>0.1,])>0){
    slice.table[slice.table$pvalue>0.1,]$point.colour<-'white'
  }
  
  plot(0,type='n',axes=FALSE,ann=FALSE)
  symbols(seq(from=1,to=21,by=4),slice.table$slope,circles=slice.table$point.size,inches=0.2,bg=slice.table$point.colour,ylim=c(-0.4,0.4),main=paste(name,'~',parameter,sep=''),xaxt='n',ylab='pgls.slope',xlab='time.interval(myr)',cex.axis=.7,las=1,xaxs='i')
  axis(1,at=seq(from=0,to=24,by=4),labels=seq(from=0,to=24,by=4),cex.axis=.7)
  legend('topleft',c('< 50','50 - 75','> 75'),pch=21,col='black',pt.cex=c(5/3,5/2,5),bty='n',cex=1.5)
  legend('topright',c('< 0.05','0.05 - 0.1','> 0.1'),pch=19,col=c('black','dark grey','white'),pt.cex=5/2,bty='n',cex=1.5)
}



#for 4ma intervals
cladesizefile<-'./clades_size_4Ma.txt'
pgls.files<-list.files('./clade_analyses/RPANDA/',pattern='_0.3_5size_0.7.tropthreshold_pgls_table_new_GBIFsampling_log_modelave.txt')
pgls.files.4Ma<-sapply(c('0_4','4_8','8_12','12_16','16_20','20_24'),function(x)pgls.files[grep(x,pgls.files)])
pgls.files.4Ma.table<-lapply(pgls.files.4Ma,function(x)read.table(file=paste('./clade_analyses/RPANDA/',x,sep=''),sep='\t',header=T,stringsAsFactors=F))
pgls.files.4Ma.table.prop.trop.lambda<-do.call('rbind',lapply(pgls.files.4Ma.table,function(x)x[x$analysis=='pgls.prop.trop.log10.lambda',c('slope','pvalue')]))
pgls.files.4Ma.table.prop.trop.turnover<-do.call('rbind',lapply(pgls.files.4Ma.table,function(x)x[x$analysis=='pgls.prop.trop.log10.turnover',c('slope','pvalue')]))
pgls.files.4Ma.table.prop.temp.lambda<-do.call('rbind',lapply(pgls.files.4Ma.table,function(x)x[x$analysis=='pgls.prop.temp.log10.lambda',c('slope','pvalue')]))
pgls.files.4Ma.table.prop.temp.turnover<-do.call('rbind',lapply(pgls.files.4Ma.table,function(x)x[x$analysis=='pgls.prop.temp.log10.turnover',c('slope','pvalue')]))

pdf('./plots/timeslices4Ma_lambdaturnover.pdf',paper='a4r')
par(mfrow=c(2,2))
plot_correlations_table_4Ma(pgls.table = pgls.files.4Ma.table.prop.trop.lambda,cladesizefile ='./clades_size_4Ma.txt',name='prop.trop',parameter='log10.lambda')
plot_correlations_table_4Ma(pgls.table = pgls.files.4Ma.table.prop.trop.turnover,cladesizefile ='./clades_size_4Ma.txt',name='prop.trop',parameter='log10.turnover')

plot_correlations_table_4Ma(pgls.table = pgls.files.4Ma.table.prop.temp.lambda,cladesizefile ='./clades_size_4Ma.txt',name='prop.temp',parameter='log10.lambda')
plot_correlations_table_4Ma(pgls.table = pgls.files.4Ma.table.prop.temp.turnover,cladesizefile ='./clades_size_4Ma.txt',name='prop.temp',parameter='log10.lambda')
dev.off()

#for 4ma intervals weighted by proportion of 
cladesizefile<-'./clades_size_4Ma_GBOTB.txt'
pgls.files<-list.files('./clade_analyses/RPANDA/',pattern='_0.3_5size_0.7.tropthreshold_pgls_table_new_GBIFsampling_log_modelave_weightedlat.txt')
pgls.files.4Ma<-sapply(c('0_4','4_8','8_12','12_16','16_20','20_24'),function(x)pgls.files[grep(x,pgls.files)])
pgls.files.4Ma.table<-lapply(pgls.files.4Ma,function(x)read.table(file=paste('./clade_analyses/RPANDA/',x,sep=''),sep='\t',header=T,stringsAsFactors=F))
pgls.files.4Ma.table.prop.trop.lambda<-do.call('rbind',lapply(pgls.files.4Ma.table,function(x)x[x$analysis=='pgls.prop.trop.log10.lambda',c('slope','pvalue')]))
pgls.files.4Ma.table.prop.trop.turnover<-do.call('rbind',lapply(pgls.files.4Ma.table,function(x)x[x$analysis=='pgls.prop.trop.log10.turnover',c('slope','pvalue')]))
pgls.files.4Ma.table.prop.temp.lambda<-do.call('rbind',lapply(pgls.files.4Ma.table,function(x)x[x$analysis=='pgls.prop.temp.log10.lambda',c('slope','pvalue')]))
pgls.files.4Ma.table.prop.temp.turnover<-do.call('rbind',lapply(pgls.files.4Ma.table,function(x)x[x$analysis=='pgls.prop.temp.log10.turnover',c('slope','pvalue')]))

pdf('./plots/timeslices4Ma_lambdaturnover_weightedlat_GBOTB.pdf',paper='a4r')
par(mfrow=c(2,2))
plot_correlations_table_4Ma(pgls.table = pgls.files.4Ma.table.prop.trop.lambda,cladesizefile ='./clades_size_4Ma_GBOTB.txt',name='prop.trop',parameter='log10.lambda')
plot_correlations_table_4Ma(pgls.table = pgls.files.4Ma.table.prop.trop.turnover,cladesizefile ='./clades_size_4Ma_GBOTB.txt',name='prop.trop',parameter='log10.turnover')

plot_correlations_table_4Ma(pgls.table = pgls.files.4Ma.table.prop.temp.lambda,cladesizefile ='./clades_size_4Ma_GBOTB.txt',name='prop.temp',parameter='log10.lambda')
plot_correlations_table_4Ma(pgls.table = pgls.files.4Ma.table.prop.temp.turnover,cladesizefile ='./clades_size_4Ma_GBOTB.txt',name='prop.temp',parameter='log10.lambda')
dev.off()
