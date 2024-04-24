library(dplyr)

library(ggplot2)

# dir with source files (data tables in .csv format)
adir<-'D:\\wrk\\Praca_13.10.2023-full\\Praca\\XScans\\VirScan\\03.23-analysis\\analysis-1.6\\'
# preparation of data tables for calculating probability of accidental results in comparison to a control
calc.p.sample<-function(i, wdir){
  #reading a file representing one sample
  dt<-read.csv(paste0(adir, i), stringsAsFactors = FALSE, header = TRUE)
  sum<-sum(dt$count) #for later calculation of signal
  #various dta for later calculations of probability
  dt$signal<-as.double(dt$count/sum)
  if (nrow(dt[dt$count==0,])>0){#setting signal for minimal value for left censored data to avoid false positives
    dt[dt$count==0,]$signal<-min(dt[dt$count>0,]$signal)
  }
  dt$log.signal<-log(dt$signal, 10)
  dt$control.signal<-10^(dt$mean.control.log.signal)
  dt$relative.signal<-dt$signal/dt$control.signal
  dt$sum<-sum

  return(dt)
}
#for calculation of probability
fd.tr<-function(x){
x[23]<-pnorm(q=as.double(x[6]), mean=as.double(x[3]), sd=as.double(x[4]), log=FALSE, lower.tail=FALSE)
x[24]<-pnorm(q=as.double(x[15]), mean=as.double(x[3]), sd=as.double(x[4]), log=FALSE, lower.tail=FALSE)
x[25]<-pnorm(q=as.double(x[21]), mean=as.double(x[3]), sd=as.double(x[4]), log=FALSE, lower.tail=FALSE)
return(x)
}
#joining technical replicates
calc.tech.reps<-function(trep1, trep2, sample){
    mdt<-merge(x=trep1, y=trep2, all.x=TRUE, all.y=TRUE, by='olgnname')
    mdt$mean.signal<-(mdt$signal.x+mdt$signal.y)/2
    mdt$relative.signal<-mdt$mean.signal/mdt$control.signal.x
    mdt$mean.log.signal<-log(mdt$mean.signal, 10)
    mdt$mean.sd.signal<-((mdt$log.signal.x-mdt$mean.log.signal)^2+(mdt$log.signal.y-mdt$mean.log.signal)^2)^0.5
    mdt$p.value.x<-as.double(1.0)
    mdt$p.value.y<-as.double(1.0)
    mdt$p.value<-as.double(1.0)
    nzmdt<-mdt[mdt$count.x>0&mdt$count.y>0,]
    zmdt<-mdt[!(mdt$olgnname %in% nzmdt$olgnname),]
    amdt<-apply(nzmdt, MARGIN=1, fd.tr)
    tamdt<-t(amdt)
    tamdt<-as.data.frame(tamdt)
    tamdt[,2:25]<-sapply(tamdt[,2:25], as.double)
    tamdt$detected.adj.p<-as.double(1.0)
    tamdt$detected.adj.p<-p.adjust(tamdt$p.value, method='fdr')
    zmdt$detected.adj.p<-as.double(1.0)
    mdt<-rbind(tamdt, zmdt)
    omdt<-mdt[order(mdt$p.value, decreasing = FALSE),]
    omdt$all.adj.p<-as.double(1.0)
    omdt$all.adj.p<-p.adjust(omdt$p.value, method='fdr')
    return(omdt)
}
fd<-function(x){
  x[8]<-pnorm(q=as.double(x[7]), mean=as.double(x[3]), sd=as.double(x[4]), log=FALSE, lower.tail=FALSE)
  return(x)
}
nm1<-'2.all-inf.all.tr1.csv'
nm2<-'2.all-inf.all.tr2.csv'
trep1<-calc.p.sample(i=nm1, wdir=wdir)
trep2<-calc.p.sample(i=nm2, wdir=wdir)
fdt<-calc.tech.reps(trep1=trep1, trep2=trep2, sample = nm)
fdt<-subset(fdt, select=c(olgnname, count.x, mean.control.log.signal.x, control.sd.log.signal.x, log.signal.x, relative.signal.x, sum.x, count.y,
                          log.signal.y, relative.signal.y, sum.y, mean.signal, relative.signal, mean.log.signal, mean.sd.signal, p.value.x, p.value.y, p.value, detected.adj.p, all.adj.p, sample))
write.csv(fdt, paste0(adir, '2.all-inf.all.csv'), row.names = FALSE)
write.csv(fdt[fdt$count.x>0&fdt$count.y>0,], paste0(adir, '2.all-inf.all-detected.csv'), row.names = FALSE)
write.csv(fdt[fdt$p.value<0.001&fdt$p.value.x<0.001&fdt$p.value.y<0.001,], paste0(adir, '2.all-inf.p0.001.csv'), row.names = FALSE)
