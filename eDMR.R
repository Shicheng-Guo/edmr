fast.read <- function(filename, header=TRUE, nrow=5, ...){
  top5=read.table(filename, header=header, stringsAsFactors=FALSE,nrow=nrow,...)
  return(read.table(filename, header=header, stringsAsFactors=FALSE,colClasses=sapply(top5, class), ...))
}

# cutoff optimization
get.dist.myDiff <- function(myDiff, cores=1){
  library(doMC, quietly=T)
  registerDoMC(cores=cores)
  x=myDiff
  chr.list=unique(as.character(x$chr))
  chr.list=chr.list[chr.list!="chrM"]
  y=foreach(chr = chr.list,.combine=c) %dopar% {diff(sort(x$start[x$chr==chr]))}
}

dist_to_mixmdl <- function(dist)
{
  library(mixtools,quietly=T)
  log2.distance=log2(dist[dist!=1])
  mixmdl=normalmixEM(log2.distance)
}

# get the optimized parameter dist value from mixmdl object
get.break_point=function(mixmdl){
  if(mixmdl$lambda[1]> mixmdl$lambda[2]) {i=1;k=2}
  else {i=2; k=1}
  break_point = optimize(
    f = function(x){
      p1=pnorm(x, mixmdl$mu[i],mixmdl$sigma[i], lower.tail = F)
      p2=pnorm(x, mixmdl$mu[k],mixmdl$sigma[k], lower.tail = T)
      cost=sum(mixmdl$lambda[i]*p1,mixmdl$lambda[k]*p2)
      return(cost)
    }
    , interval = c(2,max(mixmdl$x))
    , maximum=F
  )$minimum
  
  break_point
}

# plot the fitted bimodal normal distribution for CpGs distances distribution
plotMdl1=function(mixmdl, subtitle="", cex.sub=1,...){
  xlim=c(0,ceiling(max(mixmdl$x)))
  if(mixmdl$lambda[1] < mixmdl$lambda[2]) { 
    mixmdl$mu=mixmdl$mu[2:1]
    mixmdl$sigma=mixmdl$sigma[2:1]
    mixmdl$lambda=mixmdl$lambda[2:1]
  }
  plot(mixmdl,which=2, xlim=xlim,  breaks=seq(0,max(xlim),by=1), ...)
  lines(density(mixmdl$x, n=50), lty=2, lwd=2)
  legend("topright", c("First model","Second model","data density"), lty=c(1,1,2), pch="", col=c("red","green","black"))
  mtext(subtitle, cex=cex.sub)
}

myDiff.to.mixmdl=function(myDiff, cores=1,plot=F, main=""){
  dist=get.dist.myDiff(myDiff, cores=cores)
  mixmdl=dist_to_mixmdl(dist)
  if(plot){
    plotMdl1(mixmdl,main)
  }
  print(2^get.break_point(mixmdl))
  mixmdl
}

get.dist.cutoff <- function(mixmdl){
  dist=round(2^get.break_point(mixmdl))
  dist
}

plotCost=function(mixmdl, ...){
  xlim=c(0,ceiling(max(mixmdl$x)))
  f = function(mixmdl,x){
    if(mixmdl$lambda[1]> mixmdl$lambda[2]) {i=1;k=2}
    else {i=2; k=1}
    p1=pnorm(x, mixmdl$mu[i],mixmdl$sigma[i], lower.tail = F)
    p2=pnorm(x, mixmdl$mu[k],mixmdl$sigma[k], lower.tail = T)
    cost=sum(mixmdl$lambda[i]*p1,mixmdl$lambda[k]*p2)
    return(cost)
  }
  x=seq(min(xlim), max(xlim), by=0.1)
  y=unlist(lapply(x, function(i){f(mixmdl,i)}))
  plot(x,y, type="n",xlab="log 2 distance", ylab="weighted sum of penelty", ...)
  lines(x, y, col="blue", lwd=1.5)
  abline(v=get.break_point(mixmdl), col="red")
}

# DMR
sigDMR=function(DMR, meanmethdiff=20, num_probes=3, num_DMCs=1, r.sl.qadj=0.05){
  idx=which(abs(DMR$meanmethdiff) >= meanmethdiff & DMR$num_probes >=num_probes & DMR$num_DMCs >=num_DMCs & DMR$r.sl.qadj <= r.sl.qadj)
  DMR[idx,]
}

DMR.to.bedfile=function(DMR, bedfile){
  write.table(DMR, file=bedfile,row.names=F, col.names=F, quote=F, sep="\t", append=F)
}

gb.theme= function(){
  library(ggplot2, quietly=TRUE)
  theme(panel.background = element_rect(fill = "white",colour = "white"), 
        plot.background = element_rect(fill = "white",colour = "white"), 
        legend.position="top", 
        axis.text.x=element_text(angle=-90, size=16),
        axis.text.y=element_text(size=16),
        strip.text.y=element_text(size=16, colour="blue"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        text=element_text(size=20, colour="brown"))
} 
count.sigDMR.list=function(sigDMR.list, main="DMR count"){
  library(doMC, quietly=T)
  library(ggplot2, quietly=T)
  dat=foreach(i = 1:length(sigDMR.list), .combine=rbind) %do%  {
    x=sigDMR.list[[i]]; 
    count=c(length(which(x$mean.meth.diff>0)), length(which(x$mean.meth.diff<0))); 
    type=c("hyper","hypo"); 
    data.frame(sample=names(sigDMR.list)[i], type=type, count=count)
  }
  ggplot(dat, aes(y=count,x=sample,fill=type)) + geom_bar(stat="identity") + ggtitle(main) + gb.theme()
}
length.sigDMR.list=function(sigDMR.list, main="DMR length", ...){
  dat=foreach(i = 1: length(sigDMR.list),.combine=rbind) %do% {x=(sigDMR.list[[i]]); data.frame(sample=names(sigDMR.list)[i], length=as.integer(x$end)-as.integer(x$start))}
  boxplot(length ~ sample, data=dat, main=main, ylab="bases",xlab="sample", ...)
}

genebody.sigDMR.list=function(files, ids, main="DMR gene body distribution"){
  library(doMC, quietly=T)
  library(ggplot2, quietly=T)
  dat=foreach( i = 1:length(files),.combine=rbind) %do% {x=read.table(files[i], header=F, stringsAsFactors=F); data.frame(sample=ids[i], type=x[,1], count=x[,2])}
  ggplot(dat, aes(y=count,x=sample,fill=type)) + geom_bar(stat="identity") + ggtitle(main) + gb.theme
}
# DMR calling
.readTableFast<-function(filename,header=T,skip=0,sep="", ...)
{
  tab5rows <- read.table(filename, header = header,skip=skip,sep=sep, nrows = 100, ...)
  classes  <- sapply(tab5rows, class)
  return( read.table(filename, header = header,skip=skip,sep=sep, colClasses = classes, ...)  )
}
getCorr=function(lag.idx,myDiff,i){
  res=c()
  lag.idx.x=c()
  lag.idx.y=c()
  for(k in 1:i){
    idx=which(lag.idx==k)
    lag.idx.x=c(lag.idx.x,myDiff$ppvalue[idx])
    lag.idx.y=c(lag.idx.y,myDiff$ppvalue[idx+1])
    corr=cor(lag.idx.x,lag.idx.y)
    if(corr==0) corr=9e-16
    res=c(res,abs(corr))
  }
  res
}
ACF=function(dist,step, myDiff){
  probes.dist=myDiff$pend[-1]-myDiff$pend[-nrow(myDiff)]
  
  dist.div=probes.dist/step
  lag.idx=ceiling(dist.div)
  dt=data.table(cbind(dist=1:dist,key=ceiling(1:dist/step)))
  acf.dist=dt[,list(start=min(dist), end=max(dist)), by=key]
  acf.dist$end[nrow(acf.dist)]=acf.dist$start[nrow(acf.dist)]+step-1
  acf=getCorr(lag.idx,myDiff,nrow(acf.dist)); 
  acfTable=cbind(V1=acf.dist$start,V2=acf.dist$end,V3=acf)
  data.frame(acfTable)
}

SLcombine=function(pvals,sigma) pnorm(sum(solve(chol(sigma))%*%qnorm(pvals))/sqrt(length(pvals)))

getSLp2=function(pvals,end,acf=acf, step=100){
  n=length(pvals)
  if(n>1){
    comb=combn(n,2)
    a=matrix(1, ncol=n, nrow=n)
    for(cc in 1:ncol(comb)){
      i=comb[1,cc]
      j=comb[2,cc]
      dist=end[j] - end[i]
      if(dist<1) dist=1
      corr=acf$V3[ceiling(dist/step)]
      a[i,j]=a[j,i]=corr
    }
    sigma=a
    pvals[pvals>=1]=1.0-9e-16
    SLcombine(pvals,sigma)      
  } else {
    pvals
  }
}

regionsToDMR=function(region_pvals, myDiffbed, step){
  regions_info=.readTableFast(region_pvals, header=F, stringsAsFactors=F)
  colnames(regions_info)=c("rchr","rstart","rend","rpval","num_probes","pchr","pstart","pend","ppvalue","pqvalue","pmethdiff")
  
  dist=max(regions_info$rend-regions_info$rstart)
  acf=ACF(dist,step, myDiffbed)
  
  dt=data.table(regions_info[regions_info$num_probes>=2,])
  res=dt[, list(medianmethdiff=median(pmethdiff), 
                meanmethdiff=mean(pmethdiff), 
                num_probes=length(ppvalue), 
                r.sl.pval=getSLp(ppvalue, pstart, pend, acf)), 
         by=list(rchr,rstart,rend)]
  DMR=cbind(res, padj=p.adjust(res$r.sl.pval,'fdr'))
  sig.DMR=DMR[DMR$padj<0.05 & abs(DMR$meanmethdiff)>20,]
  DMR
}
getPeaks2=function(allMyDiff, pcutoff=0.1, dist=100){
  print("raw myDiff:")
  print(dim(allMyDiff))

  myDiff=allMyDiff[allMyDiff$ppvalue<=pcutoff,]
  probes.dist=diff(myDiff$pend)
  print(paste("max cpgs dist for region definition:", dist))
  bpoints=which(probes.dist>=dist | probes.dist<0)
  first.peak=c(as.character(myDiff[1,1]), myDiff[1,2],end=myDiff[bpoints[1],3])
  mid.peaks=cbind(myDiff[(bpoints[-length(bpoints)]+1),1:2],end=myDiff[bpoints[-1],3])
  last.idx=bpoints[length(bpoints)]+1
  last.peak=c(as.character(myDiff[last.idx,1]), myDiff[last.idx,2], end=myDiff[nrow(myDiff),3])
  peaks=as.data.frame(rbind(first.peak, mid.peaks, last.peak))
  colnames(peaks)=c("rchr","rstart","rend")
  peaks
}
getDMR=function(peaks, allMyDiff, pcutoff=0.1,step=100, DMC.qvalue=0.01, DMC.methdiff=25, num.DMCs=1, num.CpGs=3, DMR.methdiff=20){
  library(GenomicRanges,quietly =TRUE)
  library(data.table,quietly =TRUE)
  myDiff=allMyDiff[allMyDiff$ppvalue<=pcutoff,]
  myDiff.gr=GRanges(seqnames=Rle(myDiff$pchr), IRanges(start=as.integer(myDiff$pstart), end=as.integer(myDiff$pend)), strand=myDiff$pstrand, pval=myDiff$ppvalue, padj=myDiff$pqvalue, methDiff=myDiff$pmethdiff)
  peaks.gr=GRanges(seqnames=Rle(peaks$rchr), IRanges(start=as.integer(peaks$rstart), end=as.integer(peaks$rend)))
  overlap.idx=findOverlaps(myDiff.gr, peaks.gr)
  dt.pk.myD=data.table(cbind(peaks[overlap.idx@subjectHits,], myDiff[overlap.idx@queryHits,]))
  refine.pk.myD=dt.pk.myD[, list(medianmethdiff=median(pmethdiff), 
                                 meanmethdiff=mean(pmethdiff), 
                                 num_probes=length(ppvalue), 
                                 num_DMCs=length(which(pqvalue<=DMC.qvalue & abs(pmethdiff)>=DMC.methdiff))),
                          by=list(rchr,rstart,rend)]
  refine.idx=which(refine.pk.myD$num_DMCs >= num.DMCs & refine.pk.myD$num_probes >= num.CpGs & abs(refine.pk.myD$meanmethdiff)>=DMR.methdiff)
  refine.pk.gr=refine.pk.myD[refine.idx,]
  refine.pk.idx=which(paste(dt.pk.myD$rchr,dt.pk.myD$rstart, dt.pk.myD$rend, sep="_") %in% paste(refine.pk.gr$rchr,refine.pk.gr$rstart, refine.pk.gr$rend, sep="_"))
  maxdist=max(dt.pk.myD[refine.pk.idx, as.integer(rend)-as.integer(rstart)])
  print(paste("auto correlation calculation. step:", step, "dist:", maxdist))
  acf=ACF(maxdist,step, allMyDiff)
  res.pk.myD=dt.pk.myD[refine.pk.idx, list(medianmethdiff=median(pmethdiff), 
                                           meanmethdiff=mean(pmethdiff), 
                                           num_probes=length(ppvalue), 
                                           num_DMCs=length(which(pqvalue <= DMC.qvalue & abs(pmethdiff) >= DMC.methdiff)),
                                           methdiffs=paste(pmethdiff, collapse=","),
                                           pvalues=paste(ppvalue,collapse=","),
										   r.sl.qval=getSLp2(pqvalue, pend, acf),
                                           r.sl.pval=getSLp2(ppvalue, pend, acf)),
                       by=list(rchr,rstart,rend)]
  r.sl.padj=p.adjust(res.pk.myD$r.sl.pval,'fdr')
  r.sl.qadj=p.adjust(res.pk.myD$r.sl.qval,'fdr')
  DMR=cbind(res.pk.myD, r.sl.qadj, r.sl.padj) 
  res=as.data.frame(DMR)
  colnames(res)=c("chr","start","end","median.meth.diff","mean.meth.diff","num.CpGs","num.DMCs","meth.diffs","pvalues","DMR.q.pvalue","DMR.p.pvalue","DMR.q.qvalue","DMR.p.qvalue")
  res
}
myDiffToDMR=function(myDiff, dist=100, step=100, DMC.qvalue=0.01, DMC.methdiff=25, num.DMCs=1, num.CpGs=3, DMR.methdiff=20){
  library(methylKit, quietly =TRUE)
  if (class(myDiff)=="methylDiff") input=getData(myDiff)
  else if(class(chr22.data)=="data.frame") input=myDiff
  else print("Input object myDiff has too be methylDiff class from methylKit or data.frame class")
  
  # prepare myDiff for DMR calling
  idx=with(input,order(chr, start))
  myDiff.new=cbind(input[idx,2], input[idx,3]-1, input[idx,4:8])
  colnames(myDiff.new)=c("pchr","pstart","pend","pstrand","ppvalue","pqvalue","pmethdiff")

  peaks=getPeaks2(myDiff.new, pcutoff=1, dist=dist)
  getDMR(peaks, myDiff.new, pcutoff=1, step=step, DMC.qvalue=DMC.qvalue, DMC.methdiff=DMC.methdiff, num.DMCs=num.DMCs, num.CpGs=num.CpGs, DMR.methdiff=DMR.methdiff)  
}
# main eDMR function
# myDiff 
eDMR=function(myDiff, step=100, dist="none", DMC.qvalue=0.01, DMC.methdiff=25, num.DMCs=1, num.CpGs=3, DMR.methdiff=20, granges=TRUE, plot=FALSE, main="", direction="both"){
  if(direction=="both"){
    print("DMR analysis for all detected CpGs...")
  } else if (direction=="hyper") {
    print("DMR analysis for hyper methylated CpGs...")
    idx=which(myDiff$meth.diff>0)
    myDiff=myDiff[idx,]
  } else if (direction=="hypo") {
    print("DMR analysis for hypo methylated CpGs...")
    idx=which(myDiff$meth.diff<0)
    myDiff=myDiff[idx,]
  } else {
    print ("parameter direction has too be both, hyper or hyper")
  }
  if(dist=="none"){
    mixmdl=myDiff.to.mixmdl(myDiff, plot=plot, main=main)
    dist=get.dist.cutoff(mixmdl)    
  }
  DMR=myDiffToDMR(myDiff, dist=dist, step=step, DMC.qvalue=DMC.qvalue, DMC.methdiff=DMC.methdiff, num.DMCs=num.DMCs, num.CpGs=num.CpGs, DMR.methdiff=DMR.methdiff)
  myDMR=DMR[,c(1:3,5:7,10,12)]
  if(granges) {
    myDMR.gr=with(myDMR, GRanges(chr, IRanges(as.numeric(start), as.numeric(end)), mean.meth.diff=mean.meth.diff, num.CpGs=num.CpGs, num.DMCs=num.DMCs, DMR.pvalue=DMR.q.pvalue, DMR.qvalue=DMR.q.qvalue))
    return(myDMR.gr)
  } else {
    return(myDMR)
  }
}

# significant DMRs
filter.dmr=function(myDMR, DMR.qvalue=0.05, mean.meth.diff=20, num.CpGs=3, num.DMCs=1){
  x=myDMR; 
  if(class(x)=="GRanges"){
    idx=which(values(myDMR)[,"DMR.qvalue"]<=DMR.qvalue & abs(values(myDMR)[,"mean.meth.diff"])>=mean.meth.diff & num.CpGs>=num.CpGs & num.DMCs >=num.DMCs)
    return(x[idx])
  } else if(x=="data.frame"){
    idx=which(myDMR$DMR.qvalue<=DMR.qvalue & abs(myDMR$mean.meth.diff) >= mean.meth.diff & num.CpGs>=num.CpGs & num.DMCs >=num.DMCs)  
  }
    return(x[idx,])
}
# annotation functions
splitn=function (strings, field, n)
{
  sapply(strsplit(strings, field), "[[", n)
}

## generate genebody GRangesList object
genebody.anno=function(file){
  print(paste("load", file))
  genes.obj=fast.read(file,header=F, sep="\t")
  colnames(genes.obj)=c("chr","start","end","id","score","strand","gene.id","gene.symbol")
  subj <- with(genes.obj, GRanges(chr, IRanges(start, end), id=id, gene.symbol=gene.symbol,gene.id=gene.id))
  types=unique(splitn(genes.obj$id,"_",3))
  types=c("up|utr5","utr5","cds","intron","utr3")
  subj.types=foreach(i = types) %do% {
    print(paste("process:",i))
    with(genes.obj[grep(i, splitn(genes.obj$id,"_",3)),], GRanges(chr, IRanges(start, end), id=id, gene.symbol=gene.symbol,gene.id=gene.id))
  }
  names(subj.types)=c("promoter","utr5","cds","intron","utr3")
  subj.types[["gene"]]=subj
  subj.types
}

## generate cpgi GRangesList object
cpgi.anno=function(file, shore.width=2000){
  ft=fast.read(file, header=F, sep="\t")
  cpgi.gr= with(ft, GRanges(ft[,1], IRanges(ft[,2], ft[,3])))
  up.shores.gr=flank(cpgi.gr, width=shore.width, start=T,both=F)
  dn.shores.gr=flank(cpgi.gr, width=shore.width, start=F,both=F)
  return(list(cpgis=cpgi.gr, shores=sort(c(up.shores.gr, dn.shores.gr))))
}
## plot eDMR distribution over the subject genomic ranges
plot.dmr.distr=function(myDMR, subject, ...){
  col.list=c("#E41A1C","#377EB8","#984EA3","#4DAF4A","#FF7F00","#FFFF33", "#A65628", "#8DD3C7"  )
  int=lapply(subject, function(x)intersect(myDMR,x))
  res0=sapply(int, length)
  int.gr=GRanges();for(i in length(int)) {int.gr=append(int.gr,int[[i]])}
  unannotated=length(myDMR)-length(int.gr)
  res0=c(res0, unannotated)
  names(res0)[length(res0)]="un-anno"
  print(res0)
  barplot(res0, col=col.list[1:length(res0)], las=2, ...)  
}

# get gene list based the genebody granges
get.dmr.genes=function(myDMR, subject, id.type="gene.symbol"){
  ind=findOverlaps(subject,filter.DMR(myDMR))
  unique(values(subject)[unique(ind@queryHits), id.type])
}