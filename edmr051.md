#Installation and Usage for edmr version 0.5.1

# Dependencies installation #
[Installation](https://code.google.com/p/edmr/wiki/Installation)

# Input1: methylKit methylDiff object #
```
# source eDMR
con <- url("http://edmr.googlecode.com/files/eDMR.v.0.5.1.R")
source(con)
close(con)

# load example myDiff data
# this can be calculated using calculateDiffMeth() from methylKit package 
con <- url("http://edmr.googlecode.com/files/example.myDiff.2013Nov6.rda")
print(load(con))
close(con)

# fitting the bimodal normal distribution to CpGs distribution
myMixmdl=myDiff.to.mixmdl(chr22.myDiff)

# plot the fittings
plotMdl1(myMixmdl, subtitle="example", cex.sub=1.2)

# plot cost function and the determined distance cutoff
plotCost(myMixmdl, main="cost function")

# calling DMRs based on the optimized distance cutoff
myDMR=eDMR(chr22.myDiff,  # myDiff object, required
step=100, # step for calculating auto-correlation, default: 100
dist="none", # distance cutoff to call a gap for DMR, default: "none", which will be automatically determined by the bimodal normal distribution
DMC.qvalue=0.01, # qvalue cutoff for DMC definition, default: 0.01
DMC.methdiff=25,  # methylation difference cutoff for DMC definition, default: 25
num.DMCs=1,  # cutoff of the number DMCs in each region to call DMR, default: 1
num.CpGs=3,  # cutoff of the number of CpGs, default: 3
DMR.methdiff=20, # cutoff of the DMR mean methylation difference, default=20
plot=FALSE, # plot the bimodal normal distribution fitting or not, default=FAlSE
main="example", # the title of the plot, if plot=TRUE 
mode=1, # the mode of call DMRs. 1: using all CpGs together. 2: use unidirectional CpGs to call DMRs. default: 1
ACF=TRUE # p-value combination test with (TRUE, default) or without (FALSE) dependency adjustment. 
) 

# obtain the significant DMR
myDMR.sig=filter.dmr(myDMR, DMR.qvalue=0.05, mean.meth.diff=20, num.CpGs=3, num.DMCs=1)

# plot the distribution of the width of the DMR
plot.dmr.width(myDMR.sig, main="DMR width distribution")

# get genebody annotation GRangesList object
genebody=genebody.anno(file="http://edmr.googlecode.com/files/hg19_refseq_all_types.bed")

# plot the eDMR genebody annotation
plot.dmr.distr(myDMR.sig, genebody, main="eDMR genebody annotation", xlab="DMR count")

# get CpG islands and shores annotation
cpgi=cpgi.anno(file="http://edmr.googlecode.com/files/hg19_cpgisland_all.bed")

# plot the eDMR CpG islands and shores annotation
plot.dmr.distr(myDMR.sig, cpgi, main="eDMR CpG islands and shores annotation", xlab="DMR count")

# To save time for creating the annotation files for eDMR
# save for repeated usage and load the annotation R data
# save(genebody, cpgi, file="anno.rda")
# load("anno.rda")

# prepare genes for pathway analysis with significant DMRs at its promoter regions 
dmr.genes=get.dmr.genes(myDMR=myDMR.sig, subject=genebody$promoter, id.type="gene.symbol")

```

# Input2: bed file #
```
myDiff=fast.read("myDiff.bed")
myDMR=eDMR(myDiff)
```
"myDiff.bed" should include the following columns:

chr	start	end	strand	pvalue	qvalue	meth.diff