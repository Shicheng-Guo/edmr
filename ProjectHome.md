# eDMR #

Comprehensive DMR analysis based on bimodal normal distribution model and cost function for regional methylation analysis.

## Citation ##
Li S, Garrett-Bakelman FE, Akalin A, Zumbo P, Levine R, To BL, Lewis ID, Brown AL, D'Andrea RJ, Melnick A, Mason CE. An optimized algorithm for detecting and annotating regional differential methylation. BMC Bioinformatics. 2013;14 Suppl 5:S10.

If you used methylKit please also cite:
Akalin A, Kormaksson M, Li S, Garrett-Bakelman FE, Figueroa ME, Melnick A, Mason CE. methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. Genome Biol. 2012 Oct 3;13(10):[R87](https://code.google.com/p/edmr/source/detail?r=87).

## Installation (version 0.6.2) ##
```
source("http://methylkit.googlecode.com/files/install.methylKit.R")
install.methylKit(ver="0.9.2",dependencies=TRUE)
install.packages( c("mixtools", "devtools"))
source("http://bioconductor.org/biocLite.R")
biocLite("IRanges")
library(devtools)
install_url("https://github.com/ShengLi/edmr/archive/v0.6.2.tar.gz")
```

Note for edmr version 0.5.1 and before, please refer [here](https://code.google.com/p/edmr/wiki/edmr051) for installation and usage .

## Usage ##
Note that the eDMR input should be include all the CpGs sites for the significant test, not just the DMCs.

Step 1. Load add-on packages and example data
```
library(edmr)
library(methylKit)
library(GenomicRanges)
library(mixtools)
library(data.table)
data(example.myDiff.2013Nov6)
```

Step 2. myDiff evalution and plotting
```
# fitting the bimodal normal distribution to CpGs distribution
myMixmdl=myDiff.to.mixmdl(chr22.myDiff, plot=T, main="example")

# plot cost function and the determined distance cutoff
plotCost(myMixmdl, main="cost function")
```

Step 3. Calculate DMRs
```
# calculate all DMRs candidate
mydmr=edmr(chr22.myDiff, mode=1, ACF=TRUE)

# further filtering the DMRs
mysigdmr=filter.dmr(mydmr)

## annotation
# get genebody annotation GRangesList object
genebody=genebody.anno(file="http://edmr.googlecode.com/files/hg19_refseq_all_types.bed")

# plot the eDMR genebody annotation
plot.dmr.distr(mysigdmr, genebody, main="eDMR genebody annotation", xlab="DMR count")

# get CpG islands and shores annotation
cpgi=cpgi.anno(file="http://edmr.googlecode.com/files/hg19_cpgisland_all.bed")

# plot the eDMR CpG islands and shores annotation
plot.dmr.distr(mysigdmr, cpgi, main="eDMR CpG islands and shores annotation", xlab="DMR count")

# prepare genes for pathway analysis with significant DMRs at its promoter regions 
dmr.genes=get.dmr.genes(myDMR=mysigdmr, subject=genebody$promoter, id.type="gene.symbol")
dmr.genes
```

## Annotation file ##
The annotation files for hg19 and mm9 is available in Downloads. Instructions for annotation generation can be find [here](https://code.google.com/p/edmr/wiki/GeneAnnotation)

For example, input for genebody.anno(): https://edmr.googlecode.com/files/hg19_refseq_all_types.bed

Input for cpgi.anno(): https://edmr.googlecode.com/files/hg19_cpgisland_all.bed
