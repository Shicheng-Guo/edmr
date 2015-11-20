#Requirement and Installation of eDMR

Back to [eDMR](http://code.google.com/p/edmr)

# R Version Requirment #
R >=  2.15.0

# R Packages Dependencies Installation #
```
# install methylKit, GenomicRanges, data.table
source("http://methylkit.googlecode.com/files/install.methylKit.R")
install.methylKit(ver="0.9.2",dependencies=TRUE)
```

```
# install mixtools
source("http://bioconductor.org/biocLite.R")
biocLite("mixtools")
```
```
# install doMC
biocLite("doMC")
```
```
# install ggplot2
biocLite("ggplot2")
```

# eDMR installation #
```
source("http://edmr.googlecode.com/files/eDMR.v.0.3.R")
```

Back to [eDMR](http://code.google.com/p/edmr)