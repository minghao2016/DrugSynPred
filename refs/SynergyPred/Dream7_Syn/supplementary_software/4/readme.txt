R CMD INSTALL syngen_0.1.0.tar.gz

Or from the R terminal:

> install.packages("syngen_0.1.0.tar.gz", repos=NULL, type="source")

Example run:

# Loading the syngen package
> library(syngen)

# Accessing the documentation
> help(syngen)
> help(phenoDrug)
> help(dream6h)

# Running syngen on example data and displaying results
> data(dream6h)
> tmp <- phenoDrug(toxicly3, dream6h, nn=50, cutoff=ncol(dream6h), method="complemental", score=0)
> plot(tmp)