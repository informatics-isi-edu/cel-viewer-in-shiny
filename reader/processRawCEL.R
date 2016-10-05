##
##  Rscript processRawCEL.R path-to-cel-files path-to-mapping-file dat-file-name-to-use
##
##   (Rscript processRawCEL.R usc usc_mapping.csv usc_dat.R)
##  
## will generate full_dat.R, full_dat.rds
## dat is the name of the data structure
##

library(bioDist)
library(limma)
library(DT)
library(oligo)

library(mouse4302.db)
library(annotate)

source('util.R')

#http://www.perlmonks.org/?node_id=644967
Args <- commandArgs(TRUE) 
if (length(Args) == 3) {
   cel_dir <- c(Args[1])
   mapping_file <- c(Args[2])
   dat_file <- c(Args[3])
} else {  return(1) }

## read in raw cel files from a relative directory 'usc'
celFiles <- list.celfiles(cel_dir,full.names=TRUE)
myRawData <- read.celfiles(celFiles)
## smooth it
myEset <- rma(myRawData) 

# grab probesets, like 123456_at
probeset <- rownames(exprs(myEset))
symbol <- lookUp(probeset, "mouse4302.db", "SYMBOL")
desc <- lookUp(probeset, "mouse4302.db", "GENENAME")
myExpr <- exprs(myEset)

# grab the filename to target sample name mapping
mymap <- myfetch(mapping_file)
if(is.null(mymap)) return (0)
myflist <- colnames(myExpr)
newcolnames <- mylookup(mymap, myflist)
# replace it in the expression table
colnames(myExpr) <- newcolnames

dat <-as.data.frame(myExpr)
dat$probeset <- probeset
dat$symbol <- probeset
dat$desc <- desc


save(dat, file=dat_file)
#saveRDS(dat, 'full_dat.rds')

#### old code
#write(myEset,"myEset.R")
#save.csv(myEset,"myEset.R")
#write.csv(exprs(myEset), file = “myExprs.csv”)
##exprs(myEset) ## readable form,
## dump into csv, 
#write.csv(myEset, file = "myEsetMatrix.csv")

