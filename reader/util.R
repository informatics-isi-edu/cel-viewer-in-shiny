

myfetch <- function(mfile) {
   d <- read.csv(file=mfile, header=TRUE, sep=",")
   if(!is.null(d)) {
     mymap = list()
     mymap$fname=d$cel_filename
     mymap$uname=d$user_samplename
     mymap$vname=d$viewewr_samplename
     return(mymap)
   } else return(NULL)
}

# alist is a list of charater strings ie. colnames(myExpr)
# given a list of cel file names returns
# a list of sample name preferred by viewer
mylookup <- function(mymap, alist) {
  ts <- character()
  for (t in alist) ts <- c(ts, mymap$vname[mymap$fname == t])
  return(ts)
}
