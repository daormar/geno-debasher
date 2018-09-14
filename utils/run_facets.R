# Author: Daniel Ortiz Mart\'inez
# *- r -*

# Basic libraries
library(R.utils)

## Collect arguments
args <- commandArgs(asValue=TRUE, excludeReserved=TRUE)[-1]
 
# Turn arguments into R variables
keys <- attachLocally(args)

# Check arguments

## Check --help option
if("help" %in% keys || length(keys)==0)
{    
    cat("run_facets [arguments]
 
Arguments:
-c     <string>        File containing snp-pileup counts
--help                 print this text\n")
  
    q(save="no",status=0)
}

## Check -c option
if(!"c" %in% keys)
{
    print("Error: -c option not provided")
    q(save="no",status=1)
}

# Specific libraries
library(pctGCdata)
library(facets)

# Run facets
rcmat <- readSnpMatrix(c)
xx <- preProcSample(rcmat)
my_cval=300
oo <- procSample(xx, cval = my_cval)
emcncf(oo)
