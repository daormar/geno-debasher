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
    cat("run_sequenza [arguments]
 
Arguments:
-s     <string>        seqz file
-o     <string>        output directory
--help                 print this text\n")
  
    q(save="no",status=0)
}

## Check -s option
if(!"s" %in% keys)
{
    print("Error: -s option not provided")
    q(save="no",status=1)
}

## Check -o option
if(!"o" %in% keys)
{
    print("Error: -o option not provided")
    q(save="no",status=1)
}

# Specific libraries
library(sequenza)

# Extract information from seqz file (incorporates normalization steps)
seqzinfo <- sequenza.extract(s)

# Save file
save(seqzinfo,file=paste(o,"/info.seqz",sep=""))

# Infer cellularity and ploidy
CP.example <- sequenza.fit(seqzinfo)

# Report results
sequenza.results(sequenza.extract = seqzinfo, cp.table = CP.example, sample.id = "seqzout", out.dir=o)
