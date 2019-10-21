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

# Assign model results  containing the global statistics and 
# a file containing the model segmentation 
fitted_model <- emcncf(oo)

# Report Results (previously check if purity was equal to NA)
if(is.null(fitted_model$loglik) | is.na(fitted_model$purity))
{
    ## Tbl containing global statistics
    global_model_stats <- data.frame("NULL",
                                     fitted_model$purity,
                                     fitted_model$ploidy,
                                     fitted_model$dipLogR)

    colnames(global_model_stats) <- c("loglik", "purity", "ploidy", "dipLogR")
    write.csv(global_model_stats, paste0(o, "/facets_glob_mdl_stats.csv"), row.names = F)

    warning("Purity was equal to NA!")
}
else
{
    ## Tbl containing global statistics
    global_model_stats <- data.frame(fitted_model$loglik,
                                     fitted_model$purity,
                                     fitted_model$ploidy,
                                     fitted_model$dipLogR)

    colnames(global_model_stats) <- c("loglik", "purity", "ploidy", "dipLogR")
    write.csv(global_model_stats, paste0(o, "/facets_glob_mdl_stats.csv"), row.names = F)

    ## Tbl with segmentation results
    write.csv(fitted_model$cncf, paste0(o, "/facets_segmentation.csv"), row.names = F)

    ## Create facets plot
    png(paste0(o, "/facets_plot.png"), width = 900, height = 900, units = "px", pointsize = 18)
    plotSample(x = oo, emfit = fitted_model)
    dev.off()
}
