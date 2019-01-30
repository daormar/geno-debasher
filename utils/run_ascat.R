#!/usr/bin/env Rscript

#################################################################################
### Jose Espinosa-Carrasco AGENDAS-IRB Group. Nov 2018                        ###
#################################################################################
### Adapted from:                                                             ###
### https://github.com/SciLifeLab/Sarek/blob/master/bin/run_ascat.r           ###
###                                                                           ###
### Runs ascat using BAF and logR files processed from WGS data of a normal   ###
### and tumoral sample to obtain CNVs info                                    ###
### Input:                                                                    ###
### BAF tumor file                                                            ###
### LogR tumor fileshould contain a header describing the data                ###
### BAF normal file                                                           ###
### LogR normal file                                                          ###
### tumor name                                                                ###
### GC correction file (optional)                                             ###
###                                                                           ###
### Output:                                                                   ###
### ASCAT plots of the resulting CNV segmentation and tables                  ###
#################################################################################

#####################
### VARIABLES
#Reading arguments
args <- commandArgs (TRUE) #if not it doesn't start to count correctly

## Default setting when no arguments passed
if ( length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      run_ascat.R

      Arguments:
      --tumor_baf=path_tumor_baf   - character
      --tumor_logr=path_tumor_logr - character
      --ctrl_baf=path_ctrl_baf     - character
      --ctrl_logr=ctrl_logr_baf    - character
      --tumor_name=tumor_name      - character
      --help                       - print this text

      Example:
      ./run_ascat.R --tumor_baf=\"path_tumor_baf\" --tumor_logr=\"path_tumor_logr\" --normal_baf=\"path_ctrl_baf\" --normal_logr=\"path_ctrl_logr\" --tumor_name=\"tumor_name\" --gc_correction=\"path_gc_correction\" --out_dir=\"path_out_dir\" \n")

  q (save="no")
}

# Use to parse arguments beginning by --
parseArgs <- function(x)
{
  strsplit (sub ("^--", "", x), "=")
}

#Parsing arguments
argsDF <- as.data.frame (do.call("rbind", parseArgs(args)))
argsL <- as.list (as.character(argsDF$V2))
names (argsL) <- argsDF$V1

# path to tumor baf file
{
  if (is.null (argsL$tumor_baf)) {
      stop ("[FATAL]: Path to tumor BAF file is mandatory")
  }
  else {
      tumor_baf <- argsL$tumor_baf
  }
}

# path to tumor LogR file
{
  if (is.null (argsL$tumor_logr)) {
      stop ("[FATAL]: Path to tumor LogR file is mandatory")
  }
  else {
      tumor_logr <- argsL$tumor_logr
  }
}

# path to normal baf file
{
  if (is.null (argsL$normal_baf)) {
      stop ("[FATAL]: Path to normal BAF file is mandatory")
  }
  else {
      normal_baf <- argsL$normal_baf
  }
}

# path to normal LogR file
{
  if (is.null (argsL$normal_logr)) {
      stop ("[FATAL]: Path to normal LogR file is mandatory")
  }
  else {
      normal_logr <- argsL$normal_logr
  }
}

# tumor name id
{
  if (is.null (argsL$tumor_name)) {
      stop ("[FATAL]: Path to tumor name id")
  }
  else {
      tumor_name <- argsL$tumor_name
  }
}

# gc_correction file
{
  if (is.null (argsL$gc_correction)) {
      gc_correction <- FALSE
  }
  else {
      gc_correction <- argsL$gc_correction
  }
}

# output directory
{
  if (is.null (argsL$out_dir)) {
      stop ("[FATAL]: Output directory must be provided")
  }
  else {
      out_dir <- argsL$out_dir
  }
}

library("ASCAT")

options(bitmapType='cairo')

## Load the  data
ascat.bc <- ascat.loadData(Tumor_LogR_file = tumor_logr,
                           Tumor_BAF_file = tumor_baf,
                           Germline_LogR_file = normal_logr,
                           Germline_BAF_file = normal_baf)

## command line to launch run_ascat
# ~/git/benchmark_cnv_WGS/R/run_ascat.R --tumor_baf="tumor_dani_file.BAF" --tumor_logr="tumor_dani_file.LogR" \
#               --normal_baf="ctrl_dani_file.BAF" --normal_logr="ctrl_dani_file.LogR" --tumor_name="tumor_name" --gc_correction="SnpGcCorrections_GRCh37_1000g.numbered_snps.tsv"

{
    if (gc_correction==FALSE) {
        # ascat.bc <- ascat.runAscat(ascat.bc) ## do not apply anything if file for GC corrections is not provided by the user
        write("[INFO]: File for GC correction not provided, skipped.", stderr())

    }
    else {        	
	write(paste0("[INFO]: GC file used for the correction.", gc_correction), stderr())
        ascat.bc <- ascat.GCcorrect(ascat.bc, gc_correction)
	write("[INFO]: GC correction has being performed.", stderr())
    }
}

## bash command
# ./run_ascat.R --tumor_baf="${path2files}tumor_dani_file.BAF" --tumor_logr="${path2files}tumor_dani_file.LogR" --normal_baf="${path2files}ctrl_dani_file.BAF" \
#               --normal_logr="${path2files}ctrl_dani_file.LogR" --tumor_name="dani_file" --gc_correction="${path2files}SnpGcCorrections_GRCh37_1000g.numbered_snps.tsv"

#######################################################
### This part is performed with the GC uncorrected data
## Plot the raw data
ascat.plotRawData(ascat.bc, img.dir=out_dir)

## Segment the data
ascat.bc <- ascat.aspcf(ascat.bc)

## Plot the segmented data
ascat.plotSegmentedData(ascat.bc, img.dir=out_dir)

#####################################################
### This part is performed with the GC corrected data
## Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination, and discrete copy numbers
ascat.output <- ascat.runAscat(ascat.bc)

#Write out segmented regions (including regions with one copy of each allele)
write.table(ascat.output$segments, file=paste0(out_dir, tumor_name, ".segments.txt"), sep="\t", quote=F, row.names=F)

## Write out CNVs in bed format
cnvs <- ascat.output$segments[ascat.output$segments[,"nMajor"]!=1 | ascat.output$segments[,"nMinor"]!=1,2:6]
write.table(cnvs, file=paste0(out_dir, tumor_name,".cnvs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

## Write out purity and ploidy info
summary <- tryCatch({
            matrix(c(ascat.output$aberrantcellfraction, ascat.output$ploidy), ncol=2, byrow=TRUE)}, error = function(err) {
			    # error handler picks up where error was generated
			    print(paste("Could not find optimal solution:  ",err))
			    return(matrix(c(0,0),nrow=1,ncol=2,byrow = TRUE))
           })

colnames(summary) <- c("aberrant_cell_fraction","ploidy")
write.table(summary, file=paste0(out_dir, tumor_name,".purityploidy.txt"), sep="\t", quote=F, row.names=F, col.names=T)
