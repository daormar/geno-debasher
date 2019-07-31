# *- r -*

#################################################################################
### Jose Espinosa-Carrasco AGENDAS-IRB Group. Nov 2018                        ###
#################################################################################
### Adapted from:                                                             ###
### https://github.com/SciLifeLab/Sarek/blob/master/bin/convertAlleleCounts.r ###
###                                                                           ###
### R-script for converting output from AlleleCount to BAF and LogR values.   ###
###                                                                           ###
### Input:                                                                    ###
### AlleleCounter output file for tumor and normal samples                    ###
### The first line should contain a header describing the data                ###
### The following columns and headers should be present:                      ###
### CHR    POS     Count_A Count_C Count_G Count_T Good_depth                 ###
###                                                                           ###
### Output:                                                                   ###
### BAF and LogR tables (tab delimited text files)                            ###
#################################################################################

## Read command line arguments
args = commandArgs(trailingOnly=TRUE)

## args is now a list of character vectors
## First check to see if arguments are passed.
if(length(args)<5){
    stop("No input files supplied\n\nUsage:\nRscript convert_allele_counts tumorid tumorac normalid normalac gender outdir\nWhere:\ntumorid - id of tumor sample\ntumorac - output from AlleleCount for the tumor\nnormalid - id of normal sample\nnormalac - output from AlleleCount for the normal\ngender - XX or XY\noutdir - output directory\n\n")
} else{
    tumor_id = args[1]
    tumor_ac = args[2]
    normal_id = args[3]
    normal_ac = args[4]
    gender = args[5]
    out_dir = args[6]
}

tumor_counts <- read.table(tumor_ac, header=F, skip=1, sep="\t")
normal_counts <- read.table(normal_ac, header=F, skip=1, sep="\t")

SNP_pos <- matrix(nrow = dim(normal_counts)[1], ncol = 2)

#Change rownames to "chr_pos" instead, such as 1_44552
#This does not exactly work:
#rownames(SNP_pos) <- apply(cbind(tumor_counts[,1], tumor_counts[,2]), 1, paste, collapse="_")
#This is for compatibility with gc correction file

colnames(SNP_pos) <- c("Chr", "Position")
SNP_pos[,1] <- as.vector(normal_counts[,1])
SNP_pos[,2] <- normal_counts[,2]

#Calculate BAF
tumor_BAF <- matrix(nrow = dim(normal_counts)[1], ncol = 1)
## rownames(tumor_BAF) <- rownames(SNP_pos)
colnames(tumor_BAF) <- c(tumor_id)
acgt <- tumor_counts[ ,c(3:6)]
acgts <- t(apply(acgt, 1, sort))
tumor_BAF[,1] <- acgts[,4] / (acgts[,3] + acgts[,4])
tumor_BAF[,1] <- ifelse(runif(length(tumor_BAF[,1])) < 0.5, tumor_BAF[,1], 1 - tumor_BAF[,1])
tumor_BAF[is.nan(tumor_BAF)] <- NA

germline_BAF <- matrix(nrow = dim(normal_counts)[1], ncol = 1)
colnames(germline_BAF) <- c(normal_id)
acgt <- normal_counts[,c(3:6)]
acgts <- t(apply(acgt, 1, sort))
germline_BAF[,1] <- acgts[,4] / (acgts[,3] + acgts[,4])
germline_BAF[,1] <- ifelse(runif(length(germline_BAF[,1])) < 0.5, germline_BAF[,1], 1 - germline_BAF[,1])
germline_BAF[is.nan(germline_BAF)] <- NA

tumor_LogR <- matrix(nrow = dim(normal_counts)[1], ncol = 1)
germline_LogR <- matrix(nrow = dim(normal_counts)[1], ncol = 1)
colnames(tumor_LogR) <- c(tumor_id)
colnames(germline_LogR) <- c(normal_id)
tumor_LogR[,1] <- log(tumor_counts[,7] / normal_counts[,7], 2)
germline_LogR[,1] <- 0
tumor_LogR[is.infinite(tumor_LogR)] <- NA

if (gender=="XY") {
    tumor_LogR[SNP_pos[,1] == "X", 1] <- tumor_LogR[SNP_pos[,1]=="X", 1] - 1
    germline_LogR[SNP_pos[,1] == "X", 1] <- germline_LogR[SNP_pos[,1] == "X", 1] - 1
}
tumor_LogR[,1] <- tumor_LogR[ ,1] - median(tumor_LogR[ ,1], na.rm=T)
# set regions with 0 reads in tumor and normal to a LogR of 0.
tumor_LogR[is.na(tumor_LogR[ ,1]), 1] <- 0

# limit the number of digits:
tumor_LogR <- round(tumor_LogR,4)
tumor_BAF <- round(tumor_BAF,4)
germline_LogR <- round(germline_LogR,4)
germline_BAF <- round(germline_BAF,4)

# write output to files
write.table(cbind (SNP_pos, tumor_LogR), paste(out_dir,"/", tumor_id, ".LogR", sep=""), sep="\t", row.names=F, col.names=T, quote=F)
write.table(cbind (SNP_pos, tumor_BAF), paste(out_dir, "/", tumor_id, ".BAF", sep=""), sep="\t", row.names=F, col.names=T, quote=F)
write.table(cbind (SNP_pos, germline_LogR), paste(out_dir,"/", normal_id, ".LogR", sep=""), sep="\t", row.names=F, col.names=T, quote=F)
write.table(cbind (SNP_pos, germline_BAF), paste(out_dir, "/", normal_id, ".BAF", sep=""), sep="\t", row.names=F, col.names=T, quote=F)
