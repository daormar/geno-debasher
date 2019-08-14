# *- bash -*

# INCLUDE BASH LIBRARY
. ${biopanpipe_bindir}/bam_common_lib || exit 1

#################
# CFG FUNCTIONS #
#################

########
bam_analysis_shared_dirs()
{
    define_shared_dir ${DATADIR_BASENAME}
}

########
bam_analysis_fifos()
{
    :
}

# INCLUDE BASH FILES IMPLEMENTING STEPS
. ${biopanpipe_bindir}/genref_steps || exit 1
. ${biopanpipe_bindir}/bam_download_steps || exit 1
. ${biopanpipe_bindir}/bam_manip_steps || exit 1
. ${biopanpipe_bindir}/bam_analysis_steps || exit 1
. ${biopanpipe_bindir}/cleaning_steps || exit 1
