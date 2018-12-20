# *- bash -*

#############
# CONSTANTS #
#############

DEFAULT_ASP_MAX_TRANS_RATE=100m
DEFAULT_BAMDIR="data"

######################
# ASCAT_WGS_R STEPS  #
######################

########
ascat_wgs_r_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -g option
    description="Sample gender: XX|XY (required)"
    explain_cmdline_opt "-g" "<string>" "$description"

    # -sg option
    description="SNP GC correction file (optional if present correction is performed)"
    explain_cmdline_opt "-sg" "<string>" "$description"
}

########
ascat_wgs_r_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -r option
    define_cmdline_infile_opt "$cmdline" "-r" optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -g option
    define_cmdline_opt "$cmdline" "-g" optlist || exit 1

    # -sg option
    define_cmdline_infile_opt "$cmdline" "-sg" optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_jobspec "$jobspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
ascat_wgs_r()
{
    display_begin_step_message

    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local gender=`read_opt_value_from_line "$*" "-g"`
    local snpgccorr=`read_opt_value_from_line "$*" "-sg"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating ascatngs conda environment..."
    conda activate ascatngs 2>&1 || exit 1

    # Run alleleCounter
    logmsg "* Executing ascat.pl..."
    ascat.pl -n ${normalbam} -t ${tumorbam} -r ${ref} -sg ${snpgccorr} -pr WGS -g ${gender} -gc ${malesexchr} -cpus ${cpus} -o ${step_outd} 2>&1 || exit 1

    # Run convert allele counter
    	
    # Deactivate conda environment
    logmsg "* Activating ascat_r conda environment..."
    conda activate '/home/jespinosa/conda_local_env/ascat_r' 2>&1
	
    # Run 

    # Signal that step execution was completed
    signal_step_completion ${step_outd}

    display_end_step_message
}

