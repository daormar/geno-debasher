# *- bash -*

#################
# CFG FUNCTIONS #
#################

########
bam_ascatngs_shared_dirs()
{
    define_shared_dir ${DATADIR_BASENAME}
}

########
bam_ascatngs_fifos()
{
    :
}

######################
# BAM ASCATNGS STEPS #
######################

########
ascatngs_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_req_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    

    # -g option
    description="Sample gender: XX|XY"
    explain_cmdline_req_opt "-g" "<string>" "$description"    

    # -sg option
    description="SNP GC correction file"
    explain_cmdline_req_opt "-sg" "<string>" "$description"    

    # -mc option
    description="Name of male sex chromosome"
    explain_cmdline_req_opt "-mc" "<string>" "$description"    
}

########
ascatngs_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

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

    # -mc option
    define_cmdline_opt "$cmdline" "-mc" optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
ascatngs()
{
    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local gender=`read_opt_value_from_line "$*" "-g"`
    local malesexchr=`read_opt_value_from_line "$*" "-mc"`
    local snpgccorr=`read_opt_value_from_line "$*" "-sg"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    
    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate ascatngs 2>&1 || exit 1

    # Run ascat
    logmsg "* Executing ascat.pl..."
    ascat.pl -n ${normalbam} -t ${tumorbam} -r ${ref} -sg ${snpgccorr} -pr WGS -g ${gender} -gc ${malesexchr} -cpus ${cpus} -o ${step_outd} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}
