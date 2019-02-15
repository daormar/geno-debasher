# *- bash -*

#############################
# SNP-PILEUP + FACETS STEPS #
#############################

########
snp_pileup_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -sv option
    description="SNP vcf file required by snp pileup"
    explain_cmdline_opt "-sv" "<string>" "$description"
}

########
snp_pileup_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -sv option
    define_cmdline_infile_opt "$cmdline" "-sv" optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
snp_pileup()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local snpvcf=`read_opt_value_from_line "$*" "-sv"`

    # Activate conda environment if needed
    if [ -z "${FACETS_HOME_DIR}" ]; then
        logmsg "* Activating conda environment..."
        conda activate facets 2>&1 || exit 1
    fi

    # Execute snp-pileup
    logmsg "* Executing snp-pileup..."
    if [ -z "${FACETS_HOME_DIR}" ]; then
        snp-pileup ${snpvcf} ${step_outd}/snp-pileup-counts.csv ${normalbam} ${tumorbam} 2>&1 || exit 1
    else
        ${FACETS_HOME_DIR}/inst/extcode/snp-pileup ${snpvcf} ${step_outd}/snp-pileup-counts.csv ${normalbam} ${tumorbam} 2>&1 || exit 1
    fi

    # Deactivate conda environment if needed
    if [ -z "${FACETS_HOME_DIR}" ]; then
        logmsg "* Dectivating conda environment..."
        conda deactivate 2>&1
    fi

    display_end_step_message
}

#######
facets_explain_cmdline_opts()
{
    # -pileup-counts option
    description="SNP pileup file (required if pileup step has not been performed)"
    explain_cmdline_opt "-sp" "<string>" "$description"
}

########
facets_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    local pileup_dep=`find_dependency_for_step "${jobspec}" snp_pileup`
    if [ ${pileup_dep} != ${DEP_NOT_FOUND} ]; then
        local pileup_outd=`get_default_outd_for_dep ${outd} "${pileup_dep}"`
        local pileup_counts_file=${pileup_outd}/snp-pileup-counts.csv
        define_opt "-pileup-counts" ${pileup_counts_file} optlist || exit 1
    else
        define_cmdline_infile_opt "${cmdline}" "-sp" optlist || exit 1
    fi

    # Save option list
    save_opt_list optlist
}

########
facets()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local pileup_counts=`read_opt_value_from_line "$*" "-pileup-counts"`

    # Define --snpPileupCounts option if output from snp-pileup is available
    # snp_pìleup_opt=`get_snp_pileup_opt "${pileup_outd}"`

    # Activate conda environment if needed
    if [ -z "${FACETS_HOME_DIR}" ]; then
        logmsg "* Activating conda environment..."
        conda activate facets 2>&1 || exit 1
    fi

    # Execute facets
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing facets..."
    Rscript ${biopanpipe_bindir}/run_facets -c ${pileup_counts} -o ${step_outd} 2>&1 || exit 1

    # Deactivate conda environment if needed
    if [ -z "${FACETS_HOME_DIR}" ]; then
        logmsg "* Dectivating conda environment..."
        conda deactivate 2>&1
    fi

    display_end_step_message
}
