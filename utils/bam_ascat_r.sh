# *- bash -*

#################
# CFG FUNCTIONS #
#################

########
bam_ascat_r_shared_dirs()
{
    define_shared_dir ${DATADIR_BASENAME}
}

########
bam_ascat_r_fifos()
{
    :
}

######################
# ASCAT_R STEPS      #
######################

########
allele_counter_norm_explain_cmdline_opts()
{
    # -l option
    description="Loci (SNP position) file"
    explain_cmdline_req_opt "-l" "<string>" "$description"

     # -r option
    description="Reference genome file"
    explain_cmdline_req_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"
}

########
allele_counter_norm_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -l option
    define_cmdline_infile_opt "$cmdline" "-l" optlist || exit 1

    # -r option
    define_cmdline_infile_opt "$cmdline" "-r" optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
allele_counter_norm()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local locis=`read_opt_value_from_line "$*" "-l"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`

    # Activate conda environment if needed    
    logmsg "* Activating conda environment..."
    conda activate ascatngs 2>&1 || exit 1

    # Execute alleleCounter
    logmsg "* Executing alleleCounter..."
    alleleCounter -l ${locis} -r ${ref} -b ${normalbam} -o ${step_outd}/allele-counter-norm.csv 2>&1 || exit 1


    # Deactivate conda environment if needed
    logmsg "* Dectivating conda environment..."
    conda deactivate 2>&1
}

########
allele_counter_tumor_explain_cmdline_opts()
{
    # -l option
    description="Loci (SNP position) file"
    explain_cmdline_req_opt "-l" "<string>" "$description"

     # -r option
    description="Reference genome file"
    explain_cmdline_req_opt "-r" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"
}

########
allele_counter_tumor_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -l option
    define_cmdline_infile_opt "$cmdline" "-l" optlist || exit 1

    # -r option
    define_cmdline_infile_opt "$cmdline" "-r" optlist || exit 1

    # -tumorbam option
    local normalbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
allele_counter_tumor()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local locis=`read_opt_value_from_line "$*" "-l"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`


    # Activate conda environment if needed
    logmsg "* Activating conda environment..."
    conda activate ascatngs 2>&1 || exit 1
    
    # Execute alleleCounter
    logmsg "* Executing alleleCounter..."
    alleleCounter -l ${locis} -r ${ref} -b ${tumorbam} -o ${step_outd}/allele-counter-tumor.csv 2>&1 || exit 1


    # Deactivate conda environment if needed
    logmsg "* Dectivating conda environment..."
    conda deactivate 2>&1
}

########
ascatr_explain_cmdline_opts()
{
    # -acn option
    description="AlleleCounter normal file"
    explain_cmdline_opt "-acn" "<string>" "$description"

    # -act option
    description="AlleleCounter tumor file"
    explain_cmdline_opt "-act" "<string>" "$description"

    # -g option
    description="Sample gender: XX|XY"
    explain_cmdline_req_opt "-g" "<string>" "$description"

    # -sg option
    description="SNP GC correction file"
    explain_cmdline_req_opt "-sg" "<string>" "$description"
}

########
ascatr_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    local optlist=""
    
    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # Define alleleCounter-normal option or retrieve dependency
    local allelecountnorm_dep=`find_dependency_for_step "${jobspec}" allele_counter_norm`
    if [ ${allelecountnorm_dep} != ${DEP_NOT_FOUND} ]; then
        local acnorm_outd=`get_default_outd_for_dep ${outd} "${allelecountnorm_dep}"`
        local allelecount_norm_file=${acnorm_outd}/allele-counter-norm.csv
        define_opt "-alleleCounter-normal" ${allelecount_norm_file} optlist || exit 1
    else
        define_cmdline_infile_opt "${cmdline}" "-acn"  optlist || exit 1
    fi
    
    # Define alleleCounter-tumor option or retrieve dependency
    local allelecounttumor_dep=`find_dependency_for_step "${jobspec}" allele_counter_tumor`
    if [ ${allelecounttumor_dep} != ${DEP_NOT_FOUND} ]; then
        local actumor_outd=`get_default_outd_for_dep ${outd} "${allelecounttumor_dep}"`
        local allelecount_tumor_file=${acnorm_outd}/allele-counter-tumor.csv
        define_opt "-alleleCounter-tumor" ${allelecount_norm_file} optlist || exit 1
    else
        define_cmdline_infile_opt "${cmdline}" "-act"  optlist || exit 1
    fi
    
    # -g option
    define_cmdline_opt "$cmdline" "-g" optlist || exit 1

    # -sg option
    define_cmdline_infile_opt "$cmdline" "-sg" optlist || exit 1
 
    # Save option list
    save_opt_list optlist
}

########
ascatr()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local allelecounternormal=`read_opt_value_from_line "$*" "-alleleCounter-normal"`
    local allelecountertumor=`read_opt_value_from_line "$*" "-alleleCounter-tumor"`
    local gender=`read_opt_value_from_line "$*" "-g"`
    local snpgccorr=`read_opt_value_from_line "$*" "-sg"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate /home/jespinosa/conda_local_env/ascat_r 2>&1 || exit 1
    
    # Run convert allele count
    logmsg "* Executing convert_allele_counts..."
    echo "****************************${allelecountertumor}"
    echo "****************************${allelecounternormal}"
    echo "step out directory****************************${step_outd}"

    Rscript ${biopanpipe_bindir}/convert_allele_counts "tumor" ${allelecountertumor} "normal" ${allelecounternormal} ${gender} ${step_outd}
    
    echo "convert_allele_counts finished****************************"
    echo "gc_correction file********************************$snpgccorr" 

    # Run ascatr 
    logmsg "* Executing run_ascat..."
    Rscript ${biopanpipe_bindir}/run_ascat --tumor_baf="${step_outd}/tumor.BAF" --tumor_logr="${step_outd}/tumor.LogR" --normal_baf="${step_outd}/normal.BAF" --normal_logr="${step_outd}/tumor.LogR" --tumor_name="tumor" --gc_correction=${snpgccorr} --out_dir="${step_outd}/" 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}
