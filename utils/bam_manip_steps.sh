# *- bash -*

#############
# CONSTANTS #
#############

DEFAULT_MAX_RECORDS_IN_RAM_GATK=1000000

##########################
# BAM MANIPULATION STEPS #
##########################

########
index_norm_bam_explain_cmdline_opts()
{
    :
}

########
index_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
index_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Remove previous index if one was created
    if [ -f ${normalbam}.bai ]; then
        rm ${normalbam}.bai || exit 1
    fi
        
    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || exit 1

    # Execute samtools
    logmsg "* Executing samtools index..."
    samtools index ${normalbam} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
index_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
index_tum_bam_explain_cmdline_opts()
{
    :
}

########
index_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
index_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Remove previous index if one was created
    if [ -f ${tumorbam}.bai ]; then
        rm ${tumorbam}.bai || exit 1
    fi

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || exit 1
    
    # Execute samtools
    logmsg "* Executing samtools index..."
    samtools index ${tumorbam} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
index_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
sort_norm_bam_explain_cmdline_opts()
{
    :
}

########
sort_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`        
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
sort_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || exit 1

    # Verify if bam file is already sorted
    local bam_is_sorted=`samtools view -H ${normalbam} | $GREP SO:coordinate | wc -l` || exit 1
    if [ ${bam_is_sorted} -eq 1 ]; then
        echo "Warning: bam file is already sorted"
    else
        # Execute samtools
        logmsg "* Executing samtools sort..."
        samtools sort -T ${step_outd} -o ${step_outd}/sorted.bam -m 2G -@ ${cpus} ${normalbam} 2>&1 || exit 1
        # NOTE: -m option is used here to increase the maximum memory per
        # thread. One lateral efect of this is that the number of tmp files
        # generated is decreased. This constitutes one possible way to avoid
        # the "Too many open files" error reported by samtools

        # Replace initial bam file by the sorted one
        mv ${step_outd}/sorted.bam ${normalbam} 2>&1 || exit 1
    fi
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
sort_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
sort_tum_bam_explain_cmdline_opts()
{
    :
}

########
sort_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
sort_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    
    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || exit 1

    # Verify if bam file is already sorted
    local bam_is_sorted=`samtools view -H ${tumorbam} | $GREP SO:coordinate | wc -l` || exit 1
    if [ ${bam_is_sorted} -eq 1 ]; then
        echo "Warning: bam file is already sorted"
    else
        # Execute samtools
        logmsg "* Executing samtools sort..."
        samtools sort -T ${step_outd} -o ${step_outd}/sorted.bam -m 2G -@ ${cpus} ${tumorbam} 2>&1 || exit 1
        # NOTE: -m option is used here to increase the maximum memory per
        # thread. One lateral efect of this is that the number of tmp files
        # generated is decreased. This constitutes one possible way to avoid
        # the "Too many open files" error reported by samtools

        # Replace initial bam file by the sorted one
        mv ${step_outd}/sorted.bam ${tumorbam} 2>&1 || exit 1
    fi
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
sort_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
sambamba_mpileup_norm_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"
}

########
sambamba_mpileup_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
get_sambamba_mpileup_l_opt()
{
    local mpbfile=$1

    if [ "${mpbfile}" = ${OPT_NOT_FOUND} ]; then
        echo ""
    else
        echo "-L ${mpbfile}"
    fi
}

########
sambamba_mpileup_norm_bam()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Obtain sambamba mpileup -L opt
    local smp_l_opt=`get_sambamba_mpileup_l_opt ${mpbfile}`

    # Generate pileup file
    logmsg "* Generating pileup file..."
    sambamba mpileup -t ${cpus} ${smp_l_opt} --tmpdir ${step_outd} -o ${step_outd}/normal.pileup $normalbam --samtools "-f ${ref}" || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/normal.pileup || exit 1
}

########
sambamba_mpileup_norm_bam_conda_envs()
{
    define_conda_env sambamba sambamba.yml
}

########
sambamba_mpileup_tum_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"
}

########
sambamba_mpileup_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
sambamba_mpileup_tum_bam()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Obtain sambamba mpileup -L opt
    local smp_l_opt=`get_sambamba_mpileup_l_opt ${mpbfile}`
    
    # Generate pileup file
    logmsg "* Generating pileup file..."
    sambamba mpileup -t ${cpus} ${smp_l_opt} --tmpdir ${step_outd} -o ${step_outd}/tumor.pileup $tumorbam --samtools "-f ${ref}" || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/tumor.pileup || exit 1
}

########
sambamba_mpileup_tum_bam_conda_envs()
{
    define_conda_env sambamba sambamba.yml
}

########
parallel_sambamba_mpileup_norm_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_sambamba_mpileup_norm_bam_define_opts()
{ 
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -datadir option    
    abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        normalbam=${abs_datadir}/normal_${contig}.bam
        define_opt "-normalbam" ${normalbam} specific_optlist || exit 1
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_sambamba_mpileup_norm_bam()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Obtain sambamba mpileup -L opt
    local smp_l_opt=`get_sambamba_mpileup_l_opt ${mpbfile}`

    # Reset tmp directory
    tmpdir=${step_outd}/tmp_${contig}
    if [ -d $tmpdir ]; then
        rm -rf $tmpdir/*
    else
        mkdir $tmpdir || exit 1
    fi

    # Generate pileup file
    logmsg "* Generating pileup file (contig $contig)..."
    sambamba mpileup -t ${cpus} ${smp_l_opt} --tmpdir ${tmpdir} -o ${step_outd}/normal_${contig}.pileup $normalbam --samtools "-f ${ref}" || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/normal_${contig}.pileup || exit 1
}

########
parallel_sambamba_mpileup_norm_bam_reset_outdir()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Remove files
    logmsg "* Resetting output directory..."
    rm -f ${step_outd}/normal_${contig}.*
}

########
parallel_sambamba_mpileup_norm_bam_conda_envs()
{
    define_conda_env sambamba sambamba.yml
}

########
parallel_sambamba_mpileup_tum_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_sambamba_mpileup_tum_bam_define_opts()
{ 
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -datadir option    
    abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        tumorbam=${abs_datadir}/tumor_${contig}.bam
        define_opt "-tumorbam" ${tumorbam} specific_optlist || exit 1
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_sambamba_mpileup_tum_bam()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Obtain sambamba mpileup -L opt
    local smp_l_opt=`get_sambamba_mpileup_l_opt ${mpbfile}`

    # Reset tmp directory
    tmpdir=${step_outd}/tmp_${contig}
    if [ -d $tmpdir ]; then
        rm -rf $tmpdir/*
    else
        mkdir $tmpdir || exit 1
    fi

    # Generate pileup file
    logmsg "* Generating pileup file (contig $contig)..."
    sambamba mpileup -t ${cpus} ${smp_l_opt} --tmpdir ${tmpdir} -o ${step_outd}/tumor_${contig}.pileup $tumorbam --samtools "-f ${ref}" || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/tumor_${contig}.pileup || exit 1
}

########
parallel_sambamba_mpileup_tum_bam_reset_outdir()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Remove files
    logmsg "* Resetting output directory..."
    rm -f ${step_outd}/tumor_${contig}.*
}

########
parallel_sambamba_mpileup_tum_bam_conda_envs()
{
    define_conda_env sambamba sambamba.yml
}

########
samtools_mpileup_norm_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"
}

########
samtools_mpileup_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # Save option list
    save_opt_list optlist    
}

########
get_samtools_mpileup_l_opt()
{
    local mpbfile=$1

    if [ "${mpbfile}" = ${OPT_NOT_FOUND} ]; then
        echo ""
    else
        echo "-l ${mpbfile}"
    fi
}

########
samtools_mpileup_norm_bam()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || exit 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=`get_samtools_mpileup_l_opt ${mpbfile}`

    # Generate pileup file
    logmsg "* Generating pileup file..."
    samtools mpileup ${smp_l_opt} -f ${ref} -o ${step_outd}/normal.pileup $normalbam || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/normal.pileup || exit 1
}

########
samtools_mpileup_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
samtools_mpileup_tum_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"
}

########
samtools_mpileup_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # Save option list
    save_opt_list optlist    
}

########
samtools_mpileup_tum_bam()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || exit 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=`get_samtools_mpileup_l_opt ${mpbfile}`
    
    # Generate pileup file
    logmsg "* Generating pileup file..."
    samtools mpileup ${smp_l_opt} -f ${ref} -o ${step_outd}/tumor.pileup $tumorbam || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/tumor.pileup || exit 1
}

########
samtools_mpileup_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
parallel_samtools_mpileup_norm_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_samtools_mpileup_norm_bam_define_opts()
{ 
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -datadir option    
    abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        normalbam=${abs_datadir}/normal_${contig}.bam
        define_opt "-normalbam" ${normalbam} specific_optlist || exit 1
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_samtools_mpileup_norm_bam()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || exit 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=`get_samtools_mpileup_l_opt ${mpbfile}`

    # Generate pileup file
    logmsg "* Generating pileup file (contig $contig)..."
    samtools mpileup ${smp_l_opt} -f ${ref} -o ${step_outd}/normal_${contig}.pileup $normalbam || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/normal_${contig}.pileup || exit 1
}

########
parallel_samtools_mpileup_norm_bam_reset_outdir()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Remove files
    logmsg "* Resetting output directory..."
    rm -f ${step_outd}/normal_${contig}.*
}

########
parallel_samtools_mpileup_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
parallel_samtools_mpileup_tum_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_samtools_mpileup_tum_bam_define_opts()
{ 
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -datadir option    
    abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        tumorbam=${abs_datadir}/tumor_${contig}.bam
        define_opt "-tumorbam" ${tumorbam} specific_optlist || exit 1
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_samtools_mpileup_tum_bam()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || exit 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=`get_samtools_mpileup_l_opt ${mpbfile}`

    # Generate pileup file
    logmsg "* Generating pileup file (contig $contig)..."
    samtools mpileup ${smp_l_opt} -f ${ref} -o ${step_outd}/tumor_${contig}.pileup $tumorbam || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/tumor_${contig}.pileup || exit 1
}

########
parallel_samtools_mpileup_tum_bam_reset_outdir()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Remove files
    logmsg "* Resetting output directory..."
    rm -f ${step_outd}/tumor_${contig}.*
}

########
parallel_samtools_mpileup_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
parallel_split_norm_bam_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_split_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -datadir option    
    abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    define_opt "-datadir" ${abs_datadir} optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_split_norm_bam()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local abs_datadir=`read_opt_value_from_line "$*" "-datadir"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    
    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || exit 1

    # Extract contig
    logmsg "* Extracting contig (contig $contig)..."
    normalcont=${step_outd}/normal_${contig}.bam
    filter_bam_contig_samtools $normalbam $contig $normalcont || exit 1

    # Index contig
    logmsg "* Indexing contig..."
    samtools index ${normalcont} || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Move bam and index file to datadir
    mv ${normalcont}* ${abs_datadir} || exit 1
}

########
parallel_split_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
parallel_split_tum_bam_explain_cmdline_opts()
{
    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"   
}

########
parallel_split_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -datadir option
    abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    define_opt "-datadir" ${abs_datadir} optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_split_tum_bam()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local abs_datadir=`read_opt_value_from_line "$*" "-datadir"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || exit 1

    # Extract contig
    logmsg "* Extracting contig (contig $contig)..."
    tumorcont=${step_outd}/tumor_${contig}.bam
    filter_bam_contig_samtools $tumorbam $contig $tumorcont || exit 1

    # Index contig
    logmsg "* Indexing contig..."
    samtools index ${tumorcont} || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
        
    # Move bam and index file to datadir
    mv ${tumorcont}* ${abs_datadir} || exit 1
}

########
parallel_split_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
bedtools_genomecov_norm_bam_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"
}

########
bedtools_genomecov_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
bedtools_genomecov_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
        
    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate bedtools 2>&1 || exit 1

    # Execute samtools
    logmsg "* Executing bedtools coverage..."
    bedtools genomecov -ibam ${normalbam} > ${step_outd}/genomecov.tsv || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
bedtools_genomecov_norm_bam_conda_envs()
{
    define_conda_env bedtools bedtools.yml
}

########
bedtools_genomecov_tum_bam_explain_cmdline_opts()
{
    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"
}

########
bedtools_genomecov_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
bedtools_genomecov_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
        
    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate bedtools 2>&1 || exit 1

    # Execute samtools
    logmsg "* Executing bedtools coverage..."
    bedtools genomecov -ibam ${tumorbam} > ${step_outd}/genomecov.tsv || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
bedtools_genomecov_tum_bam_conda_envs()
{
    define_conda_env bedtools bedtools.yml
}

########
norm_bam_to_ubam_explain_cmdline_opts()
{
    # -n option
    description="Normal unmapped bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_cmdline_opt "-mrec" "<int>" "$description"
}

########
norm_bam_to_ubam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
norm_bam_to_ubam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local max_records=`read_opt_value_from_line "$*" "-mrec"`

    # Create tmpdir for gatk
    tmpdir=${step_outd}/tmp
    mkdir ${tmpdir} || exit 1

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gatk4 2>&1 || exit 1

    # Execute gatk RevertSam
    logmsg "* Executing gatk RevertSam..."
    gatk --java-options "-Xmx4G" RevertSam --INPUT ${normalbam} --OUTPUT ${step_outd}/unmapped.bam --SANITIZE true --SORT_ORDER queryname --TMP_DIR ${tmpdir} --MAX_RECORDS_IN_RAM ${max_records} || exit 1

    # Replace initial bam file by the mapped one
    mv ${step_outd}/unmapped.bam ${normalbam} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
norm_bam_to_ubam_conda_envs()
{
    define_conda_env gatk4 gatk4.yml
}

########
tum_bam_to_ubam_explain_cmdline_opts()
{
    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_cmdline_opt "-mrec" "<int>" "$description"
}

########
tum_bam_to_ubam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
tum_bam_to_ubam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local max_records=`read_opt_value_from_line "$*" "-mrec"`

    # Create tmpdir for gatk
    tmpdir=${step_outd}/tmp
    mkdir ${tmpdir} || exit 1

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gatk4 2>&1 || exit 1

    # Execute gatk RevertSam
    logmsg "* Executing gatk RevertSam..."
    gatk --java-options "-Xmx4G" RevertSam --INPUT ${tumorbam} --OUTPUT ${step_outd}/unmapped.bam --SANITIZE true --SORT_ORDER queryname --TMP_DIR ${tmpdir} --MAX_RECORDS_IN_RAM ${max_records} || exit 1

    # Replace initial bam file by the mapped one
    mv ${step_outd}/unmapped.bam ${tumorbam} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
tum_bam_to_ubam_conda_envs()
{
    define_conda_env gatk4 gatk4.yml
}

########
align_norm_ubam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal unmapped bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_cmdline_opt "-mrec" "<int>" "$description"
}

########
align_norm_ubam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
align_norm_ubam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local max_records=`read_opt_value_from_line "$*" "-mrec"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Create tmpdir for gatk
    tmpdir=${step_outd}/tmp
    mkdir ${tmpdir} || exit 1

    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda activate gatk4 2>&1 || exit 1

    # Execute gatk SamToFastq
    logmsg "* Executing gatk SamToFastq..."
    gatk --java-options "-Xmx4G" SamToFastq --INPUT ${normalbam} --FASTQ ${step_outd}/reads_r1.fastq.gz --SECOND_END_FASTQ ${step_outd}/reads_r2.fastq.gz --TMP_DIR ${tmpdir} --MAX_RECORDS_IN_RAM ${max_records} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bwa)..."
    conda activate gatk4 2>&1 || exit 1

    # Execute bwa
    logmsg "* Executing bwa mem..."
    bwa mem -t ${cpus} ${ref} <(${GZIP} -d -c ${step_outd}/reads_r1.fastq.gz) <(${GZIP} -d -c ${step_outd}/reads_r2.fastq.gz) | ${GZIP} > ${step_outd}/aln.sam.gz ; pipe_fail || exit 1

    # Remove fastq files
    rm ${step_outd}/reads_r1.fastq.gz ${step_outd}/reads_r2.fastq.gz || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda activate gatk4 2>&1 || exit 1

    # Create dictionary
    if [ ! -f ${ref}.dict ]; then
        logmsg "* Creating dictionary for reference..."
        gatk --java-options "-Xmx4G" CreateSequenceDictionary -I ${ref}
    fi

    # Execute gatk
    logmsg "* Executing gatk CreateSequenceDictionary..."
    gatk --java-options "-Xmx4G" MergeBamAlignment --REFERENCE_SEQUENCE ${ref} --UNMAPPED_BAM ${normalbam} --ALIGNED_BAM ${step_outd}/aln.sam.gz --OUTPUT ${step_outd}/merged.bam --SORT_ORDER coordinate --TMP_DIR ${tmpdir} --MAX_RECORDS_IN_RAM ${max_records} || exit 1

    # Replace initial unmapped bam file by the mapped one
    mv ${step_outd}/merged.bam ${normalbam} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
align_norm_ubam_conda_envs()
{
    define_conda_env bwa bwa.yml
    define_conda_env gatk4 gatk4.yml
}

########
align_tum_ubam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_cmdline_opt "-mrec" "<int>" "$description"
}

########
align_tum_ubam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
align_tum_ubam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local max_records=`read_opt_value_from_line "$*" "-mrec"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Create tmpdir for gatk
    tmpdir=${step_outd}/tmp
    mkdir ${tmpdir} || exit 1
    
    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda activate gatk4 2>&1 || exit 1

    # Execute gatk SamToFastq
    logmsg "* Executing gatk SamToFastq..."    
    gatk --java-options "-Xmx4G" SamToFastq --INPUT ${tumorbam} --FASTQ ${step_outd}/reads_r1.fastq.gz --SECOND_END_FASTQ ${step_outd}/reads_r2.fastq.gz --TMP_DIR ${tmpdir} --MAX_RECORDS_IN_RAM ${max_records} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bwa)..."
    conda activate gatk4 2>&1 || exit 1

    # Execute bwa
    logmsg "* Executing bwa mem..."
    bwa mem -t ${cpus} ${ref} <(${GZIP} -d -c ${step_outd}/reads_r1.fastq.gz) <(${GZIP} -d -c ${step_outd}/reads_r2.fastq.gz) | ${GZIP} > ${step_outd}/aln.sam.gz ; pipe_fail || exit 1

    # Remove fastq files
    rm ${step_outd}/reads_r1.fastq.gz ${step_outd}/reads_r2.fastq.gz || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda activate gatk4 2>&1 || exit 1

    # Create dictionary
    if [ ! -f ${ref}.dict ]; then
        logmsg "* Executing gatk CreateSequenceDictionary..."
        gatk --java-options "-Xmx4G" CreateSequenceDictionary -I ${ref}
    fi

    # Execute gatk
    logmsg "* Executing gatk MergeBamAlignment..."
    gatk --java-options "-Xmx4G" MergeBamAlignment --REFERENCE_SEQUENCE ${ref} --UNMAPPED_BAM ${tumorbam} --ALIGNED_BAM ${step_outd}/aln.sam.gz --OUTPUT ${step_outd}/merged.bam --SORT_ORDER coordinate --TMP_DIR ${tmpdir} --MAX_RECORDS_IN_RAM ${max_records} || exit 1

    # Replace initial unmapped bam file by the mapped one
    mv ${step_outd}/merged.bam ${tumorbam} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
align_tum_ubam_conda_envs()
{
    define_conda_env bwa bwa.yml
    define_conda_env gatk4 gatk4.yml
}
