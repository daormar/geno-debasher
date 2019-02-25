# *- bash -*

#############
# CONSTANTS #
#############

DEFAULT_NUMBER_OF_DOWNLOAD_TRIES=5
DEFAULT_NUMBER_OF_EGA_DOWNLOAD_STREAMS=50
DEFAULT_ASP_MAX_TRANS_RATE=100m

#################
# CFG FUNCTIONS #
#################

########
bam_analysis_shared_dirs()
{
    define_shared_dir "data"
}

########
bam_analysis_fifos()
{
    :
}

######################
# BAM ANALYSIS STEPS #
######################

########
get_contig_list_from_file()
{
    local file=$1
    cat $file
}

########
get_ref_contig_list()
{
    local ref=$1

    if [ ! -f ${ref}.fai ]; then
        conda activate samtools 2>&1 || exit 1
        samtools faidx ${ref}
        conda deactivate
    fi
    
    $AWK '{printf " %s",$1}' ${ref}.fai
}

########
get_ref_filename()
{
    local cmdline=$1
    local given=0
    local ref
    ref=`read_opt_value_from_line "$cmdline" "-r"` && given=1
    if [ $given -eq 1 ]; then
        # -r option was given
        file_exists $ref || { errmsg "file $ref does not exist" ; return 1; }
        echo $ref
    else
        errmsg "-r option should be given"
        echo $ref
        return 1
    fi
}

########
filter_bam_stats()
{
    ${AWK} '{if($3>0 || $4>0) printf" %s",$1}'
}

########
get_bam_contig_list()
{
    local bam=$1

    conda activate samtools 2>&1 || exit 1
    samtools idxstats $bam | filter_bam_stats
    conda deactivate
}

########
get_normal_bam_filename()
{
    local cmdline=$1
    local given=0
    local normalbam
    normalbam=`read_opt_value_from_line "$cmdline" "-n"` && given=1
    if [ $given -eq 1 ]; then
        # -n option was given
        file_exists $normalbam || { errmsg "file $normalbam does not exist" ; return 1; }
        echo $normalbam
    else
        # Check -extn option
        check_opt_given "$cmdline" "-extn" || { errmsg "-n or -extn option should be given" ; return 1; }
        local abs_bamdir=`get_absolute_shdirname "data"`
        normalbam=${abs_bamdir}/normal.bam
        echo $normalbam
    fi
}

########
get_tumor_bam_filename()
{
    local cmdline=$1
    local given=0
    local tumorbam
    tumorbam=`read_opt_value_from_line "$cmdline" "-t"` && given=1
    if [ $given -eq 1 ]; then
        # -t option was given
        file_exists $tumorbam || { errmsg "file $tumorbam does not exist" ; return 1; }
        echo $tumorbam
    else
        # Check -extt option
        check_opt_given "$cmdline" "-extt" || { errmsg "-t or -extt option should be given" ; return 1; }
        local abs_bamdir=`get_absolute_shdirname "data"`
        tumorbam=${abs_bamdir}/tumor.bam
        echo $tumorbam
    fi
}

########
manta_germline_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call for Manta and Strelka (optional)"
    explain_cmdline_opt "-cr" "<string>" "$description"
}

########
manta_germline_define_opts()
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

    # -cr option
    define_cmdline_nonmandatory_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
get_callreg_opt()
{
    local callregf=$1
    local basecallregf=`$BASENAME ${callregf}`

    if [ "${basecallregf}" = ${NOFILE} -o "${callregf}" = "" ]; then
        echo ""
    else
        echo "--callRegions ${callregf}"
    fi
}

########
manta_germline()
{
    display_begin_step_message

    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local callregf=`read_opt_value_from_line "$*" "-cr"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate manta 2>&1 || exit 1

    # Configure Manta
    logmsg "* Executing configManta.py..."
    configManta.py --bam ${normalbam} --referenceFasta ${ref} ${call_reg_opt} --runDir ${step_outd} 2>&1 || exit 1

    # Execute Manta
    logmsg "* Executing runWorkflow.py..."
    ${step_outd}/runWorkflow.py -m local -j ${cpus} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
manta_germline_conda_envs()
{
    define_conda_env manta manta.yml
}

########
cnvkit_explain_cmdline_opts()
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
}

########
cnvkit_define_opts()
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

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
cnvkit()
{
    display_begin_step_message

    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    
    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate cnvkit 2>&1 || exit 1

    # Run cnvkit
    logmsg "* Executing cnvkit.py..."
    cnvkit.py batch ${tumorbam} -n ${normalbam} -m wgs -f ${ref}  -d ${step_outd} -p ${cpus} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
cnvkit_conda_envs()
{
    define_conda_env cnvkit cnvkit.yml
}

########
manta_somatic_explain_cmdline_opts()
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

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call for Manta and Strelka (optional)"
    explain_cmdline_opt "-cr" "<string>" "$description"
}

########
manta_somatic_define_opts()
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

    # -cr option
    define_cmdline_nonmandatory_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
manta_somatic()
{
    display_begin_step_message

    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local callregf=`read_opt_value_from_line "$*" "-cr"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate manta 2>&1 || exit 1
    
    # Configure Manta
    logmsg "* Executing configManta.py..."
    configManta.py --normalBam ${normalbam} --tumorBam ${tumorbam} --referenceFasta ${ref} ${call_reg_opt} --runDir ${step_outd} 2>&1 || exit 1

    # Execute Manta
    logmsg "* Executing runWorkflow.py..."
    ${step_outd}/runWorkflow.py -m local -j ${cpus} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
manta_somatic_conda_envs()
{
    define_conda_env manta manta.yml
}

########
strelka_germline_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call for Manta and Strelka (optional)"
    explain_cmdline_opt "-cr" "<string>" "$description"
}

########
strelka_germline_define_opts()
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

    # -cr option
    define_cmdline_nonmandatory_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
strelka_germline()
{
    display_begin_step_message

    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local callregf=`read_opt_value_from_line "$*" "-cr"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate strelka 2>&1 || exit 1

    # Configure Strelka
    logmsg "* Executing configureStrelkaGermlineWorkflow.py..."
    configureStrelkaGermlineWorkflow.py --bam ${normalbam} --referenceFasta ${ref} ${call_reg_opt} --runDir ${step_outd} 2>&1 || exit 1

    # Execute Strelka
    logmsg "* Executing runWorkflow.py..."
    ${step_outd}/runWorkflow.py -m local -j ${cpus} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
strelka_germline_conda_envs()
{
    define_conda_env strelka strelka.yml
}

########
strelka_somatic_explain_cmdline_opts()
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

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call for Manta and Strelka (optional)"
    explain_cmdline_opt "-cr" "<string>" "$description"
}

########
strelka_somatic_define_opts()
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

    # -manta-outd option
    local manta_dep=`find_dependency_for_step "${stepspec}" manta_somatic`
    if [ ${manta_dep} != ${DEP_NOT_FOUND} ]; then
        local manta_outd=`get_outd_for_dep "${manta_dep}"`
        define_opt "-manta-outd" ${manta_outd} optlist || exit 1
    fi
    
    # -cr option
    define_cmdline_nonmandatory_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
get_indel_cand_opt()
{
    local manta_outd=$1

    if [ -z "${manta_outd}" ]; then
        echo ""
    else
        manta_indel_file="${manta_outd}/results/variants/candidateSmallIndels.vcf.gz"
        if [ -f ${manta_indel_file} ]; then
            echo "--indelCandidates ${manta_indel_file}"
        else
            echo "WARNING: Manta indel file for Strelka not found! (${manta_indel_file})" >&2
            echo ""
        fi
    fi
}

########
strelka_somatic()
{
    display_begin_step_message

    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local manta_outd=`read_opt_value_from_line "$*" "-manta-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local callregf=`read_opt_value_from_line "$*" "-cr"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Define --indelCandidates option if output from Manta is available
    indel_cand_opt=`get_indel_cand_opt "${manta_outd}"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate strelka 2>&1 || exit 1

    # Configure Strelka
    logmsg "* Executing configureStrelkaSomaticWorkflow.py..."
    configureStrelkaSomaticWorkflow.py --normalBam ${normalbam} --tumorBam ${tumorbam} --referenceFasta ${ref} ${indel_cand_opt} ${call_reg_opt} --runDir ${step_outd} 2>&1 || exit 1

    # Execute Strelka
    logmsg "* Executing runWorkflow.py..."
    ${step_outd}/runWorkflow.py -m local -j ${cpus} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    display_end_step_message
}

########
strelka_somatic_conda_envs()
{
    define_conda_env strelka strelka.yml
}

########
platypus_germline_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"    
}

########
platypus_germline_define_opts()
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

    # Save option list
    save_opt_list optlist    
}

########
platypus_germline_conda()
{
    display_begin_step_message
    
    # Initialize variables
    local ref=$1
    local normalbam=$2
    local step_outd=$3

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate platypus 2>&1 || exit 1

    # Run Platypus
    logmsg "* Executing Platypus.py..."
    Platypus.py callVariants --bamFiles=${normalbam} --refFile=${ref} --output=${step_outd}/output.vcf --verbosity=1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
platypus_germline_local()
{
    display_begin_step_message

    # Initialize variables
    local ref=$1
    local normalbam=$2
    local step_outd=$3

    # Run Platypus
    logmsg "* Executing Platypus.py..."
    python ${PLATYPUS_HOME_DIR}/bin/Platypus.py callVariants --bamFiles=${normalbam} --refFile=${ref} --output=${step_outd}/output.vcf --verbosity=1 2>&1 || exit 1

    display_end_step_message
}

########
platypus_germline()
{
    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`

    if [ -z "${PLATYPUS_HOME_DIR}" ]; then
        platypus_germline_conda ${ref} ${normalbam} ${step_outd}
    else
        platypus_germline_local ${ref} ${normalbam} ${step_outd}
    fi
}

########
platypus_conda_envs()
{
    define_conda_env platypus platypus.yml
}

########
msisensor_explain_cmdline_opts()
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
}

########
msisensor_define_opts()
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

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
msisensor()
{
    display_begin_step_message

    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate msisensor 2>&1 || exit 1

    # Create homopolymer and microsatellites file
    logmsg "* Executing msisensor scan..."
    command msisensor scan -d ${ref} -o ${step_outd}/msisensor.list 2>&1 || exit 1

    # Run MSIsensor analysis
    logmsg "* Executing msisensor msi..."
    command msisensor msi -d ${step_outd}/msisensor.list -n ${normalbam} -t ${tumorbam} -o ${step_outd}/output -l 1 -q 1 -b ${cpus} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
    
    display_end_step_message
}

########
msisensor_conda_envs()
{
    define_conda_env msisensor msisensor.yml
}

########
wisecondorx_explain_cmdline_opts()
{
    # -wcr option
    description="Reference file in npz format for WisecondorX"
    explain_cmdline_opt "-wcr" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    
}

########
wisecondorx_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -wcr option
    define_cmdline_infile_opt "$cmdline" "-wcr" optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
wisecondorx()
{
    display_begin_step_message

    # Initialize variables
    local wcref=`read_opt_value_from_line "$*" "-wcr"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    
    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate wisecondorx > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Convert tumor bam file into npz
    logmsg "* Executing WisecondorX convert..."
    local BINSIZE=5000
    WisecondorX convert ${tumorbam} ${step_outd}/tumor.npz --binsize $BINSIZE 2>&1 || exit 1
    
    # Use WisecondorX for prediction
    logmsg "* Executing WisecondorX predict..."
    WisecondorX predict ${step_outd}/tumor.npz ${wcref} ${step_outd}/out 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
wisecondorx_conda_envs()
{
    define_conda_env wisecondorx wisecondorx.yml
}

########
snp_pileup_plus_facets_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"        

    # -sv option
    description="SNP vcf file required by Facets"
    explain_cmdline_opt "-sv" "<string>" "$description"        
}

########
snp_pileup_plus_facets_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
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
snp_pileup_plus_facets()
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
    
    # Execute facets
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing facets..."
    Rscript ${biopanpipe_bindir}/run_facets -c ${step_outd}/snp-pileup-counts.csv -o ${step_outd} 2>&1 || exit 1

    # Deactivate conda environment if needed
    if [ -z "${FACETS_HOME_DIR}" ]; then
        logmsg "* Deactivating conda environment..."
        conda deactivate 2>&1
    fi

    display_end_step_message
}

########
snp_pileup_plus_facets_conda_envs()
{
    define_conda_env facets facets.yml
}

########
mpileup_plus_sequenza_explain_cmdline_opts()
{
    # -gcc option
    description="GC content wiggle file for sequenza (required)"
    explain_cmdline_opt "-gcc" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    
}

########
mpileup_plus_sequenza_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -gcc option
    define_cmdline_infile_opt "$cmdline" "-gcc" optlist || exit 1

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
mpileup_plus_sequenza()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local gccont=`read_opt_value_from_line "$*" "-gcc"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || exit 1

    # Generate pileup files
    logmsg "* Generating pileup files..."
    samtools mpileup -f $ref $normalbam | ${GZIP} > ${step_outd}/normal.pileup.gz ; pipe_fail || exit 1
    samtools mpileup -f $ref $tumorbam | ${GZIP} > ${step_outd}/tumor.pileup.gz ; pipe_fail || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || exit 1
    
    # Generate seqz file
    logmsg "* Generating seqz file..."
    sequenza-utils bam2seqz --pileup -gc ${gccont} -n ${step_outd}/normal.pileup.gz -t ${step_outd}/tumor.pileup.gz | ${GZIP} > ${step_outd}/seqz.gz ; pipe_fail || exit 1

    # Execute sequenza
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing sequenza..."
    Rscript ${biopanpipe_bindir}/run_sequenza -s ${step_outd}/seqz.gz -o ${step_outd} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
mpileup_plus_sequenza_conda_envs()
{
    define_conda_env samtools samtools.yml
    define_conda_env sequenza sequenza.yml
}

########
sequenza_explain_cmdline_opts()
{
    # -gcc option
    description="GC content wiggle file for sequenza (required)"
    explain_cmdline_opt "-gcc" "<string>" "$description"
}

########
sequenza_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -gcc option
    define_cmdline_infile_opt "$cmdline" "-gcc" optlist || exit 1

    # Get normal pileup file
    npileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" sambamba_mpileup_norm_bam` || { errmsg "Error: dependency sambamba_mpileup_norm_bam not defined for sequenza"; exit 1; }
    npileup=${npileupdir}/normal.pileup.gz
    define_opt "-npileup" ${npileup} optlist || exit 1

    # Get tumor pileup file
    tpileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" sambamba_mpileup_tum_bam` || { errmsg "Error: dependency sambamba_mpileup_tum_bam not defined for sequenza"; exit 1; }
    tpileup=${tpileupdir}/tumor.pileup.gz
    define_opt "-tpileup" ${tpileup} optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
sequenza()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local gccont=`read_opt_value_from_line "$*" "-gcc"`
    local npileup=`read_opt_value_from_line "$*" "-npileup"`
    local tpileup=`read_opt_value_from_line "$*" "-tpileup"`

    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || exit 1
    
    # Generate seqz file
    logmsg "* Generating seqz file..."
    sequenza-utils bam2seqz --pileup -gc ${gccont} -n ${npileup} -t ${tpileup} | ${GZIP} > ${step_outd}/seqz.gz ; pipe_fail || exit 1

    # Execute sequenza
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing sequenza..."
    Rscript ${biopanpipe_bindir}/run_sequenza -s ${step_outd}/seqz.gz -o ${step_outd} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
sequenza_conda_envs()
{
    define_conda_env sequenza sequenza.yml
}

########
parallel_bam2seqz_explain_cmdline_opts()
{
    # -gcc option
    description="GC content wiggle file for bam2seqz (required)"
    explain_cmdline_opt "-gcc" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"   
}

########
parallel_bam2seqz_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -gcc option
    define_cmdline_infile_opt "$cmdline" "-gcc" optlist || exit 1

    # Get normal pileup directory
    npileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" parallel_sambamba_mpileup_norm_bam` || { errmsg "Error: dependency parallel_sambamba_mpileup_norm_bam not defined for parallel_bam2seqz"; exit 1; }

    # Get tumor pileup directory
    tpileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" parallel_sambamba_mpileup_tum_bam` || { errmsg "Error: dependency parallel_sambamba_mpileup_tum_bam not defined for parallel_bam2seqz"; exit 1; }

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"`

    # Generate option lists for each contig
    local contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        npileup=${npileupdir}/normal_${contig}.pileup.gz
        define_opt "-npileup" ${npileup} specific_optlist || exit 1
        tpileup=${tpileupdir}/tumor_${contig}.pileup.gz
        define_opt "-tpileup" ${tpileup} specific_optlist || exit 1
        define_opt "-contig" $contig specific_optlist || exit 1
        
        save_opt_list specific_optlist
    done
}

########
parallel_bam2seqz()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local gccont=`read_opt_value_from_line "$*" "-gcc"`
    local npileup=`read_opt_value_from_line "$*" "-npileup"`
    local tpileup=`read_opt_value_from_line "$*" "-tpileup"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || exit 1
    
    # Generate seqz file
    logmsg "* Generating seqz file (contig $contig)..."
    sequenza-utils bam2seqz --pileup -gc ${gccont} -n ${npileup} -t ${tpileup} | ${GZIP} > ${step_outd}/${contig}_seqz.gz ; pipe_fail || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
parallel_bam2seqz_conda_envs()
{
    define_conda_env sequenza sequenza.yml
}

########
seqzmerge_plus_sequenza_explain_cmdline_opts()
{
    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"   
}

########
seqzmerge_plus_sequenza_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -gcc option
    define_cmdline_infile_opt "$cmdline" "-gcc" optlist || exit 1

    # Get seqz directory
    seqzdir=`get_outd_for_dep_given_stepspec "${stepspec}" parallel_bam2seqz` || { errmsg "Error: dependency parallel_bam2seqz not defined for seqzmerge_plus_sequenza"; exit 1; }
    define_opt "-seqzdir" ${seqzdir} optlist || exit 1

    # -lc option
    define_cmdline_infile_opt "$cmdline" "-lc" optlist || exit 1

    save_opt_list optlist
}


########
seqzmerge()
{
    local clist=$1
    local seqzdir=$2

    local contigs=`get_contig_list_from_file $clist` || exit 1
    local filenames=""
    local contig
    for contig in ${contigs}; do
        local seqzfname=${seqzdir}/${contig}_seqz.gz
        if [ "$filenames" = "" ]; then
            filenames=${seqzfname}
        else
            filenames="${filenames} ${seqzfname}"
        fi
    done

    ${ZCAT} ${filenames} | $AWK '{if (NR!=1 && $1 != "chromosome") {print $0}}' | ${GZIP}
}
 
########
seqzmerge_plus_sequenza()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local gccont=`read_opt_value_from_line "$*" "-gcc"`
    local seqzdir=`read_opt_value_from_line "$*" "-seqzdir"`
    local clist=`read_opt_value_from_line "$*" "-lc"`

    # Merge seqz files
    logmsg "* Merging seqz files..."
    seqzmerge ${clist} ${seqzdir}  > ${step_outd}/merged_seqz.gz

    # Activate conda environment
    logmsg "* Activating conda environment (tabix)..."
    conda activate tabix 2>&1 || exit 1

    logmsg "* Applying tabix over merged seqz file..."
    tabix -f -s 1 -b 2 -e 2 -S 1 ${step_outd}/merged_seqz.gz
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
        
    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || exit 1
    
    # Execute sequenza
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing sequenza..."
    Rscript ${biopanpipe_bindir}/run_sequenza -s ${step_outd}/merged_seqz.gz -o ${step_outd} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
seqzmerge_plus_sequenza_conda_envs()
{
    define_conda_env sequenza sequenza.yml
    define_conda_env sequenza tabix.yml
}

########
lumpy_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    
}

########
lumpy_define_opts()
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

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
lumpy()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate lumpy 2>&1 || exit 1

    logmsg "* Executing lumpyexpress..."
    lumpyexpress -B ${tumorbam},${normalbam} -o ${step_outd}/out.vcf || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
manta_somatic_conda_envs()
{
    define_conda_env lumpy lumpy.yml
}

########
parallel_exclude_plus_lumpy_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"   
}

########
parallel_exclude_plus_lumpy_define_opts()
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

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"`

    # Generate option lists for each contig
    local contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
get_contig_names_from_bam()
{
    local bam=$1

    conda activate samtools 2>&1 || return 1
    samtools idxstats $bam | $AWK '{printf"%s\n",$1}' ; pipe_fail || return 1
    conda deactivate
}

########
gen_exclusion_bed_given_bam()
{
    local bam=$1
    local contig=$2

    get_contig_names_from_bam $bam | $AWK -v contig=$contig '{if($1!=contig){printf"%s\n",$1}}' ; pipe_fail || return 1
}

########
parallel_exclude_plus_lumpy()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Generate exclusion bed file
    gen_exclusion_bed_given_bam ${normalbam} ${contig} > ${step_outd}/${contig}.bed
    
    # Activate conda environment
    logmsg "* Activating conda environment... (lumpy)"
    conda activate lumpy 2>&1 || exit 1
    
    logmsg "* Executing lumpyexpress (contig $contig)..."
    lumpyexpress -B ${tumorbam},${normalbam} -T ${step_outd}/tmp_${contig} -x ${step_outd}/${contig}.bed -o ${step_outd}/out${contig}.vcf || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
    
    display_end_step_message
}

########
parallel_exclude_plus_lumpy_conda_envs()
{
    define_conda_env samtools samtools.yml
    define_conda_env lumpy lumpy.yml
}

########
parallel_split_plus_lumpy_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"   
}

########
parallel_split_plus_lumpy_define_opts()
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

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"`

    # Generate option lists for each contig
    local contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
check_contig_does_not_exist_given_log_file()
{
    local logfile=$1

    if $GREP "does not exist" ${logfile} >/dev/null 2>&1; then
        return 0
    else
        return 1
    fi
}

########
filter_bam_contig()
{
    local inbam=$1
    local contig=$2
    local outbam=$3    
    local error=0
    
    sambamba view -h -f bam $inbam $contig > ${outbam} 2> ${outbam}.log || error=1

    if [ $error -eq 1 ]; then
        if check_contig_does_not_exist_given_log_file ${outbam}.log; then
            errmsg "Warning: contig ${contig} does not exist in ${inbam} file (see ${outbam}.log)"
            return 0
        else
            errmsg "Error while filtering ${contig} in ${inbam} file (see ${outbam}.log)"
            return 1
        fi
    else
        return 0
    fi
}

########
parallel_split_plus_lumpy()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Extract contigs
    logmsg "* Extracting contigs (contig $contig)..."
    normalcont=${step_outd}/normal_${contig}.bam
    filter_bam_contig $normalbam $contig $normalcont || exit 1
    tumorcont=${step_outd}/tumor_${contig}.bam
    filter_bam_contig $tumorbam $contig $tumorcont || exit 1

    # Index contigs
    logmsg "* Indexing contigs..."
    sambamba index ${normalcont} || exit 1
    sambamba index ${tumorcont} || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
    
    # Activate conda environment
    logmsg "* Activating conda environment... (lumpy)"
    conda activate lumpy 2>&1 || exit 1
    
    logmsg "* Executing lumpyexpress (contig $contig)..."
    lumpyexpress -B ${tumorcont},${normalcont} -T ${step_outd}/tmp_${contig} -o ${step_outd}/out${contig}.vcf || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Delete extracted contigs and related files
    rm ${normalcont}* ${tumorcont}*
    
    display_end_step_message
}

########
parallel_split_plus_lumpy_clean()
{
    logmsg "Cleaning directory..."

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Delete extracted contigs and related files
    normalcont=${step_outd}/normal_${contig}.bam
    tumorcont=${step_outd}/tumor_${contig}.bam
    rm -f ${normalcont}* ${tumorcont}*

    logmsg "Cleaning finished"
}

########
parallel_split_plus_lumpy_conda_envs()
{
    define_conda_env sambamba sambamba.yml
    define_conda_env lumpy lumpy.yml
}

########
parallel_lumpy_explain_cmdline_opts()
{
    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"   
}

########
parallel_lumpy_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # Get bam directory
    local abs_bamdir=`get_absolute_shdirname "data"`

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"`

    # Generate option lists for each contig
    local contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        normalbam=${abs_bamdir}/normal_${contig}.bam
        define_opt "-normalbam" ${normalbam} specific_optlist || exit 1
        tumorbam=${abs_bamdir}/tumor_${contig}.bam
        define_opt "-tumorbam" ${tumorbam} specific_optlist || exit 1
        define_opt "-contig" ${contig} specific_optlist || exit 1
        
        save_opt_list specific_optlist
    done
}

########
parallel_lumpy()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    
    # Activate conda environment
    logmsg "* Activating conda environment... (lumpy)"
    conda activate lumpy 2>&1 || exit 1
    
    logmsg "* Executing lumpyexpress (contig $contig)..."
    lumpyexpress -B ${tumorbam},${normalbam} -T ${step_outd}/tmp_${contig} -o ${step_outd}/out${contig}.vcf || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
    
    display_end_step_message
}

########
parallel_lumpy_conda_envs()
{
    define_conda_env lumpy lumpy.yml
}

########
delly_explain_cmdline_opts()
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

    # -dx option
    description="File with regions to exclude in bed format for Delly"
    explain_cmdline_opt "-dx" "<string>" "$description"    
}

########
delly_define_opts()
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

    # -dx option
    define_cmdline_infile_opt "$cmdline" "-dx" optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
delly()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local exclude=`read_opt_value_from_line "$*" "-dx"`

    # Activate conda environment
    logmsg "* Activating conda environment... (delly)"
    conda activate delly 2>&1 || exit 1

    logmsg "* Executing delly..."
    # "command" built-in is used here to execute the "delly" program
    # instead of the "delly" function
    command delly call -g ${ref} -x ${exclude} -o ${step_outd}/out.bcf ${tumorbam} ${normalbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment... (bcftools)"
    conda activate bcftools 2>&1 || exit 1

    # Convert bcf output to vcf
    logmsg "* Converting bcf output into vcf... (bcftools)"
    bcftools view ${step_outd}/out.bcf > ${step_outd}/out.vcf

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
parallel_lumpy_conda_envs()
{
    define_conda_env bcftools bcftools.yml
    define_conda_env delly delly.yml
}

########
parallel_split_plus_delly_explain_cmdline_opts()
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

    # -dx option
    description="File with regions to exclude in bed format for Delly"
    explain_cmdline_opt "-dx" "<string>" "$description"    

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"   
}

########
parallel_split_plus_delly_define_opts()
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

    # -dx option
    define_cmdline_infile_opt "$cmdline" "-dx" optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"`

    # Generate option lists for each contig
    local contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_split_plus_delly()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local exclude=`read_opt_value_from_line "$*" "-dx"`

    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Extract contigs
    logmsg "* Extracting contigs..."
    normalcont=${step_outd}/normal_${contig}.bam
    filter_bam_contig $normalbam $contig $normalcont || exit 1
    tumorcont=${step_outd}/tumor_${contig}.bam
    filter_bam_contig $tumorbam $contig $tumorcont || exit 1

    # Index contigs
    logmsg "* Indexing contigs..."
    sambamba index ${normalcont} || exit 1
    sambamba index ${tumorcont} || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
    
    # Activate conda environment
    logmsg "* Activating conda environment... (delly)"
    conda activate delly 2>&1 || exit 1
    
    logmsg "* Executing delly (contig $contig)..."
    # "command" built-in is used here to execute the "delly" program
    # instead of the "delly" function
    command delly call -g $ref -x ${exclude} -o ${step_outd}/out${contig}.bcf ${tumorcont} ${normalcont} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment... (bcftools)"
    conda activate bcftools 2>&1 || exit 1

    # Convert bcf output to vcf
    logmsg "* Converting bcf output into vcf... (bcftools)"
    bcftools view ${step_outd}/out${contig}.bcf > ${step_outd}/out${contig}.vcf

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Delete extracted contigs and related files
    rm ${normalcont}* ${tumorcont}*
    
    display_end_step_message
}

########
parallel_split_plus_delly_clean()
{
    logmsg "Cleaning directory..."

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local exclude=`read_opt_value_from_line "$*" "-dx"`

    # Delete extracted contigs and related files
    normalcont=${step_outd}/normal_${contig}.bam
    tumorcont=${step_outd}/tumor_${contig}.bam
    rm ${normalcont}* ${tumorcont}*

    logmsg "Cleaning finished"
}

########
parallel_split_plus_delly_conda_envs()
{
    define_conda_env sambamba sambamba.yml
    define_conda_env bcftools bcftools.yml
    define_conda_env delly delly.yml
}

########
parallel_delly_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -dx option
    description="File with regions to exclude in bed format for Delly"
    explain_cmdline_opt "-dx" "<string>" "$description"    

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"   
}

########
parallel_delly_define_opts()
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

    # Get bam directory
    local abs_bamdir=`get_absolute_shdirname "data"`

    # -dx option
    define_cmdline_infile_opt "$cmdline" "-dx" optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"`

    # Generate option lists for each contig
    local contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        normalbam=${abs_bamdir}/normal_${contig}.bam
        define_opt "-normalbam" ${normalbam} specific_optlist || exit 1
        tumorbam=${abs_bamdir}/tumor_${contig}.bam
        define_opt "-tumorbam" ${tumorbam} specific_optlist || exit 1
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_delly()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local exclude=`read_opt_value_from_line "$*" "-dx"`
    
    # Activate conda environment
    logmsg "* Activating conda environment... (delly)"
    conda activate delly 2>&1 || exit 1
    
    logmsg "* Executing delly (contig $contig)..."
    # "command" built-in is used here to execute the "delly" program
    # instead of the "delly" function
    command delly call -g $ref -x ${exclude} -o ${step_outd}/out${contig}.bcf ${tumorbam} ${normalbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment... (bcftools)"
    conda activate bcftools 2>&1 || exit 1

    # Convert bcf output to vcf
    logmsg "* Converting bcf output into vcf... (bcftools)"
    bcftools view ${step_outd}/out${contig}.bcf > ${step_outd}/out${contig}.vcf

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
    
    display_end_step_message
}

########
parallel_delly_conda_envs()
{
    define_conda_env bcftools bcftools.yml
    define_conda_env delly delly.yml
}

########
parallel_svtyper_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"   
}

########
get_vcfdir_for_svtyper()
{
    local stepspec=$1

    # Check dependency with parallel_lumpy
    local parallel_lumpy_dep=`find_dependency_for_step "${stepspec}" parallel_lumpy`
    if [ ${parallel_lumpy_dep} != ${DEP_NOT_FOUND} ]; then
        local vcfdir=`get_outd_for_dep "${parallel_lumpy_dep}"`
        echo $vcfdir
        return 0
    fi

    # Check dependency with parallel_split_plus_lumpy
    local parallel_split_plus_lumpy_dep=`find_dependency_for_step "${stepspec}" parallel_split_plus_lumpy`
    if [ ${parallel_split_plus_lumpy_dep} != ${DEP_NOT_FOUND} ]; then
        local vcfdir=`get_outd_for_dep "${parallel_split_plus_lumpy_dep}"`
        echo $vcfdir
        return 0
    fi

    return 1
}

########
parallel_svtyper_define_opts()
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

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"`

    # Determine vcf directory
    vcfdir=`get_vcfdir_for_svtyper "${stepspec}"` || { errmsg "Error: vcf directory for svtyper could not be determined"; exit 1; }
    
    # Generate option lists for each contig
    local contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        define_opt "-contig" $contig specific_optlist || exit 1
        vcf=${vcfdir}/out${contig}.vcf
        define_opt "-vcf" $vcf specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_svtyper()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local vcf=`read_opt_value_from_line "$*" "-vcf"`

    # Activate conda environment
    logmsg "* Activating conda environment (svtyper)..."
    conda activate svtyper 2>&1 || exit 1

    # Execute svtyper
    logmsg "* Executing svtyper (contig $contig)..."
    svtyper -i ${vcf} -B ${tumorbam},${normalbam} > ${step_outd}/out${contig}.vcf || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
    
    display_end_step_message
}

########
parallel_svtyper_conda_envs()
{
    define_conda_env svtyper svtyper.yml
}

########
download_ega_norm_bam_explain_cmdline_opts()
{
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"

    # -extn option
    description="External database id of normal bam file to download (required)"
    explain_cmdline_opt "-extn" "<string>" "$description"

    # -egastr option
    description="Number of streams used by the EGA download client (${DEFAULT_NUMBER_OF_EGA_DOWNLOAD_STREAMS} by default)"
    explain_cmdline_opt "-egastr" "<int>" "$description"

    # -egacred option
    description="File with EGA download client credentials (required)"
    explain_cmdline_opt "-egacred" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}


########
download_ega_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -egastr option
    define_cmdline_opt "$cmdline" "-egastr" optlist || exit 1

    # -egacred option
    define_cmdline_opt "$cmdline" "-egacred" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local abs_bamdir=`get_absolute_shdirname "data"`
    local normalbam=${abs_bamdir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
ega_download_retry()
{
    # Initialize variables
    local egastr=$1
    local egacred=$2
    local egaid=$3
    local outf=$4
    local download_tries=$5
    local step_outd=`${DIRNAME} ${outf}`
    
    # Start download with multiple tries
    local ntry=1
    while [ ${ntry} -le ${download_tries} ]; do
        logmsg "Starting download try number ${ntry}..."

        # Remove previously downloaded file (if any)
        if [ -f ${outf} ]; then
            rm ${outf}
        fi

        # Download file
        pyega3 -c ${egastr} -cf ${egacred} fetch ${egaid} ${outf} 2>&1
        
        # Check if download was successful
        if [ $? -eq 0 -a -f ${outf} ]; then
            return 0
        fi

        ntry=`expr ${ntry} + 1`
    done

    logmsg "All download attempts failed!"

    return 1
}

########
download_ega_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local egaid_normalbam=`read_opt_value_from_line "$*" "-extn"`
    local egastr=`read_opt_value_from_line "$*" "-egastr"`
    local egacred=`read_opt_value_from_line "$*" "-egacred"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate pyega3 2>&1 || exit 1

    # Download file (with multiple tries)
    ega_download_retry ${egastr} ${egacred} ${egaid_normalbam} ${step_outd}/normal.bam ${download_tries} || exit 1

    # Move file
    mv ${step_outd}/normal.bam ${normalbam} || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Create file indicating that execution was finished
    touch ${step_outd}/finished

    display_end_step_message
}

########
download_ega_norm_bam_conda_envs()
{
    define_conda_env pyega3 pyega3.yml
}

########
download_ega_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of tumor bam file to download (required)"
    explain_cmdline_opt "-extt" "<string>" "$description"

    # -egastr option
    description="Number of streams used by the EGA download client (${DEFAULT_NUMBER_OF_EGA_DOWNLOAD_STREAMS} by default)"
    explain_cmdline_opt "-egastr" "<int>" "$description"

    # -egacred option
    description="File with EGA download client credentials (required)"
    explain_cmdline_opt "-egacred" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_ega_tum_bam_define_opts()
{
    echo "$tumorbam ${extid_tumorbam} $egastr $egacred ${download_tries} ${step_outd}"
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -egastr option
    define_cmdline_opt "$cmdline" "-egastr" optlist || exit 1

    # -egacred option
    define_cmdline_opt "$cmdline" "-egacred" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local abs_bamdir=`get_absolute_shdirname "data"`
    local tumorbam=${abs_bamdir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_ega_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local egaid_normalbam=`read_opt_value_from_line "$*" "-extn"`
    local egastr=`read_opt_value_from_line "$*" "-egastr"`
    local egacred=`read_opt_value_from_line "$*" "-egacred"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate pyega3 > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Download file (with multiple tries)
    ega_download_retry ${egastr} ${egacred} ${egaid_tumorbam} ${step_outd}/tumor.bam ${download_tries} || exit 1

    # Move file
    mv ${step_outd}/tumor.bam ${tumorbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    display_end_step_message
}

########
download_ega_tum_bam_conda_envs()
{
    define_conda_env pyega3 pyega3.yml
}

########
download_aws_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="External database id of normal bam file to download (required)"
    explain_cmdline_opt "-extn" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_aws_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local abs_bamdir=`get_absolute_shdirname "data"`
    local normalbam=${abs_bamdir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
find_bam_filename()
{
    local step_outd=$1
    local result=""
    
    for f in ${step_outd}/*.bam; do
        if [ -f $f ]; then
            result=$f
        fi
    done

    echo ${result}
}

########
download_aws_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local icgcid_normalbam=`read_opt_value_from_line "$*" "-extn"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing icgc-storage-client..."
    ${ICGCSTOR_HOME_DIR}/bin/icgc-storage-client --profile aws download --object-id ${icgcid_normalbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1

    display_end_step_message
}

########
download_aws_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of tumor bam file to download (required)"
    explain_cmdline_opt "-extt" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_aws_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -tumorbam option
    local abs_bamdir=`get_absolute_shdirname "data"`
    local tumorbam=${abs_bamdir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_aws_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local icgcid_tumorbam=`read_opt_value_from_line "$*" "-extt"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing icgc-storage-client..."
    ${ICGCSTOR_HOME_DIR}/bin/icgc-storage-client --profile aws download --object-id ${icgcid_tumorbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1

    display_end_step_message
}

########
download_collab_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="External database id of normal bam file to download (required)"
    explain_cmdline_opt "-extn" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_collab_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local abs_bamdir=`get_absolute_shdirname "data"`
    local normalbam=${abs_bamdir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_collab_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local icgcid_normalbam=`read_opt_value_from_line "$*" "-extn"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing icgc-storage-client..."
    ${ICGCSTOR_HOME_DIR}/bin/icgc-storage-client --profile collab download --object-id ${icgcid_normalbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1

    display_end_step_message
}

########
download_collab_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of tumor bam file to download (required)"
    explain_cmdline_opt "-extt" "<string>" "$description"

    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"    
}

########
download_collab_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -tumorbam option
    local abs_bamdir=`get_absolute_shdirname "data"`
    local tumorbam=${abs_bamdir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_collab_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local icgcid_tumorbam=`read_opt_value_from_line "$*" "-extt"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing icgc-storage-client..."
    ${ICGCSTOR_HOME_DIR}/bin/icgc-storage-client --profile collab download --object-id ${icgcid_tumorbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1

    display_end_step_message
}

########
download_ega_asp_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="External database id of normal bam file to download (required)"
    explain_cmdline_opt "-extn" "<string>" "$description"

    # -asperausr option
    description="Username for Aspera server (required)"
    explain_cmdline_opt "-asperausr" "<string>" "$description"

    # -asperapwd option
    description="Password for Aspera server (required)"
    explain_cmdline_opt "-asperapwd" "<string>" "$description"

    # -asperaserv option
    description="Name of Aspera server (required)"
    explain_cmdline_opt "-asperaserv" "<string>" "$description"

    # -egadecrpwd option
    description="File with EGA decryptor password (required)"
    explain_cmdline_opt "-egadecrpwd" "<string>" "$description"
    
    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_ega_asp_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -asperausr option
    define_cmdline_opt "$cmdline" "-asperausr" optlist || exit 1

    # -asperapwd option
    define_cmdline_opt "$cmdline" "-asperapwd" optlist || exit 1

    # -asperaserv option
    define_cmdline_opt "$cmdline" "-asperaserv" optlist || exit 1

    # -egadecrpwd option
    define_cmdline_opt "$cmdline" "-egadecrpwd" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local abs_bamdir=`get_absolute_shdirname "data"`
    local normalbam=${abs_bamdir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
download_ega_asp_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local normalbam_file=`read_opt_value_from_line "$*" "-extn"`
    local aspera_user=`read_opt_value_from_line "$*" "-asperausr"`
    local aspera_passwd=`read_opt_value_from_line "$*" "-asperapwd"`
    local aspera_server=`read_opt_value_from_line "$*" "-asperaserv"`
    local egadecrypt_pwd=`read_opt_value_from_line "$*" "-egadecrpwd"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local max_trans_rate=${DEFAULT_ASP_MAX_TRANS_RATE}
    
    # Download file
    logmsg "* Executing ascp (${normalbam_file})..."
    ASPERA_SCP_PASS=${aspera_passwd} ${ASPERA_HOME_DIR}/bin/ascp --ignore-host-key -QTl ${max_trans_rate} ${aspera_user}@${aspera_server}:${normalbam_file} ${step_outd}/normal.bam.crypt 2>&1 || exit 1

    # Decrypt file
    logmsg "* Executing decryptor.jar..."
    $JAVA -jar ${EGADECRYPT_HOME_DIR}/decryptor.jar ${egadecrypt_pwd} ${step_outd}/normal.bam.crypt 2>&1 || exit 1
    
    # Obtain file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1

    # Remove encrypted file
    rm ${step_outd}/normal.bam.crypt || exit 1
    
    display_end_step_message
}

########
download_ega_asp_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of normal bam file to download (required)"
    explain_cmdline_opt "-extt" "<string>" "$description"

    # -asperausr option
    description="Username for Aspera server (required)"
    explain_cmdline_opt "-asperausr" "<string>" "$description"

    # -asperapwd option
    description="Password for Aspera server (required)"
    explain_cmdline_opt "-asperapwd" "<string>" "$description"

    # -asperaserv option
    description="Name of Aspera server (required)"
    explain_cmdline_opt "-asperaserv" "<string>" "$description"

    # -egadecrpwd option
    description="File with EGA decryptor password (required)"
    explain_cmdline_opt "-egadecrpwd" "<string>" "$description"
    
    # -nt option
    description="Number of download tries per file (${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} by default)"
    explain_cmdline_opt "-nt" "<int>" "$description"
}

########
download_ega_asp_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -asperausr option
    define_cmdline_opt "$cmdline" "-asperausr" optlist || exit 1

    # -asperapwd option
    define_cmdline_opt "$cmdline" "-asperapwd" optlist || exit 1

    # -asperaserv option
    define_cmdline_opt "$cmdline" "-asperaserv" optlist || exit 1

    # -egadecrpwd option
    define_cmdline_infile_opt "$cmdline" "-egadecrpwd" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -tumorbam option
    local abs_bamdir=`get_absolute_shdirname "data"`
    local tumorbam=${abs_bamdir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
download_ega_asp_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local tumorbam_file=`read_opt_value_from_line "$*" "-extt"`
    local aspera_user=`read_opt_value_from_line "$*" "-asperausr"`
    local aspera_passwd=`read_opt_value_from_line "$*" "-asperapwd"`
    local aspera_server=`read_opt_value_from_line "$*" "-asperaserv"`
    local egadecrypt_pwd=`read_opt_value_from_line "$*" "-egadecrpwd"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local max_trans_rate=${DEFAULT_ASP_MAX_TRANS_RATE}

    # Download file
    logmsg "* Executing ascp (${tumorbam_file})..."
    ASPERA_SCP_PASS=${aspera_passwd} ${ASPERA_HOME_DIR}/bin/ascp --ignore-host-key -QTl ${max_trans_rate} ${aspera_user}@${aspera_server}:${tumorbam_file} ${step_outd}/tumor.bam.crypt 2>&1 || exit 1

    # Decrypt file
    logmsg "* Executing decryptor.jar..."
    $JAVA -jar ${EGADECRYPT_HOME_DIR}/decryptor.jar ${egadecrypt_pwd} ${step_outd}/tumor.bam.crypt 2>&1 || exit 1

    # Obtain file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1

    # Remove encrypted file
    rm ${step_outd}/tumor.bam.crypt || exit 1

    display_end_step_message
}

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
    local abs_bamdir=`get_absolute_shdirname "data"`
    local normalbam=${abs_bamdir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
index_norm_bam()
{
    display_begin_step_message

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
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    display_end_step_message
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
    local abs_bamdir=`get_absolute_shdirname "data"`
    local tumorbam=${abs_bamdir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
index_tum_bam()
{
    display_begin_step_message

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

    display_end_step_message
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
    local abs_bamdir=`get_absolute_shdirname "data"`        
    local normalbam=${abs_bamdir}/normal.bam
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
    display_begin_step_message

    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools > ${step_outd}/conda_activate.log 2>&1 || exit 1

    # Verify if bam file is already sorted
    local bam_is_sorted=`samtools view -H ${normalbam} | $GREP SO:coordinate | wc -l` || exit 1
    if [ ${bam_is_sorted} -eq 1 ]; then
        echo "Warning: bam file is already sorted"
    else
        # Execute samtools
        logmsg "* Executing samtools sort..."
        samtools sort -T ${step_outd} -o ${step_outd}/sorted.bam -m 2G -@ ${cpus} ${normalbam} >  ${step_outd}/samtools.log 2>&1 || exit 1
        # NOTE: -m option is used here to increase the maximum memory per
        # thread. One lateral efect of this is that the number of tmp files
        # generated is decreased. This constitutes one possible way to avoid
        # the "Too many open files" error reported by samtools

        # Replace initial bam file by the sorted one
        mv ${step_outd}/sorted.bam ${normalbam} 2> ${step_outd}/mv.log || exit 1
    fi
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    display_end_step_message
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
    local abs_bamdir=`get_absolute_shdirname "data"`
    local tumorbam=${abs_bamdir}/tumor.bam
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
    display_begin_step_message

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
        samtools sort -T ${step_outd} -o ${step_outd}/sorted.bam -m 2G -@ ${cpus} ${tumorbam} >  ${step_outd}/samtools.log 2>&1 || exit 1
        # NOTE: -m option is used here to increase the maximum memory per
        # thread. One lateral efect of this is that the number of tmp files
        # generated is decreased. This constitutes one possible way to avoid
        # the "Too many open files" error reported by samtools

        # Replace initial bam file by the sorted one
        mv ${step_outd}/sorted.bam ${tumorbam} 2> ${step_outd}/mv.log || exit 1
    fi
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    display_end_step_message
}

########
sort_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
filter_norm_bam_contigs_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"
}

########
filter_norm_bam_contigs_define_opts()
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
    local abs_bamdir=`get_absolute_shdirname "data"`
    local normalbam=${abs_bamdir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
get_ref_contigs()
{
    local faifile=$1
    
    ${AWK} '{print $1}' ${faifile}
}

########
remove_line_breaks_from_file()
{
    local file=$1
    echo `cat ${file}`
}

########
filter_norm_bam_contigs()
{
    display_begin_step_message

    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || exit 1

    # Obtain contigs given in reference
    get_ref_contigs ${ref}.fai > ${step_outd}/refcontigs

    # Obtain new header
    
    ## Extract sam header
    samtools view -H ${step_outd}/filtered.bam > ${step_outd}/original_header || exit 1
    
    ## Generate new sam header
    ${biopanpipe_bindir}/get_filtered_sam_header -h ${step_outd}/original_header -l ${step_outd}/refcontigs > ${step_outd}/new_header || exit 1

    # Generate filtered bam
    contigs=`remove_line_breaks_from_file ${step_outd}/refcontigs`
    samtools view ${normalbam} ${contigs} | samtools view -bo ${step_outd}/filtered.bam -t ${step_outd}/new_header -

    # Move bam file
    mv ${step_outd}/filtered.bam ${normalbam}
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
filter_norm_bam_contigs_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
filter_tum_bam_contigs_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"
}

########
filter_tum_bam_contigs_define_opts()
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

    # -tumorbam option
    local abs_bamdir=`get_absolute_shdirname "data"`
    local tumorbam=${abs_bamdir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
filter_tum_bam_contigs()
{
    display_begin_step_message

    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || exit 1

    # Obtain contigs given in reference
    get_ref_contigs ${ref}.fai > ${step_outd}/refcontigs

    # Obtain new header
    
    ## Extract sam header
    samtools view -H ${step_outd}/filtered.bam > ${step_outd}/original_header || exit 1
    
    ## Generate new sam header
    ${biopanpipe_bindir}/get_filtered_sam_header -h ${step_outd}/original_header -l ${step_outd}/refcontigs > ${step_outd}/new_header || exit 1

    # Generate filtered bam
    contigs=`remove_line_breaks_from_file ${step_outd}/refcontigs`
    samtools view ${tumorbam} ${contigs} | samtools view -bo ${step_outd}/filtered.bam -t ${step_outd}/new_header -

    # Move bam file
    mv ${step_outd}/filtered.bam ${tumorbam}

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
filter_tum_bam_contigs_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
sambamba_mpileup_norm_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup (optional)"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || exit 1

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
    local mbpfile=$1

    if [ "${mbpfile}" = ${OPT_NOT_FOUND} ]; then
        echo ""
    else
        echo "-L ${mbpfile}"
    fi
}

########
sambamba_mpileup_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local mbpfile=`read_opt_value_from_line "$*" "-mpb"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Obtain sambamba mpileup -L opt
    local smp_l_opt=`get_sambamba_mpileup_l_opt ${mbpfile}`

    # Generate pileup file
    logmsg "* Generating pileup file..."
    sambamba mpileup -t ${cpus} ${smp_l_opt} --tmpdir ${step_outd} -o ${step_outd}/normal.pileup $normalbam --samtools "-f ${ref}" || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/normal.pileup

    display_end_step_message
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
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup (optional)"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || exit 1

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
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local mbpfile=`read_opt_value_from_line "$*" "-mpb"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Obtain sambamba mpileup -L opt
    local smp_l_opt=`get_sambamba_mpileup_l_opt ${mbpfile}`
    
    # Generate pileup file
    logmsg "* Generating pileup file..."
    sambamba mpileup -t ${cpus} ${smp_l_opt} --tmpdir ${step_outd} -o ${step_outd}/tumor.pileup $tumorbam --samtools "-f ${ref}" || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/tumor.pileup

    display_end_step_message
}

########
sambamba_mpileup_tum_bam_conda_envs()
{
    define_conda_env sambamba sambamba.yml
}

########
parallel_sambamba_mpileup_norm_bam_explain_cmdline_opts()
{
    # -mpb option
    description="BED file for mpileup (optional)"
    explain_cmdline_opt "-mpb" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || exit 1

    # -bamdir option    
    abs_bamdir=`get_absolute_shdirname "data"`

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"`

    # Generate option lists for each contig
    local contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        normalbam=${abs_bamdir}/normal_${contig}.bam
        define_opt "-normalbam" ${normalbam} specific_optlist || exit 1
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_sambamba_mpileup_norm_bam()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local mbpfile=`read_opt_value_from_line "$*" "-mpb"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Obtain sambamba mpileup -L opt
    local smp_l_opt=`get_sambamba_mpileup_l_opt ${mbpfile}`

    # Generate pileup file
    logmsg "* Generating pileup file (contig $contig)..."
    sambamba mpileup -t ${cpus} ${smp_l_opt} --tmpdir ${step_outd} -o ${step_outd}/normal_${contig}.pileup $normalbam --samtools "-f ${ref}" || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/normal_${contig}.pileup

    display_end_step_message
}

########
parallel_sambamba_mpileup_norm_bam_conda_envs()
{
    define_conda_env sambamba sambamba.yml
}

########
parallel_sambamba_mpileup_tum_bam_explain_cmdline_opts()
{
    # -mpb option
    description="BED file for mpileup (optional)"
    explain_cmdline_opt "-mpb" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"
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
    define_cmdline_infile_opt "$cmdline" "-r" optlist || exit 1

    # -bamdir option    
    abs_bamdir=`get_absolute_shdirname "data"`

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"`

    # Generate option lists for each contig
    local contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        tumorbam=${abs_bamdir}/tumor_${contig}.bam
        define_opt "-tumorbam" ${tumorbam} specific_optlist || exit 1
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_sambamba_mpileup_tum_bam()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local mbpfile=`read_opt_value_from_line "$*" "-mpb"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Obtain sambamba mpileup -L opt
    local smp_l_opt=`get_sambamba_mpileup_l_opt ${mbpfile}`

    # Generate pileup file
    logmsg "* Generating pileup file (contig $contig)..."
    sambamba mpileup -t ${cpus} ${smp_l_opt} --tmpdir ${step_outd} -o ${step_outd}/tumor_${contig}.pileup $tumorbam --samtools "-f ${ref}" || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    ${GZIP} ${step_outd}/tumor_${contig}.pileup

    display_end_step_message
}

########
parallel_sambamba_mpileup_tum_bam_conda_envs()
{
    define_conda_env sambamba sambamba.yml
}

########
parallel_split_norm_bam_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"
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

    # -bamdir option    
    abs_bamdir=`get_absolute_shdirname "data"`
    define_opt "-bamdir" ${abs_bamdir} optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"`

    # Generate option lists for each contig
    local contigs=`get_contig_list_from_file $clist` || exit 1
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
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local abs_bamdir=`read_opt_value_from_line "$*" "-bamdir"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    
    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Extract contig
    logmsg "* Extracting contig (contig $contig)..."
    normalcont=${step_outd}/normal_${contig}.bam
    filter_bam_contig $normalbam $contig $normalcont || exit 1

    # Index contig
    logmsg "* Indexing contig..."
    sambamba index ${normalcont} || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Move bam and index file to bamdir
    mv ${normalcont}* ${abs_bamdir} || exit 1

    display_end_step_message
}

########
parallel_split_norm_bam_conda_envs()
{
    define_conda_env sambamba sambamba.yml
}

########
parallel_split_tum_bam_explain_cmdline_opts()
{
    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_opt "-lc" "<string>" "$description"   
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

    # -bamdir option
    abs_bamdir=`get_absolute_shdirname "data"`
    define_opt "-bamdir" ${abs_bamdir} optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"`

    # Generate option lists for each contig
    local contigs=`get_contig_list_from_file $clist` || exit 1
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
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local abs_bamdir=`read_opt_value_from_line "$*" "-bamdir"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Activate conda environment
    logmsg "* Activating conda environment (sambamba)..."
    conda activate sambamba 2>&1 || exit 1

    # Extract contig
    logmsg "* Extracting contig (contig $contig)..."
    tumorcont=${step_outd}/tumor_${contig}.bam
    filter_bam_contig $tumorbam $contig $tumorcont || exit 1

    # Index contig
    logmsg "* Indexing contig..."
    sambamba index ${tumorcont} || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
        
    # Move bam and index file to bamdir
    mv ${tumorcont}* ${abs_bamdir} || exit 1

    display_end_step_message
}

########
parallel_split_tum_bam_conda_envs()
{
    define_conda_env sambamba sambamba.yml
}

########
delete_bam_files_explain_cmdline_opts()
{
    :
}

########
delete_bam_files_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -bamdir option
    abs_bamdir=`get_absolute_shdirname "data"`
    define_opt "-bamdir" ${abs_bamdir} optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
delete_bam_files()
{
    display_begin_step_message

    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local abs_bamdir=`read_opt_value_from_line "$*" "-bamdir"`

    # Delete bam files
    rm -f ${abs_bamdir}/*.bam

    display_end_step_message
}
