# *- bash -*

#############
# CONSTANTS #
#############

DEFAULT_NUMBER_OF_DOWNLOAD_TRIES=5
DEFAULT_NUMBER_OF_EGA_DOWNLOAD_STREAMS=50
DEFAULT_ASP_MAX_TRANS_RATE=100m
DEFAULT_BAMDIR="data"

######################
# BAM ANALYSIS STEPS #
######################

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
        local bamdir_fullname
        bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || { return 1; }
        normalbam=${bamdir_fullname}/normal.bam
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
        local bamdir_fullname
        bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || { return 1; }
        tumorbam=${bamdir_fullname}/tumor.bam
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

    # -callregf option
    description="bgzipped and tabixed bed file to specify regions to call for Manta and Strelka"
    explain_cmdline_opt "-cr" "<string>" "$description"
}

########
manta_germline_define_opts()
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

    # -callregf option
    define_cmdline_infile_opt "$cmdline" "-cr" optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_jobspec "$jobspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
get_callreg_opt()
{
    local callregf=$1

    if [ ${callregf} = ${NOFILE} ]; then
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

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_jobspec "$jobspec"` || exit 1
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
}

########
manta_somatic_define_opts()
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

    # -callregf option
    define_cmdline_infile_opt "$cmdline" "-cr" optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_jobspec "$jobspec"` || exit 1
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
strelka_germline_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"
}

########
strelka_germline_define_opts()
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
    
    # -cpus option
    local cpus
    cpus=`extract_cpus_from_jobspec "$jobspec"` || exit 1
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
}

########
strelka_somatic_define_opts()
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

    # -manta-outd option
    local manta_dep=`find_dependency_for_step "${jobspec}" manta_somatic`
    if [ ${manta_dep} != ${DEP_NOT_FOUND} ]; then
        local manta_outd=`get_default_outd_for_dep "${cmdline}" "${manta_dep}"`
        define_opt "-manta-outd" ${manta_outd} optlist || exit 1
    fi
    
    # -callregf option
    define_cmdline_infile_opt "$cmdline" "-cr" optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_jobspec "$jobspec"` || exit 1
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

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_jobspec "$jobspec"` || exit 1
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
    msisensor scan -d ${ref} -o ${step_outd}/msisensor.list 2>&1 || exit 1

    # Run MSIsensor analysis
    logmsg "* Executing msisensor msi..."
    msisensor msi -d ${step_outd}/msisensor.list -n ${normalbam} -t ${tumorbam} -o ${step_outd}/output -l 1 -q 1 -b ${cpus} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Dectivating conda environment..."
    conda deactivate 2>&1
    
    display_end_step_message
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
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -wcr option
    define_cmdline_infile_opt "$cmdline" "-wcr" optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_jobspec "$jobspec"` || exit 1
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
    logmsg "* Dectivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
facets_explain_cmdline_opts()
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
facets_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

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
facets()
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
    Rscript ${bindir}/run_facets -c ${step_outd}/snp-pileup-counts.csv 2>&1 || exit 1

    # Deactivate conda environment if needed
    if [ -z "${FACETS_HOME_DIR}" ]; then
        logmsg "* Dectivating conda environment..."
        conda deactivate 2>&1
    fi

    display_end_step_message
}

########
ascatngs_explain_cmdline_opts()
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
    description="SNP GC correction file required by AscatNGS"
    explain_cmdline_opt "-sg" "<string>" "$description"    

    # -mc option
    description="Name of male sex chromosome required by AscatNGS"
    explain_cmdline_opt "-mc" "<string>" "$description"    
}

########
ascatngs_define_opts()
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

    # -mc option
    define_cmdline_opt "$cmdline" "-mc" optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_jobspec "$jobspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
ascatngs()
{
    display_begin_step_message

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

    display_end_step_message
}

########
sequenza_explain_cmdline_opts()
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
sequenza_define_opts()
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

    # Save option list
    save_opt_list optlist    
}

########
sequenza()
{
    display_begin_step_message

    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || exit 1

    # Generate pileup files
    samtools mpileup -f $ref $normalbam | gzip > ${step_outd}/normal.pileup.gz || exit 1
    samtools mpileup -f $ref $tumorbam | gzip > ${step_outd}/tumor.pileup.gz || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || exit 1
    
    # Generate GC content file
    sequenza-utils.py GC-windows -w 50 $ref | gzip > ${step_outd}/ref.gc50Base.txt.gz || exit 1

    # Generate seqz file
    sequenza-utils.py pileup2seqz -gc ${step_outd}/ref.gc50Base.txt.gz -n ${step_outd}/normal.pileup.gz -t ${step_outd}/tumor.pileup.gz | gzip > ${step_outd}/seqz.gz || exit 1

    # Execute sequenza
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    Rscript ${bindir}/run_sequenza -s ${step_outd}/seqz.gz -o ${step_outd} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
lumpy_explain_cmdline_opts()
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
lumpy_define_opts()
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

    # Save option list
    save_opt_list optlist    
}

########
lumpy()
{
    display_begin_step_message

    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate lumpy 2>&1 || exit 1

    lumpyexpress -B ${normalbam},${tumorbam} -o ${step_outd}/out.vcf

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
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
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1

    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -egastr option
    define_cmdline_opt "$cmdline" "-egastr" optlist || exit 1

    # -egacred option
    define_cmdline_opt "$cmdline" "-egacred" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local normalbam=${bamdir_fullname}/normal.bam
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
download_ega_tum_bam_explain_cmdline_opts()
{
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"

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
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1
    
    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -egastr option
    define_cmdline_opt "$cmdline" "-egastr" optlist || exit 1

    # -egacred option
    define_cmdline_opt "$cmdline" "-egacred" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local tumorbam=${bamdir_fullname}/tumor.bam
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
download_aws_norm_bam_explain_cmdline_opts()
{
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"

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
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1
    
    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local normalbam=${bamdir_fullname}/normal.bam
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
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"

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
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1

    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -tumorbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local tumorbam=${bamdir_fullname}/tumor.bam
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
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"

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
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1
    
    # -extn option
    define_cmdline_opt "$cmdline" "-extn" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -normalbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local normalbam=${bamdir_fullname}/normal.bam
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
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"

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
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1
    
    # -extt option
    define_cmdline_opt "$cmdline" "-extt" optlist || exit 1

    # -nt option
    define_cmdline_nonmandatory_opt "$cmdline" "-nt" ${DEFAULT_NUMBER_OF_DOWNLOAD_TRIES} optlist || exit 1

    # -tumorbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local tumorbam=${bamdir_fullname}/tumor.bam
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
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"

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
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1
    
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
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local normalbam=${bamdir_fullname}/normal.bam
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
    logmsg "* Executing ascp..."
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
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"

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
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1
    
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
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local tumorbam=${bamdir_fullname}/tumor.bam
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
    local tumorbam_file=`read_opt_value_from_line "$*" "-extn"`
    local aspera_user=`read_opt_value_from_line "$*" "-asperausr"`
    local aspera_passwd=`read_opt_value_from_line "$*" "-asperapwd"`
    local aspera_server=`read_opt_value_from_line "$*" "-asperaserv"`
    local egadecrypt_pwd=`read_opt_value_from_line "$*" "-egadecrpwd"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local max_trans_rate=${DEFAULT_ASP_MAX_TRANS_RATE}

    # Download file
    logmsg "* Executing ascp..."
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
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"
}

########
index_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1

    # -normalbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local normalbam=${bamdir_fullname}/normal.bam
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
index_tum_bam_explain_cmdline_opts()
{
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"
}

########
index_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1

    # -tumorbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local tumorbam=${bamdir_fullname}/tumor.bam
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
sort_norm_bam_explain_cmdline_opts()
{
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"
}

########
sort_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1

    # -normalbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local normalbam=${bamdir_fullname}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_jobspec "$jobspec"` || exit 1
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
sort_tum_bam_explain_cmdline_opts()
{
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"
}

########
sort_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1

    # -tumorbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local tumorbam=${bamdir_fullname}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_jobspec "$jobspec"` || exit 1
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
    logmsg "* Dectivating conda environment..."
    conda deactivate > ${step_outd}/conda_deactivate.log 2>&1

    display_end_step_message
}

########
filter_norm_bam_contigs_explain_cmdline_opts()
{
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"

    # -r option
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"
}

########
filter_norm_bam_contigs_define_opts()
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

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1

    # -normalbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local normalbam=${bamdir_fullname}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist
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

    # Generate bed file for genome reference
    logmsg "* Executing gen_bed_for_genome..."
    ${bindir}/gen_bed_for_genome -r ${ref} -o ${step_outd}/genref
    
    # Filter normal bam file
    logmsg "* Executing samtools view..."
    samtools view -b -L ${step_outd}/genref.bed ${normalbam} > ${step_outd}/filtered.bam || exit 1

    # Replace initial bam file by the filtered one
    mv ${step_outd}/filtered.bam ${normalbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
filter_tum_bam_contigs_explain_cmdline_opts()
{
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"

    # -r option
    description="Reference genome file (required)"
    explain_cmdline_opt "-r" "<string>" "$description"
}

########
filter_tum_bam_contigs_define_opts()
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

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1

    # -tumorbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local tumorbam=${bamdir_fullname}/tumor.bam
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

    # Generate bed file for genome reference
    logmsg "* Executing gen_bed_for_genome..."
    ${bindir}/gen_bed_for_genome -r ${ref} -o ${step_outd}/genref
    
    # Filter tumor bam file
    logmsg "* Executing samtools view..."
    samtools view -b -L ${step_outd}/genref.bed ${tumorbam} > ${step_outd}/filtered.bam || exit 1

    # Replace initial bam file by the filtered one
    mv ${step_outd}/filtered.bam ${tumorbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    display_end_step_message
}

########
delete_bam_files_explain_cmdline_opts()
{    
    # -bamdir option
    description="Name of shared directory (without path) to perform operations on bam files (${DEFAULT_BAMDIR} by default)"
    explain_cmdline_opt "-bamdir" "<string>" "$description"
}

########
delete_bam_files_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    optlist=""

    # Define the -step-outd option, the output directory for the step,
    # which will have the same name of the step
    define_default_step_outd_opt "$cmdline" "$jobspec" optlist || exit 1

    # -bamdir option
    define_cmdline_nonmandatory_opt_shdir "$cmdline" "-bamdir" ${DEFAULT_BAMDIR} optlist || exit 1

    # -normalbam option
    local bamdir_fullname
    bamdir_fullname=`get_default_nonmandatory_opt_shdirname "${cmdline}" "-bamdir" ${DEFAULT_BAMDIR}` || exit 1
    local normalbam=${bamdir_fullname}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam=${bamdir_fullname}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
delete_bam_files()
{
    display_begin_step_message

    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Delete normal bam file
    logmsg "Removing normal bam file..."
    rm ${normalbam} 2>&1 || exit 1
    
    # Delete tumor bam file
    logmsg "Removing tumor bam file..."
    rm ${tumorbam} 2>&1 || exit 1

    display_end_step_message
}
