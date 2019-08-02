# *- bash -*

# INCLUDE BASH LIBRARY
. ${biopanpipe_bindir}/bam_common_lib || exit 1

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
    define_shared_dir ${DATADIR_BASENAME}
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
create_genref_for_bam_document()
{
    step_description "Creates a genome reference file for a given \`bam\` file. For this purpose, the step starts from a basic genome reference file, removing those contigs not present in the \`bam\` file and downloading or copying missing ones from the Internet or from previously existing files."
}

########
create_genref_for_bam_explain_cmdline_opts()
{
    # -br option
    description="Base reference genome file"
    explain_cmdline_req_opt "-br" "<string>" "$description"

    # -bam option
    description="bam file (required if no downloading steps or paths of normal or tumor bam files have been defined)"
    explain_cmdline_opt "-bam" "<string>" "$description"

    # -cm option
    description="File containing a mapping between contig names and accession numbers"
    explain_cmdline_opt "-cm" "<string>" "$description"

    # -fbr option
    description="Name of fallback genome reference file. If creation process fails, this file is copied as reference output file instead"
    explain_cmdline_opt "-fbr" "<string>" "$description"
}

########
get_bam_filename()
{
    local cmdline=$1
    local given=0

    # Check -bam option
    local bam
    bam=`read_opt_value_from_line "$cmdline" "-bam"` && given=1
    if [ $given -eq 1 ]; then
        # -bam option was given
        file_exists $bam || { errmsg "file $bam does not exist" ; return 1; }
        echo $bam
        return 0
    fi
    
    # Check -n option
    local normalbam
    normalbam=`read_opt_value_from_line "$cmdline" "-n"` && given=1
    if [ $given -eq 1 ]; then
        file_exists $normalbam || { errmsg "file $normalbam does not exist" ; return 1; }
        echo $normalbam
        return 0
    fi
    
    # Check -extn option
    if check_opt_given "$cmdline" "-extn"; then
        local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
        normalbam=${abs_datadir}/normal.bam
        echo $normalbam
        return 0
    fi

    # Check -t option
    local tumorbam
    tumorbam=`read_opt_value_from_line "$cmdline" "-t"` && given=1
    if [ $given -eq 1 ]; then
        file_exists $tumorbam || { errmsg "file $tumorbam does not exist" ; return 1; }
        echo $tumorbam
        return 0
    fi

    # Check -extt option
    if check_opt_given "$cmdline" "-extt"; then
        local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
        tumorbam=${abs_datadir}/tumor.bam
        echo $tumorbam
        return 0
    fi

    errmsg "-bam, -n, -extn, -t or -extt options should be given"
    return 1
}

########
create_genref_for_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""
    
    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -br option
    define_cmdline_infile_opt "$cmdline" "-br" optlist || exit 1

    # -bam option
    local bam
    bam=`get_bam_filename "$cmdline"` || exit 1
    define_opt "-bam" $bam optlist || exit 1

    # -cm option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cm" ${NOFILE} optlist || exit 1

    # -fbr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-fbr" ${NOFILE} optlist || exit 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`

    # -outfile option
    define_opt "-outfile" ${abs_datadir}/genref.fa optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
get_create_genref_for_bam_cm_opt()
{
    local value=$1

    if [ "${value}" = ${NOFILE} ]; then
        echo ""
    else
        echo "-cm ${value}"
    fi
}

########
create_genref_for_bam()
{
    # Initialize variables
    local baseref=`read_opt_value_from_line "$*" "-br"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local bam=`read_opt_value_from_line "$*" "-bam"`
    local contig_mapping=`read_opt_value_from_line "$*" "-cm"`
    local fallback_genref=`read_opt_value_from_line "$*" "-fbr"`
    local outfile=`read_opt_value_from_line "$*" "-outfile"`

    # Create genome reference
    local cm_opt=`get_create_genref_for_bam_cm_opt ${contig_mapping}`
    if ${biopanpipe_bindir}/create_genref_for_bam -r ${baseref} -b ${bam} ${cm_opt} -o ${step_outd}; then
        # Move resulting files
        mv ${step_outd}/genref_for_bam.fa ${outfile}
        mv ${step_outd}/genref_for_bam.fa.fai ${outfile}.fai
    else
        # Genome reference creation failed, check if fallback file was
        # provided
        if [ "${fallback_genref}" = ${NOFILE} ]; then
            exit 1
        else
            logmsg "Genome reference creation failed but fallback file was provided"
            # Copy fallback file
            logmsg "* Copying fallback file (${fallback_genref})..."
            cp ${fallback_genref} ${outfile} || exit 1
            # Index fallback reference
            logmsg "* Indexing fallback reference..."
            conda activate samtools 2>&1 || exit 1
            samtools faidx ${outfile} || exit 1
            conda deactivate
        fi
    fi
}

########
create_genref_for_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
get_contig_list_from_file()
{
    local file=$1
    file_exists $file || { errmsg "file $file containing contig list does not exist" ; return 1; }
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
        # Check -br option
        if check_opt_given "$cmdline" "-br"; then
            local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
            ref=${abs_datadir}/genref.fa
            echo $ref
            return 0
        fi

        errmsg "-r or -br options should be given"
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
manta_germline_document()
{
    step_description "Analyzes a normal \`bam\` file using Manta."
}

########
manta_germline_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call"
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
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -cr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

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
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
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
}

########
manta_germline_conda_envs()
{
    define_conda_env manta manta.yml
}

########
manta_somatic_document()
{
    step_description "Analyzes a pair of normal and tumor \`bam\` files using Manta."
}

########
manta_somatic_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call"
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
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -cr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

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
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
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
}

########
manta_somatic_conda_envs()
{
    define_conda_env manta manta.yml
}

########
cnvkit_document()
{
    step_description "Analyzes a pair of normal and tumor \`bam\` files using CNVkit."
}

########
cnvkit_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
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
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

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
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
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
}

########
cnvkit_conda_envs()
{
    define_conda_env cnvkit cnvkit.yml
}

########
strelka_germline_document()
{
    step_description "Analyzes a normal \`bam\` files using Strelka."
}

########
strelka_germline_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call"
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
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -cr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

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
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
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
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call"
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
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

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
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

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
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
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
    conda deactivate 2>&1
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
    description="Reference genome file"
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
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

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
    # Initialize variables
    local ref=$1
    local normalbam=$2
    local step_outd=$3

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate platypus 2>&1 || exit 1

    # Run Platypus
    logmsg "* Executing Platypus.py..."
    Platypus.py callVariants --bamFiles=${normalbam} --refFile=${ref} --output=${step_outd}/output.vcf --logFileName=${step_outd}/platypus.log --verbosity=1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
platypus_germline_local()
{
    # Initialize variables
    local ref=$1
    local normalbam=$2
    local step_outd=$3

    # Run Platypus
    logmsg "* Executing Platypus.py..."
    python ${PLATYPUS_HOME_DIR}/bin/Platypus.py callVariants --bamFiles=${normalbam} --refFile=${ref} --output=${step_outd}/output.vcf --verbosity=1 2>&1 || exit 1
}

########
platypus_germline()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
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
    description="Reference genome file"
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
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

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
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
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
    explain_cmdline_req_opt "-wcr" "<string>" "$description"

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
    # Initialize variables
    local wcref=`read_opt_value_from_line "$*" "-wcr"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    
    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate wisecondorx 2>&1 || exit 1

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
    description="SNP vcf file"
    explain_cmdline_req_opt "-sv" "<string>" "$description"        
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
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local snpvcf=`read_opt_value_from_line "$*" "-sv"`

    # Activate conda environment
    logmsg "* Activating conda environment (snp-pileup)..."
    conda activate snp-pileup 2>&1 || exit 1

    # Execute snp-pileup
    logmsg "* Executing snp-pileup..."
    snp-pileup ${snpvcf} ${step_outd}/snp-pileup-counts.csv ${normalbam} ${tumorbam} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment if needed
    logmsg "* Activating conda environment (facets)..."
    conda activate facets 2>&1 || exit 1
            
    # Execute facets
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing facets..."
    Rscript ${biopanpipe_bindir}/run_facets -c ${step_outd}/snp-pileup-counts.csv -o ${step_outd} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress snp-pileup file
    logmsg "* Compressing snp-pileup file..."
    ${GZIP} ${step_outd}/snp-pileup-counts.csv || exit 1
}

########
snp_pileup_plus_facets_conda_envs()
{
    define_conda_env facets snp-pileup.yml
    define_conda_env facets facets.yml
}

########
gen_sequenza_gcc_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"
}

########
gen_sequenza_gcc_define_opts()
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

    # Get data directory
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`

    # -outfile option
    define_opt "-outfile" ${abs_datadir}/sequenza_gccfile.txt.gz optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
gen_sequenza_gcc()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local outfile=`read_opt_value_from_line "$*" "-outfile"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate sequenza 2>&1 || exit 1

    # Generate GC content file
    logmsg "* Generating GC content file..."    
    sequenza-utils gc_wiggle -w 50 -f $ref -o ${step_outd}/sequenza_gccfile.txt.gz || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Move result file to final location
    mv ${step_outd}/sequenza_gccfile.txt.gz $outfile || exit 1
}

########
sequenza_explain_cmdline_opts()
{
    # -gcc option
    description="GC content wiggle file for sequenza"
    explain_cmdline_opt "-gcc" "<string>" "$description"
}

########
get_gcc_filename()
{
    local cmdline=$1
    local stepspec=$2
    local given=0

    gccfile=`read_opt_value_from_line "$cmdline" "-gcc"` && given=1
    if [ $given -eq 1 ]; then
        # -gcc option was given
        file_exists $gccfile || { errmsg "file $gccfile does not exist" ; return 1; }
        echo $gccfile
    else
        # Check if gen_sequenza_gcc step dependency was defined
        local gen_sequenza_gcc_dep=`find_dependency_for_step "${stepspec}" gen_sequenza_gcc`
        if [ ${gen_sequenza_gcc_dep} != ${DEP_NOT_FOUND} ]; then
            local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
            gccfile=${abs_datadir}/sequenza_gccfile.txt.gz
            echo $gccfile
            return 0            
        else
            errmsg "-gcc or dependency with gen_sequenza_gcc_dep step should be given"
            return 1
        fi            
    fi
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
    local gccfile
    gccfile=`get_gcc_filename "$cmdline" "$stepspec"` || exit 1
    define_opt "-gcc" $gccfile optlist || exit 1

    # Get normal pileup file
    local npileupdir
    npileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" samtools_mpileup_norm_bam` || { errmsg "Error: dependency samtools_mpileup_norm_bam not defined for sequenza"; exit 1; }
    local npileup=${npileupdir}/normal.pileup.gz
    define_opt "-npileup" ${npileup} optlist || exit 1

    # Get tumor pileup file
    tpileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" samtools_mpileup_tum_bam` || { errmsg "Error: dependency samtools_mpileup_tum_bam not defined for sequenza"; exit 1; }
    tpileup=${tpileupdir}/tumor.pileup.gz
    define_opt "-tpileup" ${tpileup} optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
sequenza()
{
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
    description="GC content wiggle file for bam2seqz"
    explain_cmdline_opt "-gcc" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"   
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
    local gccfile
    gccfile=`get_gcc_filename "$cmdline" "$stepspec"` || exit 1
    define_opt "-gcc" $gccfile optlist || exit 1

    # Get normal pileup directory
    local npileupdir
    npileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" parallel_samtools_mpileup_norm_bam` || { errmsg "Error: dependency parallel_samtools_mpileup_norm_bam not defined for parallel_bam2seqz"; exit 1; }

    # Get tumor pileup directory
    local tpileupdir
    tpileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" parallel_samtools_mpileup_tum_bam` || { errmsg "Error: dependency parallel_samtools_mpileup_tum_bam not defined for parallel_bam2seqz"; exit 1; }

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
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
    explain_cmdline_req_opt "-lc" "<string>" "$description"   
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

    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
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

    # NOTE: the use of the bgzip command requires to have the sequenza
    # environment activated (since it has tabix installed)
    ${ZCAT} ${filenames} | $AWK '{if (NR!=1 && $1 != "chromosome") {print $0}}' | bgzip ; pipe_fail || exit 1
}
 
########
seqzmerge_plus_sequenza()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local seqzdir=`read_opt_value_from_line "$*" "-seqzdir"`
    local clist=`read_opt_value_from_line "$*" "-lc"`

    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || exit 1

    # Merge seqz files
    logmsg "* Merging seqz files..."
    seqzmerge ${clist} ${seqzdir}  > ${step_outd}/merged_seqz.gz || exit 1

    logmsg "* Applying tabix over merged seqz file..."
    tabix -f -s 1 -b 2 -e 2 -S 1 ${step_outd}/merged_seqz.gz || exit 1
                
    # Execute sequenza
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing sequenza..."
    Rscript ${biopanpipe_bindir}/run_sequenza -s ${step_outd}/merged_seqz.gz -o ${step_outd} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
seqzmerge_plus_sequenza_conda_envs()
{
    define_conda_env sequenza sequenza.yml
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

    # -lx option
    description="File with regions to exclude in bed format"
    explain_cmdline_opt "-lx" "<string>" "$description"    
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

    # -lx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-lx" ${NOFILE} optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
get_lumpyexpress_x_opt()
{
    local value=$1

    if [ "${value}" = ${NOFILE} ]; then
        echo ""
    else
        echo "-x ${value}"
    fi
}

########
lumpy()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local exclude=`read_opt_value_from_line "$*" "-lx"`

    if [ -z "${LUMPY_HOME_DIR}" ]; then
        # Activate conda environment
        logmsg "* Activating conda environment..."
        conda activate lumpy 2>&1 || exit 1

        logmsg "* Executing lumpyexpress..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        lumpyexpress -B ${tumorbam},${normalbam} ${x_opt} -o ${step_outd}/out.vcf || exit 1
        
        # Deactivate conda environment
        logmsg "* Deactivating conda environment..."
        conda deactivate 2>&1
    else
        logmsg "* Executing lumpyexpress..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        ${LUMPY_HOME_DIR}/bin/lumpyexpress -B ${tumorbam},${normalbam} ${x_opt} -o ${step_outd}/out.vcf || exit 1
    fi
}

########
lumpy_conda_envs()
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
    explain_cmdline_req_opt "-lc" "<string>" "$description"   
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
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Generate exclusion bed file
    gen_exclusion_bed_given_bam ${normalbam} ${contig} > ${step_outd}/${contig}.bed || exit 1

    if [ -z "${LUMPY_HOME_DIR}" ]; then
        # Activate conda environment
        logmsg "* Activating conda environment (lumpy)..."
        conda activate lumpy 2>&1 || exit 1
    
        logmsg "* Executing lumpyexpress (contig $contig)..."
        lumpyexpress -B ${tumorbam},${normalbam} -T ${step_outd}/tmp_${contig} -x ${step_outd}/${contig}.bed -o ${step_outd}/out${contig}.vcf || exit 1

        # Deactivate conda environment
        logmsg "* Deactivating conda environment..."
        conda deactivate 2>&1
    else
        logmsg "* Executing lumpyexpress (contig $contig)..."
        ${LUMPY_HOME_DIR}/bin/lumpyexpress -B ${tumorbam},${normalbam} -T ${step_outd}/tmp_${contig} -x ${step_outd}/${contig}.bed -o ${step_outd}/out${contig}.vcf || exit 1
    fi
}

########
parallel_exclude_plus_lumpy_conda_envs()
{
    define_conda_env samtools samtools.yml
    define_conda_env lumpy lumpy.yml
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
filter_bam_contig_samtools()
{
    local inbam=$1
    local contig=$2
    local outbam=$3    
    local error=0
    
    samtools view -h -O BAM $inbam $contig > ${outbam} 2> ${outbam}.log || error=1

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
filter_bam_contig_sambamba()
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
parallel_lumpy_explain_cmdline_opts()
{
    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"

    # -lx option
    description="File with regions to exclude in bed format"
    explain_cmdline_opt "-lx" "<string>" "$description"    
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

    # -lx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-lx" ${NOFILE} optlist || exit 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`

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
        tumorbam=${abs_datadir}/tumor_${contig}.bam
        define_opt "-tumorbam" ${tumorbam} specific_optlist || exit 1
        define_opt "-contig" ${contig} specific_optlist || exit 1
        
        save_opt_list specific_optlist
    done
}

########
parallel_lumpy()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local exclude=`read_opt_value_from_line "$*" "-lx"`

    if [ -z "${LUMPY_HOME_DIR}" ]; then
        # Activate conda environment
        logmsg "* Activating conda environment (lumpy)..."
        conda activate lumpy 2>&1 || exit 1
        
        logmsg "* Executing lumpyexpress (contig $contig)..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        lumpyexpress -B ${tumorbam},${normalbam} ${x_opt} -T ${step_outd}/tmp_${contig} -o ${step_outd}/out${contig}.vcf || exit 1
        
        # Deactivate conda environment
        logmsg "* Deactivating conda environment..."
        conda deactivate 2>&1
    else
        logmsg "* Executing lumpyexpress (contig $contig)..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        ${LUMPY_HOME_DIR}/bin/lumpyexpress -B ${tumorbam},${normalbam} ${x_opt} -T ${step_outd}/tmp_${contig} -o ${step_outd}/out${contig}.vcf || exit 1
    fi
}

########
parallel_lumpy_conda_envs()
{
    define_conda_env lumpy lumpy.yml
}

########
smoove_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -lx option
    description="File with regions to exclude in bed format"
    explain_cmdline_opt "-lx" "<string>" "$description"    
}

########
smoove_define_opts()
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

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -lx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-lx" ${NOFILE} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
get_smoove_exclude_opt()
{
    local value=$1

    if [ "${value}" = ${NOFILE} ]; then
        echo ""
    else
        echo "--exclude ${value}"
    fi
}

########
smoove()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local exclude=`read_opt_value_from_line "$*" "-lx"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate smoove 2>&1 || exit 1

    logmsg "* Executing smoove..."
    export TMPDIR=${step_outd}
    local exclude_opt=`get_smoove_exclude_opt ${exclude}`
    local project_name="smoove"
    command smoove call --outdir ${step_outd} ${exclude_opt} --name ${project_name} --fasta ${ref} -p ${cpus} --genotype ${normalbam} ${tumorbam} || exit 1
        
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
smoove_conda_envs()
{
    define_conda_env smoove smoove.yml
}

########
delly_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    

    # -dx option
    description="File with regions to exclude in bed format"
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
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -dx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-dx" ${NOFILE} optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
get_delly_x_opt()
{
    local value=$1

    if [ "${value}" = ${NOFILE} ]; then
        echo ""
    else
        echo "-x ${value}"
    fi
}

########
delly()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local exclude=`read_opt_value_from_line "$*" "-dx"`

    # Activate conda environment
    logmsg "* Activating conda environment (delly)..."
    conda activate delly 2>&1 || exit 1

    logmsg "* Executing delly..."
    # "command" built-in is used here to execute the "delly" program
    # instead of the "delly" function
    local x_opt=`get_delly_x_opt ${exclude}`
    command delly call -g ${ref} ${x_opt} -o ${step_outd}/out.bcf ${tumorbam} ${normalbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bcftools)..."
    conda activate bcftools 2>&1 || exit 1

    # Convert bcf output to vcf
    logmsg "* Converting bcf output into vcf..."
    bcftools view ${step_outd}/out.bcf > ${step_outd}/out.vcf || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
parallel_lumpy_conda_envs()
{
    define_conda_env bcftools bcftools.yml
    define_conda_env delly delly.yml
}

########
parallel_delly_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -dx option
    description="File with regions to exclude in bed format"
    explain_cmdline_opt "-dx" "<string>" "$description"    

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"   
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
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`

    # -dx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-dx" ${NOFILE} optlist || exit 1

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
        tumorbam=${abs_datadir}/tumor_${contig}.bam
        define_opt "-tumorbam" ${tumorbam} specific_optlist || exit 1
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_delly()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local exclude=`read_opt_value_from_line "$*" "-dx"`
    
    # Activate conda environment
    logmsg "* Activating conda environment (delly)..."
    conda activate delly 2>&1 || exit 1
    
    logmsg "* Executing delly (contig $contig)..."
    # "command" built-in is used here to execute the "delly" program
    # instead of the "delly" function
    local x_opt=`get_delly_x_opt ${exclude}`
    command delly call -g $ref ${x_opt} -o ${step_outd}/out${contig}.bcf ${tumorbam} ${normalbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bcftools)..."
    conda activate bcftools 2>&1 || exit 1

    # Convert bcf output to vcf
    logmsg "* Converting bcf output into vcf..."
    bcftools view ${step_outd}/out${contig}.bcf > ${step_outd}/out${contig}.vcf || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
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
    explain_cmdline_req_opt "-lc" "<string>" "$description"   
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
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Determine vcf directory
    vcfdir=`get_vcfdir_for_svtyper "${stepspec}"` || { errmsg "Error: vcf directory for svtyper could not be determined"; exit 1; }
    
    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
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
}

########
parallel_svtyper_conda_envs()
{
    define_conda_env svtyper svtyper.yml
}

########
download_ega_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="External database id of normal bam file to download"
    explain_cmdline_opt "-extn" "<string>" "$description"

    # -egastr option
    description="Number of streams used by the EGA download client (${DEFAULT_NUMBER_OF_EGA_DOWNLOAD_STREAMS} by default)"
    explain_cmdline_opt "-egastr" "<int>" "$description"

    # -egacred option
    description="File with EGA download client credentials"
    explain_cmdline_req_opt "-egacred" "<string>" "$description"

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
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
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

        ntry=$((ntry+1))
    done

    logmsg "All download attempts failed!"

    return 1
}

########
download_ega_norm_bam()
{
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
    description="External database id of tumor bam file to download"
    explain_cmdline_req_opt "-extt" "<string>" "$description"

    # -egastr option
    description="Number of streams used by the EGA download client (${DEFAULT_NUMBER_OF_EGA_DOWNLOAD_STREAMS} by default)"
    explain_cmdline_opt "-egastr" "<int>" "$description"

    # -egacred option
    description="File with EGA download client credentials"
    explain_cmdline_req_opt "-egacred" "<string>" "$description"

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
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_ega_tum_bam()
{
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
    ega_download_retry ${egastr} ${egacred} ${egaid_tumorbam} ${step_outd}/tumor.bam ${download_tries} || exit 1

    # Move file
    mv ${step_outd}/tumor.bam ${tumorbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
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
    description="External database id of normal bam file to download"
    explain_cmdline_req_opt "-extn" "<string>" "$description"

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
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
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
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local icgcid_normalbam=`read_opt_value_from_line "$*" "-extn"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing score-client..."
    ${ICGCSTOR_HOME_DIR}/bin/score-client --profile aws download --object-id ${icgcid_normalbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1
}

########
download_aws_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of tumor bam file to download"
    explain_cmdline_req_opt "-extt" "<string>" "$description"

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
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_aws_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local icgcid_tumorbam=`read_opt_value_from_line "$*" "-extt"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing score-client..."
    ${ICGCSTOR_HOME_DIR}/bin/score-client --profile aws download --object-id ${icgcid_tumorbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1
}

########
download_collab_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="External database id of normal bam file to download"
    explain_cmdline_req_opt "-extn" "<string>" "$description"

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
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_collab_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local icgcid_normalbam=`read_opt_value_from_line "$*" "-extn"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing score-client..."
    ${ICGCSTOR_HOME_DIR}/bin/score-client --profile collab download --object-id ${icgcid_normalbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${normalbam} || exit 1
}

########
download_collab_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of tumor bam file to download"
    explain_cmdline_req_opt "-extt" "<string>" "$description"

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
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
download_collab_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local icgcid_tumorbam=`read_opt_value_from_line "$*" "-extt"`
    local download_tries=`read_opt_value_from_line "$*" "-nt"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`

    # Download file
    logmsg "* Executing score-client..."
    ${ICGCSTOR_HOME_DIR}/bin/score-client --profile collab download --object-id ${icgcid_tumorbam} --output-dir ${step_outd} 2>&1 || exit 1

    # Find bam file name
    local bam_file_name=`find_bam_filename ${step_outd}`
    
    if [ -z "${bam_file_name}" ]; then
        logmsg "Error: bam file not found after download process was completed"
        exit 1
    fi

    # Move file
    mv ${bam_file_name} ${tumorbam} || exit 1
}

########
download_ega_asp_norm_bam_explain_cmdline_opts()
{
    # -extn option
    description="External database id of normal bam file to download"
    explain_cmdline_req_opt "-extn" "<string>" "$description"

    # -asperausr option
    description="Username for Aspera server"
    explain_cmdline_req_opt "-asperausr" "<string>" "$description"

    # -asperapwd option
    description="Password for Aspera server"
    explain_cmdline_req_opt "-asperapwd" "<string>" "$description"

    # -asperaserv option
    description="Name of Aspera server"
    explain_cmdline_req_opt "-asperaserv" "<string>" "$description"

    # -egadecrpwd option
    description="File with EGA decryptor password"
    explain_cmdline_req_opt "-egadecrpwd" "<string>" "$description"
    
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
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam=${abs_datadir}/normal.bam
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
download_ega_asp_norm_bam()
{
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
}

########
download_ega_asp_tum_bam_explain_cmdline_opts()
{
    # -extt option
    description="External database id of normal bam file to download"
    explain_cmdline_req_opt "-extt" "<string>" "$description"

    # -asperausr option
    description="Username for Aspera server"
    explain_cmdline_req_opt "-asperausr" "<string>" "$description"

    # -asperapwd option
    description="Password for Aspera server"
    explain_cmdline_req_opt "-asperapwd" "<string>" "$description"

    # -asperaserv option
    description="Name of Aspera server"
    explain_cmdline_req_opt "-asperaserv" "<string>" "$description"

    # -egadecrpwd option
    description="File with EGA decryptor password"
    explain_cmdline_req_opt "-egadecrpwd" "<string>" "$description"
    
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
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local tumorbam=${abs_datadir}/tumor.bam
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
download_ega_asp_tum_bam()
{
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

    # -datadir option
    abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    define_opt "-datadir" ${abs_datadir} optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
delete_bam_files()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local abs_datadir=`read_opt_value_from_line "$*" "-datadir"`

    # Delete bam files
    rm -f ${abs_datadir}/*.bam || exit 1
}

########
clear_datadir_explain_cmdline_opts()
{
    :
}

########
clear_datadir_define_opts()
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

    # Save option list
    save_opt_list optlist
}

########
clear_datadir()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local abs_datadir=`read_opt_value_from_line "$*" "-datadir"`

    # Delete bam files
    rm -rf ${abs_datadir}/* || exit 1

    # Print README.txt file
    echo "NOTE: This directory was cleared by means of the 'clear_datadir' step" > ${abs_datadir}/README.txt || exit 1
}
