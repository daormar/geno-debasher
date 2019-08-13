# *- bash -*

##########################
# GENOME REFERENCE STEPS #
##########################

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
