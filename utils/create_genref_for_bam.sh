# *- bash -*

# INCLUDE BASH LIBRARY
. ${PANPIPE_HOME_DIR}/bin/panpipe_lib || exit 1

########
print_desc()
{
    echo "create_genref_for_bam create genome reference specific for bam file"
    echo "type \"create_genref_for_bam --help\" to get usage information"
}

########
usage()
{
    echo "create_genref_for_bam  -r <string> -b <string> [-c2a <string>]"
    echo "                       -o <string> [--help]"
    echo ""
    echo "-r <string>            File with reference genome"
    echo "-b <string>            bam file"
    echo "-c2a <string>          File containing a mapping between contig names and"
    echo "                       accession numbers"
    echo "-o <string>            Output directory"
    echo "--help                 Display this help and exit"
}

########
read_pars()
{
    r_given=0
    b_given=0
    c2a_given=0
    contig_to_acc=${NOFILE}
    o_given=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "-r") shift
                  if [ $# -ne 0 ]; then
                      baseref=$1
                      r_given=1
                  fi
                  ;;
            "-b") shift
                  if [ $# -ne 0 ]; then
                      bam=$1
                      b_given=1
                  fi
                  ;;
            "-c2a") shift
                  if [ $# -ne 0 ]; then
                      contig_to_acc=$1
                      c2a_given=1
                  fi
                  ;;
            "-o") shift
                  if [ $# -ne 0 ]; then
                      outd=$1
                      o_given=1
                  fi
                  ;;
        esac
        shift
    done   
}

########
check_pars()
{
    if [ ${r_given} -eq 0 ]; then   
        echo "Error! -r parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${baseref} ]; then
            echo "Error! file ${baseref} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${b_given} -eq 0 ]; then   
        echo "Error! -b parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${bam} ]; then
            echo "Error! file ${bam} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${c2a_given} -eq 1 ]; then   
        if [ ! -f ${contig_to_acc} ]; then
            echo "Error! file ${contig_to_acc} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${o_given} -eq 0 ]; then   
        echo "Error! -o parameter not given!" >&2
        exit 1
    else
        if [ ! -d ${outd} ]; then
            echo "Error! directory ${outd} does not exist" >&2
            exit 1
        fi
    fi
}

########
print_pars()
{
    if [ ${r_given} -eq 1 ]; then
        echo "-r is ${baseref}" >&2
    fi

    if [ ${b_given} -eq 1 ]; then
        echo "-b is ${bam}" >&2
    fi

    if [ ${c2a_given} -eq 1 ]; then
        echo "-c2a is ${contig_to_acc}" >&2
    fi

    if [ ${o_given} -eq 1 ]; then
        echo "-o is ${outd}" >&2
    fi
}

########
contig_in_list()
{
    local contig=$1
    local clist=$2

    while read c; do
        if [ "$contig" = "$c" ]; then
            return 0
        fi
    done < ${clist}

    return 1
}

########
get_ref_contigs()
{
    local ref=$1

    samtools faidx ${baseref} || return 1
    $AWK '{printf "%s\n",$1}' ${baseref}.fai
}

########
get_bam_contigs()
{
    local bam=$1

    samtools view -H $bam | $AWK '{if($1=="@SQ") print substr($2,4)}'
}

########
get_missing_contig_names()
{
    local refcontigs=$1
    local bamcontigs=$2
            
    while read bamcontigname; do
        if ! contig_in_list $bamcontigname $refcontigs; then
            echo $bamcontigname
        fi
    done < $bamcontigs    
}

########
get_ref_contigs_to_keep()
{
    local refcontigs=$1
    local bamcontigs=$2
            
    while read refcontigname; do
        if contig_in_list $refcontigname $bamcontigs; then
            echo $refcontigname
        fi
    done < $refcontigs    
}

########
contig_is_accession()
{
    local contig=$1

    if [[ $contig == *"."* ]]; then
        return 0
    else
        return 1
    fi
}

########
map_contig_to_acc_using_file()
{
    local contig_to_acc=$1
    local contig=$2

    while read entry; do
        local fields=($entry)
        local num_fields=${#fields[@]}
        if [ ${num_fields} -eq 2 ]; then
            if [ ${fields[0]} = $contig ]; then
                echo ${fields[1]}
                break
            fi
        fi
    done < ${contig_to_acc}
}

########
map_contig_to_accession()
{
    local contig_to_acc=$1
    local contig=$2

    if contig_is_accession ${contig}; then
        echo ${contig}
    else
        if [ ${contig_to_acc} != "${NOFILE}" ]; then
            map_contig_to_acc_using_file ${contig_to_acc} ${contig} || return 1
        fi
    fi
}

########
get_contigs()
{
    local contig_to_acc=$1
    local contiglist=$2

    while read contig; do
        local accession=`map_contig_to_accession ${contig_to_acc} ${contig}` || return 1
        if [ "$accession" = "" ]; then
            echo "Error: contig $contig is not a valid accession nor there were mappings for it" >&2
            return 1
        else
            echo "Getting data for contig ${contig} (mapped to accession $accession)..." >&2
            ${biopanpipe_bindir}/get_entrez_fasta -a ${accession} | ${SED} "s/${accession}/${contig}/"; pipe_fail || return 1
        fi
    done < ${contiglist}
}

########
process_pars()
{
    outfile=$outd/enriched_genref.fa
    
    # Activate conda environment
    echo "* Activating conda environment (samtools)..." >&2
    conda activate samtools || exit 1

    # Get reference contigs
    echo "* Obtaining list of current reference contigs..." >&2
    get_ref_contigs $baseref > ${outd}/refcontigs

    # Get bam contigs
    echo "* Obtaining list of bam contigs..." >&2
    get_bam_contigs $bam > ${outd}/bamcontigs
    
    # Obtain list of contigs to keep in the reference file
    echo "* Obtaining list of reference contigs to keep..." >&2
    get_ref_contigs_to_keep ${outd}/refcontigs ${outd}/bamcontigs > ${outd}/ref_contigs_to_keep.txt || exit 1

    # Copy base genome reference without extra contigs
    echo "* Copying base genome reference without extra contigs..." >&2
    ${biopanpipe_bindir}/filter_contig_from_genref -g $baseref -l ${outd}/ref_contigs_to_keep.txt > $outfile

    # Obtain list of missing contigs
    echo "* Obtaining list of missing contigs..." >&2
    get_missing_contig_names ${outd}/refcontigs_to_keep ${outd}/bamcontigs > ${outd}/missing_contigs.txt || exit 1

    # Enrich base reference
    echo "* Enriching base reference..." >&2
    get_contigs ${contig_to_acc} ${outd}/missing_contigs.txt >> $outfile || { echo "Error during FASTA data downloading" >&2; exit 1; }

    # Index enriched reference
    echo "* Indexing created reference..." >&2
    samtools faidx ${outfile} || exit 1
    
    # Deactivate conda environment
    echo "* Deactivating conda environment..." >&2
    conda deactivate
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

print_pars || exit 1

process_pars || exit 1
