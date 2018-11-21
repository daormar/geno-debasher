# *- bash -*

# INCLUDE BASH LIBRARY
. ${bindir}/bam_utils_lib

########
print_desc()
{
    echo "check_ref_and_bam_contigs check contig equivalence between bam and reference genome"
    echo "type \"check_ref_and_bam_contigs --help\" to get usage information"
}

########
usage()
{
    echo "check_ref_and_bam_contigs  -r <string> -b <string> -o <string>"
    echo "                           [--help]"
    echo ""
    echo "-r <string>                File with reference genome"
    echo "-b <string>                bam file"
    echo "-o <string>                Prefix of output files"
    echo "--help                     Display this help and exit"
}

########
read_pars()
{
    r_given=0
    b_given=0
    o_given=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "--version") version
                         exit 1
                         ;;
            "-r") shift
                  if [ $# -ne 0 ]; then
                      ref=$1
                      r_given=1
                  fi
                  ;;
            "-b") shift
                  if [ $# -ne 0 ]; then
                      bam=$1
                      b_given=1
                  fi
                  ;;
            "-o") shift
                  if [ $# -ne 0 ]; then
                      outpref=$1
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
        if [ ! -f ${ref} ]; then
            echo "Error! file ${ref} does not exist" >&2
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

    if [ ${o_given} -eq 0 ]; then
        echo "Error! -o parameter not given!" >&2
        exit 1
    fi
}

########
print_pars()
{
    if [ ${r_given} -eq 1 ]; then
        echo "-r is ${ref}" >&2
    fi

    if [ ${b_given} -eq 1 ]; then
        echo "-b is ${bam}" >&2
    fi

    if [ ${o_given} -eq 1 ]; then
        echo "-o is ${outpref}" >&2
    fi
}

########
filter_bam_stats()
{
    ${AWK} '{if($3>0 || $4>0) printf"%s %d\n",$1,$2}'
}

########
process_pars()
{
    conda activate base || exit 1

    # Obtain reference contigs
    samtools faidx ${ref}
    $AWK '{printf "%s %d\n",$1,$2}' ${ref}.fai | ${SORT} > ${outpref}.refcontigs

    # Obtain bam contigs
    samtools idxstats $bam > ${outpref}.bamstats || exit 1
    cat ${outpref}.bamstats | filter_bam_stats | ${SORT} > ${outpref}.bamcontigs || exit 1
    
    conda deactivate

    # Execute diff command
    $DIFF ${outpref}.refcontigs ${outpref}.bamcontigs > ${outpref}.contigdiffs || { echo "Contigs in reference and bam files differ!" >&2 ; exit 1; }
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
