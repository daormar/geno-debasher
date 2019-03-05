# *- bash -*

########
print_desc()
{
    echo "get_bam_contigs get list of contigs contained in bam file"
    echo "type \"get_bam_contigs --help\" to get usage information"
}

########
usage()
{
    echo "get_bam_contigs  -b <string> [--help]"
    echo ""
    echo "-b <string>      bam file"
    echo "--help           Display this help and exit"
}

########
read_pars()
{
    b_given=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "-b") shift
                  if [ $# -ne 0 ]; then
                      bam=$1
                      b_given=1
                  fi
                  ;;
        esac
        shift
    done   
}

########
check_pars()
{
    if [ ${b_given} -eq 0 ]; then   
        echo "Error! -b parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${bam} ]; then
            echo "Error! file ${bam} does not exist" >&2
            exit 1
        fi
    fi
}

########
print_pars()
{
    if [ ${b_given} -eq 1 ]; then
        echo "-b is ${bam}" >&2
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
    conda activate samtools || exit 1

    # Obtain bam contigs
    samtools idxstats $bam > ${outpref}.bamstats || exit 1
    cat ${outpref}.bamstats | filter_bam_stats | ${SORT}
    
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
