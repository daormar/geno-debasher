# *- bash -*

########
print_desc()
{
    echo "convert_snppos_to_snpgcc convert snp positions file to gc correction file"
    echo "type \"convert_snppos_to_snpgcc --help\" to get usage information"
}

########
usage()
{
    echo "convert_snppos_to_snpgcc   -r <string> -s <string> -o <string>"
    echo "                           [--help]"
    echo ""
    echo "-r <string>                File with reference genome"
    echo "-s <string>                SNP positions file"
    echo "-o <string>                Output file"
    echo "--help                     Display this help and exit"
}

########
read_pars()
{
    r_given=0
    o_given=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "-r") shift
                  if [ $# -ne 0 ]; then
                      ref=$1
                      r_given=1
                  fi
                  ;;
            "-s") shift
                  if [ $# -ne 0 ]; then
                      snpposfile=$1
                      s_given=1
                  fi
                  ;;
            "-o") shift
                  if [ $# -ne 0 ]; then
                      outfile=$1
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

    if [ ${s_given} -eq 0 ]; then
        echo "Error! -s parameter not given!" >&2
        exit 1
    fi

    if [ ${o_given} -eq 0 ]; then
        echo "Error! -o parameter not given!" >&2
        exit 1
    fi
}

########
process_pars()
{
    # Initialize variables
    TMPDIR=`${MKTEMP} -d /tmp/convsnp.XXXXX`

    # Create directories
    mkdir ${TMPDIR}/splitPos ${TMPDIR}/splitGc ${TMPDIR}/splitGcLogs

    # Check variables
    if [ "${ASCAT_GCC_UTIL}" = "" ]; then
        echo "ERROR: ASCAT_GCC_UTIL shell variable with path to 'ascatSnpPanelGcCorrections.pl' tool is not defined. This tool is provided by the AscatNGS package (PERL5LIB variable may also need to be exported)" >&2
        return 1
    else
        echo "* Testing that ${ASCAT_GCC_UTIL} given in ASCAT_GCC_UTIL variable can be executed correctly..." >&2
        ${ASCAT_GCC_UTIL} > ${TMPDIR}/ascat_gcc_util_test.txt 2>&1
        if $GREP "USAGE" ${TMPDIR}/ascat_gcc_util_test.txt; then
            echo "Test of ${ASCAT_GCC_UTIL} successful"
        else
            cat ${TMPDIR}/ascat_gcc_util_test.txt
            echo "Error while testing ${ASCAT_GCC_UTIL}" >&2
            return 1
        fi
    fi

    # Split file
    echo "* Splitting input file..." >&2
    $SPLIT --number=l/10 -d ${snpposfile} ${TMPDIR}/splitPos/snpPos.
    
    # Process fragments
    echo "* Processing fragments..." >&2
    for file in `ls ${TMPDIR}/splitPos/`; do
        ${ASCAT_GCC_UTIL} ${ref} ${TMPDIR}/splitPos/${file} > ${TMPDIR}/splitGc/${file} 2> ${TMPDIR}/splitGcLogs/${file}.log &
    done

    # Wait until all processes finish
    wait
    
    # Merge solutions
    echo "* Merging solutions..." >&2
    $HEAD -n 1 ${TMPDIR}/splitGc/snpPos.00 > ${outfile}
    cat ${TMPDIR}/splitGc/snpPos.* | $GREP -vP 'Chr\tPosition' >> ${outfile}
    
    # Remove temporary directory
    rm -rf ${TMPDIR}
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

process_pars || exit 1
