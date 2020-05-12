# Bio-PanPipe package
# Copyright (C)  Daniel Ortiz-Mart\'inez
#  
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#  
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.
  
# *- bash -*

########
print_desc()
{
    echo "convert_snppos_to_snpgcc converts snp positions file to gc correction file"
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
    s_given=0
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

    # Create temporary directories
    mkdir ${TMPDIR}/splitPos ${TMPDIR}/splitGc ${TMPDIR}/splitGcLogs

    # Remove temporary directories on exit
    trap "rm -rf ${TMPDIR} 2>/dev/null" EXIT
    
    # Check variables
    if [ "${ASCAT_GCC_UTIL}" = "" ]; then
        echo "ERROR: ASCAT_GCC_UTIL shell variable with path to 'ascatSnpPanelGcCorrections.pl' tool is not defined. This tool is provided by the AscatNGS package (PERL5LIB variable may also need to be exported)" >&2
        return 1
    else
        echo "* Testing that ${ASCAT_GCC_UTIL} given in ASCAT_GCC_UTIL variable can be executed correctly..." >&2
        ${ASCAT_GCC_UTIL} > ${TMPDIR}/ascat_gcc_util_test.txt 2>&1
        if $GREP "USAGE" ${TMPDIR}/ascat_gcc_util_test.txt > /dev/null; then
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
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

process_pars || exit 1
