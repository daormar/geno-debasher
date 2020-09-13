# Bio-PanPipe package
# Copyright (C) 2019,2020 Daniel Ortiz-Mart\'inez
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
    echo "reorder_fa_seqs_lowmem reorder sequences in fasta file"
    echo "type \"reorder_fa_seqs_lowmem --help\" to get usage information"
}

########
usage()
{
    echo "reorder_fa_seqs_lowmem     -f <string> -l <string>"
    echo "                           [--help]"
    echo ""
    echo "-f <string>                Fasta file"
    echo "-l <string>                File with list of sequence names"
    echo "--help                     Display this help and exit"
}

########
read_pars()
{
    f_given=0
    l_given=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "-f") shift
                  if [ $# -ne 0 ]; then
                      fasta=$1
                      f_given=1
                  fi
                  ;;
            "-l") shift
                  if [ $# -ne 0 ]; then
                      listseq=$1
                      l_given=1
                  fi
                  ;;
        esac
        shift
    done   
}

########
check_pars()
{
    if [ ${f_given} -eq 0 ]; then   
        echo "Error! -f parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${fasta} ]; then
            echo "Error! file ${fasta} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${l_given} -eq 0 ]; then   
        echo "Error! -l parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${listseq} ]; then
            echo "Error! file ${listseq} does not exist" >&2
            exit 1
        fi
    fi
}

########
process_pars()
{
    while read seqname; do

        ${biopanpipe_bindir}/filter_contig_from_genref -g ${fasta} -k ${seqname}
        
    done < $listseq
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

process_pars || exit 1
