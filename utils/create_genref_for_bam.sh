# Geno-PanPipe package
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

# INCLUDE BASH LIBRARY
. "${PANPIPE_HOME_DIR}"/panpipe_lib || exit 1

########
print_desc()
{
    echo "create_genref_for_bam creates genome reference specific for bam file"
    echo "type \"create_genref_for_bam --help\" to get usage information"
}

########
usage()
{
    echo "create_genref_for_bam  -r <string> -b <string> [-cm <string>]"
    echo "                       -o <string> [--help]"
    echo ""
    echo "-r <string>            File with base reference genome"
    echo "-b <string>            bam file"
    echo "-cm <string>           File containing a mapping from contig names to"
    echo "                       GenBank/NCBI RefSeq accession numbers or file names"
    echo "                       (when the mapping starts with a '/' character, it is"
    echo "                       considered a file, hence, absolute paths should be"
    echo "                       given)"
    echo "-o <string>            Output directory"
    echo "--help                 Display this help and exit"
}

########
read_pars()
{
    r_given=0
    b_given=0
    cm_given=0
    contig_mapping=${NOFILE}
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
            "-cm") shift
                  if [ $# -ne 0 ]; then
                      contig_mapping=$1
                      cm_given=1
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
        if [ ! -f "${baseref}" ]; then
            echo "Error! file ${baseref} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${b_given} -eq 0 ]; then   
        echo "Error! -b parameter not given!" >&2
        exit 1
    else
        if [ ! -f "${bam}" ]; then
            echo "Error! file ${bam} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${cm_given} -eq 1 ]; then   
        if [ ! -f "${contig_mapping}" ]; then
            echo "Error! file ${contig_mapping} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${o_given} -eq 0 ]; then   
        echo "Error! -o parameter not given!" >&2
        exit 1
    else
        if [ ! -d "${outd}" ]; then
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

    if [ ${cm_given} -eq 1 ]; then
        echo "-cm is ${contig_mapping}" >&2
    fi

    if [ ${o_given} -eq 1 ]; then
        echo "-o is ${outd}" >&2
    fi
}

########
contig_in_list()
{
    local contig=$1
    local contiglen=$2
    local clist=$3

    while read cname clen; do
        if [ "$contig" = "$cname" -a "$contiglen" = "$clen" ]; then
            return 0
        fi
    done < "${clist}"

    return 1
}

########
extract_contig_info_from_fai()
{
    local faifile=$1
    $AWK '{printf "%s %s\n",$1,$2}' "${faifile}"
}

########
get_ref_contig_names()
{
    local ref=$1

    samtools faidx "${baseref}" || return 1
    extract_contig_info_from_fai "${baseref}".fai
}

########
get_bam_contig_names()
{
    local bam=$1

    samtools view -H "$bam" | $AWK '{if($1=="@SQ") printf "%s %s\n",substr($2,4),substr($3,4)}'
}

########
get_missing_contig_names()
{
    local refcontigs=$1
    local bamcontigs=$2
            
    while read bamcontigname contiglen; do
        if ! contig_in_list $bamcontigname $contiglen $refcontigs; then
            echo $bamcontigname $contiglen
        fi
    done < "$bamcontigs"
}

########
get_ref_contig_names_to_keep()
{
    local refcontigs=$1
    local bamcontigs=$2
            
    while read refcontigname contiglen; do
        if contig_in_list $refcontigname $contiglen $bamcontigs; then
            echo $refcontigname $contiglen
        fi
    done < "$refcontigs"    
}

########
contig_is_accession()
{
    # Contig is classified as an accession if it is a string containing
    # a dot in the middle
    local contig=$1

    if [[ $contig == *"."* ]]; then
        return 0
    else
        return 1
    fi
}

########
map_contig_with_len_using_file()
{
    local contig_mapping=$1
    local contig=$2
    local contiglen=$3
    local contig_plus_len="${contig}_${contiglen}"
    
    while read entry; do
        local fields=($entry)
        local num_fields=${#fields[@]}
        if [ ${num_fields} -eq 2 ]; then
            if [ ${fields[0]} = ${contig_plus_len} ]; then
                echo ${fields[1]}
                break
            fi
        fi
    done < "${contig_mapping}"
}

########
map_contig_without_len_using_file()
{
    local contig_mapping=$1
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
    done < "${contig_mapping}"
}

########
map_contig_using_file()
{
    local contig_mapping=$1
    local contig=$2
    local contiglen=$3

    # Try to map contig taking into account contig length
    mapping=`map_contig_with_len_using_file "${contig_mapping}" ${contig} ${contiglen}`
    if [ "${mapping}" != "" ]; then
        echo ${mapping}
    else
        # Try to map contig without taking into account contig length
        mapping=`map_contig_without_len_using_file "${contig_mapping}" ${contig}`
        if [ "${mapping}" != "" ]; then
            echo ${mapping}
        fi
    fi    
}

########
map_contig()
{
    local contig_mapping=$1
    local contig=$2
    local contiglen=$3

    if contig_is_accession ${contig}; then
        echo ${contig}
    else
        if [ "${contig_mapping}" != "${NOFILE}" ]; then
            map_contig_using_file "${contig_mapping}" ${contig} ${contiglen} || return 1
        fi
    fi
}

########
replace_contig_name()
{
    local mapping=$1
    local contig=$2

    ${AWK} -v mapping=${mapping} -v contig=${contig} '{
             if(FNR==1)
             {
              printf ">%s",contig
              for(i=2;i<=NF;++i)
               printf" %s",$i
              printf"\n"
             }
             else print $0
            }'
}

########
get_contigs()
{
    local contig_mapping=$1
    local contiglist=$2

    while read contig contiglen; do
        local mapping=`map_contig "${contig_mapping}" ${contig} ${contiglen}` || return 1
        if [ "$mapping" = "" ]; then
            echo "Error: contig $contig is not a valid accession nor there were mappings for it" >&2
            return 1
        else
            # Determine whether the mapping is an accession number or a
            # file name (absolute file paths should be given)
            if is_absolute_path "${mapping}"; then
                echo "Getting data for contig ${contig} with length ${contiglen} (mapped to file $mapping)..." >&2
                cat "${mapping}" || return 1
            else
                echo "Getting data for contig ${contig} with length ${contiglen} (mapped to accession $mapping)..." >&2
                "${genopanpipe_bindir}"/get_entrez_fasta -a "${mapping}" | replace_contig_name "${mapping}" ${contig}; pipe_fail || return 1
            fi
        fi
    done < ${contiglist}
}

########
get_uniq_contigs()
{
    local cfile1=$1
    local cfile2=$2

    "${SORT}" "$cfile1" "$cfile2" | "${UNIQ}" -u
}

########
process_pars()
{
    outfile=$outd/genref_for_bam.fa
    
    # Activate conda environment
    echo "* Activating conda environment (samtools)..." >&2
    conda activate samtools || return 1

    # Get reference contigs
    echo "* Obtaining list of current reference contig names and their lengths..." >&2
    get_ref_contig_names "$baseref" > "${outd}"/refcontigs || return 1

    # Get bam contigs
    echo "* Obtaining list of bam contig names and their lengths..." >&2
    get_bam_contig_names "$bam" > "${outd}"/bamcontigs || return 1
    
    # Obtain list of contigs to keep in the reference file
    echo "* Obtaining list of reference contigs to keep..." >&2
    get_ref_contig_names_to_keep "${outd}"/refcontigs "${outd}"/bamcontigs > "${outd}"/refcontigs_to_keep || return 1

    # Copy base genome reference without extra contigs
    echo "* Copying base genome reference without extra contigs..." >&2
    "${genopanpipe_bindir}"/filter_contig_from_genref -g "$baseref" -l ${outd}/refcontigs_to_keep > "${outd}"/unordered_ref.fa || return 1

    # Obtain list of missing contigs
    echo "* Obtaining list of missing contigs..." >&2
    get_missing_contig_names "${outd}"/refcontigs_to_keep "${outd}"/bamcontigs > "${outd}"/missing_contigs || return 1

    # Enrich reference
    echo "* Enriching reference..." >&2
    get_contigs "${contig_mapping}" "${outd}"/missing_contigs >> "${outd}"/unordered_ref.fa || { echo "Error during reference enrichment" >&2; return 1; }

    # Reorder contigs
    echo "* Reordering reference contigs..." >&2
    "${genopanpipe_bindir}"/reorder_fa_seqs -f "${outd}"/unordered_ref.fa -l "${outd}"/bamcontigs > "$outfile" || { echo "Error during contig reordering" >&2; return 1; }
    rm "${outd}"/unordered_ref.fa || return 1
    
    # Index created reference
    echo "* Indexing created reference..." >&2
    samtools faidx "${outfile}" || return 1
    extract_contig_info_from_fai "${outfile}".fai > "${outd}"/created_ref_contigs || return 1

    # Check created reference
    echo "* Checking created reference..." >&2
    get_uniq_contigs "${outd}"/bamcontigs "${outd}"/created_ref_contigs > "${outd}"/uniq_contigs
    num_uniq_contigs=`$WC -l ${outd}/uniq_contigs | $AWK '{print $1}'`
    if [ ${num_uniq_contigs} -gt 0 ]; then
        echo "Bam file and created genome reference do not have the exact same contigs (see ${outd}/uniq_contigs file)" >&2
        return 1
    fi
    
    # Deactivate conda environment
    echo "* Deactivating conda environment..." >&2
    conda deactivate
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars "$@" || exit 1

check_pars || exit 1

print_pars || exit 1

process_pars || exit 1
