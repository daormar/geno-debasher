# Geno-DeBasher package
# Copyright (C) 2019-2024 Daniel Ortiz-Mart\'inez
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

##############################
# GENOME REFERENCE PROCESSES #
##############################

########
create_genref_for_bam_document()
{
    document_process "Creates a genome reference file for a given \`bam\` file. For this purpose, the process starts from a basic genome reference file, removing those contigs not present in the \`bam\` file and downloading or copying missing ones from the Internet or from previously existing files."
}

########
create_genref_for_bam_explain_opts()
{
    # -br option
    description="Base reference genome file"
    explain_opt "-br" "<file>" "$description"

    # -bam option
    description="bam file (required if no downloading processes or paths of normal or tumor bam files have been defined)"
    explain_opt "-bam" "<file>" "$description"

    # -cm option
    description="File containing a mapping between contig names and accession numbers"
    explain_opt "-cm" "<file>" "$description"

    # -fbr option
    description="Name of fallback genome reference file. If creation process fails, this file is copied as reference output file instead"
    explain_opt "-fbr" "<file>" "$description"
}

########
create_genref_for_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-br"
    opt_is_cmdline "-bam"
    opt_is_cmdline "-cm"
    opt_is_cmdline "-fbr"
}

########
create_genref_for_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -br option
    define_cmdline_infile_opt "$cmdline" "-br" optlist || return 1

    # -bam option
    define_cmdline_infile_opt "$cmdline" "-bam" optlist || return 1

    # -cm option
    define_cmdline_infile_opt_if_given "$cmdline" "-cm" optlist || return 1

    # -fbr option
    define_cmdline_infile_opt_if_given "$cmdline" "-fbr" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "data"`

    # -outfile option
    local outfile="${abs_datadir}"/genref.fa
    define_opt "-outfile" "$outfile" optlist || return 1

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
create_seq_dict_for_ref()
{
    local ref=$1

    conda activate gatk4 2>&1 || return 1
    gatk CreateSequenceDictionary -R "${ref}" || return 1
    conda deactivate
}

########
index_ref()
{
    local ref=$1

    conda activate samtools 2>&1 || return 1
    samtools faidx "${ref}" || return 1
    conda deactivate
}

########
create_genref_for_bam()
{
    # Initialize variables
    local baseref=`read_opt_value_from_func_args "-br" "$@"`
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local bam=`read_opt_value_from_func_args "-bam" "$@"`
    local contig_mapping=`read_opt_value_from_func_args "-cm" "$@"`
    if [ "$contig_mapping" = "${DEBASHER_OPT_NOT_FOUND}" ]; then
        contig_mapping=${NOFILE}
    fi
    local fallback_genref=`read_opt_value_from_func_args "-fbr" "$@"`
    if [ "$fallback_genref" = "${DEBASHER_OPT_NOT_FOUND}" ]; then
        fallback_genref=${NOFILE}
    fi
    local outfile=`read_opt_value_from_func_args "-outfile" "$@"`

    # Create genome reference
    local cm_opt=`get_create_genref_for_bam_cm_opt ${contig_mapping}`
    if "${genodebasher_bindir}"/genodb_create_genref_for_bam -r "${baseref}" -b "${bam}" ${cm_opt} -o "${process_outd}"; then
        # Move resulting files
        "${MV}" "${process_outd}"/genref_for_bam.fa "${outfile}"
        "${MV}" "${process_outd}"/genref_for_bam.fa.fai "${outfile}".fai

        # Create sequence dictionary for reference
        logmsg "* Creating sequence dictionary for reference..."
        create_seq_dict_for_ref "${outfile}" || return 1
    else
        # Genome reference creation failed, check if fallback file was
        # provided
        if [ "${fallback_genref}" = ${NOFILE} ]; then
            return 1
        else
            logmsg "Genome reference creation failed but fallback file was provided"
            # Copy fallback file
            logmsg "* Copying fallback file (${fallback_genref})..."
            cp "${fallback_genref}" "${outfile}" || return 1

            # Index fallback reference
            logmsg "* Indexing fallback reference..."
            index_ref "${outfile}" || return 1

            # Create sequence dictionary for fallback reference
            logmsg "* Creating sequence dictionary for fallback reference..."
            create_seq_dict_for_ref "${outfile}" || return 1
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
    [ -f "$file" ] || { errmsg "file $file containing contig list does not exist" ; return 1; }
    cat "$file"
}

########
get_ref_contig_list()
{
    local ref=$1

    if [ ! -f "${ref}".fai ]; then
        conda activate samtools 2>&1 || return 1
        samtools faidx "${ref}"
        conda deactivate
    fi

    "$AWK" '{printf " %s",$1}' "${ref}".fai
}

########
filter_bam_stats()
{
    "${AWK}" '{if($3>0 || $4>0) printf" %s",$1}'
}

########
get_bam_contig_list()
{
    local bam=$1

    conda activate samtools 2>&1 || return 1
    samtools idxstats "$bam" | filter_bam_stats
    conda deactivate
}
