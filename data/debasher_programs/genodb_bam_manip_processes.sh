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

#############
# CONSTANTS #
#############

DEFAULT_MAX_RECORDS_IN_RAM_GATK=1000000

##############################
# BAM MANIPULATION PROCESSES #
##############################

########
genodb.index_norm_bam_explain_opts()
{
    # -normalbam option
    local description="normal bam file"
    explain_opt "-normalbam" "<file>" "$description"

    # -out-nbidx option
    local description="output file"
    explain_opt "-out-nbidx" "<file>" "$description"

    # -out-nb option
    local description="indexed normal bam file"
    explain_opt "-out-nb" "<file>" "$description"
}

########
genodb.index_norm_bam_identify_cmdline_opts()
{
    :
}

########
genodb.index_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # -normalbam option
    local abs_datadir=$(get_absolute_shdirname "data")
    local normalbam="${abs_datadir}"/normal.bam
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -out-nbidx option
    define_opt "-out-nbidx" "$normalbam".bai optlist || return 1

    # -out-nb option (republished once indexed, so downstream processes
    # can connect to it and depend on indexing having completed)
    define_opt "-out-nb" "$normalbam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.index_norm_bam()
{
    # Initialize variables
    local normalbam=$(read_opt_value_from_func_args "-normalbam" "$@")
    local idx_file=$(read_opt_value_from_func_args "-out-nbidx" "$@")

    # Remove previous index if one was created
    if [ -f "${idx_file}" ]; then
        "${RM}" "${idx_file}" || return 1
    fi

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda_activate samtools 2>&1 || return 1

    # Execute samtools
    logmsg "* Executing samtools index..."
    samtools index "${normalbam}" "${idx_file}" 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.index_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
genodb.index_tum_bam_explain_opts()
{
    # -tumorbam option
    local description="tumor bam file"
    explain_opt "-tumorbam" "<file>" "$description"

    # -out-tbidx option
    local description="output file"
    explain_opt "-out-tbidx" "<file>" "$description"

    # -out-tb option
    local description="indexed tumor bam file"
    explain_opt "-out-tb" "<file>" "$description"
}

########
index_tum_bam_cmdline_opts()
{
    :
}

########
genodb.index_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # -tumorbam option
    local abs_datadir=$(get_absolute_shdirname "data")
    local tumorbam="${abs_datadir}"/tumor.bam
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -out-tbidx option
    define_opt "-out-tbidx" "$tumorbam".bai optlist || return 1

    # -out-tb option (republished once indexed, so downstream processes
    # can connect to it and depend on indexing having completed)
    define_opt "-out-tb" "$tumorbam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.index_tum_bam()
{
    # Initialize variables
    local tumorbam=$(read_opt_value_from_func_args "-tumorbam" "$@")
    local idx_file=$(read_opt_value_from_func_args "-out-tbidx" "$@")

    # Remove previous index if one was created
    if [ -f "${idx_file}" ]; then
        "${RM}" "${idx_file}" || return 1
    fi

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda_activate samtools 2>&1 || return 1

    # Execute samtools
    logmsg "* Executing samtools index..."
    samtools index "${tumorbam}" "${idx_file}" 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.index_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
genodb.sort_norm_bam_explain_opts()
{
    # -normalbam option
    local description="normal bam file"
    explain_opt "-normalbam" "<file>" "$description"

    # -process-outd option
    local description="output directory"
    explain_opt "-process-outd" "<file>" "$description"

    # -cpus option
    local description="number of cpus"
    explain_opt "-cpus" "<int>" "$description"
}

########
genodb.sort_norm_bam_identify_cmdline_opts()
{
    :
}

########
genodb.sort_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    define_opt "-process-outd" "${process_outdir}" optlist || return 1

    # -normalbam option
    local abs_datadir=$(get_absolute_shdirname "data")
    local normalbam="${abs_datadir}"/normal.bam
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.sort_norm_bam()
{
    # Initialize variables
    local normalbam=$(read_opt_value_from_func_args "-normalbam" "$@")
    local process_outd=$(read_opt_value_from_func_args "-process-outd" "$@")
    local cpus=$(read_opt_value_from_func_args "-cpus" "$@")

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda_activate samtools 2>&1 || return 1

    # Verify if bam file is already sorted
    local bam_is_sorted=$(samtools view -H "${normalbam}" | $GREP SO:coordinate | wc -l) || return 1
    if [ ${bam_is_sorted} -eq 1 ]; then
        echo "Warning: bam file is already sorted"
    else
        # Execute samtools
        logmsg "* Executing samtools sort..."
        samtools sort -T "${process_outd}" -o "${process_outd}"/sorted.bam -m 2G -@ ${cpus} "${normalbam}" 2>&1 || return 1
        # NOTE: -m option is used here to increase the maximum memory per
        # thread. One lateral efect of this is that the number of tmp files
        # generated is decreased. This constitutes one possible way to avoid
        # the "Too many open files" error reported by samtools

        # Replace initial bam file by the sorted one
        "${MV}" "${process_outd}"/sorted.bam "${normalbam}" 2>&1 || return 1
    fi

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.sort_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
genodb.sort_tum_bam_explain_opts()
{
    # -tumorbam option
    local description="tumor bam file"
    explain_opt "-tumorbam" "<file>" "$description"

    # -process-outd option
    local description="output directory"
    explain_opt "-process-outd" "<file>" "$description"

    # -cpus option
    local description="number of cpus"
    explain_opt "-cpus" "<int>" "$description"
}

########
genodb.sort_tum_bam_identify_cmdline_opts()
{
    :
}

########
genodb.sort_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    define_opt "-process-outd" "${process_outdir}" optlist || return 1

    # -tumorbam option
    local abs_datadir=$(get_absolute_shdirname "data")
    local tumorbam="${abs_datadir}"/tumor.bam
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.sort_tum_bam()
{
    # Initialize variables
    local tumorbam=$(read_opt_value_from_func_args "-tumorbam" "$@")
    local process_outd=$(read_opt_value_from_func_args "-process-outd" "$@")
    local cpus=$(read_opt_value_from_func_args "-cpus" "$@")

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda_activate samtools 2>&1 || return 1

    # Verify if bam file is already sorted
    local bam_is_sorted=$(samtools view -H "${tumorbam}" | $GREP SO:coordinate | wc -l) || return 1
    if [ ${bam_is_sorted} -eq 1 ]; then
        echo "Warning: bam file is already sorted"
    else
        # Execute samtools
        logmsg "* Executing samtools sort..."
        samtools sort -T "${process_outd}" -o "${process_outd}"/sorted.bam -m 2G -@ ${cpus} "${tumorbam}" 2>&1 || return 1
        # NOTE: -m option is used here to increase the maximum memory per
        # thread. One lateral efect of this is that the number of tmp files
        # generated is decreased. This constitutes one possible way to avoid
        # the "Too many open files" error reported by samtools

        # Replace initial bam file by the sorted one
        "${MV}" "${process_outd}"/sorted.bam "${tumorbam}" 2>&1 || return 1
    fi

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.sort_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
genodb.samtools_mpileup_norm_bam_explain_opts()
{
    # -r option
    description="Reference genome file"
    explain_opt "-r" "<file>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_opt "-mpb" "<file>" "$description"

    # -normalbam option
    local description="normal bam file"
    explain_opt "-normalbam" "<file>" "$description"

    # -process-outd option
    local description="output directory"
    explain_opt "-process-outd" "<file>" "$description"

    # -outfile option
    local description="output pileup file"
    explain_opt "-outfile" "<file>" "$description"
}

########
genodb.samtools_mpileup_norm_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
    opt_is_cmdline "-normalbam"
    opt_is_non_mandatory_cmdline "-mpb"
}

########
genodb.samtools_mpileup_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    define_opt "-process-outd" "${process_outdir}" optlist || return 1

    # -r option
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # -normalbam option
    define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # -outfile option
    local outfile="${process_outdir}"/normal.pileup.gz
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
get_samtools_mpileup_l_opt()
{
    local mpbfile=$1

    if [ "${mpbfile}" = "${DEBASHER_OPT_NOT_FOUND}" ]; then
        echo ""
    else
        echo "-l ${mpbfile}"
    fi
}

########
genodb.samtools_mpileup_norm_bam()
{
    # Initialize variables
    local process_outd=$(read_opt_value_from_func_args "-process-outd" "$@")
    local ref=$(read_opt_value_from_func_args "-r" "$@")
    local normalbam=$(read_opt_value_from_func_args "-normalbam" "$@")
    local mpbfile=$(read_opt_value_from_func_args "-mpb" "$@")

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda_activate samtools 2>&1 || return 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=$(get_samtools_mpileup_l_opt ${mpbfile})

    # Generate pileup file
    logmsg "* Generating pileup file..."
    samtools mpileup ${smp_l_opt} -f "${ref}" -o "${process_outd}"/normal.pileup "$normalbam" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    "${GZIP}" "${process_outd}"/normal.pileup || return 1
}

########
genodb.samtools_mpileup_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
genodb.samtools_mpileup_tum_bam_explain_opts()
{
    # -r option
    description="Reference genome file"
    explain_opt "-r" "<file>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_opt "-mpb" "<file>" "$description"

    # -tumorbam option
    local description="tumor bam file"
    explain_opt "-tumorbam" "<file>" "$description"

    # -process-outd option
    local description="output directory"
    explain_opt "-process-outd" "<file>" "$description"

    # -outfile option
    local description="output pileup file"
    explain_opt "-outfile" "<file>" "$description"
}

########
genodb.samtools_mpileup_tum_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
    opt_is_cmdline "-tumorbam"
    opt_is_non_mandatory_cmdline "-mpb"
}

########
genodb.samtools_mpileup_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    define_opt "-process-outd" "${process_outdir}" optlist || return 1

    # -r option
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # -tumorbam option
    define_cmdline_infile_opt "$cmdline" "-tumorbam" optlist || return 1

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # -outfile option
    local outfile="${process_outdir}"/tumor.pileup.gz
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.samtools_mpileup_tum_bam()
{
    # Initialize variables
    local process_outd=$(read_opt_value_from_func_args "-process-outd" "$@")
    local ref=$(read_opt_value_from_func_args "-r" "$@")
    local tumorbam=$(read_opt_value_from_func_args "-tumorbam" "$@")
    local mpbfile=$(read_opt_value_from_func_args "-mpb" "$@")

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda_activate samtools 2>&1 || return 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=$(get_samtools_mpileup_l_opt ${mpbfile})

    # Generate pileup file
    logmsg "* Generating pileup file..."
    samtools mpileup ${smp_l_opt} -f "${ref}" -o "${process_outd}"/tumor.pileup "$tumorbam" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    "${GZIP}" "${process_outd}"/tumor.pileup || return 1
}

########
genodb.samtools_mpileup_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
genodb.parallel_samtools_mpileup_norm_bam_explain_opts()
{
    # -r option
    description="Reference genome file"
    explain_opt "-r" "<file>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_opt "-mpb" "<file>" "$description"

    # -normalbam option
    local description="normal bam file of the contig"
    explain_opt "-normalbam" "<file>" "$description"

    # -outd option
    local description="directory of the task, below the output directory of the process"
    explain_opt "-outd" "<file>" "$description"

    # -outfile option
    local description="output file"
    explain_opt "-outfile" "<file>" "$description"

    # -contig option
    local description="contig id"
    explain_opt "-contig" "<string>" "$description"
}

########
genodb.parallel_samtools_mpileup_norm_bam_explain_task_shaping_opts()
{
    # -lc option
    local description="File with list of contig names to process, one task per contig"
    explain_task_shaping_opt "-lc" "<file>" "$description"
}

########
genodb.parallel_samtools_mpileup_norm_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
    opt_is_non_mandatory_cmdline "-mpb"
}

########
genodb.parallel_samtools_mpileup_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -outd option, a directory of the task's own below the output
        # directory of the process, which the engine creates and empties
        # before the task runs
        define_opt_from_process_outdir "-outd" optlist --subdir "${idx}" || return 1

        # -r option
        define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

        # -mpb option
        define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

        # -normalbam option
        local contig=${array[$idx]}
        define_opt_from_proc_task_out "-normalbam" "genodb.parallel_split_norm_bam" "${idx}" "-outfile" optlist || return 1

        # -contig option
        define_opt "-contig" "$contig" optlist || return 1

        # -outfile option, under a fixed name in the directory of the task
        local outfile="${process_outdir}"/${idx}/normal.pileup.gz
        define_opt "-outfile" "$outfile" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.parallel_samtools_mpileup_norm_bam()
{
    # Initialize variables
    local outd=$(read_opt_value_from_func_args "-outd" "$@")
    local ref=$(read_opt_value_from_func_args "-r" "$@")
    local normalbam=$(read_opt_value_from_func_args "-normalbam" "$@")
    local mpbfile=$(read_opt_value_from_func_args "-mpb" "$@")
    local contig=$(read_opt_value_from_func_args "-contig" "$@")
    local outfile=$(read_opt_value_from_func_args "-outfile" "$@")

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda_activate samtools 2>&1 || return 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=$(get_samtools_mpileup_l_opt ${mpbfile})

    # Generate pileup file
    logmsg "* Generating pileup file (contig $contig)..."
    samtools mpileup ${smp_l_opt} -f "${ref}" -o "${outd}"/normal.pileup "$normalbam" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    "${GZIP}" -c "${outd}"/normal.pileup > "${outfile}" || return 1
}

########
genodb.parallel_samtools_mpileup_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
genodb.parallel_samtools_mpileup_tum_bam_explain_opts()
{
    # -r option
    description="Reference genome file"
    explain_opt "-r" "<file>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_opt "-mpb" "<file>" "$description"

    # -tumorbam option
    local description="tumor bam file of the contig"
    explain_opt "-tumorbam" "<file>" "$description"

    # -outd option
    local description="directory of the task, below the output directory of the process"
    explain_opt "-outd" "<file>" "$description"

    # -outfile option
    local description="output file"
    explain_opt "-outfile" "<file>" "$description"

    # -contig option
    local description="contig id"
    explain_opt "-contig" "<string>" "$description"
}

########
genodb.parallel_samtools_mpileup_tum_bam_explain_task_shaping_opts()
{
    # -lc option
    local description="File with list of contig names to process, one task per contig"
    explain_task_shaping_opt "-lc" "<file>" "$description"
}

########
genodb.parallel_samtools_mpileup_tum_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
    opt_is_non_mandatory_cmdline "-mpb"
}

########
genodb.parallel_samtools_mpileup_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -outd option, a directory of the task's own below the output
        # directory of the process, which the engine creates and empties
        # before the task runs
        define_opt_from_process_outdir "-outd" optlist --subdir "${idx}" || return 1

        # -r option
        define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

        # -mpb option
        define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

        # -tumorbam option
        local contig=${array[$idx]}
        define_opt_from_proc_task_out "-tumorbam" "genodb.parallel_split_tum_bam" "${idx}" "-outfile" optlist || return 1

        # -contig option
        define_opt "-contig" "$contig" optlist || return 1

        # -outfile option, under a fixed name in the directory of the task
        local outfile="${process_outdir}"/${idx}/tumor.pileup.gz
        define_opt "-outfile" "$outfile" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.parallel_samtools_mpileup_tum_bam()
{
    # Initialize variables
    local outd=$(read_opt_value_from_func_args "-outd" "$@")
    local ref=$(read_opt_value_from_func_args "-r" "$@")
    local tumorbam=$(read_opt_value_from_func_args "-tumorbam" "$@")
    local mpbfile=$(read_opt_value_from_func_args "-mpb" "$@")
    local contig=$(read_opt_value_from_func_args "-contig" "$@")
    local outfile=$(read_opt_value_from_func_args "-outfile" "$@")

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda_activate samtools 2>&1 || return 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=$(get_samtools_mpileup_l_opt ${mpbfile})

    # Generate pileup file
    logmsg "* Generating pileup file (contig $contig)..."
    samtools mpileup ${smp_l_opt} -f "${ref}" -o "${outd}"/tumor.pileup "${tumorbam}" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    "${GZIP}" -c "${outd}"/tumor.pileup > "${outfile}" || return 1
}

########
genodb.parallel_samtools_mpileup_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
filter_bam_contig_samtools()
{
    # Writes into outbam the reads of inbam that map to contig, with the
    # header of inbam. A contig that the header of inbam does not list
    # gives a warning and an outbam with the header and no reads: the
    # normal and the tumor bam files share the contig list, and one of
    # them may lack a contig that the other has
    local inbam=$1
    local contig=$2
    local outbam=$3

    local header
    header=$(samtools view -H "${inbam}") || { errmsg "Error while reading the header of ${inbam}"; return 1; }
    "${GREP}" -q "^@SQ" <<< "${header}" || { errmsg "Error: the header of ${inbam} lists no contig"; return 1; }

    if "${AWK}" -v contig="${contig}" '$1 == "@SQ" { for (i = 2; i <= NF; i++) if ($i == "SN:" contig) found = 1 } END { exit !found }' <<< "${header}"; then
        samtools view -h -O BAM "${inbam}" "${contig}" > "${outbam}" || { errmsg "Error while extracting contig ${contig} from ${inbam}"; return 1; }
    else
        errmsg "Warning: contig ${contig} does not exist in ${inbam}, the output has no reads"
        samtools view -H -O BAM "${inbam}" > "${outbam}" || { errmsg "Error while writing the header of ${inbam}"; return 1; }
    fi
}

########
genodb.parallel_split_norm_bam_explain_opts()
{
    # -normalbam option
    local description="normal bam file"
    explain_opt "-normalbam" "<file>" "$description"

    # -outfile option
    local description="output file"
    explain_opt "-outfile" "<file>" "$description"

    # -contig option
    local description="contig id"
    explain_opt "-contig" "<string>" "$description"
}

########
genodb.parallel_split_norm_bam_explain_task_shaping_opts()
{
    # -lc option
    local description="File with list of contig names to process, one task per contig"
    explain_task_shaping_opt "-lc" "<file>" "$description"
}

########
genodb.parallel_split_norm_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
}

########
genodb.parallel_split_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -normalbam option
        define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

        # -contig option
        local contig=${array[$idx]}
        define_opt "-contig" "$contig" optlist || return 1

        # -outfile option, an output of the task's own in the output
        # directory of the process
        local outfile="${process_outdir}"/normal_${contig}.bam
        define_opt "-outfile" "$outfile" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.parallel_split_norm_bam()
{
    # Initialize variables
    local normalbam=$(read_opt_value_from_func_args "-normalbam" "$@")
    local contig=$(read_opt_value_from_func_args "-contig" "$@")
    local outfile=$(read_opt_value_from_func_args "-outfile" "$@")

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda_activate samtools 2>&1 || return 1

    # Extract contig
    logmsg "* Extracting contig (contig $contig)..."
    filter_bam_contig_samtools "$normalbam" $contig "$outfile" || return 1

    # Index contig
    logmsg "* Indexing contig..."
    samtools index "${outfile}" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.parallel_split_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
genodb.parallel_split_tum_bam_explain_opts()
{
    # -tumorbam option
    local description="tumor bam file"
    explain_opt "-tumorbam" "<file>" "$description"

    # -outfile option
    local description="output file"
    explain_opt "-outfile" "<file>" "$description"

    # -contig option
    local description="contig id"
    explain_opt "-contig" "<string>" "$description"
}

########
genodb.parallel_split_tum_bam_explain_task_shaping_opts()
{
    # -lc option
    local description="File with list of contig names to process, one task per contig"
    explain_task_shaping_opt "-lc" "<file>" "$description"
}

########
genodb.parallel_split_tum_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-tumorbam"
}

########
genodb.parallel_split_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4

    # Get name of contig list file
    local clist
    clist=$(read_opt_value_from_line "$cmdline" "-lc") || { errmsg "Error: -lc option not found"; return 1; }

    # Array of contigs to process, one task per contig
    array=( $(get_contig_list_from_file $clist) ) || return 1

    for idx in "${!array[@]}"; do
        local optlist=""

        # -tumorbam option
        define_cmdline_infile_opt "$cmdline" "-tumorbam" optlist || return 1

        # -contig option
        local contig=${array[$idx]}
        define_opt "-contig" "$contig" optlist || return 1

        # -outfile option, an output of the task's own in the output
        # directory of the process
        local outfile="${process_outdir}"/tumor_${contig}.bam
        define_opt "-outfile" "$outfile" optlist || return 1

        # Save option list
        save_opt_list optlist
    done
}

########
genodb.parallel_split_tum_bam()
{
    # Initialize variables
    local tumorbam=$(read_opt_value_from_func_args "-tumorbam" "$@")
    local contig=$(read_opt_value_from_func_args "-contig" "$@")
    local outfile=$(read_opt_value_from_func_args "-outfile" "$@")

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda_activate samtools 2>&1 || return 1

    # Extract contig
    logmsg "* Extracting contig (contig $contig)..."
    filter_bam_contig_samtools "$tumorbam" $contig "$outfile" || return 1

    # Index contig
    logmsg "* Indexing contig..."
    samtools index ${outfile} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.parallel_split_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
genodb.bedtools_genomecov_norm_bam_explain_opts()
{
    # -normalbam option
    local description="normal bam file"
    explain_opt "-normalbam" "<file>" "$description"

    # -out-processdir option
    local description="output directory"
    explain_opt "-out-processdir" "<file>" "$description"
}

########
genodb.bedtools_genomecov_norm_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
}

########
genodb.bedtools_genomecov_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -normalbam option
    define_cmdline_infile_opt "$cmdline" "-normalbam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.bedtools_genomecov_norm_bam()
{
    # Initialize variables
    local process_outd=$(read_opt_value_from_func_args "-out-processdir" "$@")
    local normalbam=$(read_opt_value_from_func_args "-normalbam" "$@")

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda_activate bedtools 2>&1 || return 1

    # Execute samtools
    logmsg "* Executing bedtools coverage..."
    bedtools genomecov -ibam "${normalbam}" > "${process_outd}"/genomecov.tsv || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.bedtools_genomecov_norm_bam_conda_envs()
{
    define_conda_env bedtools bedtools.yml
}

########
genodb.bedtools_genomecov_tum_bam_explain_opts()
{
    # -tumorbam option
    local description="tumor bam file"
    explain_opt "-tumorbam" "<file>" "$description"

    # -out-processdir option
    local description="output directory"
    explain_opt "-out-processdir" "<file>" "$description"
}

########
genodb.bedtools_genomecov_tum_bam_identify_cmdline_opts()
{
    opt_is_cmdline "-tumorbam"
}

########
genodb.bedtools_genomecov_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -tumorbam option
    define_cmdline_infile_opt "$cmdline" "-tumorbam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.bedtools_genomecov_tum_bam()
{
    # Initialize variables
    local process_outd=$(read_opt_value_from_func_args "-out-processdir" "$@")
    local tumorbam=$(read_opt_value_from_func_args "-tumorbam" "$@")

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda_activate bedtools 2>&1 || return 1

    # Execute samtools
    logmsg "* Executing bedtools coverage..."
    bedtools genomecov -ibam "${tumorbam}" > "${process_outd}"/genomecov.tsv || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.bedtools_genomecov_tum_bam_conda_envs()
{
    define_conda_env bedtools bedtools.yml
}

########
genodb.norm_bam_to_ubam_explain_opts()
{
    # -normalbam option
    description="Normal unmapped bam file (required if no downloading processes have been defined)"
    explain_opt "-normalbam" "<file>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_opt "-mrec" "<int>" "$description"

    # -out-processdir option
    local description="output directory"
    explain_opt "-out-processdir" "<file>" "$description"

    # -outfile option
    local description="output file"
    explain_opt "-outfile" "<file>" "$description"
}

########
genodb.norm_bam_to_ubam_identify_cmdline_opts()
{
    opt_is_cmdline "-normalbam"
    opt_is_non_mandatory_cmdline "-mrec"
}

########
genodb.norm_bam_to_ubam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # Get data directory
    local abs_datadir=$(get_absolute_shdirname "data")

    # -normalbam option
    local normalbam="${abs_datadir}"/normal.bam
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -mrec option
    define_cmdline_opt_if_given "$cmdline" "-mrec" optlist || return 1

    # -outfile option
    local outfile="${abs_datadir}"/normal_unmapped.bam
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.norm_bam_to_ubam()
{
    # Initialize variables
    local process_outd=$(read_opt_value_from_func_args "-out-processdir" "$@")
    local normalbam=$(read_opt_value_from_func_args "-normalbam" "$@")
    local max_records=$(read_opt_value_from_func_args "-mrec" "$@")
    if [ "$max_records" = "${DEBASHER_OPT_NOT_FOUND}" ]; then
        max_records=${DEFAULT_MAX_RECORDS_IN_RAM_GATK}
    fi
    local outfile=$(read_opt_value_from_func_args "-outfile" "$@")

    # Create tmpdir for gatk
    tmpdir="${process_outd}"/tmp
    mkdir "${tmpdir}" || return 1

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda_activate gatk4 2>&1 || return 1

    # Execute gatk RevertSam
    logmsg "* Executing gatk RevertSam..."
    gatk --java-options "-Xmx4G" RevertSam --INPUT "${normalbam}" --OUTPUT "${outfile}" --SANITIZE true --SORT_ORDER queryname --TMP_DIR "${tmpdir}" --MAX_RECORDS_IN_RAM ${max_records} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.norm_bam_to_ubam_conda_envs()
{
    define_conda_env gatk4 gatk4.yml
}

########
genodb.tum_bam_to_ubam_explain_opts()
{
    # -tumorbam option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_opt "-tumorbam" "<file>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_opt "-mrec" "<int>" "$description"

    # -out-processdir option
    local description="output directory"
    explain_opt "-out-processdir" "<file>" "$description"

    # -outfile option
    local description="output file"
    explain_opt "-outfile" "<file>" "$description"
}

########
genodb.tum_bam_to_ubam_identify_cmdline_opts()
{
    opt_is_cmdline "-tumorbam"
    opt_is_non_mandatory_cmdline "-mrec"
}

########
genodb.tum_bam_to_ubam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # Get data directory
    local abs_datadir=$(get_absolute_shdirname "data")

    # -tumorbam option
    local tumorbam="${abs_datadir}"/tumor.bam
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -mrec option
    define_cmdline_opt_if_given "$cmdline" "-mrec" optlist || return 1

    # -outfile option
    local outfile="${abs_datadir}"/tumor_unmapped.bam
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.tum_bam_to_ubam()
{
    # Initialize variables
    local process_outd=$(read_opt_value_from_func_args "-out-processdir" "$@")
    local tumorbam=$(read_opt_value_from_func_args "-tumorbam" "$@")
    local max_records=$(read_opt_value_from_func_args "-mrec" "$@")
    if [ "$max_records" = "${DEBASHER_OPT_NOT_FOUND}" ]; then
        max_records=${DEFAULT_MAX_RECORDS_IN_RAM_GATK}
    fi
    local outfile=$(read_opt_value_from_func_args "-outfile" "$@")

    # Create tmpdir for gatk
    tmpdir="${process_outd}"/tmp
    mkdir "${tmpdir}" || return 1

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda_activate gatk4 2>&1 || return 1

    # Execute gatk RevertSam
    logmsg "* Executing gatk RevertSam..."
    gatk --java-options "-Xmx4G" RevertSam --INPUT "${tumorbam}" --OUTPUT "${outfile}" --SANITIZE true --SORT_ORDER queryname --TMP_DIR "${tmpdir}" --MAX_RECORDS_IN_RAM ${max_records} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.tum_bam_to_ubam_conda_envs()
{
    define_conda_env gatk4 gatk4.yml
}

########
genodb.align_norm_ubam_explain_opts()
{
    # -r option
    description="Reference genome file"
    explain_opt "-r" "<file>" "$description"

    # -normalbam option
    description="Normal unmapped bam file (required if no downloading processes have been defined)"
    explain_opt "-normalbam" "<file>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_opt "-mrec" "<int>" "$description"

    # -out-processdir option
    local description="output directory"
    explain_opt "-out-processdir" "<file>" "$description"

    # -outfile option
    local description="output file"
    explain_opt "-outfile" "<file>" "$description"

    # -cpus option
    local description="number of cpus"
    explain_opt "-cpus" "<int>" "$description"
}

########
genodb.align_norm_ubam_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
    opt_is_cmdline "-normalbam"
    opt_is_non_mandatory_cmdline "-mrec"
}

########
genodb.align_norm_ubam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -r option
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # Get data directory
    local abs_datadir=$(get_absolute_shdirname "data")

    # -normalbam option
    local normalbam="${abs_datadir}"/normal_unmapped.bam
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -mrec option
    define_cmdline_opt_if_given "$cmdline" "-mrec" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # -outfile option
    local outfile="${abs_datadir}"/normal_aligned.bam
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
gatk_dict_exists()
{
    local ref=$1
    local ref_wo_ext="${ref%.*}"

    if [ -f "${ref}" -a -f "${ref_wo_ext}.dict" ]; then
        return 0
    else
        return 1
    fi
}

########
genodb.align_norm_ubam()
{
    # Initialize variables
    local process_outd=$(read_opt_value_from_func_args "-out-processdir" "$@")
    local normalbam=$(read_opt_value_from_func_args "-normalbam" "$@")
    local ref=$(read_opt_value_from_func_args "-r" "$@")
    local max_records=$(read_opt_value_from_func_args "-mrec" "$@")
    if [ "$max_records" = "${DEBASHER_OPT_NOT_FOUND}" ]; then
        max_records=${DEFAULT_MAX_RECORDS_IN_RAM_GATK}
    fi
    local cpus=$(read_opt_value_from_func_args "-cpus" "$@")
    local outfile=$(read_opt_value_from_func_args "-outfile" "$@")

    # Create tmpdir for gatk
    tmpdir="${process_outd}"/tmp
    mkdir "${tmpdir}" || return 1

    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda_activate gatk4 2>&1 || return 1

    # Execute gatk SamToFastq
    logmsg "* Executing gatk SamToFastq..."
    gatk --java-options "-Xmx4G" SamToFastq --INPUT "${normalbam}" --FASTQ "${process_outd}"/reads_r1.fastq.gz --SECOND_END_FASTQ "${process_outd}"/reads_r2.fastq.gz --TMP_DIR "${tmpdir}" --MAX_RECORDS_IN_RAM ${max_records} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bwa)..."
    conda_activate bwa 2>&1 || return 1

    # Execute bwa
    logmsg "* Executing bwa index..."
    bwa index -a bwtsw "${ref}"

    logmsg "* Executing bwa mem..."
    bwa mem -t ${cpus} "${ref}" <("${GZIP}" -d -c "${process_outd}"/reads_r1.fastq.gz) <("${GZIP}" -d -c "${process_outd}"/reads_r2.fastq.gz) | "${GZIP}" > "${process_outd}"/aln.sam.gz ; pipe_fail || return 1

    # Remove fastq files
    "${RM}" "${process_outd}"/reads_r1.fastq.gz "${process_outd}"/reads_r2.fastq.gz || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda_activate gatk4 2>&1 || return 1

    # Create dictionary
    if ! gatk_dict_exists "${ref}"; then
        logmsg "* Creating dictionary for reference..."
        gatk --java-options "-Xmx4G" CreateSequenceDictionary -I "${ref}"
    fi

    # Execute gatk
    logmsg "* Executing gatk CreateSequenceDictionary..."
    gatk --java-options "-Xmx4G" MergeBamAlignment --REFERENCE_SEQUENCE "${ref}" --UNMAPPED_BAM "${normalbam}" --ALIGNED_BAM "${process_outd}"/aln.sam.gz --OUTPUT "${outfile}" --SORT_ORDER coordinate --TMP_DIR "${tmpdir}" --MAX_RECORDS_IN_RAM ${max_records} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.align_norm_ubam_conda_envs()
{
    define_conda_env bwa bwa.yml
    define_conda_env gatk4 gatk4.yml
}

########
genodb.align_tum_ubam_explain_opts()
{
    # -r option
    description="Reference genome file"
    explain_opt "-r" "<file>" "$description"

    # -tumorbam option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_opt "-tumorbam" "<file>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_opt "-mrec" "<int>" "$description"

    # -out-processdir option
    local description="output directory"
    explain_opt "-out-processdir" "<file>" "$description"

    # -outfile option
    local description="output file"
    explain_opt "-outfile" "<file>" "$description"

    # -cpus option
    local description="number of cpus"
    explain_opt "-cpus" "<int>" "$description"
}

########
genodb.align_tum_ubam_identify_cmdline_opts()
{
    opt_is_cmdline "-r"
    opt_is_cmdline "-tumorbam"
    opt_is_non_mandatory_cmdline "-mrec"
}

########
genodb.align_tum_ubam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -r option
    define_cmdline_infile_opt "$cmdline" "-r" optlist || return 1

    # Get data directory
    local abs_datadir=$(get_absolute_shdirname "data")

    # -tumorbam option
    local tumorbam="${abs_datadir}"/tumor_unmapped.bam
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -mrec option
    define_cmdline_opt_if_given "$cmdline" "-mrec" optlist || return 1

    # -cpus option
    define_procspec_opt "${process_spec}" "-cpus" "cpus" optlist || return 1

    # -outfile option
    local outfile="${abs_datadir}"/tumor_aligned.bam
    define_opt "-outfile" "$outfile" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
genodb.align_tum_ubam()
{
    # Initialize variables
    local process_outd=$(read_opt_value_from_func_args "-out-processdir" "$@")
    local tumorbam=$(read_opt_value_from_func_args "-tumorbam" "$@")
    local ref=$(read_opt_value_from_func_args "-r" "$@")
    local max_records=$(read_opt_value_from_func_args "-mrec" "$@")
    if [ "$max_records" = "${DEBASHER_OPT_NOT_FOUND}" ]; then
        max_records=${DEFAULT_MAX_RECORDS_IN_RAM_GATK}
    fi
    local cpus=$(read_opt_value_from_func_args "-cpus" "$@")
    local outfile=$(read_opt_value_from_func_args "-outfile" "$@")

    # Create tmpdir for gatk
    tmpdir="${process_outd}"/tmp
    mkdir "${tmpdir}" || return 1

    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda_activate gatk4 2>&1 || return 1

    # Execute gatk SamToFastq
    logmsg "* Executing gatk SamToFastq..."
    gatk --java-options "-Xmx4G" SamToFastq --INPUT "${tumorbam}" --FASTQ "${process_outd}"/reads_r1.fastq.gz --SECOND_END_FASTQ "${process_outd}"/reads_r2.fastq.gz --TMP_DIR "${tmpdir}" --MAX_RECORDS_IN_RAM ${max_records} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bwa)..."
    conda_activate bwa 2>&1 || return 1

    # Execute bwa
    logmsg "* Executing bwa mem..."
    bwa mem -t ${cpus} "${ref}" <("${GZIP}" -d -c "${process_outd}"/reads_r1.fastq.gz) <("${GZIP}" -d -c "${process_outd}"/reads_r2.fastq.gz) | "${GZIP}" > "${process_outd}"/aln.sam.gz ; pipe_fail || return 1

    # Remove fastq files
    "${RM}" "${process_outd}"/reads_r1.fastq.gz "${process_outd}"/reads_r2.fastq.gz || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda_activate gatk4 2>&1 || return 1

    # Create dictionary
    if ! gatk_dict_exists "${ref}"; then
        logmsg "* Executing gatk CreateSequenceDictionary..."
        gatk --java-options "-Xmx4G" CreateSequenceDictionary -I "${ref}"
    fi

    # Execute gatk
    logmsg "* Executing gatk MergeBamAlignment..."
    gatk --java-options "-Xmx4G" MergeBamAlignment --REFERENCE_SEQUENCE "${ref}" --UNMAPPED_BAM "${tumorbam}" --ALIGNED_BAM "${process_outd}"/aln.sam.gz --OUTPUT "${outfile}" --SORT_ORDER coordinate --TMP_DIR "${tmpdir}" --MAX_RECORDS_IN_RAM ${max_records} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
genodb.align_tum_ubam_conda_envs()
{
    define_conda_env bwa bwa.yml
    define_conda_env gatk4 gatk4.yml
}
