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

#############
# CONSTANTS #
#############

DEFAULT_MAX_RECORDS_IN_RAM_GATK=1000000

##############################
# BAM MANIPULATION PROCESSES #
##############################

########
index_norm_bam_explain_cmdline_opts()
{
    :
}

########
index_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`
    local normalbam="${abs_datadir}"/normal.bam
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -out-nbidx option
    define_opt "-out-nbidx" "$normalbam".bai optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
index_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local idx_file=`read_opt_value_from_line "$*" "-out-nbidx"`

    # Remove previous index if one was created
    if [ -f "${idx_file}" ]; then
        rm "${idx_file}" || return 1
    fi

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || return 1

    # Execute samtools
    logmsg "* Executing samtools index..."
    samtools index "${normalbam}" "${idx_file}" 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
index_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
index_tum_bam_explain_cmdline_opts()
{
    :
}

########
index_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`
    local tumorbam="${abs_datadir}"/tumor.bam
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -out-tbidx option
    define_opt "-out-tbidx" "$tumorbam".bai optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
index_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local idx_file=`read_opt_value_from_line "$*" "-out-tbidx"`

    # Remove previous index if one was created
    if [ -f "${idx_file}" ]; then
        rm "${idx_file}" || return 1
    fi

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || return 1

    # Execute samtools
    logmsg "* Executing samtools index..."
    samtools index "${tumorbam}" "${idx_file}" 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
index_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
sort_norm_bam_explain_cmdline_opts()
{
    :
}

########
sort_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    local process_outd=`get_process_outdir_given_process_spec "$process_spec"`
    define_opt "-process-outd" "${process_outd}" optlist || return 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`
    local normalbam="${abs_datadir}"/normal.bam
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
sort_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local process_outd=`read_opt_value_from_line "$*" "-process-outd"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || return 1

    # Verify if bam file is already sorted
    local bam_is_sorted=`samtools view -H "${normalbam}" | $GREP SO:coordinate | wc -l` || return 1
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
        mv "${process_outd}"/sorted.bam "${normalbam}" 2>&1 || return 1
    fi

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
sort_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
sort_tum_bam_explain_cmdline_opts()
{
    :
}

########
sort_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    local process_outd=`get_process_outdir_given_process_spec "$process_spec"`
    define_opt "-process-outd" "${process_outd}" optlist || return 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`
    local tumorbam="${abs_datadir}"/tumor.bam
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
sort_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local process_outd=`read_opt_value_from_line "$*" "-process-outd"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || return 1

    # Verify if bam file is already sorted
    local bam_is_sorted=`samtools view -H "${tumorbam}" | $GREP SO:coordinate | wc -l` || return 1
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
        mv "${process_outd}"/sorted.bam "${tumorbam}" 2>&1 || return 1
    fi

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
sort_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
samtools_mpileup_norm_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"
}

########
samtools_mpileup_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    local process_outd=`get_process_outdir_given_process_spec "$process_spec"`
    define_opt "-process-outd" "${process_outd}" optlist || return 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || return 1
    define_opt "-r" "$genref" optlist || return 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # Save option list
    save_opt_list optlist
}

########
get_samtools_mpileup_l_opt()
{
    local mpbfile=$1

    if [ "${mpbfile}" = ${OPT_NOT_FOUND} ]; then
        echo ""
    else
        echo "-l ${mpbfile}"
    fi
}

########
samtools_mpileup_norm_bam()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-process-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || return 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=`get_samtools_mpileup_l_opt ${mpbfile}`

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
samtools_mpileup_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
samtools_mpileup_tum_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"
}

########
samtools_mpileup_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    local process_outd=`get_process_outdir_given_process_spec "$process_spec"`
    define_opt "-process-outd" "${process_outd}" optlist || return 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || return 1
    define_opt "-r" "$genref" optlist || return 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # Save option list
    save_opt_list optlist
}

########
samtools_mpileup_tum_bam()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-process-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || return 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=`get_samtools_mpileup_l_opt ${mpbfile}`

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
samtools_mpileup_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
parallel_samtools_mpileup_norm_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_samtools_mpileup_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # Obtain splitdir directory
    abs_splitdir=`get_absolute_shdirname "${SPLITDIR_BASENAME}"`

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || return 1
    define_opt "-r" "$genref" optlist || return 1

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; return 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || return 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        normalbam="${abs_splitdir}"/normal_${contig}.bam
        define_opt "-normalbam" "${normalbam}" specific_optlist || return 1
        define_opt "-contig" $contig specific_optlist || return 1
        define_opt "-outfile" "${process_outdir}"/normal_${contig}.pileup.gz specific_optlist || return 1
        save_opt_list specific_optlist
    done
}

########
parallel_samtools_mpileup_norm_bam()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-out-processdir"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local outfile=`read_opt_value_from_line "$*" "-outfile"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || return 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=`get_samtools_mpileup_l_opt ${mpbfile}`

    # Generate pileup file
    logmsg "* Generating pileup file (contig $contig)..."
    samtools mpileup ${smp_l_opt} -f "${ref}" -o "${process_outd}"/normal_${contig}.pileup "$normalbam" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    "${GZIP}" -c "${process_outd}"/normal_${contig}.pileup > "${outfile}" || return 1
}

########
parallel_samtools_mpileup_norm_bam_reset_outdir()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-out-processdir"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Remove files
    logmsg "* Resetting output directory..."
    rm -f "${process_outd}"/normal_${contig}.*
}

########
parallel_samtools_mpileup_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
parallel_samtools_mpileup_tum_bam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -mpb option
    description="BED file for mpileup"
    explain_cmdline_opt "-mpb" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_samtools_mpileup_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # Obtain splitdir directory
    abs_splitdir=`get_absolute_shdirname "${SPLITDIR_BASENAME}"`

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || return 1
    define_opt "-r" "$genref" optlist || return 1

    # -mpb option
    define_cmdline_opt_if_given "$cmdline" "-mpb" optlist

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; return 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || return 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        tumorbam="${abs_splitdir}"/tumor_${contig}.bam
        define_opt "-tumorbam" "${tumorbam}" specific_optlist || return 1
        define_opt "-contig" $contig specific_optlist || return 1
        define_opt "-outfile" "${process_outdir}"/tumor_${contig}.pileup.gz specific_optlist || return 1
        save_opt_list specific_optlist
    done
}

########
parallel_samtools_mpileup_tum_bam()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-out-processdir"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local mpbfile=`read_opt_value_from_line "$*" "-mpb"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local outfile=`read_opt_value_from_line "$*" "-outfile"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || return 1

    # Obtain samtools mpileup -L opt
    local smp_l_opt=`get_samtools_mpileup_l_opt ${mpbfile}`

    # Generate pileup file
    logmsg "* Generating pileup file (contig $contig)..."
    samtools mpileup ${smp_l_opt} -f "${ref}" -o "${process_outd}"/tumor_${contig}.pileup "${tumorbam}" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress pileup file
    logmsg "* Compressing pileup file..."
    "${GZIP}" -c "${process_outd}"/tumor_${contig}.pileup > "${outfile}" || return 1
}

########
parallel_samtools_mpileup_tum_bam_reset_outdir()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-out-processdir"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Remove files
    logmsg "* Resetting output directory..."
    rm -f "${process_outd}"/tumor_${contig}.*
}

########
parallel_samtools_mpileup_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
parallel_split_norm_bam_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_split_norm_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Obtain splitdir directory
    abs_splitdir=`get_absolute_shdirname "${SPLITDIR_BASENAME}"`

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; return 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || return 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        define_opt "-contig" $contig specific_optlist || return 1
        define_opt "-outfile" "${abs_splitdir}"/normal_${contig}.bam specific_optlist || return 1
        save_opt_list specific_optlist
    done
}

########
parallel_split_norm_bam()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local outfile=`read_opt_value_from_line "$*" "-outfile"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || return 1

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
parallel_split_norm_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
parallel_split_tum_bam_explain_cmdline_opts()
{
    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_split_tum_bam_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Obtain splitdir directory
    abs_splitdir=`get_absolute_shdirname "${SPLITDIR_BASENAME}"`

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; return 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || return 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        define_opt "-contig" $contig specific_optlist || return 1
        define_opt "-outfile" "${abs_splitdir}"/tumor_${contig}.bam specific_optlist || return 1
        save_opt_list specific_optlist
    done
}

########
parallel_split_tum_bam()
{
    # Initialize variables
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local outfile=`read_opt_value_from_line "$*" "-outfile"`

    # Activate conda environment
    logmsg "* Activating conda environment (samtools)..."
    conda activate samtools 2>&1 || return 1

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
parallel_split_tum_bam_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
bedtools_genomecov_norm_bam_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"
}

########
bedtools_genomecov_norm_bam_define_opts()
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
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
bedtools_genomecov_norm_bam()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-out-processdir"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate bedtools 2>&1 || return 1

    # Execute samtools
    logmsg "* Executing bedtools coverage..."
    bedtools genomecov -ibam "${normalbam}" > "${process_outd}"/genomecov.tsv || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
bedtools_genomecov_norm_bam_conda_envs()
{
    define_conda_env bedtools bedtools.yml
}

########
bedtools_genomecov_tum_bam_explain_cmdline_opts()
{
    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"
}

########
bedtools_genomecov_tum_bam_define_opts()
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
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
bedtools_genomecov_tum_bam()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-out-processdir"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate bedtools 2>&1 || return 1

    # Execute samtools
    logmsg "* Executing bedtools coverage..."
    bedtools genomecov -ibam "${tumorbam}" > "${process_outd}"/genomecov.tsv || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
bedtools_genomecov_tum_bam_conda_envs()
{
    define_conda_env bedtools bedtools.yml
}

########
norm_bam_to_ubam_explain_cmdline_opts()
{
    # -n option
    description="Normal unmapped bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_cmdline_opt "-mrec" "<int>" "$description"
}

########
norm_bam_to_ubam_define_opts()
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
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -normalbam option
    define_opt "-normalbam" "${abs_datadir}"/normal.bam optlist || return 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || return 1

    # -outfile option
    define_opt "-outfile" "${abs_datadir}"/normal_unmapped.bam optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
norm_bam_to_ubam()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-out-processdir"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local max_records=`read_opt_value_from_line "$*" "-mrec"`
    local outfile=`read_opt_value_from_line "$*" "-outfile"`

    # Create tmpdir for gatk
    tmpdir="${process_outd}"/tmp
    mkdir "${tmpdir}" || return 1

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gatk4 2>&1 || return 1

    # Execute gatk RevertSam
    logmsg "* Executing gatk RevertSam..."
    gatk --java-options "-Xmx4G" RevertSam --INPUT "${normalbam}" --OUTPUT "${outfile}" --SANITIZE true --SORT_ORDER queryname --TMP_DIR "${tmpdir}" --MAX_RECORDS_IN_RAM ${max_records} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
norm_bam_to_ubam_conda_envs()
{
    define_conda_env gatk4 gatk4.yml
}

########
tum_bam_to_ubam_explain_cmdline_opts()
{
    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_cmdline_opt "-mrec" "<int>" "$description"
}

########
tum_bam_to_ubam_define_opts()
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
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -tumorbam option
    define_opt "-tumorbam" "${abs_datadir}"/tumor.bam optlist || return 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || return 1

    # -outfile option
    define_opt "-outfile" "${abs_datadir}"/tumor_unmapped.bam optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
tum_bam_to_ubam()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-out-processdir"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local max_records=`read_opt_value_from_line "$*" "-mrec"`
    local outfile=`read_opt_value_from_line "$*" "-outfile"`

    # Create tmpdir for gatk
    tmpdir="${process_outd}"/tmp
    mkdir "${tmpdir}" || return 1

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gatk4 2>&1 || return 1

    # Execute gatk RevertSam
    logmsg "* Executing gatk RevertSam..."
    gatk --java-options "-Xmx4G" RevertSam --INPUT "${tumorbam}" --OUTPUT "${outfile}" --SANITIZE true --SORT_ORDER queryname --TMP_DIR "${tmpdir}" --MAX_RECORDS_IN_RAM ${max_records} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
tum_bam_to_ubam_conda_envs()
{
    define_conda_env gatk4 gatk4.yml
}

########
align_norm_ubam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal unmapped bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_cmdline_opt "-mrec" "<int>" "$description"
}

########
align_norm_ubam_define_opts()
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
    local genref
    genref=`get_ref_filename "$cmdline"` || return 1
    define_opt "-r" "$genref" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -normalbam option
    define_opt "-normalbam" "${abs_datadir}"/normal_unmapped.bam optlist || return 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # -outfile option
    define_opt "-outfile" "${abs_datadir}"/normal_aligned.bam optlist || return 1

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
align_norm_ubam()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-out-processdir"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local max_records=`read_opt_value_from_line "$*" "-mrec"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    local outfile=`read_opt_value_from_line "$*" "-outfile"`

    # Create tmpdir for gatk
    tmpdir="${process_outd}"/tmp
    mkdir "${tmpdir}" || return 1

    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda activate gatk4 2>&1 || return 1

    # Execute gatk SamToFastq
    logmsg "* Executing gatk SamToFastq..."
    gatk --java-options "-Xmx4G" SamToFastq --INPUT "${normalbam}" --FASTQ "${process_outd}"/reads_r1.fastq.gz --SECOND_END_FASTQ "${process_outd}"/reads_r2.fastq.gz --TMP_DIR "${tmpdir}" --MAX_RECORDS_IN_RAM ${max_records} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bwa)..."
    conda activate bwa 2>&1 || return 1

    # Execute bwa
    logmsg "* Executing bwa index..."
    bwa index -a bwtsw "${ref}"

    logmsg "* Executing bwa mem..."
    bwa mem -t ${cpus} "${ref}" <("${GZIP}" -d -c "${process_outd}"/reads_r1.fastq.gz) <("${GZIP}" -d -c "${process_outd}"/reads_r2.fastq.gz) | "${GZIP}" > "${process_outd}"/aln.sam.gz ; pipe_fail || return 1

    # Remove fastq files
    rm "${process_outd}"/reads_r1.fastq.gz "${process_outd}"/reads_r2.fastq.gz || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda activate gatk4 2>&1 || return 1

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
align_norm_ubam_conda_envs()
{
    define_conda_env bwa bwa.yml
    define_conda_env gatk4 gatk4.yml
}

########
align_tum_ubam_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -mrec option
    description="Maximum number of records stored in RAM required by GATK. The higher the number, the more RAM is required but the lower the number of files created for external sorting (${DEFAULT_MAX_RECORDS_IN_RAM_GATK} by default)"
    explain_cmdline_opt "-mrec" "<int>" "$description"
}

########
align_tum_ubam_define_opts()
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
    local genref
    genref=`get_ref_filename "$cmdline"` || return 1
    define_opt "-r" "$genref" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -tumorbam option
    define_opt "-tumorbam" "${abs_datadir}"/tumor_unmapped.bam optlist || return 1

    # -mrec option
    define_cmdline_nonmandatory_opt "$cmdline" "-mrec" ${DEFAULT_MAX_RECORDS_IN_RAM_GATK} optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # -outfile option
    define_opt "-outfile" "${abs_datadir}"/tumor_aligned.bam optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
align_tum_ubam()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-out-processdir"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local max_records=`read_opt_value_from_line "$*" "-mrec"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    local outfile=`read_opt_value_from_line "$*" "-outfile"`

    # Create tmpdir for gatk
    tmpdir="${process_outd}"/tmp
    mkdir "${tmpdir}" || return 1

    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda activate gatk4 2>&1 || return 1

    # Execute gatk SamToFastq
    logmsg "* Executing gatk SamToFastq..."
    gatk --java-options "-Xmx4G" SamToFastq --INPUT "${tumorbam}" --FASTQ "${process_outd}"/reads_r1.fastq.gz --SECOND_END_FASTQ "${process_outd}"/reads_r2.fastq.gz --TMP_DIR "${tmpdir}" --MAX_RECORDS_IN_RAM ${max_records} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bwa)..."
    conda activate bwa 2>&1 || return 1

    # Execute bwa
    logmsg "* Executing bwa mem..."
    bwa mem -t ${cpus} "${ref}" <("${GZIP}" -d -c "${process_outd}"/reads_r1.fastq.gz) <("${GZIP}" -d -c "${process_outd}"/reads_r2.fastq.gz) | "${GZIP}" > "${process_outd}"/aln.sam.gz ; pipe_fail || return 1

    # Remove fastq files
    rm "${process_outd}"/reads_r1.fastq.gz "${process_outd}"/reads_r2.fastq.gz || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (gatk)..."
    conda activate gatk4 2>&1 || return 1

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
align_tum_ubam_conda_envs()
{
    define_conda_env bwa bwa.yml
    define_conda_env gatk4 gatk4.yml
}
