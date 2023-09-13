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

# INCLUDE BASH LIBRARY
. ${biopanpipe_bindir}/bam_common_lib || exit 1

#################
# CFG FUNCTIONS #
#################

########
bam_filter_shared_dirs()
{
    define_shared_dir "${DATADIR_BASENAME}"
}

########
bam_filter_fifos()
{
    :
}

########################
# BAM FILTER PROCESSES #
########################

########
filter_norm_bam_contigs_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"
}

########
filter_norm_bam_contigs_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    local process_outd=`get_process_outdir_given_process_spec "$process_spec"`
    define_opt "-process-outd" "${process_outd}" optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" "$genref" optlist || exit 1

    # -normalbam option
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
    local normalbam="${abs_datadir}"/normal.bam
    define_opt "-normalbam" "$normalbam" optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
get_ref_contigs()
{
    local faifile=$1

    "${AWK}" '{print $1}' ${faifile}
}

########
remove_line_breaks_from_file()
{
    local file=$1
    echo `cat "${file}"`
}

########
get_contigs_from_header()
{
    local header=$1

    while read line; do
        local fields=( $line )
        local num_fields=${#fields[@]}
        if [ ${num_fields} -ge 3 ]; then
            if [ ${fields[0]} = "@SQ" ]; then
                contigfield=${fields[1]}
                contig=${contigfield:3}
                echo ${contig}
            fi
        fi
    done < ${header}
}

########
filter_norm_bam_contigs()
{
    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local process_outd=`read_opt_value_from_line "$*" "-process-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || exit 1

    # Obtain contigs given in reference
    get_ref_contigs "${ref}".fai > "${process_outd}"/refcontigs

    # Obtain new header

    ## Extract sam header
    samtools view -H "${normalbam}" > "${process_outd}"/original_header || exit 1

    ## Generate new sam header
    "${biopanpipe_bindir}"/get_filtered_sam_header -h "${process_outd}"/original_header -l "${process_outd}"/refcontigs > "${process_outd}"/new_header || exit 1

    # Generate filtered bam
    {
        # Print header
        cat "${process_outd}"/new_header

        # Print contig information
        contigs=`get_contigs_from_header "${process_outd}"/new_header`
        samtools view "${normalbam}" "${contigs}" | "${biopanpipe_bindir}"/get_filtered_sam_align -l "${process_outd}"/refcontigs
    } | samtools view -bo "${process_outd}"/filtered.bam -

    # Move bam file
    mv "${process_outd}"/filtered.bam "${normalbam}"

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
filter_norm_bam_contigs_conda_envs()
{
    define_conda_env samtools samtools.yml
}

########
filter_tum_bam_contigs_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"
}

########
filter_tum_bam_contigs_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    local process_outd=`get_process_outdir_given_process_spec "$process_spec"`
    define_opt "-process-outd" "${process_outd}" optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" "$genref" optlist || exit 1

    # -tumorbam option
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`
    local tumorbam="${abs_datadir}"/tumor.bam
    define_opt "-tumorbam" "$tumorbam" optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
filter_tum_bam_contigs()
{
    # Initialize variables
    local ref=`read_opt_value_from_line "$*" "-r"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local process_outd=`read_opt_value_from_line "$*" "-process-outd"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate samtools 2>&1 || exit 1

    # Obtain contigs given in reference
    get_ref_contigs "${ref}".fai > "${process_outd}"/refcontigs

    # Obtain new header

    ## Extract sam header
    samtools view -H "${tumorbam}" > "${process_outd}"/original_header || exit 1

    ## Generate new sam header
    "${biopanpipe_bindir}"/get_filtered_sam_header -h "${process_outd}"/original_header -l "${process_outd}"/refcontigs > "${process_outd}"/new_header || exit 1

    # Generate filtered bam
    {
        # Print header
        cat "${process_outd}"/new_header

        # Print contig information
        contigs=`get_contigs_from_header "${process_outd}"/new_header`
        samtools view "${tumorbam}" ${contigs} | "${biopanpipe_bindir}"/get_filtered_sam_align -l "${process_outd}"/refcontigs
    } | samtools view -bo "${process_outd}"/filtered.bam -

    # Move bam file
    mv "${process_outd}"/filtered.bam "${tumorbam}"

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
filter_tum_bam_contigs_conda_envs()
{
    define_conda_env samtools samtools.yml
}
