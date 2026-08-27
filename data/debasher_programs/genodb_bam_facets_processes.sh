1# Geno-DeBasher package
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

#################################
# SNP-PILEUP + FACETS PROCESSES #
#################################

########
snp_pileup_explain_opts()
{
    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_opt "-t" "<string>" "$description"

    # -sv option
    description="SNP vcf file"
    explain_opt "-sv" "<string>" "$description"

    # -normalbam option
    local description="normal bam file"
    explain_opt "-normalbam" "<string>" "$description"

    # -tumorbam option
    local description="tumor bam file"
    explain_opt "-tumorbam" "<string>" "$description"

    # -process-outd option
    local description="output directory"
    explain_opt "-process-outd" "<string>" "$description"
}

########
snp_pileup_identify_cmdline_opts()
{
    opt_is_cmdline "-n"
    opt_is_cmdline "-t"
    opt_is_cmdline "-sv"
}

########
snp_pileup_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    local process_outd=`get_process_outdir_given_process_spec "$process_spec"`
    define_opt "-process-outd" ${process_outd} optlist || exit 1

    # -sv option
    define_cmdline_infile_opt "$cmdline" "-sv" optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`genodb_bam_common::get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" "$normalbam" optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`genodb_bam_common::get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" "$tumorbam" optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
snp_pileup()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-process-outd" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local snpvcf=`read_opt_value_from_func_args "-sv" "$@"`

    # Activate conda environment if needed
    logmsg "* Activating conda environment..."
    conda activate snp-pileup 2>&1 || exit 1

    # Execute snp-pileup
    logmsg "* Executing snp-pileup..."
    snp-pileup ${snpvcf} "${process_outd}"/snp-pileup-counts.csv "${normalbam}" "${tumorbam}" 2>&1 || exit 1

    # Deactivate conda environment if needed
    logmsg "* Dectivating conda environment..."
    conda deactivate 2>&1
}

########
snp_pileup_conda_envs()
{
    define_conda_env snp-pileup snp-pileup.yml
}

#######
facets_explain_opts()
{
    # -pileup-counts option
    description="SNP pileup file (required if pileup process has not been performed)"
    explain_opt "-sp" "<string>" "$description"

    # -pileup-counts option
    local description="pile up counts file"
    explain_opt "-pileup-counts" "<string>" "$description"
}

########
facets_identify_cmdline_opts()
{
    opt_is_cmdline "-sp"
}

########
facets_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local jobspec=$2
    local optlist=""

    # Define the -process-outd option, the output directory for the process
    local process_outd=`get_process_outdir_given_process_spec "$process_spec"`
    define_opt "-process-outd" "${process_outd}" optlist || exit 1

    local pileup_dep=`find_dependency_for_process "${jobspec}" snp_pileup`
    if [ ${pileup_dep} != ${DEP_NOT_FOUND} ]; then
        local pileup_outd=`get_default_outd_for_dep ${outd} "${pileup_dep}"`
        local pileup_counts_file="${pileup_outd}"/snp-pileup-counts.csv
        define_opt "-pileup-counts" "${pileup_counts_file}" optlist || exit 1
    else
        define_cmdline_infile_opt "${cmdline}" "-sp" optlist || exit 1
    fi

    # Save option list
    save_opt_list optlist
}

########
facets()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-process-outd" "$@"`
    local pileup_counts=`read_opt_value_from_func_args "-pileup-counts" "$@"`

    # Activate conda environment if needed
    logmsg "* Activating conda environment..."
    conda activate facets 2>&1 || exit 1

    # Execute facets
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing facets..."
    Rscript "${genodebasher_libexecdir}"/genodb_run_facets -c "${pileup_counts}" -o "${process_outd}" 2>&1 || exit 1

    # Deactivate conda environment if needed
    logmsg "* Dectivating conda environment..."
    conda deactivate 2>&1
}

########
facets_conda_envs()
{
    define_conda_env facets facets.yml
}
