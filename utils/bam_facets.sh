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
. "${biopanpipe_bindir}"/bam_common_lib || exit 1

#################
# CFG FUNCTIONS #
#################

########
bam_facets_shared_dirs()
{
    define_shared_dir "${DATADIR_BASENAME}"
}

########
bam_facets_fifos()
{
    :
}

#################################
# SNP-PILEUP + FACETS PROCESSES #
#################################

########
snp_pileup_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -sv option
    description="SNP vcf file"
    explain_cmdline_req_opt "-sv" "<string>" "$description"
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
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" "$normalbam" optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" "$tumorbam" optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
snp_pileup()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_line "$*" "-process-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local snpvcf=`read_opt_value_from_line "$*" "-sv"`

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
facets_explain_cmdline_opts()
{
    # -pileup-counts option
    description="SNP pileup file (required if pileup process has not been performed)"
    explain_cmdline_opt "-sp" "<string>" "$description"
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
    local process_outd=`read_opt_value_from_line "$*" "-process-outd"`
    local pileup_counts=`read_opt_value_from_line "$*" "-pileup-counts"`

    # Define --snpPileupCounts option if output from snp-pileup is available
    # snp_pìleup_opt=`get_snp_pileup_opt "${pileup_outd}"`

    # Activate conda environment if needed
    logmsg "* Activating conda environment..."
    conda activate facets 2>&1 || exit 1

    # Execute facets
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing facets..."
    Rscript "${biopanpipe_bindir}"/run_facets -c "${pileup_counts}" -o "${process_outd}" 2>&1 || exit 1

    # Deactivate conda environment if needed
    logmsg "* Dectivating conda environment..."
    conda deactivate 2>&1
}

########
facets_conda_envs()
{
    define_conda_env facets facets.yml
}
