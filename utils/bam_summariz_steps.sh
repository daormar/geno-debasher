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

#############
# CONSTANTS #
#############

###########################
# BAM SUMMARIZATION STEPS #
###########################

########
generate_vcf_list()
{
    # Initialize variables
    local summarydir=$1

    # Generate vcf list
    for file in ${summarydir}/*.${SUMMARY_FILE_EXT}; do
        cat ${file}
    done
}

########
merge_germline_snvs_document()
{
    step_description "Merge generated germline vcfs."
}

########
merge_germline_snvs_explain_cmdline_opts()
{
    :
}

########
merge_germline_snvs_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # Get germline snvs summary directory
    local abs_sumdir=`get_absolute_shdirname ${GERM_SNVS_SUM_DIR_BASENAME}`

    # -summarydir option
    define_opt "-summarydir" ${abs_sumdir} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # -mem option
    local mem
    mem=`extract_mem_from_stepspec "$stepspec"` || exit 1
    mem=`slurm_to_java_mem_spec ${mem}` || exit 1
    define_opt "-mem" $mem optlist

    # Save option list
    save_opt_list optlist    
}

########
merge_germline_snvs()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local summarydir=`read_opt_value_from_line "$*" "-summarydir"`
    local mem=`read_opt_value_from_line "$*" "-mem"`

    # Generate list file
    generate_vcf_list ${summarydir} > ${summarydir}/variant_files.list || exit 1
    
    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gatk4 2>&1 || exit 1

    logmsg "* Executing gatk MergeVcfs..."
    gatk --java-options "-Xmx${mem}" MergeVcfs -I ${summarydir}/variant_files.list -O ${summarydir}/summarized_variants.vcf.gz || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
merge_germline_snvs_conda_envs()
{
    define_conda_env gatk4 gatk4.yml
}
