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

DEFAULT_MIN_SEQ_DEPTH_FACETS_PREPROC=35

##########################
# BAM ANALYSIS PROCESSES #
##########################

########
manta_germline_document()
{
    process_description "Analyzes a normal \`bam\` file using Manta."
}

########
manta_germline_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call"
    explain_cmdline_opt "-cr" "<string>" "$description"
}

########
manta_germline_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -cr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
get_callreg_opt()
{
    local callregf=$1
    local basecallregf=`$BASENAME ${callregf}`

    if [ "${basecallregf}" = ${NOFILE} -o "${callregf}" = "" ]; then
        echo ""
    else
        echo "--callRegions ${callregf}"
    fi
}

########
manta_germline()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local callregf=`read_opt_value_from_func_args "-cr" "$@"`
    local cpus=`read_opt_value_from_func_args "-cpus" "$@"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate manta 2>&1 || return 1

    # Configure Manta
    logmsg "* Executing configManta.py..."
    configManta.py --bam "${normalbam}" --referenceFasta "${ref}" ${call_reg_opt} --runDir "${process_outd}" 2>&1 || return 1

    # Execute Manta
    logmsg "* Executing runWorkflow.py..."
    "${process_outd}"/runWorkflow.py -m local -j ${cpus} 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
manta_germline_conda_envs()
{
    define_conda_env manta manta.yml
}

########
manta_somatic_document()
{
    process_description "Analyzes a pair of normal and tumor \`bam\` files using Manta."
}

########
manta_somatic_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call"
    explain_cmdline_opt "-cr" "<string>" "$description"
}

########
manta_somatic_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -cr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
manta_somatic()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local callregf=`read_opt_value_from_func_args "-cr" "$@"`
    local cpus=`read_opt_value_from_func_args "-cpus" "$@"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate manta 2>&1 || return 1

    # Configure Manta
    logmsg "* Executing configManta.py..."
    configManta.py --normalBam "${normalbam}" --tumorBam "${tumorbam}" --referenceFasta "${ref}" ${call_reg_opt} --runDir "${process_outd}" 2>&1 || return 1

    # Execute Manta
    logmsg "* Executing runWorkflow.py..."
    "${process_outd}"/runWorkflow.py -m local -j ${cpus} 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
manta_somatic_conda_envs()
{
    define_conda_env manta manta.yml
}

########
strelka_germline_document()
{
    process_description "Analyzes a normal \`bam\` files using Strelka."
}

########
strelka_germline_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call"
    explain_cmdline_opt "-cr" "<string>" "$description"
}

########
strelka_germline_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -cr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || return 1

    # Get germline snvs summary directory
    local abs_sumdir=`get_absolute_shdirname ${GERM_SNVS_SUM_DIR_BASENAME}`

    # -summarydir option
    define_opt "-summarydir" "${abs_sumdir}" optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
strelka_germline()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local callregf=`read_opt_value_from_func_args "-cr" "$@"`
    local summarydir=`read_opt_value_from_func_args "-summarydir" "$@"`
    local cpus=`read_opt_value_from_func_args "-cpus" "$@"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate strelka 2>&1 || return 1

    # Configure Strelka
    logmsg "* Executing configureStrelkaGermlineWorkflow.py..."
    configureStrelkaGermlineWorkflow.py --bam "${normalbam}" --referenceFasta "${ref}" ${call_reg_opt} --runDir "${process_outd}" 2>&1 || return 1

    # Execute Strelka
    logmsg "* Executing runWorkflow.py..."
    "${process_outd}"/runWorkflow.py -m local -j ${cpus} 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Create file in summary directory
    local label=strelka_germline
    vcf="${process_outd}"/results/variants/variants.vcf.gz
    create_summary_file "${summarydir}" ${label} "${vcf}"
}

########
strelka_germline_conda_envs()
{
    define_conda_env strelka strelka.yml
}

########
platypus_germline_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"
}

########
platypus_germline_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # Get germline snvs summary directory
    local abs_sumdir=`get_absolute_shdirname "${GERM_SNVS_SUM_DIR_BASENAME}"`

    # -summarydir option
    define_opt "-summarydir" "${abs_sumdir}" optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
platypus_germline_conda()
{
    # Initialize variables
    local ref=$1
    local normalbam=$2
    local process_outd=$3
    local summarydir=$4
    local cpus=$5

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate platypus 2>&1 || return 1

    # Run Platypus
    logmsg "* Executing Platypus.py..."
    platypus callVariants --bamFiles="${normalbam}" --refFile="${ref}" --output="${process_outd}"/output.vcf --nCPU=${cpus} --logFileName="${process_outd}"/platypus.log --verbosity=1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Create file in summary directory
    local label=platypus_germline
    vcf="${process_outd}"/output.vcf
    create_summary_file "${summarydir}" ${label} "${vcf}"
}

########
platypus_germline_local()
{
    # Initialize variables
    local ref=$1
    local normalbam=$2
    local process_outd=$3
    local summarydir=$4
    local cpus=$5

    # Run Platypus
    logmsg "* Executing Platypus.py..."
    python ${PLATYPUS_HOME_DIR}/bin/Platypus.py callVariants --bamFiles="${normalbam}" --refFile="${ref}" --nCPU=${cpus} --output="${process_outd}"/output.vcf --verbosity=1 2>&1 || return 1

    # Create file in summary directory
    local label=platypus_germline
    vcf="${process_outd}"/output.vcf
    create_summary_file "${summarydir}" ${label} "${vcf}"
}

########
platypus_germline()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local summarydir=`read_opt_value_from_func_args "-summarydir" "$@"`
    local cpus=`read_opt_value_from_func_args "-cpus" "$@"`

    if [ -z "${PLATYPUS_HOME_DIR}" ]; then
        platypus_germline_conda "${ref}" "${normalbam}" "${process_outd}" "${summarydir}" ${cpus}
    else
        platypus_germline_local "${ref}" "${normalbam}" "${process_outd}" "${summarydir}" ${cpus}
    fi
}

########
platypus_germline_conda_envs()
{
    define_conda_env platypus platypus.yml
}

########
gatk_haplotypecaller_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -sample-name option
    description="Sample name"
    explain_cmdline_req_opt "-sample-name" "<string>" "$description"
}

########
gatk_haplotypecaller_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -sample-name option
    define_cmdline_opt "$cmdline" "-sample-name" optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # -mem option
    local mem
    mem=`extract_mem_from_process_spec "$process_spec"` || return 1
    mem=`slurm_to_java_mem_spec ${mem}` || return 1
    define_opt "-mem" $mem optlist

    # Save option list
    save_opt_list optlist
}

########
gatk_haplotypecaller()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local sample_name=`read_opt_value_from_func_args "-sample-name" "$@"`
    local cpus=`read_opt_value_from_func_args "-cpus" "$@"`
    local mem=`read_opt_value_from_func_args "-mem" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gatk4 2>&1 || return 1

    # Run gatk HaplotypeCaller
    logmsg "* Executing gatk HaplotypeCaller..."
    gatk --java-options "-Xmx${mem}" HaplotypeCaller -R "${ref}" -I "${normalbam}" -O "${process_outd}"/output.g.vcf.gz -ERC GVCF --sample-name "${sample_name}" --native-pair-hmm-threads ${cpus} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
gatk_haplotypecaller_conda_envs()
{
    define_conda_env gatk4 gatk4.yml
}

########
strelka_somatic_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -cr option
    description="bgzipped and tabixed bed file to specify regions to call"
    explain_cmdline_opt "-cr" "<string>" "$description"
}

########
strelka_somatic_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -manta-outd option
    local manta_outd
    manta_outd=`get_process_outdir_adaptive manta_somatic`
    define_opt "-manta-outd" "${manta_outd}" optlist || return 1

    # -cr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
get_indel_cand_opt()
{
    local manta_outd=$1

    if [ -z "${manta_outd}" ]; then
        echo ""
    else
        manta_indel_file="${manta_outd}/results/variants/candidateSmallIndels.vcf.gz"
        if [ -f "${manta_indel_file}" ]; then
            echo "--indelCandidates ${manta_indel_file}"
        else
            echo "WARNING: Manta indel file for Strelka not found! (${manta_indel_file})" >&2
            echo ""
        fi
    fi
}

########
strelka_somatic()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local manta_outd=`read_opt_value_from_func_args "-manta-outd" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local callregf=`read_opt_value_from_func_args "-cr" "$@"`
    local cpus=`read_opt_value_from_func_args "-cpus" "$@"`

    # Define --indelCandidates option if output from Manta is available
    indel_cand_opt=`get_indel_cand_opt "${manta_outd}"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate strelka 2>&1 || return 1

    # Configure Strelka
    logmsg "* Executing configureStrelkaSomaticWorkflow.py..."
    configureStrelkaSomaticWorkflow.py --normalBam "${normalbam}" --tumorBam "${tumorbam}" --referenceFasta "${ref}" ${indel_cand_opt} ${call_reg_opt} --runDir "${process_outd}" 2>&1 || return 1

    # Execute Strelka
    logmsg "* Executing runWorkflow.py..."
    "${process_outd}"/runWorkflow.py -m local -j ${cpus} 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
strelka_somatic_conda_envs()
{
    define_conda_env strelka strelka.yml
}

########
mutect2_somatic_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -norm-sample-name option
    description="Normal sample name"
    explain_cmdline_req_opt "-norm-sample-name" "<string>" "$description"

    # -panel-of-normals
    description="File name with panel of normals"
    explain_cmdline_req_opt "-panel-of-normals" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"
}

########
mutect2_somatic_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -norm-sample-name option
    define_cmdline_opt "$cmdline" "-norm-sample-name" optlist || return 1

    # -panel-of-normals option
    define_cmdline_opt "$cmdline" "-panel-of-normals" optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # -mem option
    local mem
    mem=`extract_mem_from_process_spec "$process_spec"` || return 1
    mem=`slurm_to_java_mem_spec ${mem}` || return 1
    define_opt "-mem" $mem optlist

    # Save option list
    save_opt_list optlist
}

########
mutect2_somatic()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local norm_sample_name=`read_opt_value_from_func_args "-norm-sample-name" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local panel_of_normals=`read_opt_value_from_func_args "-panel-of-normals" "$@"`
    local cpus=`read_opt_value_from_func_args "-cpus" "$@"`
    local mem=`read_opt_value_from_func_args "-mem" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gatk4 2>&1 || return 1

    # Run Mutect2
    logmsg "* Executing gatk Mutect2..."
    gatk --java-options "-Xmx${mem}" Mutect2 -R "${ref}" -I "${normalbam}" -I "${tumorbam}" -O "${process_outd}"/somatic.vcf.gz -normal "${norm_sample_name}" --panel-of-normals "${panel_of_normals}" --native-pair-hmm-threads ${cpus} || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
mutect2_somatic_conda_envs()
{
    define_conda_env gatk4 gatk4.yml
}

########
lofreq_somatic_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"
}

########
lofreq_somatic_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
lofreq_somatic()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local cpus=`read_opt_value_from_func_args "-cpus" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate lofreq 2>&1 || return 1

    # Run lofreq somatic
    logmsg "* Executing lofeq somatic..."
    lofreq somatic -f "${ref}" -n "${normalbam}" -t "${tumorbam}" -o "${process_outd}"/somatic.vcf.gz --threads ${cpus} -o "${process_outd}"/out_ || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
lofreq_somatic_conda_envs()
{
    define_conda_env lofreq lofreq.yml
}

########
cnvkit_document()
{
    process_description "Analyzes a pair of normal and tumor \`bam\` files using CNVkit."
}

########
cnvkit_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"
}

########
cnvkit_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
cnvkit()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local cpus=`read_opt_value_from_func_args "-cpus" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate cnvkit 2>&1 || return 1

    # Run cnvkit
    logmsg "* Executing cnvkit.py..."
    cd "${process_outd}" # This is done since current implementation of
                    # cnvkit generates a file in current directory
                    # (genref.bed). Changing directory avoids possible
                    # racing conditions
    cnvkit.py batch "${tumorbam}" -n "${normalbam}" -m wgs -f "${ref}"  -d "${process_outd}" -p ${cpus} 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
cnvkit_conda_envs()
{
    define_conda_env cnvkit cnvkit.yml
}

########
snp_pileup_plus_facets_explain_cmdline_opts()
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

    # -md option
    description="Minimum sequencing depth to keep when preprocessing sample (${DEFAULT_MIN_SEQ_DEPTH_FACETS_PREPROC} by default)"
    explain_cmdline_opt "-md" "<int>" "$description"
}

########
snp_pileup_plus_facets_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # -sv option
    define_cmdline_infile_opt "$cmdline" "-sv" optlist || return 1

    # -md option
    define_cmdline_nonmandatory_opt "$cmdline" "-md" ${DEFAULT_MIN_SEQ_DEPTH_FACETS_PREPROC} optlist || return 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
snp_pileup_plus_facets()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local snpvcf=`read_opt_value_from_func_args "-sv" "$@"`
    local mindepth=`read_opt_value_from_func_args "-md" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment (snp-pileup)..."
    conda activate snp-pileup 2>&1 || return 1

    # Execute snp-pileup
    logmsg "* Executing snp-pileup..."
    snp-pileup "${snpvcf}" "${process_outd}"/snp-pileup-counts.csv "${normalbam}" "${tumorbam}" 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment if needed
    logmsg "* Activating conda environment (facets)..."
    conda activate facets 2>&1 || return 1

    # Execute facets
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing facets..."
    Rscript "${genodebasher_libexecdir}"/genodb_run_facets -c "${process_outd}"/snp-pileup-counts.csv -d ${mindepth} -o "${process_outd}" 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress snp-pileup file
    logmsg "* Compressing snp-pileup file..."
    "${GZIP}" "${process_outd}"/snp-pileup-counts.csv || return 1
}

########
snp_pileup_plus_facets_conda_envs()
{
    define_conda_env snp-pileup snp-pileup.yml
    define_conda_env facets facets.yml
}

########
gen_sequenza_gcc_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"
}

########
gen_sequenza_gcc_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || return 1
    define_opt "-r" "$genref" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -outfile option
    define_opt "-outfile" "${abs_datadir}"/sequenza_gccfile.txt.gz optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
gen_sequenza_gcc()
{
    # Initialize variables
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local outfile=`read_opt_value_from_func_args "-outfile" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate sequenza 2>&1 || return 1

    # Generate GC content file
    logmsg "* Generating GC content file..."
    sequenza-utils gc_wiggle -w 50 -f "$ref" -o "$outfile" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
gen_sequenza_gcc_conda_envs()
{
    define_conda_env sequenza sequenza.yml
}

########
sequenza_explain_cmdline_opts()
{
    # -gcc option
    description="GC content wiggle file for sequenza (required if no gen_sequenza_gcc process is defined)"
    explain_cmdline_opt "-gcc" "<string>" "$description"
}

########
sequenza_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_spec=$2
    local process_name=$3
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname "${DATADIR_BASENAME}"`

    # -gcc option
    local gccfile
    gccfile="${abs_datadir}"/sequenza_gccfile.txt.gz
    define_opt "-gcc" "$gccfile" optlist || return 1

    # Get normal pileup file
    local npileupdir
    npileupdir=`get_outd_for_dep_given_process_spec "${process_spec}" samtools_mpileup_norm_bam` || { errmsg "Error: dependency samtools_mpileup_norm_bam not defined for sequenza"; return 1; }
    local npileup="${npileupdir}"/normal.pileup.gz
    define_opt "-npileup" "${npileup}" optlist || return 1

    # Get tumor pileup file
    tpileupdir=`get_outd_for_dep_given_process_spec "${process_spec}" samtools_mpileup_tum_bam` || { errmsg "Error: dependency samtools_mpileup_tum_bam not defined for sequenza"; return 1; }
    tpileup="${tpileupdir}"/tumor.pileup.gz
    define_opt "-tpileup" "${tpileup}" optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
sequenza()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local gccont=`read_opt_value_from_func_args "-gcc" "$@"`
    local npileup=`read_opt_value_from_func_args "-npileup" "$@"`
    local tpileup=`read_opt_value_from_func_args "-tpileup" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || return 1

    # Generate seqz file
    logmsg "* Generating seqz file..."
    sequenza-utils bam2seqz --pileup -gc "${gccont}" -n "${npileup}" -t "${tpileup}" | "${GZIP}" > "${process_outd}"/seqz.gz ; pipe_fail || return 1

    # Execute sequenza
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing sequenza..."
    Rscript "${genodebasher_bindir}"/run_sequenza -s "${process_outd}"/seqz.gz -o "${process_outd}" 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
sequenza_conda_envs()
{
    define_conda_env sequenza sequenza.yml
}

########
parallel_bam2seqz_explain_cmdline_opts()
{
    # -gcc option
    description="GC content wiggle file for bam2seqz"
    explain_cmdline_opt "-gcc" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_bam2seqz_define_opts()
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

    # -gcc option
    local gccfile
    gccfile="${abs_datadir}"/sequenza_gccfile.txt.gz
    define_opt "-gcc" "$gccfile" optlist || return 1

    # Get normal pileup directory
    local npileupdir
    npileupdir=`get_process_outdir_adaptive parallel_samtools_mpileup_norm_bam`

    # Get tumor pileup directory
    local tpileupdir
    tpileupdir=`get_process_outdir_adaptive parallel_samtools_mpileup_tum_bam`

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; return 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || return 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        npileup="${npileupdir}"/normal_${contig}.pileup.gz
        define_opt "-npileup" "${npileup}" specific_optlist || return 1
        tpileup="${tpileupdir}"/tumor_${contig}.pileup.gz
        define_opt "-tpileup" "${tpileup}" specific_optlist || return 1
        define_opt "-contig" $contig specific_optlist || return 1
        define_opt "-outfile" "${process_outdir}"/${contig}_seqz.gz specific_optlist || return 1

        save_opt_list specific_optlist
    done
}

########
parallel_bam2seqz()
{
    # Initialize variables
    local gccont=`read_opt_value_from_func_args "-gcc" "$@"`
    local npileup=`read_opt_value_from_func_args "-npileup" "$@"`
    local tpileup=`read_opt_value_from_func_args "-tpileup" "$@"`
    local contig=`read_opt_value_from_func_args "-contig" "$@"`
    local outfile=`read_opt_value_from_func_args "-outfile" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || return 1

    # Generate seqz file
    logmsg "* Generating seqz file (contig $contig)..."
    sequenza-utils bam2seqz --pileup -gc "${gccont}" -n "${npileup}" -t "${tpileup}" | "${GZIP}" > "${outfile}" ; pipe_fail || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
parallel_bam2seqz_conda_envs()
{
    define_conda_env sequenza sequenza.yml
}

########
seqzmerge_plus_sequenza_explain_cmdline_opts()
{
    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
seqzmerge_plus_sequenza_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # Define the -out-processdir option, the output directory for the process
    define_opt "-out-processdir" "${process_outdir}" optlist || return 1

    # Get seqz directory
    local seqzdir
    seqzdir=`get_process_outdir_adaptive parallel_bam2seqz`
    define_opt "-seqzdir" "${seqzdir}" optlist || return 1

    # -lc option
    define_cmdline_infile_opt "$cmdline" "-lc" optlist || return 1

    save_opt_list optlist
}

########
seqzmerge()
{
    local clist=$1
    local seqzdir=$2

    local contigs
    contigs=`get_contig_list_from_file $clist` || return 1
    local filenames=""
    local contig
    for contig in ${contigs}; do
        local seqzfname="${seqzdir}"/${contig}_seqz.gz
        if [ "$filenames" = "" ]; then
            filenames="${seqzfname}"
        else
            filenames="${filenames} ${seqzfname}"
        fi
    done

    # NOTE: the use of the bgzip command requires to have the sequenza
    # environment activated (since it has tabix installed)
    "${ZCAT}" "${filenames}" | "$AWK" '{if (NR!=1 && $1 != "chromosome") {print $0}}' | bgzip ; pipe_fail || return 1
}

########
seqzmerge_plus_sequenza()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local seqzdir=`read_opt_value_from_func_args "-seqzdir" "$@"`
    local clist=`read_opt_value_from_func_args "-lc" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || return 1

    # Merge seqz files
    logmsg "* Merging seqz files..."
    seqzmerge "${clist}" "${seqzdir}"  > "${process_outd}"/merged_seqz.gz || return 1

    logmsg "* Applying tabix over merged seqz file..."
    tabix -f -s 1 -b 2 -e 2 -S 1 "${process_outd}"/merged_seqz.gz || return 1

    # Execute sequenza
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing sequenza..."
    Rscript "${genodebasher_bindir}"/run_sequenza -s "${process_outd}"/merged_seqz.gz -o "${process_outd}" 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
seqzmerge_plus_sequenza_conda_envs()
{
    define_conda_env sequenza sequenza.yml
}

########
lumpy_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -lx option
    description="File with regions to exclude in bed format"
    explain_cmdline_opt "-lx" "<string>" "$description"
}

########
lumpy_define_opts()
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

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -lx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-lx" ${NOFILE} optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
get_lumpyexpress_x_opt()
{
    local value=$1

    if [ "${value}" = ${NOFILE} ]; then
        echo ""
    else
        echo "-x ${value}"
    fi
}

########
lumpy()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local exclude=`read_opt_value_from_func_args "-lx" "$@"`

    if [ -z "${LUMPY_HOME_DIR}" ]; then
        # Activate conda environment
        logmsg "* Activating conda environment..."
        conda activate lumpy 2>&1 || return 1

        logmsg "* Executing lumpyexpress..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        lumpyexpress -B "${tumorbam}","${normalbam}" ${x_opt} -o "${process_outd}"/out.vcf || return 1

        # Deactivate conda environment
        logmsg "* Deactivating conda environment..."
        conda deactivate 2>&1
    else
        logmsg "* Executing lumpyexpress..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        ${LUMPY_HOME_DIR}/bin/lumpyexpress -B "${tumorbam}","${normalbam}" ${x_opt} -o "${process_outd}"/out.vcf || return 1
    fi
}

########
lumpy_conda_envs()
{
    define_conda_env lumpy lumpy.yml
}

########
parallel_lumpy_explain_cmdline_opts()
{
    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"

    # -lx option
    description="File with regions to exclude in bed format"
    explain_cmdline_opt "-lx" "<string>" "$description"
}

########
parallel_lumpy_define_opts()
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

    # -lx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-lx" ${NOFILE} optlist || return 1

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
        tumorbam="${abs_splitdir}"/tumor_${contig}.bam
        define_opt "-tumorbam" "${tumorbam}" specific_optlist || return 1
        define_opt "-contig" ${contig} specific_optlist || return 1
        define_opt "-outfile" "${process_outdir}"/out${contig}.vcf specific_optlist || return 1

        save_opt_list specific_optlist
    done
}

########
parallel_lumpy()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local contig=`read_opt_value_from_func_args "-contig" "$@"`
    local exclude=`read_opt_value_from_func_args "-lx" "$@"`
    local outfile=`read_opt_value_from_func_args "-outfile" "$@"`

    if [ -z "${LUMPY_HOME_DIR}" ]; then
        # Activate conda environment
        logmsg "* Activating conda environment (lumpy)..."
        conda activate lumpy 2>&1 || return 1

        logmsg "* Executing lumpyexpress (contig $contig)..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        lumpyexpress -B "${tumorbam}","${normalbam}" ${x_opt} -T "${process_outd}"/tmp_${contig} -o "${outfile}" || return 1

        # Deactivate conda environment
        logmsg "* Deactivating conda environment..."
        conda deactivate 2>&1
    else
        logmsg "* Executing lumpyexpress (contig $contig)..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        "${LUMPY_HOME_DIR}"/bin/lumpyexpress -B "${tumorbam}","${normalbam}" ${x_opt} -T "${process_outd}"/tmp_${contig} -o "${outfile}" || return 1
    fi
}

########
parallel_lumpy_conda_envs()
{
    define_conda_env lumpy lumpy.yml
}

########
smoove_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -lx option
    description="File with regions to exclude in bed format"
    explain_cmdline_opt "-lx" "<string>" "$description"
}

########
smoove_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -lx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-lx" ${NOFILE} optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
get_smoove_exclude_opt()
{
    local value=$1

    if [ "${value}" = ${NOFILE} ]; then
        echo ""
    else
        echo "--exclude ${value}"
    fi
}

########
smoove()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local exclude=`read_opt_value_from_func_args "-lx" "$@"`
    local cpus=`read_opt_value_from_func_args "-cpus" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate smoove 2>&1 || return 1

    logmsg "* Executing smoove..."
    export TMPDIR="${process_outd}"
    local exclude_opt=`get_smoove_exclude_opt ${exclude}`
    local project_name="smoove"
    command smoove call --outdir "${process_outd}" ${exclude_opt} --name ${project_name} --fasta "${ref}" -p ${cpus} --genotype "${normalbam}" "${tumorbam}" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
smoove_conda_envs()
{
    define_conda_env smoove smoove.yml
}

########
delly_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -dx option
    description="File with regions to exclude in bed format"
    explain_cmdline_opt "-dx" "<string>" "$description"
}

########
delly_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -dx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-dx" ${NOFILE} optlist || return 1

    # Save option list
    save_opt_list optlist
}

########
get_delly_x_opt()
{
    local value=$1

    if [ "${value}" = ${NOFILE} ]; then
        echo ""
    else
        echo "-x ${value}"
    fi
}

########
delly()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local exclude=`read_opt_value_from_func_args "-dx" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment (delly)..."
    conda activate delly 2>&1 || return 1

    logmsg "* Executing delly..."
    # "command" built-in is used here to execute the "delly" program
    # instead of the "delly" function
    local x_opt=`get_delly_x_opt ${exclude}`
    command delly call -g "${ref}" ${x_opt} -o "${process_outd}"/out.bcf "${tumorbam}" "${normalbam}" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bcftools)..."
    conda activate bcftools 2>&1 || return 1

    # Convert bcf output to vcf
    logmsg "* Converting bcf output into vcf..."
    bcftools view "${process_outd}"/out.bcf > "${process_outd}"/out.vcf || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
parallel_delly_conda_envs()
{
    define_conda_env bcftools bcftools.yml
    define_conda_env delly delly.yml
}

########
parallel_delly_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -dx option
    description="File with regions to exclude in bed format"
    explain_cmdline_opt "-dx" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_delly_define_opts()
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

    # -dx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-dx" ${NOFILE} optlist || return 1

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
        tumorbam="${abs_splitdir}"/tumor_${contig}.bam
        define_opt "-tumorbam" "${tumorbam}" specific_optlist || return 1
        define_opt "-contig" $contig specific_optlist || return 1
        save_opt_list specific_optlist
    done
}

########
parallel_delly()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local contig=`read_opt_value_from_func_args "-contig" "$@"`
    local exclude=`read_opt_value_from_func_args "-dx" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment (delly)..."
    conda activate delly 2>&1 || return 1

    logmsg "* Executing delly (contig $contig)..."
    # "command" built-in is used here to execute the "delly" program
    # instead of the "delly" function
    local x_opt=`get_delly_x_opt ${exclude}`
    command delly call -g "$ref" ${x_opt} -o "${process_outd}"/out${contig}.bcf "${tumorbam}" "${normalbam}" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bcftools)..."
    conda activate bcftools 2>&1 || return 1

    # Convert bcf output to vcf
    logmsg "* Converting bcf output into vcf..."
    bcftools view "${process_outd}"/out${contig}.bcf > "${process_outd}"/out${contig}.vcf || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
parallel_delly_conda_envs()
{
    define_conda_env bcftools bcftools.yml
    define_conda_env delly delly.yml
}

########
parallel_svtyper_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"
}

########
parallel_svtyper_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local process_spec=$2
    local process_name=$3
    local process_outdir=$4
    local optlist=""

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; return 1; }

    # Determine vcf directory
    local vcfdir
    vcfdir=`get_process_outdir_adaptive parallel_lumpy`

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || return 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        define_opt "-contig" $contig specific_optlist || return 1
        vcf="${vcfdir}"/out${contig}.vcf
        define_opt "-vcf" "$vcf" specific_optlist || return 1
        define_opt "-outfile"  "${process_outdir}"/out${contig}.vcf specific_optlist || return 1
        save_opt_list specific_optlist
    done
}

########
parallel_svtyper()
{
    # Initialize variables
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local contig=`read_opt_value_from_func_args "-contig" "$@"`
    local vcf=`read_opt_value_from_func_args "-vcf" "$@"`
    local outfile=`read_opt_value_from_func_args "-outfile" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment (svtyper)..."
    conda activate svtyper 2>&1 || return 1

    # Execute svtyper
    logmsg "* Executing svtyper (contig $contig)..."
    svtyper -i "${vcf}" -B "${tumorbam}","${normalbam}" > "${outfile}" || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
parallel_svtyper_conda_envs()
{
    define_conda_env svtyper svtyper.yml
}

########
msisensor_pro_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading processes have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"
}

########
msisensor_pro_define_opts()
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

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || return 1
    define_opt "-normalbam" "$normalbam" optlist || return 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || return 1
    define_opt "-tumorbam" "$tumorbam" optlist || return 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_process_spec "$process_spec"` || return 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
msisensor_pro()
{
    # Initialize variables
    local process_outd=`read_opt_value_from_func_args "-out-processdir" "$@"`
    local ref=`read_opt_value_from_func_args "-r" "$@"`
    local normalbam=`read_opt_value_from_func_args "-normalbam" "$@"`
    local tumorbam=`read_opt_value_from_func_args "-tumorbam" "$@"`
    local cpus=`read_opt_value_from_func_args "-cpus" "$@"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate msisensor-pro 2>&1 || return 1

    # Create homopolymer and microsatellites file
    logmsg "* Executing msisensor-pro scan..."
    msisensor-pro scan -d "${ref}" -o "${process_outd}"/msisensor_pro.list 2>&1 || return 1

    # Run Msisensor_Pro analysis
    logmsg "* Executing msisensor-pro msi..."
    msisensor-pro msi -d "${process_outd}"/msisensor_pro.list -n "${normalbam}" -t "${tumorbam}" -o "${process_outd}"/output -l 1 -q 1 -b ${cpus} 2>&1 || return 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
msisensor_pro_conda_envs()
{
    define_conda_env msisensor-pro msisensor_pro.yml
}
