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

DEFAULT_MIN_SEQ_DEPTH_FACETS_PREPROC=35

######################
# BAM ANALYSIS STEPS #
######################

########
manta_germline_document()
{
    step_description "Analyzes a normal \`bam\` file using Manta."
}

########
manta_germline_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
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
    local stepspec=$2
    local optlist=""
    
    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -cr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
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
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local callregf=`read_opt_value_from_line "$*" "-cr"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate manta 2>&1 || exit 1

    # Configure Manta
    logmsg "* Executing configManta.py..."
    configManta.py --bam ${normalbam} --referenceFasta ${ref} ${call_reg_opt} --runDir ${step_outd} 2>&1 || exit 1

    # Execute Manta
    logmsg "* Executing runWorkflow.py..."
    ${step_outd}/runWorkflow.py -m local -j ${cpus} 2>&1 || exit 1

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
    step_description "Analyzes a pair of normal and tumor \`bam\` files using Manta."
}

########
manta_somatic_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
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
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -cr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
manta_somatic()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local callregf=`read_opt_value_from_line "$*" "-cr"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate manta 2>&1 || exit 1
    
    # Configure Manta
    logmsg "* Executing configManta.py..."
    configManta.py --normalBam ${normalbam} --tumorBam ${tumorbam} --referenceFasta ${ref} ${call_reg_opt} --runDir ${step_outd} 2>&1 || exit 1

    # Execute Manta
    logmsg "* Executing runWorkflow.py..."
    ${step_outd}/runWorkflow.py -m local -j ${cpus} 2>&1 || exit 1

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
gridss_somatic_document()
{
    step_description "Analyzes a pair of normal and tumor \`bam\` files using Gridss."
}

########
gridss_somatic_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"

    # -bl option
    description="BED file containing regions to ignore."
    explain_cmdline_opt "-bl" "<string>" "$description"
}

########
gridss_somatic_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -bl option
    define_cmdline_infile_nonmand_opt "$cmdline" "-bl" ${NOFILE} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
get_blacklist_opt()
{
    local blacklistf=$1
    local baseblacklistf=`$BASENAME ${blacklistf}`

    if [ "${baseblacklistf}" = ${NOFILE} -o "${blacklistf}" = "" ]; then
        echo ""
    else
        echo "--blacklist ${blacklistf}"
    fi
}

########
gridss_somatic()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local blacklist=`read_opt_value_from_line "$*" "-bl"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Define --blacklist option
    blacklist_opt=`get_blacklist_opt "${blacklist}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gridss 2>&1 || exit 1

    # Create working directory
    mkdir ${step_outd}/workingdir

    # Execute Gridss
    logmsg "* Executing Gridss..."
    gridss --reference ${ref} --output ${step_outd}/output.vcf.gz --workingdir ${step_outd}/workingdir --assembly ${step_outd}/assembly.bam --threads ${cpus} ${blacklist_opt} ${normalbam} ${tumorbam} 2>&1 || exit 1

    # Remove working directory
    rm -rf ${step_outd}/workingdir

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
gridss_somatic_conda_envs()
{
    define_conda_env gridss gridss.yml
}

########
strelka_germline_document()
{
    step_description "Analyzes a normal \`bam\` files using Strelka."
}

########
strelka_germline_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
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
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -cr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

    # Get germline snvs summary directory
    local abs_sumdir=`get_absolute_shdirname ${GERM_SNVS_SUM_DIR_BASENAME}`

    # -summarydir option
    define_opt "-summarydir" ${abs_sumdir} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
strelka_germline()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local callregf=`read_opt_value_from_line "$*" "-cr"`
    local summarydir=`read_opt_value_from_line "$*" "-summarydir"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate strelka 2>&1 || exit 1

    # Configure Strelka
    logmsg "* Executing configureStrelkaGermlineWorkflow.py..."
    configureStrelkaGermlineWorkflow.py --bam ${normalbam} --referenceFasta ${ref} ${call_reg_opt} --runDir ${step_outd} 2>&1 || exit 1

    # Execute Strelka
    logmsg "* Executing runWorkflow.py..."
    ${step_outd}/runWorkflow.py -m local -j ${cpus} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Create file in summary directory
    local label=strelka_germline
    vcf=${step_outd}/results/variants/variants.vcf.gz
    create_summary_file ${summarydir} ${label} ${vcf}
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
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"    
}

########
platypus_germline_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # Get germline snvs summary directory
    local abs_sumdir=`get_absolute_shdirname ${GERM_SNVS_SUM_DIR_BASENAME}`

    # -summarydir option
    define_opt "-summarydir" ${abs_sumdir} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
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
    local step_outd=$3
    local summarydir=$4
    local cpus=$5

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate platypus 2>&1 || exit 1

    # Run Platypus
    logmsg "* Executing Platypus.py..."
    Platypus.py callVariants --bamFiles=${normalbam} --refFile=${ref} --output=${step_outd}/output.vcf --nCPU=${cpus} --logFileName=${step_outd}/platypus.log --verbosity=1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Create file in summary directory
    local label=platypus_germline
    vcf=${step_outd}/output.vcf
    create_summary_file ${summarydir} ${label} ${vcf}
}

########
platypus_germline_local()
{
    # Initialize variables
    local ref=$1
    local normalbam=$2
    local step_outd=$3
    local summarydir=$4
    local cpus=$5

    # Run Platypus
    logmsg "* Executing Platypus.py..."
    python ${PLATYPUS_HOME_DIR}/bin/Platypus.py callVariants --bamFiles=${normalbam} --refFile=${ref} --nCPU=${cpus} --output=${step_outd}/output.vcf --verbosity=1 2>&1 || exit 1

    # Create file in summary directory
    local label=platypus_germline
    vcf=${step_outd}/output.vcf
    create_summary_file ${summarydir} ${label} ${vcf}
}

########
platypus_germline()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local summarydir=`read_opt_value_from_line "$*" "-summarydir"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    if [ -z "${PLATYPUS_HOME_DIR}" ]; then
        platypus_germline_conda ${ref} ${normalbam} ${step_outd} ${summarydir} ${cpus}
    else
        platypus_germline_local ${ref} ${normalbam} ${step_outd} ${summarydir} ${cpus}
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
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"    
}

########
gatk_haplotypecaller_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

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
gatk_haplotypecaller()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    local mem=`read_opt_value_from_line "$*" "-mem"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gatk4 2>&1 || exit 1

    # Run gatk HaplotypeCaller
    logmsg "* Executing gatk HaplotypeCaller..."
    gatk --java-options "-Xmx${mem}" HaplotypeCaller -R ${ref} -I ${normalbam} -O ${step_outd}/output.g.vcf.gz -ERC GVCF --native-pair-hmm-threads ${cpus} --tmp-dir ${step_outd} || exit 1

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
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
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
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -manta-outd option
    local manta_dep=`find_dependency_for_step "${stepspec}" manta_somatic`
    if [ ${manta_dep} != ${DEP_NOT_FOUND} ]; then
        local manta_outd=`get_outd_for_dep "${manta_dep}"`
        define_opt "-manta-outd" ${manta_outd} optlist || exit 1
    fi
    
    # -cr option
    define_cmdline_infile_nonmand_opt "$cmdline" "-cr" ${NOFILE} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
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
        if [ -f ${manta_indel_file} ]; then
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
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local manta_outd=`read_opt_value_from_line "$*" "-manta-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local callregf=`read_opt_value_from_line "$*" "-cr"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Define --indelCandidates option if output from Manta is available
    indel_cand_opt=`get_indel_cand_opt "${manta_outd}"`

    # Define --callRegions option
    call_reg_opt=`get_callreg_opt "${callregf}"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate strelka 2>&1 || exit 1

    # Configure Strelka
    logmsg "* Executing configureStrelkaSomaticWorkflow.py..."
    configureStrelkaSomaticWorkflow.py --normalBam ${normalbam} --tumorBam ${tumorbam} --referenceFasta ${ref} ${indel_cand_opt} ${call_reg_opt} --runDir ${step_outd} 2>&1 || exit 1

    # Execute Strelka
    logmsg "* Executing runWorkflow.py..."
    ${step_outd}/runWorkflow.py -m local -j ${cpus} 2>&1 || exit 1

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
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"    

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"
}

########
mutect2_somatic_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

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
mutect2_somatic()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    local mem=`read_opt_value_from_line "$*" "-mem"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate gatk4 2>&1 || exit 1

    # Run Mutect2
    logmsg "* Executing gatk Mutect2..."
    gatk --java-options "-Xmx${mem}" Mutect2 -R ${ref} -I ${normalbam} -I ${tumorbam} -O ${step_outd}/somatic.vcf.gz --native-pair-hmm-threads ${cpus} --tmp-dir ${step_outd} || exit 1

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
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"    

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"
}

########
lofreq_somatic_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
lofreq_somatic()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate lofreq 2>&1 || exit 1

    # Run lofreq somatic
    logmsg "* Executing lofeq somatic..."
    lofreq somatic -f ${ref} -n ${normalbam} -t ${tumorbam} -o ${step_outd}/somatic.vcf.gz --threads ${cpus} -o ${step_outd}/out_ || exit 1

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
    step_description "Analyzes a pair of normal and tumor \`bam\` files using CNVkit."
}

########
cnvkit_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"
}

########
cnvkit_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist
}

########
cnvkit()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    
    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate cnvkit 2>&1 || exit 1

    # Run cnvkit
    logmsg "* Executing cnvkit.py..."
    cd ${step_outd} # This is done since current implementation of
                    # cnvkit generates a file in current directory
                    # (genref.bed). Changing directory avoids possible
                    # racing conditions
    cnvkit.py batch ${tumorbam} -n ${normalbam} -m wgs -f ${ref}  -d ${step_outd} -p ${cpus} 2>&1 || exit 1

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
wisecondorx_explain_cmdline_opts()
{
    # -wcr option
    description="Reference file in npz format for WisecondorX"
    explain_cmdline_req_opt "-wcr" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    
}

########
wisecondorx_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -wcr option
    define_cmdline_infile_opt "$cmdline" "-wcr" optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
wisecondorx()
{
    # Initialize variables
    local wcref=`read_opt_value_from_line "$*" "-wcr"`
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`
    
    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate wisecondorx 2>&1 || exit 1

    # Convert tumor bam file into npz
    logmsg "* Executing WisecondorX convert..."
    local BINSIZE=5000
    WisecondorX convert ${tumorbam} ${step_outd}/tumor.npz --binsize $BINSIZE 2>&1 || exit 1
    
    # Use WisecondorX for prediction
    logmsg "* Executing WisecondorX predict..."
    WisecondorX predict ${step_outd}/tumor.npz ${wcref} ${step_outd}/out 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
wisecondorx_conda_envs()
{
    define_conda_env wisecondorx wisecondorx.yml
}

########
snp_pileup_plus_facets_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
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
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -sv option
    define_cmdline_infile_opt "$cmdline" "-sv" optlist || exit 1

    # -md option
    define_cmdline_nonmandatory_opt "$cmdline" "-md" ${DEFAULT_MIN_SEQ_DEPTH_FACETS_PREPROC} optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
snp_pileup_plus_facets()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local snpvcf=`read_opt_value_from_line "$*" "-sv"`
    local mindepth=`read_opt_value_from_line "$*" "-md"`

    # Activate conda environment
    logmsg "* Activating conda environment (snp-pileup)..."
    conda activate snp-pileup 2>&1 || exit 1

    # Execute snp-pileup
    logmsg "* Executing snp-pileup..."
    snp-pileup ${snpvcf} ${step_outd}/snp-pileup-counts.csv ${normalbam} ${tumorbam} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment if needed
    logmsg "* Activating conda environment (facets)..."
    conda activate facets 2>&1 || exit 1
            
    # Execute facets
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing facets..."
    Rscript ${biopanpipe_bindir}/run_facets -c ${step_outd}/snp-pileup-counts.csv -d ${mindepth} -o ${step_outd} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Compress snp-pileup file
    logmsg "* Compressing snp-pileup file..."
    ${GZIP} ${step_outd}/snp-pileup-counts.csv || exit 1
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
    local stepspec=$2
    local optlist=""
    
    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1
    
    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`

    # -outfile option
    define_opt "-outfile" ${abs_datadir}/sequenza_gccfile.txt.gz optlist || exit 1

    # Save option list
    save_opt_list optlist
}

########
gen_sequenza_gcc()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local outfile=`read_opt_value_from_line "$*" "-outfile"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate sequenza 2>&1 || exit 1

    # Generate GC content file
    logmsg "* Generating GC content file..."    
    sequenza-utils gc_wiggle -w 50 -f $ref -o ${step_outd}/sequenza_gccfile.txt.gz || exit 1
    
    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Move result file to final location
    mv ${step_outd}/sequenza_gccfile.txt.gz $outfile || exit 1
}

########
sequenza_explain_cmdline_opts()
{
    # -gcc option
    description="GC content wiggle file for sequenza (required if no gen_sequenza_gcc step is defined)"
    explain_cmdline_opt "-gcc" "<string>" "$description"
}

########
get_gcc_filename()
{
    local cmdline=$1
    local stepspec=$2
    local given=0

    gccfile=`read_opt_value_from_line "$cmdline" "-gcc"` && given=1
    if [ $given -eq 1 ]; then
        # -gcc option was given
        file_exists $gccfile || { errmsg "file $gccfile does not exist" ; return 1; }
        echo $gccfile
    else
        # Check if gen_sequenza_gcc step dependency was defined
        local gen_sequenza_gcc_dep=`find_dependency_for_step "${stepspec}" gen_sequenza_gcc`
        if [ ${gen_sequenza_gcc_dep} != ${DEP_NOT_FOUND} ]; then
            local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`
            gccfile=${abs_datadir}/sequenza_gccfile.txt.gz
            echo $gccfile
            return 0            
        else
            errmsg "-gcc or dependency with gen_sequenza_gcc step should be given"
            return 1
        fi            
    fi
}

########
sequenza_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -gcc option
    local gccfile
    gccfile=`get_gcc_filename "$cmdline" "$stepspec"` || exit 1
    define_opt "-gcc" $gccfile optlist || exit 1

    # Get normal pileup file
    local npileupdir
    npileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" samtools_mpileup_norm_bam` || { errmsg "Error: dependency samtools_mpileup_norm_bam not defined for sequenza"; exit 1; }
    local npileup=${npileupdir}/normal.pileup.gz
    define_opt "-npileup" ${npileup} optlist || exit 1

    # Get tumor pileup file
    tpileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" samtools_mpileup_tum_bam` || { errmsg "Error: dependency samtools_mpileup_tum_bam not defined for sequenza"; exit 1; }
    tpileup=${tpileupdir}/tumor.pileup.gz
    define_opt "-tpileup" ${tpileup} optlist || exit 1

    # Save option list
    save_opt_list optlist    
}

########
sequenza()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local gccont=`read_opt_value_from_line "$*" "-gcc"`
    local npileup=`read_opt_value_from_line "$*" "-npileup"`
    local tpileup=`read_opt_value_from_line "$*" "-tpileup"`

    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || exit 1
    
    # Generate seqz file
    logmsg "* Generating seqz file..."
    sequenza-utils bam2seqz --pileup -gc ${gccont} -n ${npileup} -t ${tpileup} | ${GZIP} > ${step_outd}/seqz.gz ; pipe_fail || exit 1

    # Execute sequenza
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing sequenza..."
    Rscript ${biopanpipe_bindir}/run_sequenza -s ${step_outd}/seqz.gz -o ${step_outd} 2>&1 || exit 1

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
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -gcc option
    local gccfile
    gccfile=`get_gcc_filename "$cmdline" "$stepspec"` || exit 1
    define_opt "-gcc" $gccfile optlist || exit 1

    # Get normal pileup directory
    local npileupdir
    npileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" parallel_samtools_mpileup_norm_bam` || { errmsg "Error: dependency parallel_samtools_mpileup_norm_bam not defined for parallel_bam2seqz"; exit 1; }

    # Get tumor pileup directory
    local tpileupdir
    tpileupdir=`get_outd_for_dep_given_stepspec "${stepspec}" parallel_samtools_mpileup_tum_bam` || { errmsg "Error: dependency parallel_samtools_mpileup_tum_bam not defined for parallel_bam2seqz"; exit 1; }

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        npileup=${npileupdir}/normal_${contig}.pileup.gz
        define_opt "-npileup" ${npileup} specific_optlist || exit 1
        tpileup=${tpileupdir}/tumor_${contig}.pileup.gz
        define_opt "-tpileup" ${tpileup} specific_optlist || exit 1
        define_opt "-contig" $contig specific_optlist || exit 1
        
        save_opt_list specific_optlist
    done
}

########
parallel_bam2seqz()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local gccont=`read_opt_value_from_line "$*" "-gcc"`
    local npileup=`read_opt_value_from_line "$*" "-npileup"`
    local tpileup=`read_opt_value_from_line "$*" "-tpileup"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || exit 1
    
    # Generate seqz file
    logmsg "* Generating seqz file (contig $contig)..."
    sequenza-utils bam2seqz --pileup -gc ${gccont} -n ${npileup} -t ${tpileup} | ${GZIP} > ${step_outd}/${contig}_seqz.gz ; pipe_fail || exit 1

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
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # Get seqz directory
    seqzdir=`get_outd_for_dep_given_stepspec "${stepspec}" parallel_bam2seqz` || { errmsg "Error: dependency parallel_bam2seqz not defined for seqzmerge_plus_sequenza"; exit 1; }
    define_opt "-seqzdir" ${seqzdir} optlist || exit 1

    # -lc option
    define_cmdline_infile_opt "$cmdline" "-lc" optlist || exit 1

    save_opt_list optlist
}

########
seqzmerge()
{
    local clist=$1
    local seqzdir=$2

    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local filenames=""
    local contig
    for contig in ${contigs}; do
        local seqzfname=${seqzdir}/${contig}_seqz.gz
        if [ "$filenames" = "" ]; then
            filenames=${seqzfname}
        else
            filenames="${filenames} ${seqzfname}"
        fi
    done

    # NOTE: the use of the bgzip command requires to have the sequenza
    # environment activated (since it has tabix installed)
    ${ZCAT} ${filenames} | $AWK '{if (NR!=1 && $1 != "chromosome") {print $0}}' | bgzip ; pipe_fail || exit 1
}
 
########
seqzmerge_plus_sequenza()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local seqzdir=`read_opt_value_from_line "$*" "-seqzdir"`
    local clist=`read_opt_value_from_line "$*" "-lc"`

    # Activate conda environment
    logmsg "* Activating conda environment (sequenza)..."
    conda activate sequenza 2>&1 || exit 1

    # Merge seqz files
    logmsg "* Merging seqz files..."
    seqzmerge ${clist} ${seqzdir}  > ${step_outd}/merged_seqz.gz || exit 1

    logmsg "* Applying tabix over merged seqz file..."
    tabix -f -s 1 -b 2 -e 2 -S 1 ${step_outd}/merged_seqz.gz || exit 1
                
    # Execute sequenza
    # IMPORTANT NOTE: Rscript is used here to ensure that conda's R
    # installation is used (otherwise, general R installation given in
    # shebang directive would be executed)
    logmsg "* Executing sequenza..."
    Rscript ${biopanpipe_bindir}/run_sequenza -s ${step_outd}/merged_seqz.gz -o ${step_outd} 2>&1 || exit 1

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
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
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
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -lx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-lx" ${NOFILE} optlist || exit 1

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
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local exclude=`read_opt_value_from_line "$*" "-lx"`

    if [ -z "${LUMPY_HOME_DIR}" ]; then
        # Activate conda environment
        logmsg "* Activating conda environment..."
        conda activate lumpy 2>&1 || exit 1

        logmsg "* Executing lumpyexpress..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        lumpyexpress -B ${tumorbam},${normalbam} ${x_opt} -o ${step_outd}/out.vcf || exit 1
        
        # Deactivate conda environment
        logmsg "* Deactivating conda environment..."
        conda deactivate 2>&1
    else
        logmsg "* Executing lumpyexpress..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        ${LUMPY_HOME_DIR}/bin/lumpyexpress -B ${tumorbam},${normalbam} ${x_opt} -o ${step_outd}/out.vcf || exit 1
    fi
}

########
lumpy_conda_envs()
{
    define_conda_env lumpy lumpy.yml
}

########
parallel_exclude_plus_lumpy_explain_cmdline_opts()
{
    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"   
}

########
parallel_exclude_plus_lumpy_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
get_contig_names_from_bam()
{
    local bam=$1

    conda activate samtools 2>&1 || return 1
    samtools idxstats $bam | $AWK '{printf"%s\n",$1}' ; pipe_fail || return 1
    conda deactivate
}

########
gen_exclusion_bed_given_bam()
{
    local bam=$1
    local contig=$2

    get_contig_names_from_bam $bam | $AWK -v contig=$contig '{if($1!=contig){printf"%s\n",$1}}' ; pipe_fail || return 1
}

########
parallel_exclude_plus_lumpy()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`

    # Generate exclusion bed file
    gen_exclusion_bed_given_bam ${normalbam} ${contig} > ${step_outd}/${contig}.bed || exit 1

    if [ -z "${LUMPY_HOME_DIR}" ]; then
        # Activate conda environment
        logmsg "* Activating conda environment (lumpy)..."
        conda activate lumpy 2>&1 || exit 1
    
        logmsg "* Executing lumpyexpress (contig $contig)..."
        lumpyexpress -B ${tumorbam},${normalbam} -T ${step_outd}/tmp_${contig} -x ${step_outd}/${contig}.bed -o ${step_outd}/out${contig}.vcf || exit 1

        # Deactivate conda environment
        logmsg "* Deactivating conda environment..."
        conda deactivate 2>&1
    else
        logmsg "* Executing lumpyexpress (contig $contig)..."
        ${LUMPY_HOME_DIR}/bin/lumpyexpress -B ${tumorbam},${normalbam} -T ${step_outd}/tmp_${contig} -x ${step_outd}/${contig}.bed -o ${step_outd}/out${contig}.vcf || exit 1
    fi
}

########
parallel_exclude_plus_lumpy_conda_envs()
{
    define_conda_env samtools samtools.yml
    define_conda_env lumpy lumpy.yml
}

########
check_contig_does_not_exist_given_log_file()
{
    local logfile=$1

    if $GREP "does not exist" ${logfile} >/dev/null 2>&1; then
        return 0
    else
        return 1
    fi
}

########
filter_bam_contig_samtools()
{
    local inbam=$1
    local contig=$2
    local outbam=$3    
    local error=0
    
    samtools view -h -O BAM $inbam $contig > ${outbam} 2> ${outbam}.log || error=1

    if [ $error -eq 1 ]; then
        if check_contig_does_not_exist_given_log_file ${outbam}.log; then
            errmsg "Warning: contig ${contig} does not exist in ${inbam} file (see ${outbam}.log)"
            return 0
        else
            errmsg "Error while filtering ${contig} in ${inbam} file (see ${outbam}.log)"
            return 1
        fi
    else
        return 0
    fi
}

########
filter_bam_contig_sambamba()
{
    local inbam=$1
    local contig=$2
    local outbam=$3    
    local error=0
    
    sambamba view -h -f bam $inbam $contig > ${outbam} 2> ${outbam}.log || error=1

    if [ $error -eq 1 ]; then
        if check_contig_does_not_exist_given_log_file ${outbam}.log; then
            errmsg "Warning: contig ${contig} does not exist in ${inbam} file (see ${outbam}.log)"
            return 0
        else
            errmsg "Error while filtering ${contig} in ${inbam} file (see ${outbam}.log)"
            return 1
        fi
    else
        return 0
    fi
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
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -lx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-lx" ${NOFILE} optlist || exit 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        normalbam=${abs_datadir}/normal_${contig}.bam
        define_opt "-normalbam" ${normalbam} specific_optlist || exit 1
        tumorbam=${abs_datadir}/tumor_${contig}.bam
        define_opt "-tumorbam" ${tumorbam} specific_optlist || exit 1
        define_opt "-contig" ${contig} specific_optlist || exit 1
        
        save_opt_list specific_optlist
    done
}

########
parallel_lumpy()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local exclude=`read_opt_value_from_line "$*" "-lx"`

    if [ -z "${LUMPY_HOME_DIR}" ]; then
        # Activate conda environment
        logmsg "* Activating conda environment (lumpy)..."
        conda activate lumpy 2>&1 || exit 1
        
        logmsg "* Executing lumpyexpress (contig $contig)..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        lumpyexpress -B ${tumorbam},${normalbam} ${x_opt} -T ${step_outd}/tmp_${contig} -o ${step_outd}/out${contig}.vcf || exit 1
        
        # Deactivate conda environment
        logmsg "* Deactivating conda environment..."
        conda deactivate 2>&1
    else
        logmsg "* Executing lumpyexpress (contig $contig)..."
        local x_opt=`get_lumpyexpress_x_opt ${exclude}`
        ${LUMPY_HOME_DIR}/bin/lumpyexpress -B ${tumorbam},${normalbam} ${x_opt} -T ${step_outd}/tmp_${contig} -o ${step_outd}/out${contig}.vcf || exit 1
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
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
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
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -lx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-lx" ${NOFILE} optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
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
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local exclude=`read_opt_value_from_line "$*" "-lx"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate smoove 2>&1 || exit 1

    logmsg "* Executing smoove..."
    export TMPDIR=${step_outd}
    local exclude_opt=`get_smoove_exclude_opt ${exclude}`
    local project_name="smoove"
    command smoove call --outdir ${step_outd} ${exclude_opt} --name ${project_name} --fasta ${ref} -p ${cpus} --genotype ${normalbam} ${tumorbam} || exit 1
        
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
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
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
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -dx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-dx" ${NOFILE} optlist || exit 1

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
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local exclude=`read_opt_value_from_line "$*" "-dx"`

    # Activate conda environment
    logmsg "* Activating conda environment (delly)..."
    conda activate delly 2>&1 || exit 1

    logmsg "* Executing delly..."
    # "command" built-in is used here to execute the "delly" program
    # instead of the "delly" function
    local x_opt=`get_delly_x_opt ${exclude}`
    command delly call -g ${ref} ${x_opt} -o ${step_outd}/out.bcf ${tumorbam} ${normalbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bcftools)..."
    conda activate bcftools 2>&1 || exit 1

    # Convert bcf output to vcf
    logmsg "* Converting bcf output into vcf..."
    bcftools view ${step_outd}/out.bcf > ${step_outd}/out.vcf || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
parallel_lumpy_conda_envs()
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
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # Get data directory
    local abs_datadir=`get_absolute_shdirname ${DATADIR_BASENAME}`

    # -dx option
    define_cmdline_infile_nonmand_opt "$cmdline" "-dx" ${NOFILE} optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        normalbam=${abs_datadir}/normal_${contig}.bam
        define_opt "-normalbam" ${normalbam} specific_optlist || exit 1
        tumorbam=${abs_datadir}/tumor_${contig}.bam
        define_opt "-tumorbam" ${tumorbam} specific_optlist || exit 1
        define_opt "-contig" $contig specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_delly()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local exclude=`read_opt_value_from_line "$*" "-dx"`
    
    # Activate conda environment
    logmsg "* Activating conda environment (delly)..."
    conda activate delly 2>&1 || exit 1
    
    logmsg "* Executing delly (contig $contig)..."
    # "command" built-in is used here to execute the "delly" program
    # instead of the "delly" function
    local x_opt=`get_delly_x_opt ${exclude}`
    command delly call -g $ref ${x_opt} -o ${step_outd}/out${contig}.bcf ${tumorbam} ${normalbam} || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1

    # Activate conda environment
    logmsg "* Activating conda environment (bcftools)..."
    conda activate bcftools 2>&1 || exit 1

    # Convert bcf output to vcf
    logmsg "* Converting bcf output into vcf..."
    bcftools view ${step_outd}/out${contig}.bcf > ${step_outd}/out${contig}.vcf || exit 1

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
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"    

    # -lc option
    description="File with list of contig names to process"
    explain_cmdline_req_opt "-lc" "<string>" "$description"   
}

########
get_vcfdir_for_svtyper()
{
    local stepspec=$1

    # Check dependency with parallel_lumpy
    local parallel_lumpy_dep=`find_dependency_for_step "${stepspec}" parallel_lumpy`
    if [ ${parallel_lumpy_dep} != ${DEP_NOT_FOUND} ]; then
        local vcfdir=`get_outd_for_dep "${parallel_lumpy_dep}"`
        echo $vcfdir
        return 0
    fi

    return 1
}

########
parallel_svtyper_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # Get name of contig list file
    local clist
    clist=`read_opt_value_from_line "$cmdline" "-lc"` || { errmsg "Error: -lc option not found"; exit 1; }

    # Determine vcf directory
    vcfdir=`get_vcfdir_for_svtyper "${stepspec}"` || { errmsg "Error: vcf directory for svtyper could not be determined"; exit 1; }
    
    # Generate option lists for each contig
    local contigs
    contigs=`get_contig_list_from_file $clist` || exit 1
    local contig
    for contig in ${contigs}; do
        local specific_optlist=${optlist}
        define_opt "-contig" $contig specific_optlist || exit 1
        vcf=${vcfdir}/out${contig}.vcf
        define_opt "-vcf" $vcf specific_optlist || exit 1
        save_opt_list specific_optlist
    done
}

########
parallel_svtyper()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local contig=`read_opt_value_from_line "$*" "-contig"`
    local vcf=`read_opt_value_from_line "$*" "-vcf"`

    # Activate conda environment
    logmsg "* Activating conda environment (svtyper)..."
    conda activate svtyper 2>&1 || exit 1

    # Execute svtyper
    logmsg "* Executing svtyper (contig $contig)..."
    svtyper -i ${vcf} -B ${tumorbam},${normalbam} > ${step_outd}/out${contig}.vcf || exit 1
    
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
msisensor_explain_cmdline_opts()
{
    # -r option
    description="Reference genome file"
    explain_cmdline_opt "-r" "<string>" "$description"

    # -n option
    description="Normal bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-n" "<string>" "$description"

    # -t option
    description="Tumor bam file (required if no downloading steps have been defined)"
    explain_cmdline_opt "-t" "<string>" "$description"        
}

########
msisensor_define_opts()
{
    # Initialize variables
    local cmdline=$1
    local stepspec=$2
    local optlist=""

    # Define the -step-outd option, the output directory for the step
    local step_outd=`get_step_outdir_given_stepspec "$stepspec"`
    define_opt "-step-outd" ${step_outd} optlist || exit 1

    # -r option
    local genref
    genref=`get_ref_filename "$cmdline"` || exit 1
    define_opt "-r" $genref optlist || exit 1

    # -normalbam option
    local normalbam
    normalbam=`get_normal_bam_filename "$cmdline"` || exit 1
    define_opt "-normalbam" $normalbam optlist || exit 1

    # -tumorbam option
    local tumorbam
    tumorbam=`get_tumor_bam_filename "$cmdline"` || exit 1
    define_opt "-tumorbam" $tumorbam optlist || exit 1

    # -cpus option
    local cpus
    cpus=`extract_cpus_from_stepspec "$stepspec"` || exit 1
    define_opt "-cpus" $cpus optlist

    # Save option list
    save_opt_list optlist    
}

########
msisensor()
{
    # Initialize variables
    local step_outd=`read_opt_value_from_line "$*" "-step-outd"`
    local ref=`read_opt_value_from_line "$*" "-r"`
    local normalbam=`read_opt_value_from_line "$*" "-normalbam"`
    local tumorbam=`read_opt_value_from_line "$*" "-tumorbam"`
    local cpus=`read_opt_value_from_line "$*" "-cpus"`

    # Activate conda environment
    logmsg "* Activating conda environment..."
    conda activate msisensor 2>&1 || exit 1

    # Create homopolymer and microsatellites file
    logmsg "* Executing msisensor scan..."
    command msisensor scan -d ${ref} -o ${step_outd}/msisensor.list 2>&1 || exit 1

    # Run MSIsensor analysis
    logmsg "* Executing msisensor msi..."
    command msisensor msi -d ${step_outd}/msisensor.list -n ${normalbam} -t ${tumorbam} -o ${step_outd}/output -l 1 -q 1 -b ${cpus} 2>&1 || exit 1

    # Deactivate conda environment
    logmsg "* Deactivating conda environment..."
    conda deactivate 2>&1
}

########
msisensor_conda_envs()
{
    define_conda_env msisensor msisensor.yml
}
