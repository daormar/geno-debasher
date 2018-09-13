# *- bash -*

# INCLUDE BASH LIBRARY
. ${bindir}/bam_utils_lib.sh

########
print_desc()
{
    echo "submit_bam_analysis performs analyses given normal and tumor bam files"
    echo "type \"submit_bam_analysis --help\" to get usage information"
}

########
usage()
{
    echo "submit_bam_analysis  -r <string>"
    echo "                     -n <string>|-egan <string> -t <string>|-egat <string>"
    echo "                     -a <string> -g <string> -o <string>"
    echo "                     [-egastr <int>] [-egacred <string>]"
    echo "                     [-debug] [--help]"
    echo ""
    echo "-r <string>          File with reference genome"
    echo "-n <string>          Normal bam file"
    echo "-t <string>          Tumor bam file"
    echo "-egan <string>       EGA id of normal bam file to download"
    echo "-egat <string>       EGA id of tumor bam file to download"
    echo "-a <string>          File with analysis steps to be performed."
    echo "                     Expected format:"
    echo "                      <stepname> <cpus> <mem> <time> <jobdeps=stepname1:...>"
    echo "-g <string>          Sample gender (XX|XY)"
    echo "-o <string>          Output directory"
    echo "-egastr <int>        Number of streams used by the EGA download client"
    echo "                     (50 by default)"
    echo "-egacred <string>    File with EGA download client credentials"
    echo "-debug               After ending, do not delete temporary files"
    echo "                     (for debugging purposes)"
    echo "--help               Display this help and exit"
}

########
read_pars()
{
    r_given=0
    n_given=0
    t_given=0
    egan_given=0
    egat_given=0
    a_given=0
    g_given=0
    gender="XX"
    o_given=0
    egastr_given=0
    egastr=50
    egacred_given=0
    egacred="cred.json"
    debug=0
    while [ $# -ne 0 ]; do
        case $1 in
            "--help") usage
                      exit 1
                      ;;
            "--version") version
                         exit 1
                         ;;
            "-r") shift
                  if [ $# -ne 0 ]; then
                      ref=$1
                      r_given=1
                  fi
                  ;;
            "-n") shift
                  if [ $# -ne 0 ]; then
                      normalbam=$1
                      n_given=1
                  fi
                  ;;
            "-t") shift
                  if [ $# -ne 0 ]; then
                      tumorbam=$1
                      t_given=1
                  fi
                  ;;
            "-egan") shift
                  if [ $# -ne 0 ]; then
                      egaid_normalbam=$1
                      egan_given=1
                  fi
                  ;;
            "-egat") shift
                  if [ $# -ne 0 ]; then
                      egaid_tumorbam=$1
                      egat_given=1
                  fi
                  ;;
            "-a") shift
                  if [ $# -ne 0 ]; then
                      afile=$1
                      a_given=1
                  fi
                  ;;
            "-g") shift
                  if [ $# -ne 0 ]; then
                      gender=$1
                      g_given=1
                  fi
                  ;;
            "-o") shift
                  if [ $# -ne 0 ]; then
                      outd=$1
                      o_given=1
                  fi
                  ;;
            "-egastr") shift
                  if [ $# -ne 0 ]; then
                      egastr=$1
                      egastr_given=1
                  fi
                  ;;
            "-egacred") shift
                  if [ $# -ne 0 ]; then
                      egacred=$1
                      egacred_given=1
                  fi
                  ;;
            "-debug") debug=1
                      debug_opt="-debug"
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
        if [ ! -f ${ref} ]; then
            echo "Error! file ${ref} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${n_given} -eq 0 -a ${egan_given} -eq 0 ]; then
        echo "Error, -n or -egan options should be given" >&2
    fi

    if [ ${n_given} -eq 1 -a ${egan_given} -eq 1 ]; then
        echo "Error, -n and -egan options cannot be given simultaneously" >&2
    fi

    if [ ${t_given} -eq 0 -a ${egat_given} -eq 0 ]; then
        echo "Error, -t or -egat options should be given" >&2
    fi

    if [ ${t_given} -eq 1 -a ${egat_given} -eq 1 ]; then
        echo "Error, -t and -egat options cannot be given simultaneously" >&2
    fi

    if [ ${n_given} -eq 1 ]; then
        if [ ! -f ${normalbam} ]; then
            echo "Error! file ${normalbam} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${t_given} -eq 1 ]; then
        if [ -a ! -f ${tumorbam} ]; then
            echo "Error! file ${tumorbam} does not exist" >&2
            exit 1
        fi
    fi
    
    if [ ${a_given} -eq 0 ]; then   
        echo "Error! -a parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${afile} ]; then
            echo "Error! file ${afile} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${g_given} -eq 0 ]; then   
        echo "Error! -g parameter not given!" >&2
        exit 1
    fi

    if [ ${o_given} -eq 0 ]; then
        echo "Error! -o parameter not given!" >&2
        exit 1
    else
        if [ -d ${outd} ]; then
            echo "Warning! output directory does exist" >&2 
        fi
    fi
}

########
print_pars()
{
    if [ ${r_given} -eq 1 ]; then
        echo "-r is ${ref}" >&2
    fi

    if [ ${n_given} -eq 1 ]; then
        echo "-n is ${normalbam}" >&2
    fi

    if [ ${t_given} -eq 1 ]; then
        echo "-t is ${tumorbam}" >&2
    fi

    if [ ${egan_given} -eq 1 ]; then
        echo "-egan is ${egaid_normalbam}" >&2
    fi

    if [ ${egat_given} -eq 1 ]; then
        echo "-egat is ${egaid_tumorbam}" >&2
    fi

    if [ ${a_given} -eq 1 ]; then
        echo "-a is ${afile}" >&2
    fi

    if [ ${g_given} -eq 1 ]; then
        echo "-g is ${gender}" >&2
    fi

    if [ ${o_given} -eq 1 ]; then
        echo "-o is ${outd}" >&2
    fi
}

########
create_dirs()
{
    mkdir -p ${outd} || { echo "Error! cannot create output directory" >&2; return 1; }

    mkdir -p ${outd}/scripts || { echo "Error! cannot create scripts directory" >&2; return 1; }
}

########
set_bam_filenames()
{
    if [ ${egan_given} -eq 1 ]; then
        normalbam=${outd}/normal.bam
    fi

    if [ ${egat_given} -eq 1 ]; then
        tumorbam=${outd}/tumor.bam
    fi
}

########
get_pars_manta_somatic()
{
    echo "$ref $normalbam $tumorbam $outd $cpus"
}

########
execute_manta_somatic()
{
    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_tumorbam=$3
    local_outd=$4
    local_cpus=$5
    
    # Initialize variables
    MANTA_OUTD=`get_step_dirname ${local_outd} manta_somatic`

    # Activate conda environment
    conda activate manta
    
    # Configure Manta
    configManta.py --normalBam ${local_normalbam} --tumorBam ${local_tumorbam} --referenceFasta ${local_ref} --runDir ${MANTA_OUTD} > ${MANTA_OUTD}/configManta.log 2>&1 || exit 1

    # Execute Manta
    ${MANTA_OUTD}/runWorkflow.py -m local -j ${local_cpus} > ${MANTA_OUTD}/runWorkflow.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${MANTA_OUTD}/finished
}

########
get_pars_strelka_somatic()
{
    echo "$ref $normalbam $tumorbam $outd $cpus"
}

########
execute_strelka_somatic()
{
    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_tumorbam=$3
    local_outd=$4
    local_cpus=$5

    # Initialize variables
    STRELKA_OUTD=`get_step_dirname ${outd} strelka_somatic`

    # Activate conda environment
    conda activate strelka

    # Configure Strelka
    configureStrelkaSomaticWorkflow.py --normalBam ${local_normalbam} --tumorBam ${local_tumorbam} --referenceFasta ${local_ref} --runDir ${STRELKA_OUTD} > ${STRELKA_OUTD}/configureStrelkaSomaticWorkflow.log 2>&1 || exit 1

    # Execute Strelka
    ${STRELKA_OUTD}/runWorkflow.py -m local -j ${local_cpus} > ${STRELKA_OUTD}/runWorkflow.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${STRELKA_OUTD}/finished
}

########
get_pars_msisensor()
{
    echo "$ref $normalbam $tumorbam $outd $cpus"
}

########
execute_msisensor()
{
    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_tumorbam=$3
    local_outd=$4
    local_cpus=$5

    # Initialize variables
    MSISENSOR_OUTD=`get_step_dirname ${outd} msisensor`
    
    # Activate conda environment
    conda activate msisensor

    # Create homopolymer and microsatellites file
    msisensor scan -d ${local_ref} -o ${MSISENSOR_OUTD}/msisensor.list > ${MSISENSOR_OUTD}/msisensor_scan.log 2>&1 || exit 1

    # Run MSIsensor analysis
    msisensor msi -d ${MSISENSOR_OUTD}/msisensor.list -n ${local_normalbam} -t ${local_tumorbam} -o ${MSISENSOR_OUTD}/output -l 1 -q 1 -b ${local_cpus} > ${MSISENSOR_OUTD}/msisensor_msi.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate
    
    # Create file indicating that execution was finished
    touch ${MSISENSOR_OUTD}/finished
}

########
execute_platypus_germline_conda()
{
    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_outd=$3

    # Initialize variables
    PLATYPUS_OUTD=`get_step_dirname ${outd} platypus_germline_conda`

    # Activate conda environment
    conda activate platypus

    # Run Platypus
    Platypus.py callVariants --bamFiles=${local_normalbam} --refFile=${local_ref} --output=${PLATYPUS_OUTD}/output.vcf --verbosity=1 > ${PLATYPUS_OUTD}/platypus.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${PLATYPUS_OUTD}/finished        
}

########
execute_platypus_germline_local()
{
    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_outd=$3

    # Initialize variables
    PLATYPUS_OUTD=`get_step_dirname ${outd} platypus_germline_local`
    
    # Run Platypus
    python ${PLATYPUS_BUILD_DIR}/bin/Platypus.py callVariants --bamFiles=${local_normalbam} --refFile=${local_ref} --output=${PLATYPUS_OUTD}/output.vcf --verbosity=1 > ${PLATYPUS_OUTD}/platypus.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${PLATYPUS_OUTD}/finished    
}

########
get_pars_platypus_germline()
{
    echo "$ref $normalbam $outd"
}

########
execute_platypus_germline()
{
    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_outd=$3

    if [ -z "${PLATYPUS_BUILD_DIR}" ]; then
        execute_platypus_germline_conda ${local_ref} ${local_normalbam} ${local_outd}
    else
        execute_platypus_germline_local ${local_ref} ${local_normalbam} ${local_outd}
    fi
}

########
get_pars_cnvkit()
{
    echo "$ref $normalbam $tumorbam $outd $cpus"
}

########
execute_cnvkit()
{
    # Initialize variables
    local_ref=$1
    local_normalbam=$2
    local_tumorbam=$3
    local_outd=$4
    local_cpus=$5

    # Initialize variables
    CNVKIT_OUTD=`get_step_dirname ${outd} cnvkit`
    
    # Activate conda environment
    conda activate cnvkit

    # Run cnvkit
    cnvkit.py batch ${local_tumorbam} -n ${local_normalbam} -m wgs -f ${local_ref}  -d ${CNVKIT_OUTD} -p ${local_cpus} > ${CNVKIT_OUTD}/cnvkit.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${CNVKIT_OUTD}/finished
}

########
get_pars_download_ega_norm_bam()
{
    echo "$normalbam ${egaid_normalbam} $egastr $egacred $outd"
}

########
execute_download_ega_norm_bam()
{
    # Initialize variables
    local_normalbam=$1
    local_egaid_normalbam=$2
    local_egastr=$3
    local_egacred=$4
    local_outd=$5

    # Initialize variables
    DOWNLOAD_EGA_NORM_BAM_OUTD=`get_step_dirname ${local_outd} download_ega_norm_bam`

    # Download file
    ${PYEGA_BUILD_DIR}/pyega3 -c ${local_egastr} -cf ${local_egacred} fetch ${local_egaid_normalbam} ${local_normalbam} > ${DOWNLOAD_EGA_NORM_BAM_OUTD}/pyega3.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${DOWNLOAD_EGA_NORM_BAM_OUTD}/finished
}

########
get_pars_download_ega_tum_bam()
{
    echo "$tumorbam ${egaid_tumorbam} $egastr $egacred $outd"
}

########
execute_download_ega_tum_bam()
{
    # Initialize variables
    local_tumorbam=$1
    local_egaid_tumorbam=$2
    local_egastr=$3
    local_egacred=$4
    local_outd=$5

    # Initialize variables
    DOWNLOAD_EGA_TUM_BAM_OUTD=`get_step_dirname ${outd} download_ega_tum_bam`

    # Download file
    ${PYEGA_BUILD_DIR}/pyega3 -c ${local_egastr} -cf ${local_egacred} fetch ${local_egaid_tumorbam} ${local_tumorbam} > ${DOWNLOAD_EGA_TUM_BAM_OUTD}/pyega3.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${DOWNLOAD_EGA_TUM_BAM_OUTD}/finished
}

########
get_pars_index_norm_bam()
{
    echo "$normalbam $outd"
}

########
execute_index_norm_bam()
{
    # Initialize variables
    local_normalbam=$1
    local_outd=$2

    # Initialize variables
    INDEX_NORM_BAM_OUTD=`get_step_dirname ${outd} index_norm_bam`

    # Index normal bam file if necessary
    if [ ! -f ${local_normalbam}.bai ]; then

        # Activate conda environment
        conda activate base

        # Execute samtools
        samtools index ${local_normalbam} > ${INDEX_NORM_BAM_OUTD}/samtools.log 2>&1 || exit 1

        # Deactivate conda environment
        conda deactivate
    fi
    
    # Create file indicating that execution was finished
    touch ${INDEX_NORM_BAM_OUTD}/finished
}

########
get_pars_index_tum_bam()
{
    echo "$tumorbam $outd"
}

########
execute_index_tum_bam()
{
    # Initialize variables
    local_tumorbam=$1
    local_outd=$2

    # Initialize variables
    INDEX_TUM_BAM_OUTD=`get_step_dirname ${outd} index_tum_bam`

    # Index tumor bam file if necessary
    if [ ! -f ${local_tumorbam}.bai ]; then

        # Activate conda environment
        conda activate base

        # Execute samtools
        samtools index ${local_tumorbam} > ${INDEX_TUM_BAM_OUTD}/samtools.log 2>&1 || exit 1

        # Deactivate conda environment
        conda deactivate
    fi
    
    # Create file indicating that execution was finished
    touch ${INDEX_TUM_BAM_OUTD}/finished
}

########
get_pars_delete_bam_files()
{
    echo "$normalbam $tumorbam $outd"
}

########
execute_delete_bam_files()
{
    # Initialize variables
    local_normalbam=$1
    local_tumorbam=$2
    local_outd=$3

    # Initialize variables
    DELETE_BAM_FILES_OUTD=`get_step_dirname ${outd} delete_bam_files`

    # Delete normal bam file
    rm ${local_normalbam} > ${DELETE_BAM_FILES_OUTD}/rm_norm.log 2>&1 || exit 1
    
    # Delete tumor bam file
    rm ${local_tumorbam} > ${DELETE_BAM_FILES_OUTD}/rm_tum.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${DELETE_BAM_FILES_OUTD}/finished
}

########
get_jobdeps()
{
    local_jobdeps_spec=$1
    case ${local_jobdeps_spec} in
            "all") echo ${step_jids}
                   ;;
            "none") echo ""
                    ;;
            *) local_jdeps=""
               for dep in `echo ${local_jobdeps_spec} | $SED 's/:/ /g'`; do
                   dep_jid=${dep}_jid
                   local_jdeps=${local_jdeps}":"${!dep_jid}
               done
               echo ${local_jdeps}
               ;;
    esac
}

########
execute_step()
{
    # Initialize variables
    local_dirname=$1
    local_stepname=$2
    local_cpus=$3
    local_mem=$4
    local_time=$5
    local_jobdeps_spec=$6

    # Execute step
    local_script_pars=`get_pars_${local_stepname}`
    create_script ${local_dirname}/scripts/execute_${local_stepname} execute_${local_stepname} "${local_script_pars}"
    local_status=`${bindir}/get_analysis_status -d ${local_dirname} -s "${local_stepname}"`
    echo "STEP: ${local_stepname} ; STATUS: ${local_status}" >&2
    if [ "${local_status}" != "FINISHED" ]; then
        reset_outdir_for_step ${local_dirname} ${local_stepname} || exit 1
        local_jobdeps="`get_jobdeps ${local_jobdeps_spec}`"
        local_stepname_jid=${local_stepname}_jid
        launch ${local_dirname}/scripts/execute_${local_stepname} ${local_cpus} ${local_mem} ${local_time} "${local_jobdeps}" ${local_stepname_jid}
        
        # Update variables storing jids
        step_jids="${step_jids}:${!local_stepname_jid}"
    fi
}

########
execute_steps_in_afile()
{
    # Read input parameters
    local_dirname=$1
    local_afile=$2

    # step_jids will store the job ids of the analysis steps
    step_jids=""
    
    # Read information about the steps to be executed
    while read entry; do
        entry_ok=`entry_is_ok "$entry"`
        if [ ${entry_ok} = "yes" ]; then
            # Extract entry information
            local_stepname=`extract_stepname_from_entry "$entry"`
            cpus=`extract_cpus_from_entry "$entry"`
            local_mem=`extract_mem_from_entry "$entry"`
            local_time=`extract_time_from_entry "$entry"`
            local_jobdeps_spec=`extract_jobdeps_spec_from_entry "$entry"`

            # Execute step
            execute_step ${local_dirname} ${local_stepname} ${cpus} ${local_mem} ${local_time} ${local_jobdeps_spec} || exit 1
        fi
    done < ${local_afile}
}

########

if [ $# -eq 0 ]; then
    print_desc
    exit 1
fi

read_pars $@ || exit 1

check_pars || exit 1

print_pars || exit 1

create_dirs || exit 1

set_bam_filenames || exit 1

execute_steps_in_afile ${outd} ${afile} || exit 1
