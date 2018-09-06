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
    echo "submit_bam_analysis       -r <string> -n <string> -t <string> -o <string>"
    echo "                          -a <string>"
    echo "                          [-debug] [--help]"
    echo ""
    echo "-r <string>               File with reference genome."
    echo "-n <string>               Normal bam file."
    echo "-t <string>               Tumor bam file."
    echo "-o <string>               Output directory."
    echo "-a <string>               File with analysis steps to be performed."
    echo "                          Expected format:"
    echo "                          <stepname> <cpus> <mem> <time>"
    echo "-debug                    After ending, do not delete temporary files"
    echo "                          (for debugging purposes)."
    echo "--help                    Display this help and exit."
}

########
read_pars()
{
    r_given=0
    n_given=0
    t_given=0
    o_given=0
    a_given=0
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
            "-o") shift
                  if [ $# -ne 0 ]; then
                      outd=$1
                      o_given=1
                  fi
                  ;;
            "-a") shift
                  if [ $# -ne 0 ]; then
                      afile=$1
                      a_given=1
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

    if [ ${n_given} -eq 0 ]; then   
        echo "Error! -n parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${normalbam} ]; then
            echo "Error! file ${normalbam} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${t_given} -eq 0 ]; then
        echo "Error! -t parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${tumorbam} ]; then
            echo "Error! file ${tumorbam} does not exist" >&2
            exit 1
        fi
    fi

    if [ ${o_given} -eq 0 ]; then
        echo "Error! -o parameter not given!" >&2
        exit 1
    else
        if [ -d ${outd} ]; then
            echo "Warning! output directory does exist" >&2 
        fi
    fi

    if [ ${a_given} -eq 0 ]; then   
        echo "Error! -s parameter not given!" >&2
        exit 1
    else
        if [ ! -f ${afile} ]; then
            echo "Error! file ${afile} does not exist" >&2
            exit 1
        fi
    fi
}

########
create_dirs()
{
    mkdir -p ${outd} || { echo "Error! cannot create output directory" >&2; return 1; }

    mkdir -p ${outd}/scripts || { echo "Error! cannot create scripts directory" >&2; return 1; }
}

########
execute_manta_somatic()
{
    # WARNING: for this function to work successfully, it is necessary
    # to define the following global variables:
    #
    # - ref: reference genome file
    # - normalbam: bam file for normal sample
    # - tumorbam: bam file for tumor sample
    # - stepname: name of step to be executed
    # - outd: output directory
    # - cpus: number of cpus used to execute the step
    
    # Initialize variables
    MANTA_OUTD=`get_step_dirname ${outd} ${stepname}`

    # Activate conda environment
    conda activate manta
    
    # Configure Manta
    configManta.py --normalBam ${normalbam} --tumorBam ${tumorbam} --referenceFasta ${ref} --runDir ${MANTA_OUTD} > ${MANTA_OUTD}/configManta.log 2>&1 || exit 1

    # Execute Manta
    ${MANTA_OUTD}/runWorkflow.py -m local -j ${cpus} > ${MANTA_OUTD}/runWorkflow.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${MANTA_OUTD}/finished
}

########
execute_strelka_somatic()
{
    # WARNING: for this function to work successfully, it is necessary
    # to define the following global variables:
    #
    # - ref: reference genome file
    # - normalbam: bam file for normal sample
    # - tumorbam: bam file for tumor sample
    # - stepname: name of step to be executed
    # - outd: output directory
    # - cpus: number of cpus used to execute the step

    # Initialize variables
    STRELKA_OUTD=`get_step_dirname ${outd} ${stepname}`

    # Activate conda environment
    conda activate strelka

    # Configure Strelka
    configureStrelkaSomaticWorkflow.py --normalBam ${normalbam} --tumorBam ${tumorbam} --referenceFasta ${ref} --runDir ${STRELKA_OUTD} > ${STRELKA_OUTD}/configureStrelkaSomaticWorkflow.log 2>&1 || exit 1

    # Execute Strelka
    ${STRELKA_OUTD}/runWorkflow.py -m local -j ${cpus} > ${STRELKA_OUTD}/runWorkflow.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${STRELKA_OUTD}/finished
}

########
execute_msisensor()
{
    # WARNING: for this function to work successfully, it is necessary
    # to define the following global variables:
    #
    # - ref: reference genome file
    # - normalbam: bam file for normal sample
    # - tumorbam: bam file for tumor sample
    # - stepname: name of step to be executed
    # - outd: output directory
    # - cpus: number of cpus used to execute the step

    # Initialize variables
    MSISENSOR_OUTD=`get_step_dirname ${outd} ${stepname}`
    
    # Activate conda environment
    conda activate msisensor

    # Create homopolymer and microsatellites file
    msisensor scan -d ${ref} -o ${MSISENSOR_OUTD}/msisensor.list > ${MSISENSOR_OUTD}/msisensor_scan.log 2>&1 || exit 1

    # Run MSIsensor analysis
    msisensor msi -d ${MSISENSOR_OUTD}/msisensor.list -n ${normalbam} -t ${tumorbam} -o ${MSISENSOR_OUTD}/output -l 1 -q 1 -b ${cpus} > ${MSISENSOR_OUTD}/msisensor_msi.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate
    
    # Create file indicating that execution was finished
    touch ${MSISENSOR_OUTD}/finished
}

########
execute_platypus_germline_conda()
{
    # Initialize variables
    PLATYPUS_OUTD=`get_step_dirname ${outd} ${stepname}`

    # Activate conda environment
    conda activate platypus

    # Run Platypus
    Platypus.py callVariants --bamFiles=${normalbam} --refFile=${ref} --output=$3 --verbosity=1 > ${PLATYPUS_OUTD}/platypus.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${PLATYPUS_OUTD}/finished        
}

########
execute_platypus_germline_local()
{
    # Initialize variables
    PLATYPUS_OUTD=`get_step_dirname ${outd} ${stepname}`
    
    # Run Platypus
    python ${PLATYPUS_BUILD_DIR}/bin/Platypus.py callVariants --bamFiles=${normalbam} --refFile=${ref} --output=$3 --verbosity=1 > ${PLATYPUS_OUTD}/platypus.log 2>&1 || exit 1

    # Create file indicating that execution was finished
    touch ${PLATYPUS_OUTD}/finished    
}

########
execute_platypus_germline()
{
    # WARNING: for this function to work successfully, it is necessary
    # to define the following global variables:
    #
    # - ref: reference genome file
    # - normalbam: bam file for normal sample
    # - stepname: name of step to be executed
    # - outd: output directory

    if [ -z "${PLATYPUS_BUILD_DIR}" ]; then
        execute_platypus_germline_conda
    else
        execute_platypus_germline_local
    fi
}

########
execute_cnvkit()
{
    # WARNING: for this function to work successfully, it is necessary
    # to define the following global variables:
    #
    # - ref: reference genome file
    # - normalbam: bam file for normal sample
    # - tumorbam: bam file for tumor sample
    # - stepname: name of step to be executed
    # - outd: output directory
    # - cpus: number of cpus used to execute the step

    # Initialize variables
    CNVKIT_OUTD=`get_step_dirname ${outd} ${stepname}`
    
    # Activate conda environment
    conda activate cnvkit

    # Run cnvkit
    cnvkit.py batch ${tumorbam} -n ${normalbam} -m wgs -f ${ref}  -d ${CNVKIT_OUTD} -p ${cpus} > ${CNVKIT_OUTD}/cnvkit.log 2>&1 || exit 1

    # Deactivate conda environment
    conda deactivate

    # Create file indicating that execution was finished
    touch ${CNVKIT_OUTD}/finished
}

########
execute_pre_analysis_actions()
{
    # Initialize variables
    PRE_ANALYSIS_ACTIONS_OUTD=`get_step_dirname ${outd} ${stepname}`

    # TBD
    echo "TBD"
    sleep 20

    # Create file indicating that execution was finished
    touch ${PRE_ANALYSIS_ACTIONS_OUTD}/finished
}

########
execute_post_analysis_actions()
{
    # Initialize variables
    POST_ANALYSIS_ACTIONS_OUTD=`get_step_dirname ${outd} ${stepname}`

    # TBD
    echo "TBD"
    sleep 10

    # Create file indicating that execution was finished
    touch ${POST_ANALYSIS_ACTIONS_OUTD}/finished
}

########
get_local_jobdeps()
{
    local_mutex=$1
    if [ ${local_mutex} = "yes" ]; then
        echo ${step_jids}
    else
        echo ${mutex_jids}
    fi
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
    local_mutex=$6

    # Execute step
    create_script ${local_dirname}/scripts/execute_${local_stepname} execute_${local_stepname}
    status=`${bindir}/get_analysis_status -d ${local_dirname} -s "${local_stepname}"`
    echo "STEP: ${local_stepname} ; STATUS: ${status}" >&2
    if [ "$status" != "FINISHED" ]; then
        reset_outdir_for_step ${local_dirname} ${local_stepname} || exit 1
        local_jobdeps="`get_local_jobdeps ${local_mutex}`"
        launch ${local_dirname}/scripts/execute_${local_stepname} ${local_cpus} ${local_mem} ${local_time} "${local_jobdeps}" job_id
        
        # Update variables storing jids
        step_jids="${step_jids}:${job_id}"
        if [ ${local_mutex} = "yes" ]; then
            mutex_jids="${mutex_jids}:${job_id}"
        fi
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

    # mutex_jids will store the job ids of those analysis steps that are
    # mutually exclusive
    mutex_jids=""
    
    # Read information about the steps to be executed
    while read entry; do
        entry_ok=`entry_is_ok "$entry"`
        if [ ${entry_ok} = "yes" ]; then
            # Extract entry information
            stepname=`extract_stepname_from_entry "$entry"`
            cpus=`extract_cpus_from_entry "$entry"`
            mem=`extract_mem_from_entry "$entry"`
            time=`extract_time_from_entry "$entry"`
            mutex=`extract_mutex_from_entry "$entry"`

            # Execute step
            execute_step ${local_dirname} ${stepname} ${cpus} ${mem} ${time} ${mutex} || exit 1
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

create_dirs || exit 1

execute_steps_in_afile ${outd} ${afile} || exit 1
